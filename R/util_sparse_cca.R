##################################################
## Sparse Canonical Correlation Analysis (CCA)
## Fast method for 2-omics integration
## Uses PMA package (Pure R, CRAN)
## Reference: Witten et al., Biostatistics 2009
## Created: 2025-01-04
##################################################

#' Perform Sparse CCA for Two Omics Layers
#'
#' Fast integration of two omics using sparse canonical correlation analysis.
#' Identifies correlated features between omics layers (10-20 seconds).
#'
#' @param omics1 First omics matrix (samples x features)
#' @param omics2 Second omics matrix (samples x features)
#' @param n_components Number of canonical components (default: 3)
#' @param penaltyx Sparsity penalty for omics1 (0-1, default: NULL = auto-select)
#' @param penaltyz Sparsity penalty for omics2 (0-1, default: NULL = auto-select)
#' @param nperms Number of permutations for penalty selection (default: 25)
#'               Lower values = faster but less reliable (min: 5, recommended: 10-25)
#' @param center Whether to center features (default: TRUE)
#' @param scale Whether to scale features (default: TRUE)
#' @param progress_callback Optional progress callback function
#'
#' @return List containing:
#'   \item{canonical_variates}{List with u (omics1) and v (omics2) variates}
#'   \item{correlations}{Canonical correlations for each component}
#'   \item{weights}{List with omics1 and omics2 feature weights}
#'   \item{selected_features}{Features with non-zero weights (sparse solution)}
#'   \item{variance_explained}{Variance explained in each omics}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Integrate transcriptomics and proteomics
#' result <- PerformSparseCCA(
#'   omics1 = transcriptomics_data,
#'   omics2 = proteomics_data,
#'   n_components = 3,
#'   penaltyx = 0.3,  # 30% of features selected
#'   penaltyz = 0.3
#' )
#'
#' # Plot canonical variates
#' plot(result$canonical_variates$u[,1],
#'      result$canonical_variates$v[,1])
#' }
PerformSparseCCA <- function(omics1,
                             omics2,
                             n_components = 3,
                             penaltyx = NULL,
                             penaltyz = NULL,
                             nperms = 25,
                             center = TRUE,
                             scale = TRUE,
                             progress_callback = NULL) {

  # Progress helper
  update_progress <- if (!is.null(progress_callback)) {
    progress_callback
  } else {
    function(pct, msg) NULL
  }

  update_progress(5, "Validating input data")

  # Input validation
  if (nrow(omics1) != nrow(omics2)) {
    stop("omics1 and omics2 must have same number of samples (rows)")
  }

  n_samples <- nrow(omics1)

  # Check if PMA package is available
  if (!requireNamespace("PMA", quietly = TRUE)) {
    stop("Package 'PMA' is required. Install with: install.packages('PMA')")
  }

  update_progress(10, "Preprocessing data")

  # Preprocessing
  if (center || scale) {
    omics1 <- scale(omics1, center = center, scale = scale)
    omics2 <- scale(omics2, center = center, scale = scale)
  }

  # Remove zero-variance features
  var1 <- apply(omics1, 2, var, na.rm = TRUE)
  var2 <- apply(omics2, 2, var, na.rm = TRUE)

  valid1 <- !is.na(var1) & var1 > 0
  valid2 <- !is.na(var2) & var2 > 0

  omics1 <- omics1[, valid1, drop = FALSE]
  omics2 <- omics2[, valid2, drop = FALSE]

  update_progress(20, "Selecting optimal sparsity penalties")

  # Auto-select penalties using cross-validation if not provided
  if (is.null(penaltyx) || is.null(penaltyz)) {
    update_progress(20, sprintf("Auto-selecting penalties (%d permutations)", nperms))
    cv_result <- PMA::CCA.permute(
      x = omics1,
      z = omics2,
      typex = "standard",
      typez = "standard",
      nperms = nperms  # User-configurable (default: 25)
    )

    if (is.null(penaltyx)) {
      penaltyx <- cv_result$bestpenaltyx
    }
    if (is.null(penaltyz)) {
      penaltyz <- cv_result$bestpenaltyz
    }

    message(sprintf("Auto-selected penalties: omics1=%.3f, omics2=%.3f",
                   penaltyx, penaltyz))
  }

  update_progress(40, sprintf("Running sparse CCA (%d components)", n_components))

  # Perform sparse CCA
  cca_result <- PMA::CCA(
    x = omics1,
    z = omics2,
    K = n_components,
    penaltyx = penaltyx,
    penaltyz = penaltyz,
    typex = "standard",
    typez = "standard",
    trace = FALSE
  )

  update_progress(70, "Extracting canonical variates")

  # PMA::CCA returns u and v as WEIGHTS (features x components), not scores
  # We need to compute sample scores by projecting the data
  # Canonical variates (scores) = data %*% weights

  # Ensure u and v are matrices
  if (is.null(cca_result$u) || is.null(cca_result$v)) {
    stop("CCA did not return canonical weights (u/v)")
  }

  # Force conversion to matrix
  u_weights <- as.matrix(cca_result$u)
  v_weights <- as.matrix(cca_result$v)

  # Ensure correct dimensions: features x components
  if (ncol(u_weights) != n_components) {
    u_weights <- matrix(u_weights, ncol = n_components)
  }
  if (ncol(v_weights) != n_components) {
    v_weights <- matrix(v_weights, ncol = n_components)
  }

  u_scores <- omics1 %*% u_weights  # samples x components
  v_scores <- omics2 %*% v_weights  # samples x components

  canonical_variates <- list(
    u = u_scores,  # Canonical variates (scores) for omics1
    v = v_scores   # Canonical variates (scores) for omics2
  )

  # Set names
  rownames(canonical_variates$u) <- rownames(omics1)
  rownames(canonical_variates$v) <- rownames(omics2)
  colnames(canonical_variates$u) <- paste0("CC", 1:n_components, "_omics1")
  colnames(canonical_variates$v) <- paste0("CC", 1:n_components, "_omics2")

  update_progress(80, "Computing correlations and feature weights")

  # Canonical correlations
  correlations <- cca_result$cors

  # Feature weights (loadings)
  # PMA::CCA returns u/v (standardized) and ws/vs (original scale)
  # Use ws/vs if available, otherwise fall back to u/v
  if (!is.null(cca_result$ws) && !is.null(cca_result$vs)) {
    ws_mat <- as.matrix(cca_result$ws)
    vs_mat <- as.matrix(cca_result$vs)
  } else {
    # Fall back to u/v if ws/vs not available
    ws_mat <- u_weights  # Already computed above
    vs_mat <- v_weights
  }

  # Ensure correct dimensions: features x components
  if (ncol(ws_mat) != n_components) {
    ws_mat <- matrix(ws_mat, ncol = n_components)
  }
  if (ncol(vs_mat) != n_components) {
    vs_mat <- matrix(vs_mat, ncol = n_components)
  }

  weights <- list(
    omics1 = ws_mat,  # Sparse weights for omics1 features
    omics2 = vs_mat   # Sparse weights for omics2 features
  )

  # Add feature names
  if (!is.null(colnames(omics1))) {
    rownames(weights$omics1) <- colnames(omics1)
  }
  if (!is.null(colnames(omics2))) {
    rownames(weights$omics2) <- colnames(omics2)
  }

  update_progress(90, "Identifying selected features")

  # Selected features (non-zero weights)
  selected_features <- list(
    omics1 = lapply(1:n_components, function(i) {
      feat_idx <- which(weights$omics1[, i] != 0)
      if (!is.null(colnames(omics1))) {
        colnames(omics1)[feat_idx]
      } else {
        feat_idx
      }
    }),
    omics2 = lapply(1:n_components, function(i) {
      feat_idx <- which(weights$omics2[, i] != 0)
      if (!is.null(colnames(omics2))) {
        colnames(omics2)[feat_idx]
      } else {
        feat_idx
      }
    })
  )

  names(selected_features$omics1) <- paste0("CC", 1:n_components)
  names(selected_features$omics2) <- paste0("CC", 1:n_components)

  # Variance explained
  var_explained_omics1 <- apply(canonical_variates$u, 2, var) / sum(apply(omics1, 2, var))
  var_explained_omics2 <- apply(canonical_variates$v, 2, var) / sum(apply(omics2, 2, var))

  variance_explained <- list(
    omics1 = var_explained_omics1 * 100,
    omics2 = var_explained_omics2 * 100
  )

  update_progress(95, "Finalizing results")

  # Package results
  result <- list(
    canonical_variates = canonical_variates,
    correlations = correlations,
    weights = weights,
    selected_features = selected_features,
    variance_explained = variance_explained,
    penalties = list(omics1 = penaltyx, omics2 = penaltyz),
    n_components = n_components,
    method = "Sparse CCA"
  )

  class(result) <- c("SparseCCA", "list")

  update_progress(100, "Sparse CCA complete")

  return(result)
}


#' Perform Multi-Omics CCA (>2 omics)
#'
#' Sequential sparse CCA for multiple omics layers.
#' Anchors first omics and aligns others to it.
#'
#' @param omics_list Named list of omics matrices
#' @param n_components Number of canonical components
#' @param ... Additional arguments passed to PerformSparseCCA
#' @return List of pairwise CCA results
#' @export
PerformMultiOmicsCCA <- function(omics_list, n_components = 3, ...) {

  if (length(omics_list) < 2) {
    stop("Need at least 2 omics layers")
  }

  if (is.null(names(omics_list))) {
    names(omics_list) <- paste0("Omics", seq_along(omics_list))
  }

  message("Performing sequential CCA: anchoring to ", names(omics_list)[1])

  # Anchor to first omics
  base_omics <- omics_list[[1]]
  base_name <- names(omics_list)[1]

  results <- list()

  for (i in 2:length(omics_list)) {
    message(sprintf("CCA: %s vs %s", base_name, names(omics_list)[i]))

    cca_result <- PerformSparseCCA(
      omics1 = base_omics,
      omics2 = omics_list[[i]],
      n_components = n_components,
      ...
    )

    results[[names(omics_list)[i]]] <- cca_result
  }

  # Aggregate results
  result <- list(
    pairwise_results = results,
    base_omics = base_name,
    n_omics = length(omics_list),
    omics_names = names(omics_list),
    method = "Multi-Omics Sequential CCA"
  )

  class(result) <- c("MultiOmicsCCA", "list")

  return(result)
}


#' Plot Sparse CCA Results
#'
#' @param cca_result Result from PerformSparseCCA
#' @param components Which components to plot (default: c(1,2))
#' @param groups Optional grouping variable for coloring
#' @param main Plot title
#' @export
PlotSparseCCA <- function(cca_result,
                         components = c(1, 2),
                         groups = NULL,
                         main = "Sparse CCA") {

  if (!inherits(cca_result, "SparseCCA")) {
    stop("Input must be a SparseCCA object")
  }

  u <- cca_result$canonical_variates$u[, components]
  v <- cca_result$canonical_variates$v[, components]
  cors <- cca_result$correlations[components]

  # Set colors
  if (is.null(groups)) {
    colors <- rep("black", nrow(u))
  } else {
    groups <- as.factor(groups)
    colors <- rainbow(nlevels(groups))[as.numeric(groups)]
  }

  # Plot for each component
  par(mfrow = c(1, length(components)))

  for (i in seq_along(components)) {
    comp <- components[i]
    plot(u[, i], v[, i],
         col = colors,
         pch = 19,
         cex = 1.5,
         xlab = sprintf("Omics1 CC%d", comp),
         ylab = sprintf("Omics2 CC%d", comp),
         main = sprintf("%s - CC%d (r=%.3f)", main, comp, cors[i]))

    # Add diagonal reference line
    abline(a = 0, b = 1, lty = 2, col = "gray")

    if (!is.null(groups) && i == 1) {
      legend("topright",
             legend = levels(groups),
             col = rainbow(nlevels(groups)),
             pch = 19,
             cex = 0.6)
    }
  }

  par(mfrow = c(1, 1))
}


#' Plot Feature Weights Heatmap
#'
#' @param cca_result Result from PerformSparseCCA
#' @param top_n Number of top features to show (default: 20)
#' @param component Which component to visualize (default: 1)
#' @export
PlotFeatureWeights <- function(cca_result, top_n = 20, component = 1) {

  if (!inherits(cca_result, "SparseCCA")) {
    stop("Input must be a SparseCCA object")
  }

  weights1 <- cca_result$weights$omics1[, component]
  weights2 <- cca_result$weights$omics2[, component]

  # Get top features by absolute weight
  top_idx1 <- head(order(abs(weights1), decreasing = TRUE), top_n)
  top_idx2 <- head(order(abs(weights2), decreasing = TRUE), top_n)

  top_weights1 <- weights1[top_idx1]
  top_weights2 <- weights2[top_idx2]

  # Create barplots
  par(mfrow = c(2, 1), mar = c(4, 8, 3, 2))

  # Omics 1
  barplot(sort(top_weights1),
          horiz = TRUE,
          las = 1,
          col = ifelse(top_weights1 > 0, "steelblue", "coral"),
          main = sprintf("Top %d Features in Omics1 (CC%d)", top_n, component),
          xlab = "Feature Weight")

  # Omics 2
  barplot(sort(top_weights2),
          horiz = TRUE,
          las = 1,
          col = ifelse(top_weights2 > 0, "steelblue", "coral"),
          main = sprintf("Top %d Features in Omics2 (CC%d)", top_n, component),
          xlab = "Feature Weight")

  par(mfrow = c(1, 1))
}


#' Print Sparse CCA Summary
#'
#' @param x SparseCCA object
#' @param ... Additional arguments (not used)
#' @export
print.SparseCCA <- function(x, ...) {
  cat("=== Sparse Canonical Correlation Analysis ===\n\n")
  cat("Number of canonical components:", x$n_components, "\n")
  cat("Sparsity penalties: omics1=", round(x$penalties$omics1, 3),
      ", omics2=", round(x$penalties$omics2, 3), "\n\n")

  cat("Canonical correlations:\n")
  for (i in 1:x$n_components) {
    n_feat1 <- length(x$selected_features$omics1[[i]])
    n_feat2 <- length(x$selected_features$omics2[[i]])
    cat(sprintf("  CC%d: r=%.3f (features: %d + %d)\n",
                i, x$correlations[i], n_feat1, n_feat2))
  }

  cat("\nVariance explained:\n")
  cat("  Omics1:", paste(sprintf("%.2f%%", x$variance_explained$omics1), collapse=", "), "\n")
  cat("  Omics2:", paste(sprintf("%.2f%%", x$variance_explained$omics2), collapse=", "), "\n")

  invisible(x)
}
