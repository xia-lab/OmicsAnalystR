##################################################
## Consensus PCA for Fast Multi-Omics Integration
## Ultra-fast, deterministic, R-native method
## Reference: ConsensusOPLS concept (J Proteome Res 2022)
## Created: 2025-01-04
##################################################

#' Perform Consensus PCA Integration
#'
#' Fast multi-omics integration using consensus PCA approach.
#' This is the FASTEST method, ideal for web platform (5-10 seconds).
#'
#' @param omics_list Named list of data matrices (samples x features)
#'                   Each matrix must have same number of rows (samples)
#' @param n_components Number of consensus components to extract (default: 5)
#' @param scale Whether to scale features (default: TRUE)
#' @param center Whether to center features (default: TRUE)
#' @param weights Optional named vector of omics weights (default: equal weights)
#' @param progress_callback Optional function(percent, message) for progress updates
#'
#' @return List containing:
#'   \item{consensus_scores}{Matrix of consensus PC scores (samples x components)}
#'   \item{omics_scores}{List of omics-specific PC scores}
#'   \item{loadings}{List of feature loadings per omics}
#'   \item{variance_explained}{Variance explained by each consensus component}
#'   \item{omics_contribution}{Contribution of each omics to consensus}
#'   \item{procrustes_stats}{Procrustes alignment statistics}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example data
#' set.seed(123)
#' omics_list <- list(
#'   transcriptomics = matrix(rnorm(100*1000), 100, 1000),
#'   proteomics = matrix(rnorm(100*500), 100, 500),
#'   metabolomics = matrix(rnorm(100*200), 100, 200)
#' )
#'
#' # Run consensus PCA
#' result <- PerformConsensusPCA(omics_list, n_components = 5)
#'
#' # Plot consensus scores
#' plot(result$consensus_scores[,1:2], col = your_groups)
#' }
PerformConsensusPCA <- function(omics_list,
                                n_components = 5,
                                scale = TRUE,
                                center = TRUE,
                                weights = NULL,
                                progress_callback = NULL) {

  # Progress callback helper
  update_progress <- if (!is.null(progress_callback)) {
    progress_callback
  } else {
    function(pct, msg) NULL
  }

  update_progress(5, "Validating input data")

  # Input validation
  if (!is.list(omics_list) || length(omics_list) < 2) {
    stop("omics_list must be a list with at least 2 omics layers")
  }

  n_samples <- unique(sapply(omics_list, nrow))
  if (length(n_samples) > 1) {
    stop("All omics must have the same number of samples (rows)")
  }
  n_samples <- n_samples[1]

  # Set omics names if not provided
  if (is.null(names(omics_list))) {
    names(omics_list) <- paste0("Omics", seq_along(omics_list))
  }
  omics_names <- names(omics_list)

  # Set equal weights if not provided
  if (is.null(weights)) {
    weights <- rep(1, length(omics_list))
    names(weights) <- omics_names
  }

  # Normalize weights
  weights <- weights / sum(weights)

  update_progress(15, "Computing PCA for each omics layer")

  # Step 1: Perform individual PCA for each omics
  pca_list <- list()
  for (i in seq_along(omics_list)) {
    omics_name <- omics_names[i]

    # Remove zero-variance features
    omics_data <- omics_list[[i]]
    feature_vars <- apply(omics_data, 2, var, na.rm = TRUE)
    valid_features <- !is.na(feature_vars) & feature_vars > 0
     omics_data <- omics_data[, valid_features, drop = FALSE]

    if (ncol(omics_data) < n_components) {
      warning(sprintf("%s has fewer features than components. Reducing n_components.", omics_name))
      n_comp_use <- min(n_components, ncol(omics_data), nrow(omics_data) - 1)
    } else {
      n_comp_use <- min(n_components, nrow(omics_data) - 1)
    }

    # Perform PCA
    pca_result <- prcomp(omics_data,
                        center = center,
                        scale. = scale,
                        rank. = n_comp_use)

    pca_list[[omics_name]] <- list(
      scores = pca_result$x[, 1:n_comp_use, drop = FALSE],
      loadings = pca_result$rotation[, 1:n_comp_use, drop = FALSE],
      variance = pca_result$sdev[1:n_comp_use]^2,
      variance_pct = summary(pca_result)$importance[2, 1:n_comp_use] * 100
    )

    update_progress(15 + (i / length(omics_list)) * 30,
                   sprintf("PCA complete for %s", omics_name))
  }

  update_progress(50, "Aligning omics to consensus space")

  # Step 2: Procrustes alignment to consensus
  # Use first omics as initial reference, then iterate
  reference_scores <- pca_list[[1]]$scores

  aligned_scores <- list()
  procrustes_stats <- list()

  for (i in seq_along(pca_list)) {
    omics_name <- omics_names[i]

    if (i == 1) {
      # First omics is reference
      aligned_scores[[omics_name]] <- reference_scores
      procrustes_stats[[omics_name]] <- list(m12 = 0, ss = 0)
    } else {
      # Align subsequent omics to reference using Procrustes
      proc_result <- ProcrustesAlign(
        target = reference_scores,
        source = pca_list[[omics_name]]$scores
      )

      aligned_scores[[omics_name]] <- proc_result$aligned
      procrustes_stats[[omics_name]] <- proc_result$stats
    }
  }

  update_progress(70, "Computing consensus scores")

  # Step 3: Compute weighted consensus scores
  # Simple weighted average of aligned scores
  consensus_scores <- Reduce("+", lapply(seq_along(aligned_scores), function(i) {
    aligned_scores[[i]] * weights[i]
  }))

  rownames(consensus_scores) <- rownames(omics_list[[1]])
  colnames(consensus_scores) <- paste0("ConsensusPC", 1:ncol(consensus_scores))

  update_progress(85, "Computing variance explained")

  # Step 4: Compute consensus variance explained
  consensus_variance <- apply(consensus_scores, 2, var)
  total_variance <- sum(consensus_variance)
  variance_explained <- (consensus_variance / total_variance) * 100

  # Step 5: Compute omics contributions to each component
  omics_contributions <- matrix(0, nrow = length(omics_list), ncol = n_components)
  rownames(omics_contributions) <- omics_names
  colnames(omics_contributions) <- paste0("PC", 1:n_components)

  for (i in seq_along(aligned_scores)) {
    n_comp <- ncol(aligned_scores[[i]])
    for (j in 1:n_comp) {
      # Contribution measured by correlation with consensus
      omics_contributions[i, j] <- cor(aligned_scores[[i]][, j],
                                       consensus_scores[, j])^2
    }
  }

  update_progress(95, "Finalizing results")

  # Package results
  result <- list(
    consensus_scores = consensus_scores,
    omics_scores = aligned_scores,
    loadings = lapply(pca_list, function(x) x$loadings),
    variance_explained = variance_explained,
    omics_contribution = omics_contributions,
    procrustes_stats = procrustes_stats,
    individual_variance = lapply(pca_list, function(x) x$variance_pct),
    weights = weights,
    n_components = n_components,
    method = "Consensus PCA"
  )

  class(result) <- c("ConsensusPCA", "list")

  update_progress(100, "Consensus PCA complete")

  return(result)
}


#' Procrustes Alignment
#'
#' Align source matrix to target matrix using Procrustes rotation
#'
#' @param target Target matrix (n x k)
#' @param source Source matrix to align (n x k)
#' @return List with aligned matrix and statistics
#' @keywords internal
ProcrustesAlign <- function(target, source) {

  # Ensure same dimensions
  n_comp <- min(ncol(target), ncol(source))
  target <- target[, 1:n_comp, drop = FALSE]
  source <- source[, 1:n_comp, drop = FALSE]

  # Center matrices
  target_centered <- scale(target, scale = FALSE)
  source_centered <- scale(source, scale = FALSE)

  # Find optimal rotation matrix using SVD
  # target = source * R  (solve for R)
  svd_result <- svd(t(source_centered) %*% target_centered)
  rotation_matrix <- svd_result$u %*% t(svd_result$v)

  # Apply rotation
  aligned <- source_centered %*% rotation_matrix

  # Compute Procrustes statistics
  ss_total <- sum(target_centered^2)
  ss_residual <- sum((target_centered - aligned)^2)
  m12 <- sqrt(ss_residual / ss_total)  # Procrustes distance

  return(list(
    aligned = aligned,
    rotation = rotation_matrix,
    stats = list(
      m12 = m12,           # Procrustes distance (0 = perfect alignment)
      ss = ss_residual     # Sum of squared residuals
    )
  ))
}


#' Plot Consensus PCA Results
#'
#' @param consensus_pca_result Result from PerformConsensusPCA
#' @param groups Optional factor/vector of sample groups for coloring
#' @param components Which components to plot (default: c(1,2))
#' @param main Plot title
#' @export
PlotConsensusPCA <- function(consensus_pca_result,
                            groups = NULL,
                            components = c(1, 2),
                            main = "Consensus PCA") {

  if (!inherits(consensus_pca_result, "ConsensusPCA")) {
    stop("Input must be a ConsensusPCA object")
  }

  scores <- consensus_pca_result$consensus_scores[, components]
  var_exp <- consensus_pca_result$variance_explained[components]

  # Set up plotting
  if (is.null(groups)) {
    colors <- rep("black", nrow(scores))
  } else {
    groups <- as.factor(groups)
    colors <- rainbow(nlevels(groups))[as.numeric(groups)]
  }

  # Create plot
  plot(scores,
       col = colors,
       pch = 19,
       cex = 1.5,
       xlab = sprintf("Consensus PC%d (%.1f%%)", components[1], var_exp[1]),
       ylab = sprintf("Consensus PC%d (%.1f%%)", components[2], var_exp[2]),
       main = main)

  # Add legend if groups provided
  if (!is.null(groups)) {
    legend("topright",
           legend = levels(groups),
           col = rainbow(nlevels(groups)),
           pch = 19,
           cex = 0.8)
  }

  grid()
}


#' Plot Omics Contributions
#'
#' @param consensus_pca_result Result from PerformConsensusPCA
#' @param components Which components to show (default: all)
#' @export
PlotOmicsContributions <- function(consensus_pca_result, components = NULL) {

  if (!inherits(consensus_pca_result, "ConsensusPCA")) {
    stop("Input must be a ConsensusPCA object")
  }

  contrib <- consensus_pca_result$omics_contribution

  if (is.null(components)) {
    components <- 1:ncol(contrib)
  }

  contrib_subset <- contrib[, components, drop = FALSE]

  # Create barplot
  barplot(contrib_subset,
          beside = TRUE,
          col = rainbow(nrow(contrib)),
          legend.text = rownames(contrib),
          args.legend = list(x = "topright", cex = 0.8),
          xlab = "Consensus Component",
          ylab = "Contribution (R²)",
          main = "Omics Layer Contributions to Consensus PCs",
          ylim = c(0, 1))

  grid()
}


#' Print Summary of Consensus PCA
#'
#' @param x ConsensusPCA object
#' @param ... Additional arguments (not used)
#' @export
print.ConsensusPCA <- function(x, ...) {
  cat("=== Consensus PCA Integration ===\n\n")
  cat("Number of omics layers:", length(x$omics_scores), "\n")
  cat("Omics names:", paste(names(x$omics_scores), collapse = ", "), "\n")
  cat("Number of samples:", nrow(x$consensus_scores), "\n")
  cat("Number of components:", x$n_components, "\n\n")

  cat("Variance explained by consensus components:\n")
  var_df <- data.frame(
    Component = paste0("PC", 1:length(x$variance_explained)),
    Variance = round(x$variance_explained, 2)
  )
  print(var_df, row.names = FALSE)

  cat("\nCumulative variance:", round(sum(x$variance_explained), 2), "%\n\n")

  cat("Omics contribution to first component:\n")
  contrib_pc1 <- round(x$omics_contribution[, 1], 3)
  print(contrib_pc1)

  invisible(x)
}
