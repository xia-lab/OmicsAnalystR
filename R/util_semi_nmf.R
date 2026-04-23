##################################################
## Semi Non-negative Matrix Factorization (Semi-NMF)
## Handles data with negative values directly
## without pre-transformation
## Reference: Ding et al. (2010) Convex and Semi-Nonnegative Matrix Factorizations
## Created: 2025-01-04
##################################################

#' Perform Semi-NMF for Multi-Omics Integration
#'
#' Semi-NMF allows negative input values while keeping factor matrices non-negative.
#' This preserves the biological meaning of negative values (e.g., down-regulation)
#' while maintaining interpretable parts-based representations.
#'
#' @param omics_list Named list of omics matrices (samples x features)
#'                   All matrices must have same number of rows (samples)
#' @param k Number of factors/components to extract (default: 10)
#' @param method Integration strategy: "concatenate" (combine features) or
#'               "joint" (joint factorization with omics-specific weights)
#' @param normalize Whether to normalize features (default: TRUE)
#' @param max_iter Maximum iterations (default: 100)
#' @param tol Convergence tolerance (default: 1e-4)
#' @param seed Random seed for reproducibility (default: 123)
#' @param progress_callback Optional progress callback function(percent, message)
#'
#' @return List containing:
#'   \item{H}{Sample factor matrix (samples x k) - non-negative}
#'   \item{W}{Feature weight matrices (list per omics) - non-negative}
#'   \item{omics_contributions}{Contribution of each omics to factors}
#'   \item{feature_importance}{Top contributing features per factor}
#'   \item{reconstruction_error}{Final reconstruction error}
#'   \item{variance_explained}{Variance explained per factor}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Integrate multiple omics with negative values
#' omics_list <- list(
#'   transcriptomics = matrix(rnorm(100*1000), 100, 1000),  # Can be negative
#'   proteomics = matrix(rnorm(100*500), 100, 500),
#'   metabolomics = matrix(rnorm(100*200), 100, 200)
#' )
#'
#' # Run Semi-NMF
#' result <- PerformSemiNMF(omics_list, k = 10)
#'
#' # Plot factor scores
#' PlotSemiNMFFactors(result, factors = c(1,2), groups = your_groups)
#' }
PerformSemiNMF <- function(omics_list,
                           k = 10,
                           method = c("concatenate", "joint"),
                           normalize = TRUE,
                           max_iter = 100,
                           tol = 1e-4,
                           seed = 123,
                           progress_callback = NULL) {

  method <- match.arg(method)

  # Progress helper
  update_progress <- if (!is.null(progress_callback)) {
    progress_callback
  } else {
    function(pct, msg) NULL
  }

  update_progress(5, "Validating input data")

  # Input validation
  if (!is.list(omics_list) || length(omics_list) < 1) {
    stop("omics_list must be a list with at least 1 omics layer")
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

  update_progress(10, "Preprocessing omics data")

  # Preprocess each omics
  processed_omics <- list()
  omics_scales <- list()

  for (i in seq_along(omics_list)) {
    omics_name <- omics_names[i]
    omics_data <- as.matrix(omics_list[[i]])

    # Remove zero-variance features
    feature_vars <- apply(omics_data, 2, var, na.rm = TRUE)
    valid_features <- !is.na(feature_vars) & feature_vars > 0
    omics_data <- omics_data[, valid_features, drop = FALSE]

    # Handle missing values (replace with column mean)
    if (any(is.na(omics_data))) {
      col_means <- colMeans(omics_data, na.rm = TRUE)
      for (j in 1:ncol(omics_data)) {
        omics_data[is.na(omics_data[, j]), j] <- col_means[j]
      }
      warning(sprintf("%s had missing values - replaced with column means", omics_name))
    }

    # Normalize features (center and scale)
    if (normalize) {
      omics_data <- scale(omics_data, center = TRUE, scale = TRUE)
      omics_scales[[omics_name]] <- attr(omics_data, "scaled:scale")
    }

    # Semi-NMF allows negative values - no transformation needed
    processed_omics[[omics_name]] <- omics_data

    update_progress(10 + (i / length(omics_list)) * 20,
                   sprintf("Preprocessed %s", omics_name))
  }

  update_progress(35, sprintf("Running Semi-NMF (%s method)", method))

  # Set seed for reproducibility
  set.seed(seed)

  if (method == "concatenate") {
    # Concatenate all omics features
    result <- PerformConcatenatedSemiNMF(
      processed_omics = processed_omics,
      k = k,
      max_iter = max_iter,
      tol = tol,
      update_progress = update_progress
    )
  } else {
    # Joint factorization (shared H, omics-specific W)
    result <- PerformJointSemiNMF(
      processed_omics = processed_omics,
      k = k,
      max_iter = max_iter,
      tol = tol,
      update_progress = update_progress
    )
  }

  update_progress(90, "Computing feature importance")

  # Compute feature importance for each factor
  feature_importance <- ComputeSemiNMFFeatureImportance(
    W_list = result$W,
    top_n = 20
  )

  update_progress(95, "Computing omics contributions")

  # Compute omics contribution to each factor
  omics_contributions <- ComputeSemiNMFOmicsContributions(
    W_list = result$W,
    processed_omics = processed_omics
  )

  update_progress(98, "Finalizing results")

  # Package results
  final_result <- list(
    H = result$H,  # Sample factor matrix (non-negative)
    W = result$W,  # List of feature weight matrices (non-negative)
    omics_contributions = omics_contributions,
    feature_importance = feature_importance,
    reconstruction_error = result$error,
    variance_explained = result$variance_explained,
    method = paste0("Semi-NMF (", method, ")"),
    k = k,
    omics_names = omics_names,
    preprocessing = list(
      normalized = normalize,
      scales = omics_scales,
      allows_negative = TRUE
    )
  )

  class(final_result) <- c("SemiNMF", "list")

  update_progress(100, "Semi-NMF complete")

  return(final_result)
}


#' Perform Concatenated Semi-NMF
#'
#' @param processed_omics List of preprocessed omics matrices
#' @param k Number of factors
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param update_progress Progress callback function
#' @return List with H, W, error, variance_explained
#' @keywords internal
PerformConcatenatedSemiNMF <- function(processed_omics, k, max_iter, tol, update_progress) {

  # Concatenate all omics
  concat_matrix <- do.call(cbind, processed_omics)
  X <- concat_matrix  # samples x features

  # Keep track of feature boundaries
  feature_boundaries <- cumsum(c(0, sapply(processed_omics, ncol)))
  names(feature_boundaries) <- c(names(processed_omics), "end")

  update_progress(40, "Initializing Semi-NMF matrices")

  n_samples <- nrow(X)
  n_features <- ncol(X)

  # Initialize H (samples x k) - non-negative
  H <- matrix(runif(n_samples * k, 0, 1), nrow = n_samples, ncol = k)

  # Initialize W (features x k) - non-negative
  W <- matrix(runif(n_features * k, 0, 1), nrow = n_features, ncol = k)

  # Semi-NMF update rules (from Ding et al. 2010)
  prev_error <- Inf

  for (iter in 1:max_iter) {
    # Update H with non-negativity constraint
    # H^T = (W^T * W)^-1 * W^T * X^T
    WtW <- t(W) %*% W + 1e-10 * diag(k)  # Regularization
    H_new <- X %*% W %*% solve(WtW)

    # Project to non-negative
    H_new[H_new < 0] <- 0
    H <- H_new

    # Update W with non-negativity constraint
    # Split W update into positive and negative parts of X*H
    XH <- t(X) %*% H  # features x k
    HtH <- t(H) %*% H + 1e-10 * diag(k)

    # Separate positive and negative contributions
    XH_pos <- pmax(XH, 0)
    XH_neg <- pmax(-XH, 0)

    # Update W
    # W = W * sqrt((XH_pos) / (W * HtH + XH_neg + eps))
    numerator <- XH_pos
    denominator <- W %*% HtH + XH_neg + 1e-10

    W <- W * sqrt(numerator / denominator)

    # Ensure W stays non-negative
    W[W < 0] <- 0
    W[is.na(W)] <- 0
    W[is.infinite(W)] <- 0

    # Compute reconstruction error
    reconstruction <- H %*% t(W)
    error <- sqrt(mean((X - reconstruction)^2))

    # Check convergence
    if (abs(prev_error - error) < tol) {
      message(sprintf("Semi-NMF converged at iteration %d", iter))
      break
    }

    prev_error <- error

    # Progress update every 10 iterations
    if (iter %% 10 == 0) {
      progress_pct <- 40 + (iter / max_iter) * 30
      update_progress(progress_pct, sprintf("Semi-NMF iteration %d/%d", iter, max_iter))
    }
  }

  update_progress(75, "Extracting factor matrices")

  # Split W back into omics-specific matrices
  W_list <- list()
  for (i in seq_along(processed_omics)) {
    omics_name <- names(processed_omics)[i]
    start_idx <- feature_boundaries[i] + 1
    end_idx <- feature_boundaries[i + 1]

    W_omics <- W[start_idx:end_idx, , drop = FALSE]
    rownames(W_omics) <- colnames(processed_omics[[i]])
    colnames(W_omics) <- paste0("Factor", 1:k)

    W_list[[omics_name]] <- W_omics
  }

  # Set names for H
  rownames(H) <- rownames(processed_omics[[1]])
  colnames(H) <- paste0("Factor", 1:k)

  # Compute variance explained by each factor
  H_vars <- apply(H, 2, var)
  variance_explained <- (H_vars / sum(H_vars)) * 100

  return(list(
    H = H,
    W = W_list,
    error = error,
    variance_explained = variance_explained
  ))
}


#' Perform Joint Semi-NMF (Shared H, Omics-specific W)
#'
#' @param processed_omics List of preprocessed omics matrices
#' @param k Number of factors
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param update_progress Progress callback function
#' @return List with H, W, error, variance_explained
#' @keywords internal
PerformJointSemiNMF <- function(processed_omics, k, max_iter, tol, update_progress) {

  n_samples <- nrow(processed_omics[[1]])
  n_omics <- length(processed_omics)

  # Initialize H randomly (shared across omics) - non-negative
  H <- matrix(runif(n_samples * k, 0, 1), nrow = n_samples, ncol = k)

  # Initialize W for each omics - non-negative
  W_list <- list()
  for (omics_name in names(processed_omics)) {
    n_features <- ncol(processed_omics[[omics_name]])
    W_list[[omics_name]] <- matrix(runif(n_features * k, 0, 1), nrow = n_features, ncol = k)
  }

  # Iterative optimization
  prev_error <- Inf

  for (iter in 1:max_iter) {
    # Update each W given H (with non-negativity)
    for (i in seq_along(processed_omics)) {
      omics_name <- names(processed_omics)[i]
      X <- processed_omics[[omics_name]]

      # W update with non-negativity
      XH <- t(X) %*% H  # features x k
      HtH <- t(H) %*% H + 1e-10 * diag(k)

      # Separate positive and negative parts
      XH_pos <- pmax(XH, 0)
      XH_neg <- pmax(-XH, 0)

      # Update W
      W <- W_list[[omics_name]]
      numerator <- XH_pos
      denominator <- W %*% HtH + XH_neg + 1e-10

      W_new <- W * sqrt(numerator / denominator)

      # Ensure non-negativity
      W_new[W_new < 0] <- 0
      W_new[is.na(W_new)] <- 0
      W_new[is.infinite(W_new)] <- 0

      W_list[[omics_name]] <- W_new
    }

    # Update H given all W (weighted average with non-negativity)
    H_new <- matrix(0, nrow = n_samples, ncol = k)

    for (i in seq_along(processed_omics)) {
      omics_name <- names(processed_omics)[i]
      X <- processed_omics[[omics_name]]
      W <- W_list[[omics_name]]

      # H update
      XtW <- t(X) %*% W  # This is actually X %*% W since we want samples x k
      # Correction: XW not XtW
      XW <- X %*% W
      WtW <- t(W) %*% W + 1e-10 * diag(k)

      H_update <- XW %*% solve(WtW)
      H_new <- H_new + H_update
    }

    H <- H_new / n_omics

    # Ensure non-negativity
    H[H < 0] <- 0
    H[is.na(H)] <- 0
    H[is.infinite(H)] <- 0

    # Compute reconstruction error
    total_error <- 0
    for (i in seq_along(processed_omics)) {
      omics_name <- names(processed_omics)[i]
      X <- processed_omics[[omics_name]]
      W <- W_list[[omics_name]]

      reconstruction <- H %*% t(W)
      total_error <- total_error + sum((X - reconstruction)^2)
    }

    error <- sqrt(total_error / sum(sapply(processed_omics, length)))

    # Check convergence
    if (abs(prev_error - error) < tol) {
      message(sprintf("Joint Semi-NMF converged at iteration %d", iter))
      break
    }

    prev_error <- error

    # Progress update every 10 iterations
    if (iter %% 10 == 0) {
      progress_pct <- 40 + (iter / max_iter) * 30
      update_progress(progress_pct, sprintf("Joint Semi-NMF iteration %d/%d", iter, max_iter))
    }
  }

  update_progress(75, "Joint Semi-NMF optimization complete")

  # Set names
  rownames(H) <- rownames(processed_omics[[1]])
  colnames(H) <- paste0("Factor", 1:k)

  for (omics_name in names(W_list)) {
    rownames(W_list[[omics_name]]) <- colnames(processed_omics[[omics_name]])
    colnames(W_list[[omics_name]]) <- paste0("Factor", 1:k)
  }

  # Compute variance explained
  H_vars <- apply(H, 2, var)
  variance_explained <- (H_vars / sum(H_vars)) * 100

  return(list(
    H = H,
    W = W_list,
    error = error,
    variance_explained = variance_explained
  ))
}


#' Compute Feature Importance for Semi-NMF
#'
#' @param W_list List of feature weight matrices
#' @param top_n Number of top features per factor
#' @return List of top features per factor per omics
#' @keywords internal
ComputeSemiNMFFeatureImportance <- function(W_list, top_n = 20) {

  feature_importance <- list()

  for (omics_name in names(W_list)) {
    W <- W_list[[omics_name]]
    k <- ncol(W)

    omics_importance <- list()

    for (j in 1:k) {
      weights <- W[, j]
      top_idx <- head(order(weights, decreasing = TRUE), top_n)

      top_features <- data.frame(
        feature = rownames(W)[top_idx],
        weight = weights[top_idx],
        stringsAsFactors = FALSE
      )

      omics_importance[[paste0("Factor", j)]] <- top_features
    }

    feature_importance[[omics_name]] <- omics_importance
  }

  return(feature_importance)
}


#' Compute Omics Contributions to Factors for Semi-NMF
#'
#' @param W_list List of feature weight matrices
#' @param processed_omics List of preprocessed omics
#' @return Matrix of omics contributions (omics x factors)
#' @keywords internal
ComputeSemiNMFOmicsContributions <- function(W_list, processed_omics) {

  k <- ncol(W_list[[1]])
  n_omics <- length(W_list)

  contributions <- matrix(0, nrow = n_omics, ncol = k)
  rownames(contributions) <- names(W_list)
  colnames(contributions) <- paste0("Factor", 1:k)

  # Contribution = sum of weights per factor
  for (i in seq_along(W_list)) {
    W <- W_list[[i]]
    contributions[i, ] <- colSums(W)
  }

  # Normalize to percentages per factor
  contributions <- sweep(contributions, 2, colSums(contributions), "/") * 100

  return(contributions)
}


#' Plot Semi-NMF Factor Scores
#'
#' @param nmf_result Result from PerformSemiNMF
#' @param factors Which factors to plot (default: c(1,2))
#' @param groups Optional grouping variable for coloring
#' @param main Plot title
#' @export
PlotSemiNMFFactors <- function(nmf_result,
                               factors = c(1, 2),
                               groups = NULL,
                               main = "Semi-NMF Factors") {

  if (!inherits(nmf_result, "SemiNMF")) {
    stop("Input must be a SemiNMF object")
  }

  H <- nmf_result$H[, factors]
  var_exp <- nmf_result$variance_explained[factors]

  # Set colors
  if (is.null(groups)) {
    colors <- rep("black", nrow(H))
  } else {
    groups <- as.factor(groups)
    colors <- rainbow(nlevels(groups))[as.numeric(groups)]
  }

  # Create plot
  plot(H,
       col = colors,
       pch = 19,
       cex = 1.5,
       xlab = sprintf("Factor%d (%.1f%%)", factors[1], var_exp[1]),
       ylab = sprintf("Factor%d (%.1f%%)", factors[2], var_exp[2]),
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


#' Print Semi-NMF Summary
#'
#' @param x SemiNMF object
#' @param ... Additional arguments (not used)
#' @export
print.SemiNMF <- function(x, ...) {
  cat("=== Semi Non-negative Matrix Factorization ===\n\n")
  cat("Method:", x$method, "\n")
  cat("Number of factors:", x$k, "\n")
  cat("Number of omics layers:", length(x$W), "\n")
  cat("Omics names:", paste(x$omics_names, collapse = ", "), "\n")
  cat("Number of samples:", nrow(x$H), "\n")
  cat("Allows negative input values: YES\n\n")

  cat("Reconstruction error:", sprintf("%.4f", x$reconstruction_error), "\n\n")

  cat("Variance explained by factors:\n")
  var_df <- data.frame(
    Factor = paste0("Factor", 1:length(x$variance_explained)),
    Variance = round(x$variance_explained, 2)
  )
  #print(var_df, row.names = FALSE)

  cat("\nCumulative variance:", round(sum(x$variance_explained), 2), "%\n\n")

  cat("Omics contributions to Factor1:\n")
  contrib_factor1 <- round(x$omics_contributions[, 1], 1)
  print(contrib_factor1)

  cat("\nTop features in Factor1:\n")
  for (omics_name in names(x$feature_importance)) {
    cat(sprintf("  %s: %s\n",
                omics_name,
                paste(head(x$feature_importance[[omics_name]]$Factor1$feature, 3),
                      collapse = ", ")))
  }

  invisible(x)
}
