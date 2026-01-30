##################################################
## R script for OmicsAnalyst Pro
## Description: Arrow utilities for zero-copy data exchange with Java
## Author: OmicsAnalyst Team
## Part of the Rserve/qs to Apache Arrow migration
###################################################

#' Sync file to disk and verify existence (Safe-Handshake pattern)
#'
#' This function ensures that a file is fully written to disk before
#' returning control to Java. Uses normalizePath with mustWork=TRUE
#' to provide a filesystem-level guarantee that the file exists.
#'
#' @param file_path Path to the file to sync and verify
#' @param delay Delay in seconds before verification (default: 0.02 = 20ms)
#' @return The normalized (absolute) path to the file, or NULL if verification fails
#' @export
sync_file <- function(file_path, delay = 0.02) {
    if (is.null(file_path) || !nzchar(file_path)) {
        return(invisible(NULL))
    }

    # Brief delay to allow filesystem buffers to flush
    Sys.sleep(delay)

    # Use normalizePath with mustWork=TRUE as filesystem verification
    # This blocks until the OS confirms the file is accessible
    tryCatch({
        verified_path <- base::normalizePath(file_path, mustWork = TRUE)
        return(verified_path)
    }, error = function(e) {
        warning(sprintf("sync_file: File verification failed for '%s': %s",
                        file_path, e$message))
        return(NULL)
    })
}

#' Write Arrow file with safe-handshake verification
#'
#' Writes an Arrow (Feather) file and returns the verified absolute path.
#' Uses normalizePath(mustWork=TRUE) to ensure the file is fully written
#' before returning control to Java.
#'
#' @param df Data frame to write
#' @param path Path for the Arrow file
#' @param compress Compression type (default: "uncompressed" for memory-mapping)
#' @return The verified absolute path to the Arrow file
#' @export
write_arrow_safe <- function(df, path, compress = "uncompressed") {
    # Ensure factors are converted to character
    for (col in names(df)) {
        if (is.factor(df[[col]])) {
            df[[col]] <- as.character(df[[col]])
        }
    }

    # CRITICAL: Remove existing file first to prevent file-lock conflicts
    if (file.exists(path)) {
        unlink(path)
        Sys.sleep(0.01)
    }

    # Write the Arrow file
    arrow::write_feather(df, path, compression = compress)

    # Brief delay then verify with normalizePath
    Sys.sleep(0.02)

    # Return verified absolute path - blocks until file is confirmed accessible
    verified_path <- base::normalizePath(path, mustWork = TRUE)
    return(verified_path)
}

#' Safe column extraction with name-first, index-fallback strategy
#'
#' This function provides a safe way to extract columns from data frames
#' during Arrow migration. It ensures data integrity by:
#' 1. Trying named column access first (preferred)
#' 2. Falling back to index if name doesn't exist (with warning)
#' 3. Returning NA vector with error if both fail
#'
#' @param tab Data frame or matrix to extract from
#' @param name Column name to try first (character)
#' @param idx Fallback column index (1-based integer)
#' @param context Optional context string for logging
#' @return Column values as vector, or NA vector if extraction fails
#' @export
safeGetCol <- function(tab, name, idx = NULL, context = "") {
  nrows <- nrow(tab)
  if (is.null(nrows) || nrows == 0) {
    return(character(0))
  }

  # Strategy 1: Try by name first (preferred - most reliable)
  if (!is.null(name) && name %in% colnames(tab)) {
    return(tab[[name]])
  }

  # Strategy 2: Fall back to index (with warning)
  if (!is.null(idx) && is.numeric(idx) && idx > 0 && ncol(tab) >= idx) {
    actualName <- colnames(tab)[idx]
    warning(sprintf("[%s] Column '%s' not found, using index %d (actual column: '%s').",
                    context, name, idx, actualName))
    return(tab[, idx])
  }

  # Strategy 3: Both failed - return NA with error
  warning(sprintf("[%s] FAILED to extract column: name='%s', idx=%s. Available columns: %s",
                  context, name, ifelse(is.null(idx), "NULL", idx),
                  paste(colnames(tab), collapse=", ")))
  return(rep(NA, nrows))
}

#' Validate that required columns exist in a data frame
#'
#' @param tab Data frame to validate
#' @param required Character vector of required column names
#' @param context Context string for error messages
#' @return TRUE if all columns exist, FALSE otherwise (with warnings)
#' @export
validateColumns <- function(tab, required, context = "") {
  missing <- setdiff(required, colnames(tab))
  if (length(missing) > 0) {
    warning(sprintf("[%s] Missing required columns: %s. Available: %s",
                    context, paste(missing, collapse=", "),
                    paste(colnames(tab), collapse=", ")))
    return(FALSE)
  }
  return(TRUE)
}

#' Shadow save function for numeric matrices (Safe-Handshake pattern)
#'
#' Saves data as both .qs (for R) and .arrow (for Java zero-copy access).
#' ALWAYS preserves rownames as the first column (row_names_id) in Arrow files.
#' Returns the verified absolute path to the Arrow file for Java to memory-map.
#' Removes existing Arrow file first to prevent file-locking conflicts.
#'
#' @param obj The object to save (matrix or data.frame)
#' @param file The .qs file path (will auto-generate .arrow path)
#' @param compress Compression for Arrow (default: "uncompressed" for memory-mapping)
#' @return The verified absolute path to the Arrow file, or NULL if save failed
#' @export
shadow_save <- function(obj, file, compress = "uncompressed") {
    # Always save to qs for R compatibility
    qs::qsave(obj, file)

    # Generate Arrow path
    arrow_path <- sub("\\.qs$", ".arrow", file)

    tryCatch({
        if (is.matrix(obj) || is.data.frame(obj)) {
            df <- as.data.frame(obj)

            # Preserve rownames as first column
            rn <- rownames(obj)
            if (!is.null(rn) && length(rn) > 0 && !all(rn == as.character(1:nrow(df)))) {
                df <- cbind(row_names_id = as.character(rn), df)
            }

            # Ensure all character columns are properly typed
            for (col in names(df)) {
                if (is.factor(df[[col]])) {
                    df[[col]] <- as.character(df[[col]])
                }
            }

            # CRITICAL: Remove existing file first to prevent file-lock conflicts
            # Java should have closed its memory-map before R reaches here
            if (file.exists(arrow_path)) {
                unlink(arrow_path)
                Sys.sleep(0.01)
            }

            arrow::write_feather(df, arrow_path, compression = compress)

            # SAFE-HANDSHAKE: Brief delay then verify with normalizePath
            # This blocks until the OS confirms the file is accessible
            Sys.sleep(0.02)
            verified_path <- base::normalizePath(arrow_path, mustWork = TRUE)
            return(verified_path)
        }
    }, error = function(e) {
        warning(paste("Arrow shadow save failed:", e$message))
    })

    return(NULL)
}

#' Shadow save for mixed-type data frames (Safe-Handshake pattern)
#'
#' Handles factors by converting to character. Returns verified path.
#' Removes existing Arrow file first to prevent file-locking conflicts.
#'
#' @param obj Data frame or matrix to save
#' @param file Path to .qs file
#' @param compress Compression for Arrow (default: "uncompressed")
#' @return The verified absolute path to the Arrow file, or NULL if save failed
#' @export
shadow_save_mixed <- function(obj, file, compress = "uncompressed") {
    # Always save qs format for backward compatibility
    qs::qsave(obj, file)

    # Derive Arrow path
    arrow_path <- sub("\\.qs$", ".arrow", file)

    tryCatch({
        if (is.matrix(obj) || is.data.frame(obj)) {
            df <- as.data.frame(obj, stringsAsFactors = FALSE)

            # Convert factors to character for Arrow compatibility
            for (col in names(df)) {
                if (is.factor(df[[col]])) {
                    df[[col]] <- as.character(df[[col]])
                }
            }

            rn <- rownames(obj)
            if (!is.null(rn) && length(rn) > 0 && !all(rn == as.character(1:nrow(df)))) {
                # Prepend row_names_id as first column
                df <- cbind(row_names_id = as.character(rn), df)
            }

            # CRITICAL: Remove existing file first to prevent file-lock conflicts
            if (file.exists(arrow_path)) {
                unlink(arrow_path)
                Sys.sleep(0.01)
            }

            arrow::write_feather(df, arrow_path, compression = compress)

            # SAFE-HANDSHAKE: Brief delay then verify with normalizePath
            Sys.sleep(0.02)
            verified_path <- base::normalizePath(arrow_path, mustWork = TRUE)
            return(verified_path)
        }
    }, error = function(e) {
        warning(paste("Arrow save (mixed) failed:", e$message))
    })

    return(NULL)
}

#' Export statistical results matrix to Arrow format (Safe-Handshake)
#'
#' For correlation, clustering, or other analysis results. Returns verified path.
#' Removes existing file first to prevent file-locking conflicts with Java.
#'
#' @param result_mat Result matrix to export
#' @param filename Output filename (without .arrow extension)
#' @return The verified absolute path to the Arrow file, or NULL on failure
#' @export
ExportResultMatArrow <- function(result_mat, filename) {
    tryCatch({
        df <- as.data.frame(result_mat)

        # Convert factors to character
        for (col in names(df)) {
            if (is.factor(df[[col]])) {
                df[[col]] <- as.character(df[[col]])
            }
        }

        rn <- rownames(result_mat)
        if (!is.null(rn)) {
            df <- cbind(row_names_id = as.character(rn), df)
        }

        arrow_path <- paste0(filename, ".arrow")

        # CRITICAL: Remove existing file first to prevent file-lock conflicts
        # Java should have closed its memory-map before R reaches here
        if (file.exists(arrow_path)) {
            unlink(arrow_path)
            Sys.sleep(0.01)  # Brief pause after delete
        }

        arrow::write_feather(df, arrow_path, compression = "uncompressed")

        # SAFE-HANDSHAKE: Verify file is ready
        Sys.sleep(0.02)
        verified_path <- base::normalizePath(arrow_path, mustWork = TRUE)
        return(verified_path)
    }, error = function(e) {
        warning(paste("ExportResultMatArrow failed:", e$message))
        return(NULL)
    })
}

# =============================================================================
# FEATURE TABLE ARROW EXPORTS (for JSF DataTable views)
# =============================================================================
# These functions export feature detail tables as Arrow files for zero-copy
# access from Java. The tables include row IDs, symbols/labels, and data columns.
# =============================================================================

#' Export Covariate Significance Table to Arrow
#'
#' Exports dataSet$analSet$cov$sig.mat as Arrow file for Java DataTable view.
#' Includes row IDs and labels as separate columns for proper display.
#'
#' @param dataName Dataset name
#' @return The verified absolute path to the Arrow file, or NULL on failure
#' @export
ExportCovSigArrow <- function(dataName) {
    tryCatch({
        dataSet <- readDataset(dataName)
        sig.mat <- dataSet$analSet$cov$sig.mat

        if (is.null(sig.mat) || nrow(sig.mat) == 0) {
            return(NULL)
        }

        # Build export data frame with IDs, labels, and numeric data
        df <- data.frame(
            row_names_id = rownames(sig.mat),
            ids = as.character(sig.mat$ids),
            label = as.character(sig.mat$label),
            stringsAsFactors = FALSE
        )

        # Add numeric columns (exclude ids and label)
        drops <- c("ids", "label")
        numeric_cols <- sig.mat[, !(names(sig.mat) %in% drops), drop = FALSE]
        df <- cbind(df, numeric_cols)

        # Convert any factors to character
        for (col in names(df)) {
            if (is.factor(df[[col]])) {
                df[[col]] <- as.character(df[[col]])
            }
        }

        arrow_path <- paste0("covsig_", dataName, ".arrow")

        # Safe-handshake: remove existing file first
        if (file.exists(arrow_path)) {
            unlink(arrow_path)
            Sys.sleep(0.01)
        }

        arrow::write_feather(df, arrow_path, compression = "uncompressed")
        Sys.sleep(0.02)

        return(base::normalizePath(arrow_path, mustWork = TRUE))
    }, error = function(e) {
        warning(paste("ExportCovSigArrow failed:", e$message))
        return(NULL)
    })
}

#' Export Loading Table to Arrow
#'
#' Exports loading.pos.xyz for a specific omics type as Arrow file.
#' Includes row IDs, labels, and loading values.
#'
#' @param dataName Dataset name (determines omics type filter)
#' @return The verified absolute path to the Arrow file, or NULL on failure
#' @export
ExportLoadingArrow <- function(dataName) {
    tryCatch({
        reductionSet <- .get.rdt.set()
        dataSet <- readDataset(dataName)

        method <- reductionSet$reductionOpt
        loading.pos.xyz <- reductionSet[[method]]$loading.pos.xyz

        if (is.null(loading.pos.xyz) || nrow(loading.pos.xyz) == 0) {
            return(NULL)
        }

        # Filter by omics type
        omicstype <- dataSet$type
        inx <- loading.pos.xyz$type %in% omicstype
        filtered <- loading.pos.xyz[inx, , drop = FALSE]

        if (nrow(filtered) == 0) {
            return(NULL)
        }

        # Build export data frame
        df <- data.frame(
            row_names_id = rownames(filtered),
            ids = as.character(filtered$ids),
            label = as.character(filtered$label),
            type = as.character(filtered$type),
            stringsAsFactors = FALSE
        )

        # Add numeric columns (exclude ids, label, type)
        drops <- c("ids", "label", "type")
        numeric_cols <- filtered[, !(names(filtered) %in% drops), drop = FALSE]
        df <- cbind(df, numeric_cols)

        # Convert any factors to character
        for (col in names(df)) {
            if (is.factor(df[[col]])) {
                df[[col]] <- as.character(df[[col]])
            }
        }

        arrow_path <- paste0("loading_", method, "_", dataName, ".arrow")

        # Safe-handshake: remove existing file first
        if (file.exists(arrow_path)) {
            unlink(arrow_path)
            Sys.sleep(0.01)
        }

        arrow::write_feather(df, arrow_path, compression = "uncompressed")
        Sys.sleep(0.02)

        return(base::normalizePath(arrow_path, mustWork = TRUE))
    }, error = function(e) {
        warning(paste("ExportLoadingArrow failed:", e$message))
        return(NULL)
    })
}

#' Export Node Table to Arrow
#'
#' Exports node_table.csv data as Arrow file for Java DataTable view.
#' Reads from CSV and exports as Arrow with proper column types.
#'
#' @return The verified absolute path to the Arrow file, or NULL on failure
#' @export
ExportNodeTableArrow <- function() {
    tryCatch({
        df <- .readDataTable('node_table.csv')

        if (is.null(df) || nrow(df) == 0) {
            return(NULL)
        }

        # Ensure Id and Label columns exist and are character
        if ("Id" %in% names(df)) {
            df$Id <- as.character(df$Id)
        }
        if ("Label" %in% names(df)) {
            df$Label <- as.character(df$Label)
        }

        # Add row_names_id as first column (use Id if available)
        if ("Id" %in% names(df)) {
            df <- cbind(row_names_id = df$Id, df)
        } else {
            df <- cbind(row_names_id = as.character(1:nrow(df)), df)
        }

        # Convert numeric columns
        numeric_cols <- names(df)[!names(df) %in% c("row_names_id", "Id", "Label")]
        for (col in numeric_cols) {
            df[[col]] <- as.numeric(as.character(df[[col]]))
        }

        # Convert any factors to character
        for (col in names(df)) {
            if (is.factor(df[[col]])) {
                df[[col]] <- as.character(df[[col]])
            }
        }

        arrow_path <- "node_table.arrow"

        # Safe-handshake: remove existing file first
        if (file.exists(arrow_path)) {
            unlink(arrow_path)
            Sys.sleep(0.01)
        }

        arrow::write_feather(df, arrow_path, compression = "uncompressed")
        Sys.sleep(0.02)

        return(base::normalizePath(arrow_path, mustWork = TRUE))
    }, error = function(e) {
        warning(paste("ExportNodeTableArrow failed:", e$message))
        return(NULL)
    })
}
