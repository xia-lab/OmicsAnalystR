
PerformVariancePartitioning <- function(omics.type, top_n = 100, meta = NA, fileName = "variance_partition_plot", dpi = 300, format = "png") {
  library(variancePartition)
  library(limma)
  library(Cairo)
  library(car)

  rdtSet <- .get.rdt.set()
  
  # Ensure that normalized data and metadata are available
  if (is.null(rdtSet$dataSet$norm) || is.null(rdtSet$dataSet$meta.info) || is.null(rdtSet$dataSet$meta.types)) {
    stop("Normalized data, metadata, or meta.types not found.")
  }
  
  # Extract normalized gene expression data
  gene_expr <- rdtSet$dataSet$norm
  meta_data <- rdtSet$dataSet$meta.info
  meta_types <- rdtSet$dataSet$meta.types  # Continuous ("cont") or discrete ("disc") types
  
  # Align the sample names between gene expression data and metadata
  common_samples <- intersect(colnames(gene_expr), rownames(meta_data))
  
  # Subset both gene expression data and metadata to have only common samples
  gene_expr <- gene_expr[, common_samples, drop = FALSE]
  meta_data <- meta_data[common_samples, , drop = FALSE]
  
  # Check if they are correctly aligned
  if (!all(colnames(gene_expr) == rownames(meta_data))) {
    stop("Sample names between gene expression data and metadata do not match.")
  }
  meta_data <- meta_data[,-3]
  meta_types <- meta_types[-3]
  meta_data[meta_data == "NA"] <- NA
  # Check for missing values in meta_data and remove rows with missing data
  if (any(is.na(meta_data))) {
    warning("Missing values found in metadata. Removing rows with missing data.")
    meta_data <- na.omit(meta_data)
    # Align gene expression data again after removing NA rows
    gene_expr <- gene_expr[, rownames(meta_data), drop = FALSE]
  }
  
  # Dynamically create the formula based on meta.types
  fixed_effects <- names(meta_types[meta_types == "cont"])
  random_effects <- names(meta_types[meta_types == "disc"])
  
  # Create formula for variance partitioning
  formula_fixed <- paste(fixed_effects, collapse = " + ")
  formula_random <- paste(paste0("(1 | ", random_effects, ")"), collapse = " + ")
  
  if (formula_fixed != "" && formula_random != "") {
    formula <- as.formula(paste("~", formula_fixed, "+", formula_random))
  } else if (formula_fixed != "") {
    formula <- as.formula(paste("~", formula_fixed))
  } else if (formula_random != "") {
    formula <- as.formula(paste("~", formula_random))
  } else {
    stop("No valid covariates found in metadata.")
  }
  
  # Ensure factors are treated as factors and numeric covariates are numeric
  for (col in colnames(meta_data)) {
    if (meta_types[col] == "disc") {
      meta_data[[col]] <- as.factor(meta_data[[col]])
    } else {
      meta_data[[col]] <- as.numeric(meta_data[[col]])
    }
  }
  
  # Perform variance partitioning
  varPart <- fitExtractVarPartModel(gene_expr, formula, meta_data)
  
  # Sort genes by the amount of variance explained (sum of variance explained by all factors)
  varExplained <- rowSums(varPart)
  
  # Get the top 'top_n' genes by variance explained
  top_genes_idx <- order(varExplained, decreasing = TRUE)[1:top_n]
  top_genes <- rownames(gene_expr)[top_genes_idx]
  
  # Return the top genes and their variance explained
  top_gene_results <- data.frame(
    Gene = top_genes,
    VarianceExplained = varExplained[top_genes_idx]
  )
  
  # Plot the variance partitioning results and save the image
  file_path <- paste0(fileName, ".", format)
  
  Cairo::Cairo(file = file_path, type = format, dpi = dpi, width = 8, height = 6, units = "in", bg = "white")
  plotVarPart(varPart)
  dev.off()
  
  # Store the top gene results in rdtSet for future use
  rdtSet$analSet$top_genes_varpart <- top_gene_results
  
  return(.set.rdt.set(rdtSet))
}
