PlotPercentBars <- function(top_n=10, fileName="", dpi=72, format="png"){
  rdtSet <- .get.rdt.set()
  varPart <- rdtSet$analSet$varPart.df[,-1];
  vp <- sortCols(varPart);
  imgName = paste(fileName, "dpi", dpi, ".", format, sep="");
  Cairo::Cairo(file = imgName, type = format, dpi = dpi, width = 10, height = 6, units = "in", bg = "white")
  p<-plotPercentBars(vp[1:top_n, ]) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),  
      axis.title.x = element_text(size = 14),                
      axis.title.y = element_text(size = 14),                
      axis.text = element_text(size = 12),                   
      legend.text = element_text(size = 12),                 
      legend.title = element_text(size = 14)                
    )
  print(p);
  dev.off();
  return(1)
}

PerformVarPartOverview <- function(dataName,top_n = 1000, fileName = "variance_partition_plot", dpi = 300, format = "png") {
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
  gene_expr <- rdtSet$dataSet$norm;
  meta_data <- rdtSet$dataSet$meta.info
  meta_types <- rdtSet$dataSet$meta.types  # Continuous ("cont") or discrete ("disc") types

  if(dataName!="all"){
    dataSet <- readDataset(dataName);
    var.nms <- paste0(dataSet[["type"]],".",dataSet[["orig.var.nms"]])
    gene_expr <-  gene_expr[rownames( gene_expr) %in% var.nms,]
   }

  if(nrow(gene_expr)>500){
    gene_expr <-  gene_expr[1:500,]

  }
   
  
  if (exists("fixed_effects.vec") && length(fixed_effects.vec) == 0 && 
      exists("random_effects.vec") && length(random_effects.vec) == 0) {
    fixed_effects <- colnames(meta_data)
    random_effects <- "";
  }else{
    
    # Check if the fixed_effects and random_effects are empty or null
    if (exists("fixed_effects.vec") && length(fixed_effects.vec) == 0) {
      fixed_effects <- ""
    } else {
      fixed_effects <- fixed_effects.vec
    }
    
    if (exists("random_effects.vec") && length(random_effects.vec) == 0) {
      random_effects <- ""
    } else {
      random_effects <- random_effects.vec
    }
  }
  
  # Align the sample names between gene expression data and metadata
  common_samples <- intersect(colnames(gene_expr), rownames(meta_data))
  
  # Subset both gene expression data and metadata to have only common samples
  gene_expr <- gene_expr[, common_samples, drop = FALSE]
  meta_data <- meta_data[common_samples, , drop = FALSE]
  
  # Check if they are correctly aligned
  if (!all(colnames(gene_expr) == rownames(meta_data))) {
    stop("Sample names between gene expression data and metadata do not match.")
  }
  

  meta_data[meta_data == "NA"] <- NA
  # Check for missing values in meta_data and remove rows with missing data
  if (any(is.na(meta_data))) {
    warning("Missing values found in metadata. Removing rows with missing data.")
    meta_data <- na.omit(meta_data)
    # Align gene expression data again after removing NA rows
    gene_expr <- gene_expr[, rownames(meta_data), drop = FALSE]
  }
  
  if (length(fixed_effects) == 0) {
    stop("At least one fixed effect must be specified.")
  }
  
  # Create formula for variance partitioning
  formula_fixed <- paste(fixed_effects, collapse = " + ")
  if(random_effects != ""){
    formula_random <- paste(paste0("(1 | ", random_effects, ")"), collapse = " + ")
  }else{
    formula_random <- "";
  }
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
  varExplained_first_col <- varPart[, 1]
  top_genes_idx <- order(varExplained_first_col, decreasing = TRUE)
  
  # Reorder the entire varPart matrix based on the top genes
  varPart <- varPart[top_genes_idx, , drop = FALSE]
  
  # Plot the variance partitioning results and save the image
  imgName = paste(fileName, "dpi", dpi, ".", format, sep="");
  
  Cairo::Cairo(file = imgName, type = format, dpi = dpi, width = 8, height = 6, units = "in", bg = "white")
  p <- plotVarPart(varPart)
  print(p)
  dev.off()
  
  # Store the top gene results in rdtSet for future use
  rdtSet$analSet$varPart.symbols <- rownames(varPart);
  varPart<- cbind(Symbol = rdtSet$analSet$varPart.symbols, varPart);
  rdtSet$analSet$varPart.fileName <- "varPart_results.csv";
  rdtSet$analSet$varPart.df <- varPart;
 
  fast.write.csv(varPart, file="varPart_results.csv")
.set.rdt.set(rdtSet)
  return(1)
}

PrepareVarData <- function(type="NA"){
  #save.image("prepare.RData");
  rdtSet <- .get.rdt.set();
  
  data.list <- list();
  omics.vec <- vector();
  featureNms <- vector();
  uniqFeats <- vector();
  
  rdtSet$dataSet$meta.info.proc <- process_metadata(rdtSet$dataSet$meta.info);
  meta.sel.inx <- mmeta.all == 1;
  meta.sel.nms <- c();  # Assuming no metadata selection for this case
  
  if(length(meta.sel.nms) > 0) {
    for(i in 1:length(meta.sel.nms)){
      data.list[[meta.sel.nms[i]]] <- rdtSet$dataSet$meta.info.proc[,meta.sel.nms[i]]
      if(i == 1){       
        featureNms <- meta.sel.nms[i]
        omics.vec <- meta.sel.nms[i]
        uniqFeats <- meta.sel.nms[i]
      } else {
        featureNms <- c(featureNms, meta.sel.nms[i]);
        omics.vec <- c(omics.vec, meta.sel.nms[i]);
        uniqFeats <- c(uniqFeats, meta.sel.nms[i])
      }
    }
  }
  
  sel.inx <- mdata.all == 1;
  sel.nms <- names(mdata.all)[sel.inx];
  
  for(i in 1:length(sel.nms)){
    dataSet = readDataset(sel.nms[i])
    
    # Sanitize row names
    sanitized_names <- gsub("[[:cntrl:]]|[^[:ascii:]]", "_", rownames(dataSet$data.proc), perl = TRUE)
    rownames(dataSet$data.proc) <- sanitized_names;
    
    # Calculate variance for each feature
    feature_variances <- apply(dataSet$data.proc, 1, var)
    ordinx <- order(feature_variances, decreasing = TRUE)
    dataSet$data.proc <-  dataSet$data.proc[ordinx,]

    # Sort features by variance and select the top 1000 (if applicable)
    #if (nrow(dataSet$data.proc) > 1000) {
    #  top_features <- order(feature_variances, decreasing = TRUE)[1:1000]
    #  dataSet$data.proc <- dataSet$data.proc[top_features, , drop = FALSE]
    #}
    
    # Add to the list
    data.list[[dataSet$type]] <- dataSet$data.proc
    
    if(i == 1){       
      featureNms <- rownames(dataSet$data.proc);
      omics.vec <- rep(dataSet$type, nrow(dataSet$data.proc));
      uniqFeats <- paste0(rownames(dataSet$data.proc), "_", dataSet$type)
    } else {
      featureNms <- c(featureNms, rownames(dataSet$data.proc));
      omics.vec <- c(omics.vec, rep(dataSet$type, nrow(dataSet$data.proc)));
      uniqFeats <- c(uniqFeats, paste0(rownames(dataSet$data.proc), "_", dataSet$type))
    }
  }
  
  # Convert vectors to data frames if necessary
  for (i in seq_along(data.list)) {
    if (is.vector(data.list[[i]])) {
      data.list[[i]] <- data.frame(value = data.list[[i]])
    }
  }
  
  # Merge all datasets
  merged_data <- do.call(rbind, data.list)
  
  # Check if there are new samples to update `norm`
  new.inx <- is.na(rdtSet$dataSet$cls.all) | rdtSet$dataSet$cls.all == "";
  if(sum(new.inx) > 0){
    rdtSet$dataSet$new.samples <- TRUE;
    rdtSet$dataSet$new.data <- rdtSet$dataSet$norm.all[new.inx, , drop = F];
    rdtSet$dataSet$norm <- merged_data;
    rdtSet$dataSet$cls <- factor(rdtSet$dataSet$meta.info[,1])
  }else{
    rdtSet$dataSet$new.samples <- FALSE;
    rdtSet$dataSet$new.data <- NULL;
    rdtSet$dataSet$norm <- merged_data;
    rdtSet$dataSet$cls <- rdtSet$dataSet$meta.info[,1]; 
  }
  
  return(.set.rdt.set(rdtSet));
}

# Function to get the variance partition matrix (as a numeric matrix, excluding the symbol column)
GetVarMat <- function() {
  rdtSet <- .get.rdt.set()
 
  # Assuming varPart is stored in rdtSet$analSet$varPart.df
  if (is.null(rdtSet$analSet$varPart.df)) {
    stop("Variance partition matrix not found.")
  }
  
  # Extract the variance partitioning matrix without the symbol column
  varPart_matrix <- as.matrix(rdtSet$analSet$varPart.df[ , -1, drop = FALSE]) # Removing the symbol column
  
  return(varPart_matrix)
}

# Function to get the row names (gene IDs) of the variance partitioning matrix
GetVarIds <- function() {
  rdtSet <- .get.rdt.set()
  
  # Assuming varPart is stored in rdtSet$analSet$varPart.df
  if (is.null(rdtSet$analSet$varPart.df)) {
    stop("Variance partition matrix not found.")
  }
  
  # Return the row names (gene IDs)
  varPart_ids <- rownames(rdtSet$analSet$varPart.df)
  
  return(varPart_ids)
}

# Function to get the gene symbols associated with the variance partitioning matrix
GetVarSymbols <- function() {
  rdtSet <- .get.rdt.set()
  
  # Assuming varPart is stored in rdtSet$analSet$varPart.df
  if (is.null(rdtSet$analSet$varPart.df)) {
    stop("Variance partition matrix not found.")
  }
  
  # Extract the symbols column
  varPart_symbols <- rdtSet$analSet$varPart.df$Symbol
  return(varPart_symbols)
}

# Function to get the column names of the variance partitioning matrix (excluding the symbol column)
GetVarColNames <- function() {
  rdtSet <- .get.rdt.set()
  
  # Assuming varPart is stored in rdtSet$analSet$varPart.df
  if (is.null(rdtSet$analSet$varPart.df)) {
    stop("Variance partition matrix not found.")
  }
  
  # Extract the column names (excluding the symbol column)
  varPart_colnames <- colnames(rdtSet$analSet$varPart.df)[-1] # Exclude the symbol column
  
  return(varPart_colnames)
}

# Function to get the name of the CSV file containing the variance partition results
GetVarFileName <- function() {
  rdtSet <- .get.rdt.set()
  
  # Assuming the file name is stored in rdtSet$analSet$varPart.fileName
  if (is.null(rdtSet$analSet$varPart.fileName)) {
    stop("Variance partition results file name not found.")
  }
  
  return(rdtSet$analSet$varPart.fileName)
}

# Function to plot variance partitioning for an individual feature
PlotVarPartFeature <- function(feature_name, symbol, fileName = "varpart_feature_plot",  format = "png", dpi=72) {
  rdtSet <- .get.rdt.set()
  
  # Ensure that varPart is available
  if (is.null(rdtSet$analSet$varPart.df)) {
    stop("Variance partition matrix not found.")
  }
  
  varPart <- rdtSet$analSet$varPart.df
  
  # Check if the feature exists in varPart
  if (!(feature_name %in% rownames(varPart))) {
    stop(paste("Feature", feature_name, "not found in the variance partition matrix."))
  }
  
  # Extract the row corresponding to the selected feature
  feature_varPart <- varPart[feature_name, -1, drop = FALSE]  # Exclude the Symbol column
  
  # Create a data frame for plotting
  feature_varPart_df <- data.frame(
    Covariate = colnames(feature_varPart),
    VarianceExplained = as.numeric(feature_varPart)
  )
  
  # Construct the image filename
  imgName <- paste(fileName, "dpi", dpi, ".", format, sep = "")
  
  # Create and save the plot using Cairo
  Cairo::Cairo(file = imgName, type = format, dpi = dpi, width = 8, height = 6, units = "in", bg = "white")
  
  # Plot the variance partitioning for the selected feature
  p <- ggplot(feature_varPart_df, aes(x = Covariate, y = VarianceExplained)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(title = paste("Variance Partitioning for", feature_name),
         x = "Covariates",
         y = "Variance Explained (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print and save the plot
  print(p)
  
  # Close the Cairo device to save the image
  dev.off()
  
  return(imgName)  # Optionally return the image file name for reference
}
