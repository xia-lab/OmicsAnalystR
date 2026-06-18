PlotPercentBars <- function(dataNm,top_n=10, fileName="", dpi=150, format="png"){
  rdtSet <- .get.rdt.set() 
  varPart <- rdtSet$analSet$varPart.df[,-c(1:2)];
  vp <- varPart[order(varPart[, 1], decreasing = TRUE),,drop=F ]#sortCols(varPart);
  imgName = paste(fileName, "dpi", dpi, ".", format, sep="");
 
   if (top_n < 10){
      h <- 6;
    } else if (top_n < 15){
      h <- top_n/1.6;
    } else if (top_n < 20){
      h <- top_n/1.8;
    } else if (top_n < 25){
      h <- top_n/2;
    } else if (top_n < 30){
      h <- top_n/2.2;
    } else if (top_n < 40){
      h <- top_n/2.5;
    } else {
      h <- top_n/4.5;
    };

  Cairo::Cairo(file = imgName, type = format, dpi = dpi, width = 10, height = h, units = "in", bg = "white")
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
   
  dataSet <- readDataset(dataNm);
  dataSet$varPar$topNum <- top_n; 
  RegisterData(dataSet)
  return(1)
}

PerformVarPartOverview <- function(selectedData,selMeta, top_n = 500, fileName = "variance_partition_plot", dpi = 300, format = "png",color) {
  library(variancePartition)
  #library(limma)
  library(Cairo)
  library(car)
   rdtSet <- .get.rdt.set() 
  # Ensure that normalized data and metadata are available
 
#  sel.inx <- mdata.all==1;
 # sel.nms <- names(mdata.all)[sel.inx];
 
dataSet <- readDataset(selectedData);
  df =  dataSet[["data.proc"]]
 sanitized_names <- gsub("[[:cntrl:]]|[^[:ascii:]]", "_", rownames(df), perl = TRUE)  
  sanitized_names <- names(dataSet[["enrich_ids"]])[match(sanitized_names,dataSet$enrich_ids)]
  uniqFeats <- paste0( dataSet$type,".",sanitized_names)
  rownames(df) <- uniqFeats;
 

# types <- unlist(lapply(dataSetList, function(x) return(x$type)))
# types <- paste(paste0("^",types,"."),collapse  = "|")

  # Extract normalized gene expression data
  gene_expr <- df
  feature_variances <- apply( gene_expr, 1, var)
  ordinx <- order(feature_variances, decreasing = TRUE)
  gene_expr <-   gene_expr[ordinx,]
  
  meta_data <- rdtSet$dataSet$meta.info
  meta_types <- rdtSet$dataSet$meta.types  # Continuous ("cont") or discrete ("disc") types
  
  if(nrow(gene_expr)>top_n){
    gene_expr <-  gene_expr[1:top_n,]
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
    AddErrMsg("Sample names between gene expression data and metadata do not match.");
    return(0);
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
    AddErrMsg("At least one fixed effect must be specified.");
    return(0);
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
    AddErrMsg("No valid covariates found in metadata."); return(0);
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


  if(ncol(varPart)>1){
  varPart <- varPart[, c(selMeta, setdiff(names(varPart),selMeta))]
  }
 
  # Plot the variance partitioning results and save the image
  imgName = paste(fileName, "dpi", dpi, ".", format, sep="");
  
  Cairo::Cairo(file = imgName, type = format, dpi = dpi, width = 8, height = 6, units = "in", bg = "white")
  p <- plotVarPart(varPart)
  print(p)
  dev.off()

  symbols <- gsub(paste0(dataSet$type,"."),"",rownames(varPart))
  ids <- unname(dataSet$enrich_ids[match(symbols,names(dataSet$enrich_ids))])
  # Store the top gene results in rdtSet for future use
  rdtSet$analSet$varPart.symbols <- symbols;
  rdtSet$analSet$varPart.ids <-  ids;
  varPart<-cbind(Symbol = rdtSet$analSet$varPart.symbols,  ID=rdtSet$analSet$varPart.ids, varPart);
  rdtSet$analSet$varPart.fileName <- "varPart_results.csv";
  rdtSet$analSet$varPart.df <- varPart; 
  dataSet$varPar$symbols <-  symbols;
  dataSet$varPar$ids <- ids;
  dataSet$varPar$topNum <- top_n;
  dataSet$varPar$varPart.df <- varPart;
  fast.write.csv(varPart, file="varPart_results.csv")
 
    RegisterData(dataSet)
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
  print(c(mmeta.all ,"mmeta.all "))
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
    AddErrMsg("Variance partition matrix not found."); return(0);
  }
  
  # Extract the variance partitioning matrix without the symbol column
  varPart_matrix <- as.matrix(rdtSet$analSet$varPart.df[ , -c(1:2), drop = FALSE]) # Removing the symbol column
  
  return(varPart_matrix)
}

# Function to get the row names (gene IDs) of the variance partitioning matrix
GetVarIds <- function() {
  rdtSet <- .get.rdt.set()
  
  # Assuming varPart is stored in rdtSet$analSet$varPart.df
  if (is.null(rdtSet$analSet$varPart.df)) {
    AddErrMsg("Variance partition matrix not found."); return(0);
  }
 
  # Return the row names (gene IDs)
   varPart_ids <-  rdtSet$analSet$varPart.ids
 #  varPart_ids <- gsub(rdtSet$dataSet$types,"", varPart_ids )
  return(varPart_ids)
}

# Function to get the gene symbols associated with the variance partitioning matrix
GetVarSymbols <- function() {
  rdtSet <- .get.rdt.set()
  
  # Assuming varPart is stored in rdtSet$analSet$varPart.df
  if (is.null(rdtSet$analSet$varPart.df)) {
    AddErrMsg("Variance partition matrix not found."); return(0);
  }
  
  # Extract the symbols column
  varPart_symbols <- rdtSet$analSet$varPart.df$ID
  return(varPart_symbols)
}

GetVarLables <- function() {
  rdtSet <- .get.rdt.set()
  
  # Assuming varPart is stored in rdtSet$analSet$varPart.df
  if (is.null(rdtSet$analSet$varPart.df)) {
    AddErrMsg("Variance partition matrix not found."); return(0);
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
    AddErrMsg("Variance partition matrix not found."); return(0);
  }
  
  # Extract the column names (excluding the symbol column)
  varPart_colnames <- colnames(rdtSet$analSet$varPart.df)[-c(1:2)] # Exclude the symbol column
  
  return(varPart_colnames)
}

# Function to get the name of the CSV file containing the variance partition results
GetVarFileName <- function() {
  rdtSet <- .get.rdt.set()
  
  # Assuming the file name is stored in rdtSet$analSet$varPart.fileName
  if (is.null(rdtSet$analSet$varPart.fileName)) {
    AddErrMsg("Variance partition results file name not found."); return(0);
  }
  
  return(rdtSet$analSet$varPart.fileName)
}

# Function to plot variance partitioning for an individual feature
PlotVarPartFeature <- function(feature_name, symbol, fileName = "varpart_feature_plot",  format = "png", dpi=150) {
  rdtSet <- .get.rdt.set()
  
  # Ensure that varPart is available
  if (is.null(rdtSet$analSet$varPart.df)) {
    AddErrMsg("Variance partition matrix not found."); return(0);
  }
  
  varPart <- rdtSet$analSet$varPart.df
  
  # Check if the feature exists in varPart
  if (!(feature_name %in% rownames(varPart))) {
    AddErrMsg(paste("Feature", feature_name, "not found in the variance partition matrix.")); return(0);
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

# ─── Multi-omics VPA (vegan::varpart) helpers ────────────────────────────────

PrepareVpaInput <- function(x1Name, x2Name, traitCols, pcaReduce, nPcsX1, nPcsX2,
                            fileName, dpi, format) {
  rdtSet <- .get.rdt.set()

  ds1 <- readDataset(x1Name)
  ds2 <- readDataset(x2Name)
  if (is.null(ds1) || is.null(ds1$data.proc)) {
    AddErrMsg(paste("Dataset not found:", x1Name)); return(0)
  }
  if (is.null(ds2) || is.null(ds2$data.proc)) {
    AddErrMsg(paste("Dataset not found:", x2Name)); return(0)
  }

  # data.proc is features × samples; use colnames directly as sample IDs before transposing
  samp1 <- colnames(ds1$data.proc)
  samp2 <- colnames(ds2$data.proc)
  if (length(samp1) == 0) {
    AddErrMsg(paste("Dataset", x1Name, "has no column names — was it processed correctly?")); return(0)
  }
  if (length(samp2) == 0) {
    AddErrMsg(paste("Dataset", x2Name, "has no column names — was it processed correctly?")); return(0)
  }

  meta <- rdtSet$dataSet$meta.info
  meta_samp <- rownames(meta)

  common_samples <- Reduce(intersect, list(samp1, samp2, meta_samp))
  cat(sprintf(">>> VPA PrepareInput: X1 n=%d, X2 n=%d, meta n=%d, common=%d\n",
              length(samp1), length(samp2), length(meta_samp), length(common_samples)))
  if (length(common_samples) < 5) {
    AddErrMsg(paste0(
      "Too few common samples (", length(common_samples), "). ",
      "Check that sample names match across datasets and metadata. ",
      "X1 samples (first 3): ", paste(head(samp1, 3), collapse=", "), "; ",
      "meta samples (first 3): ", paste(head(meta_samp, 3), collapse=", ")
    ))
    return(0)
  }

  # Diagnostic: log names at both samp-detection and subset points
  cat(sprintf(">>> VPA: common_samples[1:3]=%s\n", paste(head(common_samples, 3), collapse=",")))
  cat(sprintf(">>> VPA: ds1 dim=%s col[1:3]=%s row[1:3]=%s\n",
              paste(dim(ds1$data.proc), collapse="x"),
              paste(head(colnames(ds1$data.proc), 3), collapse=","),
              paste(head(rownames(ds1$data.proc), 3), collapse=",")))

  # Orientation-agnostic subset: try columns first (features×samples), then rows (samples×features)
  orient_and_subset <- function(dp, samples) {
    col_sel <- intersect(colnames(dp), samples)
    if (length(col_sel) > 0) {
      t(as.matrix(dp[, col_sel, drop = FALSE]))  # features×samples → samples×features
    } else {
      row_sel <- intersect(rownames(dp), samples)
      if (length(row_sel) == 0) stop("Sample names not found in either rows or columns of data.proc")
      cat(sprintf(">>> VPA: using ROW orientation (samples×features), %d rows match\n", length(row_sel)))
      as.matrix(dp[row_sel, , drop = FALSE])     # already samples×features
    }
  }

  X1_mat <- orient_and_subset(ds1$data.proc, common_samples)
  X2_mat <- orient_and_subset(ds2$data.proc, common_samples)
  cat(sprintf(">>> VPA: matrices — X1 %s, X2 %s\n",
              paste(dim(X1_mat), collapse="x"), paste(dim(X2_mat), collapse="x")))

  trait_cols <- intersect(traitCols, colnames(meta))
  if (length(trait_cols) == 0) {
    AddErrMsg("No matching trait columns found in metadata."); return(0)
  }

  # Build Y_mat: numeric columns kept; binary categorical → 0/1 dummy; multi-class → dropped
  Y_raw <- meta[common_samples, trait_cols, drop = FALSE]
  # Authoritative per-column typing from OmicsAnalyst (meta.types: "cont"/"disc",
  # set at upload and finalized by SanityCheckMeta), with value inference as fallback.
  mtypes <- rdtSet$dataSet$meta.types
  Y_list <- lapply(trait_cols, function(col) {
    v <- Y_raw[, col]
    mtype <- if (!is.null(mtypes) && col %in% names(mtypes)) unname(mtypes[[col]]) else NA_character_
    v_char <- as.character(v)
    nonblank <- !is.na(v_char) & nzchar(trimws(v_char))
    v_num <- suppressWarnings(as.numeric(v_char))
    parses_num <- any(nonblank) && !any(is.na(v_num[nonblank]))
    is_cont <- identical(mtype, "cont") || (is.na(mtype) && (is.numeric(v) || parses_num))
    if (is_cont) {
      cat(sprintf(">>> VPA: column '%s' → numeric (%s)\n", col, if (is.na(mtype)) "inferred" else "cont"))
      return(v_num)
    }
    lvls <- unique(na.omit(v_char[nonblank]))
    if (length(lvls) == 2) {
      cat(sprintf(">>> VPA: column '%s' → 0/1 (%s=0, %s=1)\n", col, lvls[1], lvls[2]))
      return(ifelse(nonblank, ifelse(v_char == lvls[1], 0, 1), NA_real_))
    }
    cat(sprintf(">>> VPA: dropping column '%s' — %d-level categorical (disc)\n", col, length(lvls)))
    return(NULL)
  })
  names(Y_list) <- trait_cols
  Y_list <- Filter(Negate(is.null), Y_list)
  if (length(Y_list) == 0) {
    AddErrMsg("No numeric or binary trait columns remain for VPA."); return(0)
  }
  Y_mat <- do.call(cbind, Y_list)
  if (!is.matrix(Y_mat)) Y_mat <- matrix(Y_mat, ncol = 1, dimnames = list(common_samples, names(Y_list)))
  rownames(Y_mat) <- common_samples
  trait_cols <- names(Y_list)

  # Drop samples with missing trait values
  ok <- complete.cases(Y_mat)
  if (!all(ok)) {
    n_removed <- sum(!ok)
    Y_mat   <- Y_mat[ok,   , drop = FALSE]
    X1_mat  <- X1_mat[ok,  , drop = FALSE]
    X2_mat  <- X2_mat[ok,  , drop = FALSE]
    warning(paste("Removed", n_removed, "samples with missing trait values."))
  }

  x1_name <- x1Name; x2_name <- x2Name
  pca_reduce <- as.logical(pcaReduce)
  n_pcs_x1 <- as.integer(nPcsX1); n_pcs_x2 <- as.integer(nPcsX2)

  save(X1_mat, X2_mat, Y_mat, x1_name, x2_name, pca_reduce, n_pcs_x1, n_pcs_x2,
       trait_cols, fileName, dpi, format, file = "vpa_input.RData")
  return(1)
}

LoadVpaResults <- function() {
  if (!file.exists("vpa_output.RData")) {
    AddErrMsg("VPA output not found — the analysis may have failed."); return(0)
  }
  e <- new.env(parent = emptyenv())
  load("vpa_output.RData", envir = e)
  rdtSet <- .get.rdt.set()
  rdtSet$analSet$vpa <- list(
    results  = e$results_df,
    x1_name  = e$x1_name,
    x2_name  = e$x2_name,
    a_val    = e$a_val,
    b_val    = e$b_val,
    c_val    = e$c_val,
    d_val    = e$d_val,
    performed = TRUE
  )
  .set.rdt.set(rdtSet)
  return(1)
}

# ─── RunVpa: core multi-omics variance partitioning ──────────────────────────
# Decomposes variation in phenotypic traits (Y) across two omics blocks (X1, X2)
# using vegan::varpart() with partial-RDA permutation tests (n=999).
# High-dimensional blocks are compressed to top PCs before partitioning.
# Computation runs in an isolated RSclient subprocess (vegan/Cairo isolated).
# Results stored in rdtSet$analSet$vpa. Returns 1L on success, 0L on failure.

# Return the comma-separated names of the supplied traits that are MULTI-CLASS
# categorical (>2 levels) under OmicsAnalyst's meta.types typing — i.e. the ones
# VPA would one-hot encode (with the Euclidean-response caveat). Continuous and
# binary traits are NOT returned. Used by the UI to warn before running. Empty
# string when none.
GetVpaCategoricalTraits <- function(traitVec) {
  rdtSet <- tryCatch(.get.rdt.set(), error = function(e) NULL)
  if (is.null(rdtSet)) return("")
  meta <- rdtSet$dataSet$meta.info
  if (is.null(meta) || nrow(meta) == 0L) return("")
  mtypes <- rdtSet$dataSet$meta.types
  cats <- character(0)
  for (col in traitVec) {
    if (!col %in% colnames(meta)) next
    v   <- meta[[col]]
    mtype <- if (!is.null(mtypes) && col %in% names(mtypes)) unname(mtypes[[col]]) else NA_character_
    v_c <- as.character(v)
    nonblank <- !is.na(v_c) & nzchar(trimws(v_c))
    v_num <- suppressWarnings(as.numeric(v_c))
    parses_num <- any(nonblank) && !any(is.na(v_num[nonblank]))
    is_cont <- identical(mtype, "cont") || (is.na(mtype) && (is.numeric(v) || parses_num))
    if (is_cont) next
    lvls <- unique(na.omit(v_c[nonblank]))
    if (length(lvls) > 2L) cats <- c(cats, col)
  }
  paste(cats, collapse = ",")
}

RunVpa <- function(x1Name, x2Name, traitCols = NULL,
                   pcaReduce = TRUE, nPcsX1 = 10L, nPcsX2 = 10L,
                   fileName = "vpa_venn_0_", dpi = 150L, format = "png",
                   includeCat = FALSE,
                   barName = NULL, rdaName = NULL, corrName = NULL,
                   colorBy = NULL, showLabels = FALSE) {

  rdtSet <- .get.rdt.set()
  meta   <- rdtSet$dataSet$meta.info
  if (is.null(meta) || nrow(meta) == 0L) {
    AddErrMsg("No metadata (rdtSet$dataSet$meta.info) for VPA."); return(0L)
  }

  # ── 1. Load data.proc from both datasets ───────────────────────────────────
  ds1 <- readDataset(x1Name); ds2 <- readDataset(x2Name)
  if (is.null(ds1$data.proc)) { AddErrMsg(paste("data.proc not found for:", x1Name)); return(0L) }
  if (is.null(ds2$data.proc)) { AddErrMsg(paste("data.proc not found for:", x2Name)); return(0L) }
  dp1 <- ds1$data.proc; dp2 <- ds2$data.proc

  # ── 2. Orientation-agnostic subset (data.proc is features×samples) ─────────
  orient_and_subset <- function(dp, samples) {
    col_sel <- intersect(colnames(dp), samples)
    if (length(col_sel) > 0L) return(t(as.matrix(dp[, col_sel, drop=FALSE])))
    row_sel <- intersect(rownames(dp), samples)
    if (length(row_sel) == 0L) stop("sample names not found in data.proc rows or columns")
    as.matrix(dp[row_sel, , drop=FALSE])
  }

  # ── 3. Common samples ────────────────────────────────────────────────────────
  # data.proc is features × samples and features >> samples, so the sample-ID
  # axis is whichever dimension overlaps the metadata samples — NOT ncol>=nrow,
  # which picks the FEATURE axis for omics data (thousands of features, dozens of
  # samples) and collapses the intersection to 0.
  meta_samp <- rownames(meta)
  pick_samples <- function(dp) {
    if (length(intersect(colnames(dp), meta_samp)) >= length(intersect(rownames(dp), meta_samp)))
      colnames(dp) else rownames(dp)
  }
  samp1 <- pick_samples(dp1)
  samp2 <- pick_samples(dp2)
  common <- Reduce(intersect, list(samp1, samp2, meta_samp))
  cat(sprintf(">>> RunVpa: X1 n=%d, X2 n=%d, meta n=%d, common=%d\n",
              length(samp1), length(samp2), nrow(meta), length(common)))
  if (length(common) < 5L) {
    AddErrMsg(sprintf(
      paste0("Variance partitioning needs >= 5 samples shared by both omics layers and the ",
             "metadata, but only %d match. %s vs %s share these sample IDs (column headers): ",
             "%s | metadata sample IDs (row names): %s. Make the sample IDs identical across all ",
             "tables and the metadata."),
      length(common), x1Name, x2Name,
      paste(utils::head(samp1, 4), collapse = ", "),
      paste(utils::head(meta_samp, 4), collapse = ", ")))
    return(0L)
  }

  X1_mat <- tryCatch(orient_and_subset(dp1, common),
                     error = function(e) { AddErrMsg(paste("X1 orient error:", e$message)); NULL })
  X2_mat <- tryCatch(orient_and_subset(dp2, common),
                     error = function(e) { AddErrMsg(paste("X2 orient error:", e$message)); NULL })
  if (is.null(X1_mat) || is.null(X2_mat)) return(0L)
  cat(sprintf(">>> RunVpa: X1 %dx%d, X2 %dx%d\n",
              nrow(X1_mat), ncol(X1_mat), nrow(X2_mat), ncol(X2_mat)))

  # ── 4. Trait (Y) matrix — binary categorical → 0/1; multi-class → drop ──────
  if (is.null(traitCols) || length(traitCols) == 0L)
    traitCols <- colnames(meta)
  traitCols <- intersect(traitCols, colnames(meta))
  if (length(traitCols) == 0L) {
    AddErrMsg("No matching trait columns in metadata."); return(0L)
  }

  Y_raw   <- meta[common, traitCols, drop=FALSE]
  # Drive trait typing off OmicsAnalyst's authoritative per-column classification
  # (rdtSet$dataSet$meta.types: "cont" vs "disc", computed by GetNumbericalInx /
  # GetDiscreteInx at upload). This keeps a continuous trait stored as character/
  # factor numeric, and keeps a numeric-CODED categorical (e.g. a subject id)
  # categorical. Falls back to value inference when no classification is present.
  mtypes <- rdtSet$dataSet$meta.types
  Y_cols  <- lapply(traitCols, function(col) {
    v <- Y_raw[[col]]
    mtype <- if (!is.null(mtypes) && col %in% names(mtypes)) unname(mtypes[[col]]) else NA_character_
    v_c  <- as.character(v)
    nonblank <- !is.na(v_c) & nzchar(trimws(v_c))
    v_num <- suppressWarnings(as.numeric(v_c))
    parses_num <- any(nonblank) && !any(is.na(v_num[nonblank]))
    is_cont <- identical(mtype, "cont") || (is.na(mtype) && (is.numeric(v) || parses_num))
    if (is_cont) {
      cat(sprintf(">>> RunVpa: column '%s' → numeric (%s)\n", col, if (is.na(mtype)) "inferred" else "cont"))
      return(v_num)
    }
    lvls <- unique(na.omit(v_c[nonblank]))
    if (length(lvls) == 2L) {
      cat(sprintf(">>> RunVpa: column '%s' → 0/1 (%s=0)\n", col, lvls[1L]))
      return(ifelse(nonblank, ifelse(v_c == lvls[1L], 0, 1), NA_real_))
    }
    if (length(lvls) > 2L && isTRUE(includeCat)) {
      # One-hot a multi-level categorical into k-1 reference-coded 0/1 dummy columns.
      # CAVEAT: RDA treats the response as Euclidean, so a nominal multi-class trait is
      # an approximation (it models group-membership indicators); for true group
      # separation DIABLO is more appropriate. Enabled only when the user opts in.
      ref <- lvls[1L]
      dummies <- lapply(lvls[-1L], function(L)
        ifelse(nonblank, ifelse(v_c == L, 1, 0), NA_real_))
      mm <- do.call(cbind, dummies)
      colnames(mm) <- paste0(col, "=", lvls[-1L])
      cat(sprintf(">>> RunVpa: column '%s' → one-hot (%d dummies, ref=%s)\n",
                  col, length(lvls) - 1L, ref))
      return(mm)
    }
    cat(sprintf(">>> RunVpa: dropping '%s' — %d-level categorical (disc; enable 'include categorical' to one-hot)\n",
                col, length(lvls)))
    NULL
  })
  names(Y_cols) <- traitCols
  Y_cols <- Filter(Negate(is.null), Y_cols)
  if (length(Y_cols) == 0L) {
    AddErrMsg("No numeric or binary trait columns for VPA after encoding."); return(0L)
  }
  Y_mat <- do.call(cbind, Y_cols)
  if (!is.matrix(Y_mat)) Y_mat <- matrix(Y_mat, ncol=1L)
  rownames(Y_mat) <- common

  # ── 5. Complete cases ─────────────────────────────────────────────────────────
  keep <- complete.cases(X1_mat, X2_mat, Y_mat)
  if (sum(keep) < 5L) {
    AddErrMsg(sprintf("Only %d complete cases — too few to partition variance.", sum(keep)))
    return(0L)
  }
  X1_mat <- X1_mat[keep, , drop=FALSE]
  X2_mat <- X2_mat[keep, , drop=FALSE]
  Y_mat  <- Y_mat[keep,  , drop=FALSE]
  cat(sprintf(">>> RunVpa: n=%d complete cases, Y=%d trait(s)\n", nrow(X1_mat), ncol(Y_mat)))

  # Primary grouping for the ordination plot: the first DISCRETE trait (2..12
  # levels) among the requested columns, aligned to the kept samples. NULL when
  # none (the ordination then uses a single colour).
  grp_vec <- NULL
  kept_samples <- rownames(Y_mat)
  for (col in intersect(traitCols, colnames(meta))) {
    mtype <- if (!is.null(mtypes) && col %in% names(mtypes)) unname(mtypes[[col]]) else NA_character_
    vc <- as.character(meta[kept_samples, col])
    lv <- unique(vc[!is.na(vc) & nzchar(trimws(vc))])
    if (!identical(mtype, "cont") && length(lv) >= 2L && length(lv) <= 12L) { grp_vec <- vc; break }
  }

  # Ordination colour vector: an explicit colorBy metadata column overrides the
  # auto-picked grp_vec. Resolve values aligned to the kept samples + the type.
  cby_vals <- NULL; cby_type <- NA_character_; cby_label <- "Group"
  cby_col <- if (!is.null(colorBy) && nzchar(colorBy) && colorBy %in% colnames(meta)) colorBy else {
    # fall back to the first discrete trait already picked into grp_vec, else first meta col
    if (!is.null(grp_vec)) NA_character_ else colnames(meta)[1]
  }
  if (!is.null(grp_vec) && (is.null(colorBy) || !nzchar(colorBy))) {
    cby_vals <- grp_vec; cby_type <- "disc"
    cby_label <- "Group"
  } else if (!is.na(cby_col) && cby_col %in% colnames(meta)) {
    cby_vals <- as.character(meta[kept_samples, cby_col])
    mtp <- if (!is.null(mtypes) && cby_col %in% names(mtypes)) unname(mtypes[[cby_col]]) else NA_character_
    nb  <- !is.na(cby_vals) & nzchar(trimws(cby_vals))
    vnum <- suppressWarnings(as.numeric(cby_vals))
    is_cont <- identical(mtp,"cont") || (is.na(mtp) && any(nb) && !any(is.na(vnum[nb])) && length(unique(cby_vals[nb]))>12L)
    cby_type <- if (is_cont) "cont" else "disc"
    cby_label <- cby_col
  }

  # ── 6. Output paths ───────────────────────────────────────────────────────────
  homeDir <- tryCatch(rdtSet$paramSet$home.path, error=function(e) getwd())
  if (is.null(homeDir) || !nzchar(homeDir)) homeDir <- getwd()
  # Filename must match SessionBean1.getCurrentImage(): "<getNewImage>dpi150.png".
  # Sibling figures reuse the Venn's numbering by swapping the "venn" key for the
  # bar / rda / corr keys, so the whole figure set shares one counter. These are
  # surfaced via the manifest (PerformOaVpa), NOT hard-wired into any view.
  derive_nm <- function(base, key)
    if (grepl("venn", base, fixed=TRUE)) sub("venn", key, base, fixed=TRUE) else paste0(base, key, "_")
  if (is.null(barName))  barName  <- derive_nm(fileName, "bar")
  if (is.null(rdaName))  rdaName  <- derive_nm(fileName, "rda")
  if (is.null(corrName)) corrName <- derive_nm(fileName, "corr")
  imgPath  <- file.path(homeDir, paste0(fileName, "dpi", dpi, ".", format))
  barPath  <- file.path(homeDir, paste0(barName,  "dpi", dpi, ".", format))
  rdaPath  <- file.path(homeDir, paste0(rdaName,  "dpi", dpi, ".", format))
  corrPath <- file.path(homeDir, paste0(corrName, "dpi", dpi, ".", format))
  csvPath  <- file.path(homeDir, "vpa_results.csv")

  # ── 7. Isolated subprocess: vegan + Cairo ─────────────────────────────────────
  vpa_res <- tryCatch(
    rsclient_isolated_exec(
      func_body = function(inp) {
        X1_mat    <- inp$X1_mat;    X2_mat    <- inp$X2_mat;    Y_mat  <- inp$Y_mat
        grp_vec   <- inp$grp_vec
        pcaReduce <- inp$pcaReduce; nPcsX1    <- inp$nPcsX1;    nPcsX2 <- inp$nPcsX2
        x1_name   <- inp$x1_name;  x2_name   <- inp$x2_name
        imgPath   <- inp$imgPath;  csvPath   <- inp$csvPath
        barPath   <- inp$barPath;  rdaPath   <- inp$rdaPath;   corrPath <- inp$corrPath
        dpi       <- inp$dpi;      fmt        <- inp$format

        impute_cm <- function(m) {
          for (j in seq_len(ncol(m))) {
            nas <- is.na(m[, j])
            if (any(nas)) m[nas, j] <- mean(m[!nas, j], na.rm=TRUE)
          }
          m
        }
        nzv_rm <- function(m) {
          vv <- apply(m, 2L, var, na.rm=TRUE)
          m[, !is.na(vv) & vv > 0, drop=FALSE]
        }
        X1_mat <- nzv_rm(impute_cm(X1_mat))
        X2_mat <- nzv_rm(impute_cm(X2_mat))
        if (ncol(X1_mat) == 0L) stop("X1: no variable features after NZV filter")
        if (ncol(X2_mat) == 0L) stop("X2: no variable features after NZV filter")

        yv <- apply(Y_mat, 2L, var, na.rm=TRUE)
        Y_mat <- Y_mat[, !is.na(yv) & yv > 0, drop=FALSE]
        if (ncol(Y_mat) == 0L) stop("all Y trait columns are constant after imputation")

        do_pca <- function(m, k) {
          k <- min(as.integer(k), ncol(m), nrow(m) - 1L)
          prcomp(m, center=TRUE, scale.=FALSE, rank.=k)$x[, seq_len(k), drop=FALSE]
        }
        if (isTRUE(pcaReduce)) {
          X1_use <- do_pca(X1_mat, nPcsX1)
          X2_use <- do_pca(X2_mat, nPcsX2)
          cat(sprintf(">>> RunVpa: PCA X1→%dPC, X2→%dPC\n", ncol(X1_use), ncol(X2_use)))
        } else {
          X1_use <- X1_mat; X2_use <- X2_mat
        }

        library(vegan)
        X1_df <- as.data.frame(X1_use); X2_df <- as.data.frame(X2_use)
        vp    <- vegan::varpart(Y_mat, X1_df, X2_df)
        cat(">>> RunVpa: varpart done\n")

        safe_p <- function(expr) tryCatch(
          anova(expr, permutations=999L)$`Pr(>F)`[1L], error=function(e) NA_real_)
        p1  <- safe_p(vegan::rda(Y_mat ~ . + Condition(as.matrix(X2_df)), data=X1_df))
        p2  <- safe_p(vegan::rda(Y_mat ~ . + Condition(as.matrix(X1_df)), data=X2_df))
        p12 <- safe_p(vegan::rda(Y_mat ~ ., data=as.data.frame(cbind(X1_df, X2_df))))

        frac  <- vp$part$indfract
        # vegan's adjusted-R2 column name is version-dependent ("Adj.R.squared" in
        # current vegan, "Adj.R.square" in older); match it flexibly and force the
        # fractions to numeric — the column can come back as a factor/character,
        # which breaks the downstream arithmetic (max/round/*100).
        adj_col <- grep("Adj", colnames(frac), ignore.case = TRUE, value = TRUE)[1L]
        if (is.na(adj_col)) adj_col <- grep("R.?squared?", colnames(frac), value = TRUE)[1L]
        adj_vals <- suppressWarnings(as.numeric(as.character(frac[, adj_col])))
        a_val <- adj_vals[1L]; b_val <- adj_vals[2L]
        c_val <- adj_vals[3L]; d_val <- adj_vals[4L]

        fmt_p    <- function(p) if (is.na(p)) "N/A" else if (p<0.001) "<0.001" else sprintf("%.3f",p)
        sig_star <- function(p) if (is.na(p)) "" else if (p<0.001) "***" else if (p<0.01) "**" else if (p<0.05) "*" else "ns"

        library(Cairo)
        # Short, extension-free block names so the Venn corner labels don't clip.
        short_nm <- function(s) {
          s <- sub("\\.(csv|tsv|txt|tabular)$", "", s, ignore.case = TRUE)
          if (nchar(s) > 18L) paste0(substr(s, 1L, 17L), "…") else s
        }
        xn <- c(short_nm(x1_name), short_nm(x2_name))

        # ── Figure 1: variance-partition Venn ───────────────────────────────────
        Cairo::Cairo(file=imgPath, type=fmt, dpi=dpi, width=8.5, height=7, units="in", bg="white")
        par(mar=c(5, 3, 3, 3), xpd=NA)
        plot(vp, digits=2, Xnames=xn, bg=c("steelblue","tomato"), alpha=80,
             main="Variance Partitioning", cex=0.95)
        mtext(sprintf("[a] Unique %s: %.3f   |   [b] Shared: %.3f   |   [c] Unique %s: %.3f   |   Residual: %.3f",
                      xn[1], max(0,a_val), max(0,b_val), xn[2], max(0,c_val), max(0,d_val)),
              side=1L, line=2.3, cex=0.78)
        mtext(sprintf("Permutation P (n=999):  [a] %s%s   |   [c] %s%s   |   joint %s%s",
                      fmt_p(p1), sig_star(p1), fmt_p(p2), sig_star(p2),
                      fmt_p(p12), sig_star(p12)),
              side=1L, line=3.5, cex=0.72)
        dev.off()

        # ── Figure 2: fraction-importance bar (adj. R2 %, with significance) ─────
        bar_ok <- tryCatch({
          fr  <- pmax(0, c(a_val, b_val, c_val, d_val)) * 100
          sig <- c(sig_star(p1), "", sig_star(p2), "")
          labs<- c(paste0("Unique: ", xn[1]), "Shared", paste0("Unique: ", xn[2]), "Residual")
          cols<- c("steelblue", "mediumpurple3", "tomato", "grey80")
          Cairo::Cairo(file=barPath, type=fmt, dpi=dpi, width=8.5, height=5, units="in", bg="white")
          par(mar=c(4.5, 13, 3, 3))
          bp <- barplot(rev(fr), horiz=TRUE, names.arg=rev(labs), col=rev(cols),
                        las=1, border=NA, xlim=c(0, max(fr)*1.2 + 2),
                        xlab="Variance explained (adj. R², %)",
                        main="VPA fraction importance", cex.names=0.82)
          text(rev(fr), bp,
               labels=paste0(sprintf("%.1f%%", rev(fr)),
                             ifelse(rev(sig)=="", "", paste0("  ", rev(sig)))),
               pos=4, cex=0.9, xpd=NA)
          dev.off(); TRUE
        }, error=function(e){ cat(">>> RunVpa: bar plot skipped:", conditionMessage(e), "\n"); FALSE })

        # ── Figure 3: sample ordination — MetaboAnalyst PCA score-plot style (ggplot2) ──
        # Filled group-coloured points (shape 21) with shaded 95% confidence ellipses
        # (categorical colorBy) or a viridis gradient (continuous colorBy), fitted trait
        # vectors (envfit), optional sample labels — mirrors PlotPCA2DScore / WfPlotSamplePCA.
        rda_ok <- tryCatch({
          library(ggplot2)
          ord <- vegan::rda(cbind(X1_df, X2_df))
          evv <- ord$CA$eig; evv <- evv[evv > 0]
          p1v <- if (length(evv) >= 1L) 100*evv[1]/sum(evv) else NA_real_
          p2v <- if (length(evv) >= 2L) 100*evv[2]/sum(evv) else NA_real_
          si  <- as.data.frame(vegan::scores(ord, display="sites", choices=1:2))
          colnames(si) <- c("PC1", "PC2")
          si$OVlab <- if (!is.null(inp$samp_names) && length(inp$samp_names) == nrow(si))
                        inp$samp_names else rownames(si)
          cby <- inp$cby_vals; cbt <- inp$cby_type
          clab <- if (!is.null(inp$cby_label)) inp$cby_label else "Group"
          has_color <- !is.null(cby) && length(cby) == nrow(si) && !all(is.na(cby))
          mode <- if (has_color && identical(cbt, "cont")) "cont" else if (has_color) "disc" else "none"
          if (mode == "cont") si$OVcol <- suppressWarnings(as.numeric(as.character(cby)))
          if (mode == "disc") si$OVgrp <- factor(as.character(cby))

          axlab <- function(n, pct) if (is.na(pct)) n else sprintf("%s (%.1f%%)", n, pct)
          p <- ggplot(si, aes(x = PC1, y = PC2)) +
            geom_hline(yintercept = 0, linetype = "dashed", colour = "grey85") +
            geom_vline(xintercept = 0, linetype = "dashed", colour = "grey85")
          if (mode == "cont") {
            p <- p + geom_point(aes(fill = OVcol), shape = 21, size = 3.4, colour = "grey30", stroke = 0.4) +
                 scale_fill_viridis_c(name = clab)
          } else if (mode == "disc") {
            pal <- grDevices::hcl.colors(max(2L, nlevels(si$OVgrp)), "Dark 3")
            if (nlevels(si$OVgrp) >= 2L)
              p <- p + stat_ellipse(aes(fill = OVgrp), geom = "polygon", alpha = 0.18,
                                    colour = NA, level = 0.95, show.legend = FALSE, na.rm = TRUE)
            p <- p + geom_point(aes(fill = OVgrp), shape = 21, size = 3.4, colour = "grey25", stroke = 0.4) +
                 scale_fill_manual(name = clab, values = pal)
          } else {
            p <- p + geom_point(shape = 21, size = 3.4, fill = "grey70", colour = "grey30", stroke = 0.4)
          }
          ft <- tryCatch(vegan::envfit(ord, as.data.frame(Y_mat), choices = 1:2, permutations = 0),
                         error = function(e) NULL)
          if (!is.null(ft)) tryCatch({
            mul <- tryCatch(vegan::ordiArrowMul(ft), error = function(e) 1)
            arr <- as.data.frame(vegan::scores(ft, display = "vectors")) * mul
            colnames(arr) <- c("PC1", "PC2"); arr$OVlab <- rownames(arr)
            p <- p +
              geom_segment(data = arr, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                           arrow = grid::arrow(length = grid::unit(0.18, "cm")),
                           colour = "grey20", linewidth = 0.4, inherit.aes = FALSE) +
              geom_text(data = arr, aes(x = PC1 * 1.08, y = PC2 * 1.08, label = OVlab),
                        colour = "grey10", size = 3.2, fontface = "italic", inherit.aes = FALSE)
          }, error = function(e) NULL)
          if (isTRUE(inp$show_labels))
            p <- p + geom_text(aes(label = OVlab), size = 2.2, colour = "grey45", vjust = -0.9, show.legend = FALSE)
          p <- p +
            labs(x = axlab("PC1", p1v), y = axlab("PC2", p2v),
                 title = "Sample ordination (omics PC space)") +
            theme_bw(base_size = 13) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                  legend.position = "right")
          Cairo::Cairo(file = rdaPath, type = fmt, dpi = dpi, width = 8, height = 6.2, units = "in", bg = "white")
          print(p); dev.off(); TRUE
        }, error=function(e){ cat(">>> RunVpa: ordination skipped:", conditionMessage(e), "\n"); FALSE })

        # ── Figure 4: predictor collinearity (leading omics PCs + traits) ───────
        # Match the Data Overview metadata-correlation / RV heatmap style: upper-triangle
        # tiles, reversed-RdYlBu fill on [-1,1], in-cell values, y-axis on the right,
        # vertical legend on the left, coord_fixed + theme_minimal, content-proportional size.
        corr_ok <- tryCatch({
          cmat <- cbind(X1_df[, 1, drop=FALSE], X2_df[, 1, drop=FALSE], Y_mat)
          colnames(cmat) <- c(paste0(xn[1], ".PC1"), paste0(xn[2], ".PC1"), colnames(Y_mat))
          cc <- round(cor(cmat, use="pairwise.complete.obs"), 3)
          k  <- ncol(cc)
          ccu <- cc; ccu[lower.tri(ccu)] <- NA          # upper triangle incl. diagonal
          df <- reshape2::melt(ccu, na.rm = TRUE)
          df$value <- signif(df$value, 3)
          p <- ggplot2::ggplot(df, ggplot2::aes(Var2, Var1, fill = value)) +
            ggplot2::geom_tile(color = "white") +
            ggplot2::scale_y_discrete("Var1", position = "right") +
            ggplot2::scale_fill_gradientn(
              colors = rev(RColorBrewer::brewer.pal(10, "RdYlBu")),
              limits = c(-1, 1), name = "r") +
            ggplot2::geom_text(ggplot2::aes(Var2, Var1, label = value),
                               color = "black", size = 4) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
                           axis.title.x = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.text.y.right = ggplot2::element_text(),
                           legend.direction = "vertical", legend.position = "left") +
            ggplot2::coord_fixed()
          side <- max(3.2, 1.6 + 0.8 * k)
          Cairo::Cairo(file=corrPath, type=fmt, dpi=dpi, width=side, height=side, units="in", bg="white")
          print(p); dev.off(); TRUE
        }, error=function(e){ cat(">>> RunVpa: corr heatmap skipped:", conditionMessage(e), "\n"); FALSE })

        results_df <- data.frame(
          Fraction = c(sprintf("Unique %s [a]",x1_name), "Shared [b]",
                       sprintf("Unique %s [c]",x2_name), "Unexplained [d]"),
          Symbol   = c("[a]","[b]","[c]","[d]"),
          Adj_R2   = round(c(a_val,b_val,c_val,d_val), 4L),
          Percent  = round(c(a_val,b_val,c_val,d_val)*100, 1L),
          P_permutation = c(fmt_p(p1), "—", fmt_p(p2), "—"),
          Significance  = c(sig_star(p1), "", sig_star(p2), ""),
          Note = c(
            sprintf("Partial RDA: unique to %s (other block as covariate)", x1_name),
            if (!is.na(b_val) && b_val < 0) "Negative — collinear blocks (adj. R2 artefact)" else "Shared signal",
            sprintf("Partial RDA: unique to %s (other block as covariate)", x2_name),
            "Residual / unexplained variance"),
          stringsAsFactors = FALSE
        )
        utils::write.csv(results_df, csvPath, row.names=FALSE)

        list(results_df=results_df, vp_summary=capture.output(print(vp)),
             a_val=a_val, b_val=b_val, c_val=c_val, d_val=d_val,
             p1=p1, p2=p2, p12=p12, n_samples=nrow(X1_mat),
             trait_cols=colnames(Y_mat), imgPath=imgPath, csvPath=csvPath,
             barPath=if (isTRUE(bar_ok)) barPath else NA_character_,
             rdaPath=if (isTRUE(rda_ok)) rdaPath else NA_character_,
             corrPath=if (isTRUE(corr_ok)) corrPath else NA_character_,
             x1_name=x1_name, x2_name=x2_name, performed=TRUE)
      },
      input_data = list(
        X1_mat=X1_mat, X2_mat=X2_mat, Y_mat=Y_mat, grp_vec=grp_vec,
        cby_vals=cby_vals, cby_type=cby_type, cby_label=cby_label,
        show_labels=isTRUE(showLabels), samp_names=rownames(Y_mat),
        pcaReduce=isTRUE(pcaReduce), nPcsX1=as.integer(nPcsX1), nPcsX2=as.integer(nPcsX2),
        x1_name=x1Name, x2_name=x2Name,
        imgPath=imgPath, barPath=barPath, rdaPath=rdaPath, corrPath=corrPath,
        csvPath=csvPath, dpi=as.integer(dpi), format=format
      ),
      packages = c("vegan", "Cairo", "qs", "ggplot2"),
      timeout  = 300L
    ),
    error = function(e) { AddErrMsg(paste("VPA subprocess error:", e$message)); NULL }
  )

  if (is.null(vpa_res) || isFALSE(vpa_res$success)) {
    msg <- if (!is.null(vpa_res$message)) vpa_res$message else "VPA subprocess failed"
    AddErrMsg(paste("RunVpa FAILED:", msg)); return(0L)
  }

  rdtSet$analSet$vpa <- vpa_res
  .set.rdt.set(rdtSet)
  cat(sprintf(">>> RunVpa: done — a=%.3f b=%.3f c=%.3f n=%d\n",
              vpa_res$a_val, vpa_res$b_val, vpa_res$c_val, vpa_res$n_samples))
  return(1L)
}
