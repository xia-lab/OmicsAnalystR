##################################################
## R script for OmicsAnalyst
## Description: Compute upset diagram
## Authors: 
## G. Zhou, guangyan.zhou@mail.mcgill.ca
###################################################


#'Prepare data for Upset diagram
#'@param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#'@param fileNm file name of the json file output 
#'@export
PrepareUpsetData <- function(fileNm, metaCol) {
  # Fetch the dataset from the current session
  rdtSet <- get.rdt.set()
  paramSet <- rdtSet$paramSet
  analSet <- rdtSet$analSet
  
  # Extract normalized data and metadata
  dataSet <- rdtSet$dataSet
  gene.mat <- dataSet$norm  # Assumed normalized gene expression matrix
  meta_data <- dataSet$meta.info  # Metadata containing the group information

  # Get the groups based on the provided metadata column
  groups <- meta_data[, metaCol]
  
  if (is.null(groups)) {
    stop(paste("The provided metadata column", metaCol, "does not exist in the dataset."))
  }

  # Ensure it's a factor for pairwise comparison
  groups <- as.factor(groups)
  
  # Initialize list for storing DE gene lists for each comparison
  sel.dats <- list()
  
  # Perform pairwise comparisons between each group
  for (i in 1:(length(levels(groups)) - 1)) {
    for (j in (i + 1):length(levels(groups))) {
      group1 <- levels(groups)[i]
      group2 <- levels(groups)[j]
      comp_label <- paste(group1, "vs", group2, sep = "_")

      # Subset the data for the two groups
      idx <- groups %in% c(group1, group2)
      comp_data <- gene.mat[, idx]
      comp_groups <- groups[idx]

      # Fit a linear model or apply a DE analysis method (example: limma)
      # For example, using limma
      design <- model.matrix(~ comp_groups)
      fit <- lmFit(comp_data, design)
      fit <- eBayes(fit)
      
      # Get DE genes for this comparison
      topTableRes <- topTable(fit, coef = 2, number = Inf, sort.by = "p")
      sig_genes <- rownames(topTableRes[topTableRes$adj.P.Val < 0.05, ])
      
      # Store the significant genes for this comparison
      sel.dats[[comp_label]] <- sig_genes
    }
  }

  # If no significant genes found, stop
  if (length(sel.dats) == 0) {
    stop("No significant features for any dataset!")
  }

  # Prepare data for upset plot
  sums <- unlist(lapply(sel.dats, length))
  names <- unlist(lapply(sel.dats, paste, collapse = ", "))
  
  require(reshape2)
  df <- reshape::melt(sel.dats, value.name = "id")
  colnames(df) <- c("name", 'set')
  uniq.nms <- unique(df$name)
  new.df <- dcast(df, name ~ set, value.var = 'set', fill = 0)
  rownames(new.df) <- new.df[, 1]
  new.df <- new.df[, -1, drop = FALSE]

  # Gene mapping
  gene.map <- queryGeneDB("entrez", paramSet$data.org)
  gene.map[] <- lapply(gene.map, as.character)

  # Prepare the json data for the upset plot
  json.list <- list()
  for (i in 1:nrow(new.df)) {
    json.list[[i]] <- list()
    json.list[[i]][["sets"]] <- new.df[i,][new.df[i,] != 0]
    entrez.vec <- rownames(new.df)[i]
    hit.inx <- match(entrez.vec, gene.map[, "gene_id"])
    symbols <- gene.map[hit.inx, "symbol"]

    # If not gene symbol, use id by itself
    na.inx <- is.na(symbols)
    symbols[na.inx] <- entrez.vec[na.inx]
    json.list[[i]][["name"]] <- symbols
    json.list[[i]][["entrez"]] <- entrez.vec
  }

  # Set colors for the upset plot
  col.vec <- gg_color_hue(length(sel.dats))

  # Save JSON data for the upset plot
  jsonNm <- paste0(fileNm, ".json")
  json.mat <- RJSONIO::toJSON(list(json.list, col.vec))
  sink(jsonNm)
  cat(json.mat)
  sink()

  return(1)
}

#Record upset intersection mode for report
SetUpsetMode <- function(mode){
      paramSet <- readSet(paramSet, "paramSet");
      paramSet$upsetMode <- mode;
      saveSet(paramSet);
}