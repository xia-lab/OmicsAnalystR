##################################################
## R scripts for OmicsAnalyst
## Description: Related to linear modeling
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#' CovariateScatter.Anal
#' @param imgName image name
#' @param imgFormat image format
#' @param analysis.var variable of analysis
#' @param ref reference group
#' @param block block name
#' @param thresh threshold
#' @param contrast.cls contrast group
#' @export
CovariateScatter.Anal <- function(dataName, 
                                  imgName="NA", 
                                  imgFormat="png", 
                                  analysis.var, 
                                  ref = NULL, 
                                  block = "NA", 
                                  thresh=0.05,
                                  contrast.cls = "anova"){

  dataSet <- readDataset(dataName);
  rdtSet <- .get.rdt.set();
  msg.lm <- ""
  # load libraries
  library(limma)
  library(dplyr)
  
  # get inputs
  if(!exists('adj.vec')){
    adj.bool = F;
    vars <- analysis.var;
    covariates.vec <- "NA" #for report generation purpose only
  }else{
    if(length(adj.vec) > 0){
      adj.bool = T;
      vars <- c(analysis.var, adj.vec)
      covariates.vec <- adj.vec;
      if(length(adj.vec) == 1){
        if(adj.vec == ""){
            adj.bool = F;
        }
      }
    }else{
      adj.bool = F;
      vars <- analysis.var;
      covariates.vec <- "NA"
    }
  }

  covariates <- rdtSet$dataSet$meta.info
  covariates <- droplevels(covariates)
  var.types <- rdtSet$dataSet[["meta.types"]]
  feature_table <- dataSet$data.proc;
  covariates <- covariates[which(rownames(covariates) %in% colnames(feature_table)),,drop=F]

  # process inputs
  thresh <- as.numeric(thresh)
  ref <- make.names(ref)
  analysis.type <- unname(rdtSet$dataSet$meta.types[analysis.var]);
  
  # process metadata table (covariates)
   for(i in c(1:length(var.types))){ # ensure all columns are the right type
    if(var.types[i] == "disc"){
      if(class(covariates[,i]) !="factor"){
        covariates[,i] <- covariates[,i] %>% make.names() %>% factor()
      }
    } else {
      if("NA" %in%covariates[,i] | any(is.na( covariates[,i])) ){
        covariates[!is.na(covariates[,i]) & covariates[,i]!="NA" ,i] <- covariates[!is.na(covariates[,i]) & covariates[,i]!="NA",i] %>% as.character() %>% as.numeric()
      }else{
        covariates[,i] <- covariates[,i] %>% as.character() %>% as.numeric()
      }
      
    }
  }
  #subset to samples contained in dataset
  covariates <- covariates[match(colnames(feature_table), rownames(covariates)),,drop=F]
  if (block != "NA"){    
    if(rdtSet$dataSet$meta.types[block] == "cont"){
      AddMsg("Blocking factor can not be continuous data type.")
      return(c(-1,-1));
    }
    # recent update: remove check for unbalanced design. Limma can handle.
  }
  
  sig.num <- 0;

  if(analysis.type == "disc"){
    # build design and contrast matrix
    #covariates[, analysis.var] <- covariates[, analysis.var] %>% make.names() %>% factor();
    grp.nms <- levels(covariates[, analysis.var]);
    design <- model.matrix(formula(paste0("~ 0", paste0(" + ", vars, collapse = ""))), data = covariates);
    colnames(design)[1:length(grp.nms)] <- grp.nms;
    myargs <- list();
    
    # perform specified contrast
    if(contrast.cls == "anova"){
      cntr.cls <- grp.nms[grp.nms != ref];
      myargs <- as.list(paste(cntr.cls, "-", ref, sep = ""));
    } else {
      myargs <- as.list(paste(contrast.cls, "-", ref, sep = ""));
    }
    myargs[["levels"]] <- design;
    contrast.matrix <- do.call(makeContrasts, myargs);
    feature_table <- feature_table[,which(colnames(feature_table) %in% rownames(design)) ]
    # handle blocking factor
    if (block == "NA") {
      fit <- tryCatch({
        lmFit(feature_table, design)
        }, error=function(e){
           msg.lm <- c(msg.lm,e)
        }, warning=function(w){
          msg.lm <- c(msg.lm,w)
        })
    } else {
      block.vec <- covariates[,block];
      corfit <- duplicateCorrelation(feature_table, design, block = block.vec)
      fit <- tryCatch({
        lmFit(feature_table, design, block = block.vec, correlation = corfit$consensus)
      }, error=function(e){
        msg.lm <- c(msg.lm,e)
      }, warning=function(w){
        msg.lm <- c(msg.lm,w)
      })
  }
    
    fit <-  tryCatch({
      contrasts.fit(fit, contrast.matrix);
    }, error=function(e){
      msg.lm <- c(msg.lm,e)
    }, warning=function(w){
      msg.lm <- c(msg.lm,w)
    })
  
    fit <-  tryCatch({
      eBayes(fit);
    }, error=function(e){
      msg.lm <- c(msg.lm,e)
    }, warning=function(w){
      msg.lm <- c(msg.lm,w)
    })
    rest <- topTable(fit, number = Inf);
    
    ### get results with no adjustment
    design <- model.matrix(formula(paste0("~ 0", paste0(" + ", analysis.var, collapse = ""))), data = covariates);
    colnames(design)[1:length(grp.nms)] <- grp.nms;
    myargs[["levels"]] <- design;
    contrast.matrix <- do.call(makeContrasts, myargs);
    fit <- lmFit(feature_table, design)
    fit <- contrasts.fit(fit, contrast.matrix);
    fit <- eBayes(fit);
    res.noadj <- topTable(fit, number = Inf);
    
  } else { 
    covariates[, analysis.var] <- covariates[, analysis.var] %>% as.numeric();
    types <- unname(rdtSet$dataSet$meta.types[vars])
    if(sum(types == "cont") == length(vars)){ #in case of single, continuous variable, must use different intercept or limma will give unreasonable results
      design <- model.matrix(formula(paste0("~", paste0(" + ", vars, collapse = ""))), data = covariates);
    } else {
      design <- model.matrix(formula(paste0("~ 0", paste0(" + ", vars, collapse = ""))), data = covariates);
    }
    
    feature_table <- feature_table[,which(colnames(feature_table) %in% rownames(design)) ]
    # recent update: enable blocking factor for continuous primary metadata
    if (block == "NA") {
      fit <-  tryCatch({
        lmFit(feature_table, design)
      }, error=function(e){
        msg.lm <- c(msg.lm,e)
      }, warning=function(w){
        msg.lm <- c(msg.lm,w)
      })
    } else {
      block.vec <- covariates[,block];
      corfit <- duplicateCorrelation(feature_table, design, block = block.vec)
      fit <- tryCatch({
        lmFit(feature_table, design, block = block.vec, correlation = corfit$consensus)
      }, error=function(e){
         msg.lm <- c(msg.lm,e)
      }, warning=function(w){
        msg.lm <- c(msg.lm,w)
      })
   
 }
    fit <-   tryCatch({
      eBayes(fit);
    }, error=function(e){
      msg.lm <- c(msg.lm,e)
    }, warning=function(w){
      msg.lm <- c(msg.lm,w)
    })
    
    rest <- topTable(fit, number = Inf, coef = analysis.var);
    colnames(rest)[1] <- analysis.var;
    
    ### get results with no adjustment
    design <- model.matrix(formula(paste0("~", analysis.var)), data = covariates);
    
    fit <- eBayes(lmFit(feature_table, design));
    res.noadj <- topTable(fit, number = Inf);
  }

  if(msg.lm!=""){
  AddMsg(msg.lm)
  return(c(-1,-1));
  }

  # make visualization
  adj.mat <- rest[, c("P.Value", "adj.P.Val")]
  noadj.mat <- res.noadj[, c("P.Value", "adj.P.Val")]
  
  colnames(adj.mat) <- c("pval.adj", "fdr.adj")
  colnames(noadj.mat) <- c("pval.no", "fdr.no")
  
  both.mat <- merge(adj.mat, noadj.mat, by = "row.names")
  both.mat$pval.adj <- -log10(both.mat$pval.adj)
  both.mat$fdr.adj <- -log10(both.mat$fdr.adj)
  both.mat$pval.no <- -log10(both.mat$pval.no)
  both.mat$fdr.no <- -log10(both.mat$fdr.no)
  both.mat$label <- invert_named_vector(dataSet$enrich_ids)[as.character(rownames(both.mat))];  

  # make plot
  if( "F" %in% colnames(rest)){
    fstat <- rest[, "F"];
  }else{
    fstat <- rest[, "t"];
  }  

  p.value <- rest[,"P.Value"];
  ord.inx <- order(rest[,"P.Value"], decreasing = FALSE);
  rest <- rest[ord.inx,,drop=F];
  colnames(rest)[1] <- "coefficient"; 
  rest$ids <- rownames(rest);

  names(fstat) <- names(p.value) <- rownames(dataSet$data.proc);
  fdr.p <- rest[,"adj.P.Val"];
  inx.imp <- p.value <= thresh;
  inx.imp <- ifelse(is.na(inx.imp), FALSE, inx.imp);
  sig.num <- length(which(inx.imp == TRUE))
  
  if(sig.num > 0){ 
    sig.p <- p.value[inx.imp];
    sig.mat <- rest[inx.imp,];
    sig.mat[,-ncol(sig.mat)] <- sapply(sig.mat[,-ncol(sig.mat)], function(x) signif(x, 5));
    rownames(sig.mat) <- make.names(rownames(rest)[inx.imp])
    # order the result simultaneously
  }else{
    return(c(0, 0));
  }

  AddMsg(paste(c("A total of", length(which(inx.imp == TRUE)), "significant features were found."), collapse=" "));
  rownames(both.mat) = both.mat[,1]
  both.mat <- both.mat[rownames(rest),]

  rest$label <- invert_named_vector(dataSet$enrich_ids)[as.character(rest$ids)];
  dataSet$comp.res <- rest;
  sig.mat$label <-  invert_named_vector(dataSet$enrich_ids)[as.character(sig.mat$ids)];
  rownames(sig.mat) <- sig.mat$ids;
  dataSet$sig.mat <- sig.mat

  if(sig.num> 0){
    res <- 1;
    fileName <- paste0("covariate_result_",sub("^(.*)[.].*", "\\1", dataSet$name), ".csv")
    fast.write.csv(rest,file=fileName);
    cov<-list (
      sig.num = sig.num,
      sig.nm = fileName,
      raw.thresh = thresh,
      thresh = -log10(thresh), # only used for plot threshold line
      p.value = p.value,
      p.value.no = both.mat$pval.no,
      p.log = -log10(p.value),
      inx.imp = inx.imp,
      sig.mat = sig.mat
    );
  }else{
    res <- 0;
    cov<-list (
      sig.num = sig.num,
      raw.thresh = thresh,
      thresh = -log10(thresh), # only used for plot threshold line
      p.value = p.value,
      p.value.no = both.mat$pval.no,
      p.log = -log10(p.value),
      inx.imp = inx.imp
    );
  }
  
  dataSet$design <- design;
  dataSet$contrast.type <- analysis.type;
  dataSet$comp.res <- rest;
  dataSet$de.method <- "limma"
  dataSet$comp.type <- "default"
  dataSet$fit.obj <- fit;

  dataSet$pval <- thresh;
  dataSet$fc.val <- 1;
  dataSet$analysis.var <- analysis.var;
  dataSet$de.adj <- covariates.vec;

  # for detail table
  dataSet$analSet$cov <- cov; 
  # for plotting adjp vs p
  dataSet$analSet$cov.mat <- both.mat; 
  both.list <- apply(both.mat, 2, function(x){unname(as.list(x))})

  both.list$thresh <- thresh;
  jsonNm <- gsub(paste0(".", imgFormat), ".json", imgName);
  jsonObj <- rjson::toJSON(both.list);
  sink(jsonNm);
  cat(jsonObj);
  sink();
    
  nonSig <- nrow(dataSet$comp.res) - sig.num;

  comp_res_path <- paste0(names(dataSets)[i], "_data/", "comp_res.csv");
  fast.write.csv(rest, file=paste0(dataName, "_data/", "comp_res.csv"))
  RegisterData(dataSet)
  return(c(sig.num, nonSig));
}

# Define function to invert named vector
invert_named_vector <- function(input_named_vec) {
  # Get names and values of input named vector
  input_names <- names(input_named_vec)
  input_values <- unname(input_named_vec)
  
  # Invert the named vector
  output_named_vec <- setNames(input_names, input_values)
  
  # Return output named vector
  return(output_named_vec)
}


PlotCovariateMap <- function(dataName, theme="default", imgName="NA", format="png", dpi=72){
  dataSet <- readDataset(dataName);
  both.mat <- dataSet$cov.mat
  both.mat <- both.mat[order(-both.mat[,"pval.adj"]),]
  logp_val <- dataSet$cov$thresh
  load_ggplot();
  library(ggrepel);
  topFeature <- 5;
  if(nrow(both.mat) < topFeature){
    topFeature <- nrow(both.mat);
  }
  if(theme == "default"){
    p <- ggplot(both.mat, mapping = aes(x = pval.no, y = pval.adj, label = Row.names)) +
      geom_rect(mapping = aes(xmin = logp_val, xmax = Inf, 
                              ymin = logp_val, ymax = Inf),
                fill = "#6699CC") +
      geom_rect(mapping = aes(xmin = -Inf, xmax = logp_val, 
                              ymin = -Inf, ymax = logp_val),
                fill = "grey") +
      geom_rect(mapping = aes(xmin = logp_val, xmax = Inf, 
                              ymin = -Inf, ymax = logp_val),
                fill = "#E2808A") +
      geom_rect(mapping = aes(xmin = -Inf, xmax = logp_val, 
                              ymin = logp_val, ymax = Inf),
                fill = "#94C973") +
      guides(size="none") +
      #annotate("text", x = 0.8, y = 0, label = "Never significant", size = 3) +
      #annotate("text", x = 2, y = 0, label = "Significant without adjustment", size = 3) +
      #annotate("text", x = 0.4, y = 1.5, label = "Significant with adjustment", size = 3) +
      #annotate("text", x = 2.25, y = 1.5, label = "Always significant", size = 3) +
      geom_point(aes(size=pval.adj), alpha=0.5) +
      geom_abline(slope=1, intercept = 0, linetype="dashed", color = "red", size = 1) +
      xlab("-log10(P-value): no covariate adjustment") +
      ylab("-log10(P-value): adjusted") +
      geom_text_repel(data = both.mat[c(1:topFeature),], 
                  aes(x=pval.no,y=pval.adj,label=Row.names)) +
      theme_bw();
  }else{
    p <- ggplot(both.mat, mapping = aes(x = pval.no, y = pval.adj, label = Row.names)) +
      guides(size="none") +
      geom_point(aes(size=pval.adj), alpha=0.5) +
      geom_abline(slope=1, intercept = 0, linetype="dashed", color = "red", size = 1) +
      geom_vline(xintercept = logp_val) +
      geom_hline(yintercept = logp_val) +
      xlab("-log10(P-value): no covariate adjustment") +
      ylab("-log10(P-value): adjusted") +
      geom_text_repel(data = both.mat[c(1:topFeature),], 
                  aes(x=pval.no,y=pval.adj,label=Row.names))
  }
  
  dataSet$covAdj <- imgName;

  width <- 8;
  height <- 8.18;
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=width, height=height, type=format, bg="white");
  print(p)
  dev.off()
  
  return(RegisterData(dataSet));
}


# general message only print when running local
AddMsg <- function(msg){
  if(!exists("msg.vec")){
    msg.vec <<- "";
  }
  msg.vec <<- c(msg.vec, msg);
  if(!.on.public.web){
    print(msg);
  }
}

PlotMultiFacCmpdSummary <- function(dataName, imgName, name, id, meta, meta2 = NA, densityBool = F, version, format = "png", dpi = 72, plotType = "boxplot") {
  rdtSet <- .get.rdt.set()

  # Set the dataset based on dataName
  if (dataName == "varPart") {
    # For variance partitioning, load the normalized dataset
    dat <- rdtSet$dataSet$norm
  } else {
    # Load data for normal conditions (non-varPart)
    dataSet <- readDataset(dataName)
    dat <- dataSet$data.proc
  }

  if (.on.public.web) {
    load_ggplot()
  }

  w <- 7.5
  
  meta.info <- rdtSet$dataSet$meta.info
  
  # Select the primary metadata
  sel.cls <- meta.info[which(rownames(meta.info) %in% colnames(dat)), meta]
  
  # Fix potential issues with empty or misaligned data
  if (length(sel.cls) == 0) {
    stop("Meta information could not be found for the selected columns.")
  }

  cls.type <- unname(rdtSet$dataSet$meta.types[meta])

  xlab = meta
  h <- 6
  imgName <- paste(imgName, "dpi", dpi, ".", format, sep = "")
  
  inx <- which(rownames(dat) == id)
  
  if (length(inx) == 0) {
    stop("ID not found in the dataset.")
  }

  sel.cls2 <- meta.info[which(rownames(meta.info) %in% colnames(dat)), meta2]

  if (!is.na(meta2) && meta2 != "" && meta2 != "NA" && length(sel.cls2) != 0) {
    cls.type2 <- unname(rdtSet$dataSet$meta.types[meta2])
    
    # Prepare the data frame based on meta and meta2
    df.norm <- data.frame(value = as.vector(t(dat)[, inx]), 
                          meta_val = if (cls.type == "cont") as.numeric(as.character(sel.cls)) else sel.cls,
                          meta2_val = if (cls.type2 == "cont") as.numeric(as.character(sel.cls2)) else sel.cls2)
    
    # Check if any NA values are introduced
    if (any(is.na(df.norm))) {
      warning("NA values were introduced in the data. Check the conversion of metadata types.")
    }

    # Get color schema for both meta and meta2
    col <- unique(GetColorSchema(sel.cls))
    col2 <- unique(GetColorSchema(sel.cls2))

    Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, width = w, height = h, type = format, bg = "white")
    
    if (cls.type == "disc" && cls.type2 == "disc") {
      # Both meta and meta2 are discrete
      p <- ggplot2::ggplot(df.norm, aes(x = meta_val, y = value, fill = meta2_val));
       # Add either boxplot or violin plot based on plotType argument
      if (plotType == "boxplot") {
        p <- p + geom_boxplot(outlier.shape = NA, outlier.colour = NA)
      } else if (plotType == "violin") {
        p <- p + geom_violin(trim = FALSE)
      }

       p <- p+theme_bw() + geom_jitter(size = 1) +
        scale_fill_manual(values = col2) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(fill = meta2) +
        ggtitle(name) +
        theme(plot.title = element_text(size = 11, hjust = 0.5, face = "bold")) +
        ylab("Abundance") + xlab(paste(meta, "vs", meta2)) +
        theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
        theme(plot.margin = margin(t = 0.15, r = 0.25, b = 0.15, l = 0.25, "cm"), axis.text = element_text(size = 10))

    } else if (cls.type == "cont" && cls.type2 == "cont") {
      # Both meta and meta2 are continuous
      p <- ggplot2::ggplot(df.norm, aes(x = meta_val, y = value, color = meta2_val)) +
        geom_point(size = 2) + theme_bw() + geom_smooth(method = lm, se = TRUE) +
        scale_color_gradient(low = "lightblue", high = "darkblue") +
        labs(color = meta2) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(name) + theme(plot.title = element_text(size = 11, hjust = 0.5, face = "bold")) +
        ylab("Abundance") + xlab(paste(meta, "vs", meta2)) +
        theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
        theme(plot.margin = margin(t = 0.15, r = 0.25, b = 0.15, l = 0.25, "cm"), axis.text = element_text(size = 10))

    } else if (cls.type == "disc" && cls.type2 == "cont") {
        p <- ggplot2::ggplot(df.norm, aes(x = meta_val, y = value, fill = meta_val))  # Fill violin/boxplot with metadata colors

        # Add either boxplot or violin plot based on plotType argument
        if (plotType == "boxplot") {
          p <- p + geom_boxplot(outlier.shape = NA, outlier.colour = NA)  # Boxplot with fill colors from metadata
        } else if (plotType == "violin") {
          p <- p + geom_violin(trim = FALSE)  # Violin plot with fill colors from metadata
        }

        p <- p + theme_bw() +
          geom_jitter(aes(color = meta2_val), size = 2) +  # Jitter points colored by continuous metadata, grey to black
          scale_color_gradient(low = "grey", high = "black") +  # Grey to black gradient for jitter points
          scale_fill_manual(values = col) +  # Fill the violin/boxplot with metadata colors
          labs(color = meta2, fill = meta) +  # Labels for color and fill
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          ggtitle(name) + theme(plot.title = element_text(size = 11, hjust = 0.5, face = "bold")) +
          ylab("Abundance") + xlab(meta) +
          theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
          theme(plot.margin = margin(t = 0.15, r = 0.25, b = 0.15, l = 0.25, "cm"), axis.text = element_text(size = 10))



    } else if (cls.type == "cont" && cls.type2 == "disc") {
      # meta is continuous, meta2 is discrete -> Line plot colored by meta2
      p <- ggplot2::ggplot(df.norm, aes(x = meta_val, y = value)) +
        geom_point(aes(color = meta2_val), size = 2) + theme_bw() +
        geom_smooth(method = lm, se = TRUE, color = "black") +
        scale_color_manual(values = col2) +
        labs(color = meta2) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(name) + theme(plot.title = element_text(size = 11, hjust = 0.5, face = "bold")) +
        ylab("Abundance") + xlab(meta) +
        theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
        theme(plot.margin = margin(t = 0.15, r = 0.25, b = 0.15, l = 0.25, "cm"), axis.text = element_text(size = 10))
    }

  } else {
    # Original approach (when meta2 is not provided)
    df.norm <- data.frame(value = as.vector(t(dat)[, inx]), 
                          name = if (cls.type == "cont") as.numeric(as.character(sel.cls)) else sel.cls)
    
    col <- unique(GetColorSchema(sel.cls))

    Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, width = w, height = h, type = format, bg = "white")

    if (cls.type == "disc") {
      p <- ggplot2::ggplot(df.norm, aes(x = name, y = value, fill = name));
      # Add either boxplot or violin plot based on plotType argument
      if (plotType == "boxplot") {
        p <- p + geom_boxplot(outlier.shape = NA, outlier.colour = NA)
      } else if (plotType == "violin") {
        p <- p + geom_violin(trim = FALSE)
      } 
      p <- p + theme_bw() + geom_jitter(size = 1) +
        scale_fill_manual(values = col) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(name) +
        theme(plot.title = element_text(size = 11, hjust = 0.5, face = "bold")) +
        ylab("Abundance") + xlab(meta) +
        theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
        theme(plot.margin = margin(t = 0.15, r = 0.25, b = 0.15, l = 0.25, "cm"), axis.text = element_text(size = 10))

    } else {
      p <- ggplot2::ggplot(df.norm, aes(x = name, y = value)) +
        geom_point(size = 2) + theme_bw() + geom_smooth(method = lm, se = TRUE, linewidth = 1) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(size = "none") +
        ggtitle(name) + theme(plot.title = element_text(size = 11, hjust = 0.5, face = "bold")) +
        ylab("Abundance") + xlab(meta) +
        theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
        theme(plot.margin = margin(t = 0.15, r = 0.25, b = 0.15, l = 0.25, "cm"), axis.text = element_text(size = 10))
    }
  }

  print(p)
  dev.off()

  if (.on.public.web) {
    return(imgName)
  } else {
    return(.set.rdt.set(rdtSet))
  }
}
