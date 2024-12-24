#'Prepare data for ROC analysis
#'@description Prepare data for ROC analysis
#'@param rdtSet Input the name of the created rdtSet (see InitDataObjects)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PrepareROCData <- function(sel.meta="NA",factor1,factor2){
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
    
    # Sort features by variance and select the top 1000 (if applicable)
    if (nrow(dataSet$data.proc) > 1000) {
      top_features <- order(feature_variances, decreasing = TRUE)[1:1000]
      dataSet$data.proc <- dataSet$data.proc[top_features, , drop = FALSE]
    }
    
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

   meta.info = rdtSet$dataSet$meta.info
  if(factor2!="NA" &factor2!="all" ){
  meta.info = rdtSet$dataSet$meta.info
    sample_include = rownames(meta.info[which(meta.info[[sel.meta]] %in% c(factor1,factor2)),])
 merged_data <- merged_data[,sample_include]
meta.info <- meta.info[rownames(meta.info) %in% sample_include,,drop=F]
   meta.info[[sel.meta]] <- droplevels(meta.info[[sel.meta]])
  }
if(factor2=="all"){
meta.info[[sel.meta]] <- as.character(meta.info[[sel.meta]])
idx = which(meta.info[[sel.meta]]!=factor1)
meta.info[[sel.meta]][idx] <- "Others"
meta.info[[sel.meta]] <- factor(meta.info[[sel.meta]],levels=c(factor1,"Others"))
}
   
  # Check if there are new samples to update `norm`
  new.inx <- is.na(rdtSet$dataSet$cls.all) | rdtSet$dataSet$cls.all == "";
  if(sum(new.inx) > 0){
    rdtSet$dataSet$new.samples <- TRUE;
    rdtSet$dataSet$new.data <- rdtSet$dataSet$norm.all[new.inx, , drop = F];
    rdtSet$dataSet$roc.norm <- merged_data;
    rdtSet$dataSet$roc.cls.orig <- factor(rdtSet$dataSet$meta.info[,sel.meta])
    rdtSet$dataSet$roc.cls <- meta.info[[sel.meta]]
  }else{
    rdtSet$dataSet$new.samples <- FALSE;
    rdtSet$dataSet$new.data <- NULL;
    rdtSet$dataSet$roc.norm <- merged_data;
    rdtSet$dataSet$roc.cls.orig <- factor(rdtSet$dataSet$meta.info[,sel.meta])
    rdtSet$dataSet$roc.cls <- meta.info[[sel.meta]]
  }
  
  return(.set.rdt.set(rdtSet));
}



#'Calculates feature importance
#'@description Perform calculation of feature importance (AUC, p value, fold change)
#'@usage CalculateFeatureRanking(rdtSet=NA, clust.num=5)
#'@param rdtSet Input the name of the created rdtSet (see InitDataObjects)
#'@param clust.num Numeric, input the number of clusters for cluster-analysis
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
CalculateFeatureRanking <- function(clust.num=5){
  #save.image("calculate.RData");
  
  rdtSet <- .get.rdt.set()
  LRConverged <<- "FALSE"; 
  
  x <- t(rdtSet$dataSet$roc.norm);
  y <- rdtSet$dataSet$roc.cls;
 
  # Check if multiclass or binary
  if(length(levels(y)) > 2){
    # Multiclass case: calculate AUC using pROC::multiclass.roc
    aucs <- sapply(1:ncol(x), function(i) {
      # Compute multiclass AUC for each feature
      roc.obj <- suppressMessages(pROC::multiclass.roc(y, x[, i]))
      return(roc.obj$auc)
    })
    
    # Macro-average AUC: take the mean across all features
    auc <- sapply(aucs, mean)
    
  } else {
    # Binary case: calculate AUC directly
    auc <- caTools::colAUC(x, y, plotROC=F)[1,]
  }
 
  # Perform univariate tests
  if(.on.public.web & RequireFastUnivTests(rdtSet)){
    res <- PerformFastUnivTests(t(x), y);
    ttp <- res[,2];  # P-values
  }else{
    ttp <- GetROCTtestP(x, y);  # Assuming GetROCTtestP also needs adaptation for multiclass
  }
  
  # Fold change computation (pairwise or aggregated for multiclass)
  data <- x;
  class_means <- sapply(levels(y), function(cls) {
    colMeans(data[which(rdtSet$dataSet$roc.cls == cls), , drop=FALSE])
  })
  
  # Compute pairwise fold changes
  ratio <- apply(class_means, 1, function(means) {
    means / means[1]  # Compare all to first class
  })
  ratio[ratio < 0] <- 0;
  fc <- signif(log2(ratio), 5);
  fc[is.infinite(fc) & fc < 0] <- -99;
  fc[is.infinite(fc) & fc > 0] <- 99;
  
  # Calculate the average fold change across classes
  if(length(levels(y)) > 2){
    fc <- colMeans(fc)
  }
  
  rdtSet$dataSet$roc_cols <- 2;
  if(rdtSet$dataSet$roc_cols > 1){
    # Perform k-means clustering for feature similarities
    kms <- ComputeKmeanClusters(t(x), clust.num);
    feat.rank.mat <- data.frame(AUC=auc, Pval=ttp, Avg_FC=fc, clusters = kms);
    rownames(feat.rank.mat) <- colnames(x);
    
    ord.inx <- order(feat.rank.mat$AUC, decreasing=T);
    feat.rank.mat  <- data.matrix(feat.rank.mat[ord.inx, , drop=FALSE]);
  }else{
    feat.rank.mat <- data.matrix(data.frame(AUC=auc, Pval=ttp, Avg_FC=fc, clusters=1))
  }
  
  # Format and return the result
  feat.rank.mat <<- signif(feat.rank.mat, digits = 5);
  
  if(rdtSet$analSet$mode == "univ"){
    fast.write.csv(feat.rank.mat, file="metaboanalyst_roc_univ.csv");
  }
  return(.set.rdt.set(rdtSet));  
}


RequireFastUnivTests <- function(rdtSet){
  if(nrow(rdtSet$dataSet$roc.norm) < 1000){
    return(FALSE);
  }else{
    return(TRUE);
  }
}

PerformFastUnivTests <- function(data, cls, var.equal=TRUE){
  if(!exists("mem.univ")){
    require("memoise");
    mem.univ <<- memoise(.perform.fast.univ.tests);
  }
  return(mem.univ(data, cls, var.equal));
}

.perform.fast.univ.tests <- function(data, cls, var.equal=TRUE){
  
  print("Performing fast univariate tests ....");
  # Convert data to matrix if not already
  data <- as.matrix(data);
  
  if(length(levels(cls)) > 2){
    # Perform one-way ANOVA for multiclass
    res <- try(rowcolFt(data, cls, var.equal = var.equal));
  }else{
    # Perform t-test for binary classification
    res <- try(rowcoltt(data, cls, FALSE, 1L, FALSE));
  }  
  
  if(class(res) == "try-error") {
    res <- cbind(NA, NA);
  }else{
    # Keep statistic and p-value
    res <- res[, c("statistic", "p.value")];
  }
  
  return(res);
}


#'Get p-values for ROC
#'@description ROC p-vaues, used in higher function
#'@param data Input data
#'@param cls Input class labels
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)

GetROCTtestP <- function(data, cls){
  if(nrow(data) < 1000){
    inx1 <- which(cls==levels(cls)[1]);
    inx2 <- which(cls==levels(cls)[2]);
    p.value <- apply(as.matrix(data), 2, function(x) {
      tmp <- try(t.test(x[inx1], x[inx2], paired = F, var.equal = T));
      if(class(tmp) == "try-error") {
        return(NA);
      }else{
        return(tmp$p.value);
      }
    })
  }else{ # use fast version
    p.value <- try(genefilter::rowttests(t(as.matrix(data)), cls)$p.value);
    if(class(p.value) == "try-error") {
      p.value <- NA;
    }
  }
  return(p.value);
}

# return a vector contain the cluster index 
ComputeKmeanClusters <- function(data, clust.num){
  set.seed(28051968);
  # if feature number is too low, change clust.num
  if(ncol(data) <= clust.num){
    clust.num <- nrow(data)-1;
  }
  if(clust.num < 2){
    clust.num <- 1;
  }
  kmeans.res <- kmeans(data, clust.num, nstart=100);
  return(kmeans.res$cluster);
}

SetAnalysisMode <- function(rdtSet=NA, mode){
  rdtSet <- .get.rdt.set()
  rdtSet$analSet$mode <- mode;
  return(.set.rdt.set(rdtSet));
}

GetUnivRankedFeatureNames <- function(){
  rownames(feat.rank.mat);
}

GetFeatureRankingMat <- function(){
  feat.rank.mat;
}

Perform.UnivROC <- function(feat.nm, 
                            version, format="png", 
                            dpi=72, isAUC, isOpt, 
                            optMethod, isPartial, measure, cutoff){
  
  #save.image("univroc.RData");
  rdtSet <- .get.rdt.set();
  
  imgName <- rdtSet$dataSet$url.var.nms[feat.nm];
  imgName = paste("roc_univ_", imgName, "_", version, "_dpi", dpi, ".", format, sep="");
  
  data_ori_norm <- rdtSet$dataSet$roc.norm
  if(!is.null(rdtSet$dataSet$roc.norm.orig)){
    data_ori_norm <- rdtSet$dataSet$roc.norm.orig
  }
  
  x <- unname(unlist(data_ori_norm[feat.nm, ]));
  y <- rdtSet$dataSet$roc.cls;
  opt.thresh = NA;
  
  # Check if y has more than two levels (multiclass)
   if (length(unique(y)) > 2) {
    # Multiclass ROC
    roc.obj <- pROC::multiclass.roc(y, x)
    class_labels <- unique(y)  # Get the class labels
    
    # Handle image creation and plotting
    w <- h <- 6; 
    rocs <- roc.obj$rocs
    colors <- rainbow(length(rocs))  # Assign different colors for each curve
    Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    par(oma=c(0,0,1,0));
    par(mar=c(4,4,4,4)+.1);
    
    # Loop through each ROC and plot
    for (i in 1:length(rocs)) {
      # Plot each ROC curve
      pROC::plot.roc(rocs[[i]], col=colors[i], legacy.axes=TRUE, grid=TRUE,
                     xlab="False positive rate", ylab="True positive rate", main=paste("Multiclass ROC: ", feat.nm),
                     add = ifelse(i > 1, TRUE, FALSE))  # Add only for subsequent plots
      
      
      # If isOpt is TRUE, calculate and display the optimal threshold for each class
      if (isOpt) {
        opt.ps <- data.frame(pROC::coords(rocs[[i]], "best", best.method=optMethod, transpose = TRUE))
        opt.thresh <- opt.ps["threshold",]
        points(opt.ps["specificity",], opt.ps["sensitivity",], pch=19, col=colors[i]);
        lbls <- paste(signif(opt.ps["threshold",], 3), " (", round(opt.ps["specificity",], 1), ", ", round(opt.ps["sensitivity",], 1), ")", sep="")
        text(opt.ps["specificity",], opt.ps["sensitivity",], adj=c(-.05, 1.25), label=lbls, col=colors[i])
      }
      
      # Optionally, you can add confidence intervals for each ROC curve (similar to binary case)
      #ci.obj <- pROC::ci.se(rocs[[i]], specificities=seq(0, 1, 0.05), boot.n=200, progress="none")
      #ROCR::plot(ci.obj, type="shape", col=adjustcolor(colors[i], alpha.f = 0.2))
    }
    
    # Plot the macro-average AUC if isAUC is TRUE
    if (isAUC) {
      auc.lbl <- paste("Multiclass AUC: ", round(roc.obj$auc,3), sep="")
      text(0.5, 0.5, auc.lbl, adj=c(0,1), col="navy")
      auc.ci <- pROC::ci.auc(roc.obj, method="bootstrap", boot.n=500, progress="none")
      auc.lbl <- paste("Multiclass AUC: ", round(auc.ci[2],3), sep="")
      ci.lbl <- paste("(", round(auc.ci[1],3), "-", round(auc.ci[3],3), ")", sep="")
      text(0.5, 0.4, paste(auc.lbl, "\n", ci.lbl, sep=""), adj=c(0,1), col="blue")
    }
    
    
    # Add a legend for class labels
    legend("bottomright", legend=class_labels, col=colors, lty=1, lwd=2, title="Classes")
    
    dev.off()
    
   }  else {
    # Binary ROC
    if(isPartial){
      if(measure == "se"){
        cutoff = cutoff;
      }else{
        cutoff = 1-cutoff;
      }
      roc.obj <- pROC::roc(y, x, partial.auc=c(1.0, cutoff), ci=TRUE, partial.auc.focus=measure, boot.n=50, percent = F, progress="none");
    }else{
      roc.obj <- pROC::roc(y, x, percent = F);
    }
    
    w <- h <- 6; 
    
    if(length(rdtSet$imgSet$roc.univ.name)==0){
      rdtSet$imgSet$roc.univ.name <- feat.nm;
      rdtSet$imgSet$roc.univ.plot <- imgName;
    } else {
      if(feat.nm %in% rdtSet$imgSet$roc.univ.name){
        idx <- which(feat.nm %in% rdtSet$imgSet$roc.univ.name);
        rdtSet$imgSet$roc.univ.name[idx] <- feat.nm;
        rdtSet$imgSet$roc.univ.plot[idx] <- imgName;
      } else {
        rdtSet$imgSet$roc.univ.name <- c(rdtSet$imgSet$roc.univ.name, feat.nm);
        rdtSet$imgSet$roc.univ.plot <- c(rdtSet$imgSet$roc.univ.plot, imgName);
      }   
    }
    
    Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    par(oma=c(0,0,1,0));
    par(mar=c(4,4,4,4)+.1);
    
    opt.thresh = NA;
    
    # Plot binary ROC curve
    if(isAUC){
      pROC::plot.roc(roc.obj, print.auc=F, legacy.axes=TRUE, col="navy", grid=T,
                     xlab = "False positive rate", ylab="True positive rate", main=feat.nm);
      ci.obj <- pROC::ci.se(roc.obj, specificities=seq(0, 1, 0.05), boot.n=200, progress="none");
      ROCR::plot(ci.obj,type="shape",col="#0000ff22");
    }else{
      pROC::plot.roc(roc.obj, print.auc=F, legacy.axes=TRUE, col="navy", grid=T,
                     xlab = "False positive rate", ylab="True positive rate",
                     auc.polygon=TRUE, auc.polygon.col="#0000ff22", main=feat.nm);
    }
    
    auc.ci <- pROC::ci.auc(roc.obj, method="bootstrap", boot.n=500, progress="none");
    roc.obj$ci <- auc.ci;
    auc.lbl <- paste("AUC: ", round(roc.obj$ci[2],3), sep="");
    ci.lbl <- paste("(", round(roc.obj$ci[1],3), "-", round(roc.obj$ci[3],3), ")", sep="");
    text(0.5, 0.5, paste(auc.lbl, "\n", ci.lbl, sep=""), adj=c(0,1), col="navy");
    
    if(isOpt){
      par(xpd=T);
      opt.ps <- data.frame(pROC::coords(roc.obj, "best", best.method=optMethod, transpose = TRUE));
      opt.thresh <- opt.ps["threshold",]
      points(opt.ps["specificity",], opt.ps["sensitivity",], pch=19, col="red");
      lbls=paste(signif(opt.ps["threshold",],3), "(", round(opt.ps["specificity",],1), ", ", round(opt.ps["sensitivity",],1), ")", sep="");
      text(opt.ps["specificity",], opt.ps["sensitivity",], adj=c(-.05,1.25), label=lbls);
    }
    
    dev.off();
  }
  
  rdtSet$analSet$opt.thresh <- opt.thresh
  
  if(.on.public.web){
    .set.rdt.set(rdtSet);
    return(imgName);
  }else{
    return(.set.rdt.set(rdtSet));
  }
}

PlotRocUnivBoxPlot <- function(feat.nm, version, format="png", dpi=72, isOpt, isQuery){
  
  rdtSet <- .get.rdt.set();
  
  if(.on.public.web){
    load_ggplot()
  }
  
  data_ori_norm <- rdtSet$dataSet$roc.norm
  if(!is.null(rdtSet$dataSet$roc.norm.orig)){
    data_ori_norm <- rdtSet$dataSet$roc.norm.orig
  }
  
  imgName <- rdtSet$dataSet$url.var.nms[feat.nm];
  imgName = paste("roc_boxplot_", imgName, "_", version, "_dpi", dpi, ".", format, sep="");
  
  x <- unname(unlist(data_ori_norm[feat.nm,]));
  y <- rdtSet$dataSet$roc.cls;
  scale <- dpi/72;
  w <- 200*scale;
  h <- 400*scale; 
  
  # Handle multiclass color schema
  col <- GetColorSchema(y)  # Ensure this supports multiclass
  
  if(length(rdtSet$imgSet$roc.univ.boxplot)==0){
    rdtSet$imgSet$roc.univ.boxplot <- imgName;
    rdtSet$imgSet$roc.univ.name2 <- feat.nm;
  } else {
    if(feat.nm %in% rdtSet$imgSet$roc.univ.name2){
      idx <- which(feat.nm %in% rdtSet$imgSet$roc.univ.name);
      rdtSet$imgSet$roc.univ.boxplot[idx] <- imgName;
    } else {
      rdtSet$imgSet$roc.univ.name2 <- c(rdtSet$imgSet$roc.univ.name2, feat.nm);
      rdtSet$imgSet$roc.univ.boxplot <- c(rdtSet$imgSet$roc.univ.boxplot, imgName);
    }   
  }
  
  Cairo::Cairo(file=imgName, width=w, height=h, type=format, bg="white", dpi=dpi);
  
  df <- data.frame(conc = x, class = y)
  
  # Generate boxplot using ggplot2
  p <- ggplot2::ggplot(df, aes(x=class, y=conc, fill=class)) + 
    geom_boxplot(notch=FALSE, outlier.shape = NA, outlier.colour=NA) + 
    theme_bw() + geom_jitter(size=1)
  
  p <- p + theme(axis.title.x = element_blank(), 
                 axis.title.y = element_blank(), 
                 legend.position = "none")
  
  p <- p + stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE)
  p <- p + theme(text = element_text(size=15), 
                 plot.margin = margin(t=0.45, r=0.25, b=1.5, l=0.25, "cm"), 
                 axis.text = element_text(size=10))
  
  # Apply color schema
  p <- p + scale_fill_manual(values=col)
  
  # Plot optimal threshold line only if `y` is binary and `isOpt` is TRUE
  if(isOpt && length(unique(y)) == 2){
    opt.thresh <- rdtSet$analSet$opt.thresh
    p <- p + geom_hline(aes(yintercept=opt.thresh), colour="red")
  }
  
  # Query threshold line if relevant (binary case)
  if(isQuery && length(unique(y)) == 2){
    thresh <- as.numeric(rdtSet$analSet$roc.obj$thresh)
    p <- p + geom_hline(aes(yintercept=thresh), colour="red")
  }
  
  # Print and save plot
  print(p)
  dev.off()
  
  # Return the image name or update the set
  if(.on.public.web){
    .set.rdt.set(rdtSet);
    return(imgName);
  }else{
    return(.set.rdt.set(rdtSet));
  }
}


GetROC.coords <- function(fld.nm, val, plot=TRUE, imgNm, classLabel=NULL){
  #save.image("coords.RData");
  rdtSet <- .get.rdt.set();
  roc.obj <- rdtSet$analSet$roc.obj
  
  if (length(unique(rdtSet$dataSet$roc.cls)) > 2) {
    # Multiclass case
    if (is.null(classLabel) || is.na(classLabel)) {
      stop("Please provide a class label to get ROC coordinates for a specific class.")
    }
    
    # Find the index of the class label
    class_labels <- unique(rdtSet$dataSet$roc.cls)
    class_index <- which(class_labels == classLabel)
    
    if (length(class_index) == 0) {
      stop(paste("Class label", classLabel, "not found in the data."))
    }
    
    # Get the ROC curve for the specified class (one-vs-rest)
    roc_class <- roc.obj$rocs[[class_index]]
    res <- pROC::coords(roc_class, val, input=fld.nm, transpose=TRUE)
    
    # Specificity and sensitivity
    sp <- res[2]
    se <- res[3]
    res <- round(res, 3)
    
    # If querying threshold, calculate confidence intervals
    if (fld.nm == "threshold") {
      ci.s <- pROC::ci.thresholds(roc_class, boot.n=100, thresholds=val, progress="none")
      specs <- round(ci.s$specificity, 3)
      sens <- round(ci.s$sensitivity, 3)
      res[2] <- paste(res[2], "(", specs[1], "-", specs[3], ")", sep="")
      res[3] <- paste(res[3], "(", sens[1], "-", sens[3], ")", sep="")
    }
    
    # Plot if requested
    if (plot) {
      PlotDetailROC(rdtSet, imgNm, res[1], sp, se)
    }
    
    # Return the result for the specified class
    return(res)
    
  } else {
    # Binary case (existing logic)
    res <- pROC::coords(roc.obj, val, input=fld.nm, transpose=TRUE)
    
    sp <- res[2]
    se <- res[3]
    res <- round(res, 3)
    
    rdtSet$analSet$thresh.obj <- NULL
    if (fld.nm == "threshold") {
      ci.s <- pROC::ci.thresholds(roc.obj, boot.n=100, thresholds=val, progress="none")
      specs <- round(ci.s$specificity, 3)
      sens <- round(ci.s$sensitivity, 3)
      res[2] <- paste(res[2], "(", specs[1], "-", specs[3], ")", sep="")
      res[3] <- paste(res[3], "(", sens[1], "-", sens[3], ")", sep="")
    }
    
    # Plot if requested
    if (plot) {
      PlotDetailROC(rdtSet, imgNm, res[1], sp, se)
    }
    
    return(res)
  }
}



PlotDetailROC <- function(imgName, thresh, sp, se, dpi=72, format="png"){
  
  rdtSet <- .get.rdt.set();
  
  imgName = paste(imgName, "_dpi", dpi, ".", format, sep="");
  
  roc.obj <- rdtSet$analSet$roc.obj;
  
  w <- h <- 6;
  rdtSet$imgSet$roc.univ <- imgName;
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  par(oma = c(0,0,1,0));
  par(mar=c(4,4,4,4) + .1);
  
  pROC::plot.roc(roc.obj, print.auc=F, legacy.axes=TRUE, col="navy", grid=T,
                 xlab = "False positive rate", ylab="True positive rate",
                 auc.polygon=TRUE, auc.polygon.col="#0000ff22", main=current.feat.nm);
  
  points(sp, se, cex=1.8, col="red");
  
  dev.off();
  return(.set.rdt.set(rdtSet));
  
}

PrepareROCDetails <- function(feat.nm){
  
  rdtSet <- .get.rdt.set();
  
  x <- unname(unlist(rdtSet$dataSet$roc.norm[feat.nm, ]));
  y <- rdtSet$dataSet$roc.cls;
  
  if (length(unique(y)) > 2) {
    # Multiclass case: Perform one-vs-rest ROC for each class
    classes <- unique(y)
    
    # List to store the ROC details for each class
    roc.details <- list()
    i <- 1;
    for (class in classes) {
      # Create a binary response for this class (1 for the current class, 0 for others)
      y_bin <- ifelse(y == class, 1, 0)
      
      # Check if there are both control (0) and case (1) observations
      if (sum(y_bin == 0) > 0 && sum(y_bin == 1) > 0) {
        # Compute ROC for this class (one-vs-rest)
        roc.res <- pROC::roc(y_bin, x, ci = TRUE, of = "auc")
        
        # Create a matrix for each class
        roc.mat <- as.matrix(data.frame(
          "Class" = i,
          "Cut.Offs" = roc.res$thresholds,
          "Sensitivity" = roc.res$sensitivities,
          "Specificity" = roc.res$specificities,
          "Sens.Spec." = roc.res$sensitivities + roc.res$specificities,
          "LRP" = roc.res$sensitivities / (1 - roc.res$specificities),
          "LRN" = (1 - roc.res$sensitivities) / roc.res$specificities
        ))
        
        # Store the result for this class
        roc.details[[class]] <- roc.mat
      } else {
        # If no control observations, skip this class
        message(paste("Skipping class", class, "due to lack of control observations."))
      }
      i <- i +1;
    }
    
    # Combine the ROC data for all classes that had valid control observations
    roc.mat.combined <- do.call(rbind, roc.details)
    
  } else {
    # Binary case: Perform standard ROC
    roc.res <- pROC::roc(y, x, ci = TRUE, of = "auc")
    
    # Create a matrix for binary case
    roc.mat.combined <- as.matrix(data.frame(
      "Class" = 1,
      "Cut.Offs" = roc.res$thresholds,
      "Sensitivity" = roc.res$sensitivities,
      "Specificity" = roc.res$specificities,
      "Sens.Spec." = roc.res$sensitivities + roc.res$specificities,
      "LRP" = roc.res$sensitivities / (1 - roc.res$specificities),
      "LRN" = (1 - roc.res$sensitivities) / roc.res$specificities
    ))
  }
  
  # Clean up any NA/Inf/-Inf values
  roc.mat.combined[!is.finite(roc.mat.combined)] <- NA
  
  # Write the ROC details to a CSV file
  filename <- paste(rdtSet$dataSet$url.var.nms[feat.nm], "_roc.csv", sep="")
  fast.write.csv(signif(roc.mat.combined, 4), file=filename, row.names=FALSE)
  
  # Save ROC data in the rdtSet for future use
  rdtSet$analSet$roc.mat <- signif(roc.mat.combined, 6)
  rdtSet$analSet$roc.obj <- roc.res
  
  current.feat.nm <<- feat.nm
  
  return(.set.rdt.set(rdtSet))
}

PerformCV.explore <- function(cls.method, rank.method = "auroc", lvNum = 2, propTraining = 2/3, target_group = NULL) {
rdtSet <- .get.rdt.set();
  
  rdtSet$analSet$exp.method <- cls.method
  rdtSet$analSet$rank.method <- rank.method
  rdtSet$analSet$exp.lvNum <- lvNum

  data <- t(rdtSet$dataSet$roc.norm)
  cls <- rdtSet$dataSet$roc.cls
  nClasses <- length(unique(cls))
  isBinary <- nClasses == 2
 
target_group <- levels(cls)[1]
other_group <- levels(cls)[2]
  # Convert to one-against-all if target_group is specified and data is multiclass
 
  # Number of subsampling runs to produce smooth curve
  nRuns <- if (nrow(data) > 500) 10 else if (nrow(data) > 200) 20 else if (nrow(data) > 100) 30 else 50

  nFeatures <- GetFeatureNumbers(ncol(data))
  
  feat.outp <- actualCls <- vector(length = nRuns, mode = "list")
  perf.outp <- vector(length = length(nFeatures), mode = "list")
  perf.outp <- lapply(perf.outp, function(x) { x <- vector(length = nRuns, mode = "list"); return(x) })
  auc.mat <- accu.mat <- matrix(nrow = nRuns, ncol = length(nFeatures))
  
  splitMat <- GetTrainTestSplitMat(cls, propTraining, nRuns)
  trainRuns <- splitMat$training.mat
  testRuns <- splitMat$testing.mat
  
  for (irun in 1:nRuns) {
 
    trainingSampleRun <- trainRuns[irun, ]
    testSampleRun <- testRuns[irun, ]
    x.in <- data[trainingSampleRun, ]
    y.train <- cls[trainingSampleRun]
    actualCls[[irun]] <- y.test <- cls[testSampleRun]
 
    # Check for missing values in y.train
    if (any(is.na(y.train))) {
      stop("y.train contains missing values. Please ensure labels are correctly assigned before training.")
    }
    
    # Feature ranking
    imp.vec <- RankFeatures(x.in, y.train, rank.method, lvNum)
    feat.outp[[irun]] <- imp.vec
    ord.inx <- order(imp.vec, decreasing = TRUE)
    ordData <- data[, ord.inx]
    
    # Build classifiers for each subset of features and test on the test data
    for (inum in seq(along = nFeatures)) {
      x.train <- ordData[trainingSampleRun, 1:nFeatures[inum]]
      x.test <- ordData[testSampleRun, 1:nFeatures[inum]]
      prob.out <- Predict.class(x.train, y.train, x.test, cls.method, lvNum)
      
      # Calculate AUC for binary classification
      pred <- ROCR::prediction(prob.out, y.test)
      auc.mat[irun, inum] <- slot(ROCR::performance(pred, "auc"), "y.values")[[1]]
      
      perf.outp[[inum]][[irun]] <- prob.out
         if (!all(levels(y.test) ==levels(factor(pred@labels[[1]]))) ) {
      auc.mat[irun, inum] <- 1-      auc.mat[irun, inum]
  perf.outp[[inum]][[irun]] <- 1-prob.out;
}
     
      pred.out <- as.factor(ifelse(prob.out > 0.5, target_group, other_group))
      accu.mat[irun, inum] <- Get.Accuracy(table(pred.out, y.test))
    }
  }
 
  # Prepare results for plotting
  preds <- vector(length = length(nFeatures), mode = "list")
  act.vec <- unlist(actualCls)
  for (m in 1:length(nFeatures)) {
    prob.vec <- unlist(perf.outp[[m]])
    pred <- ROCR::prediction(prob.vec, act.vec)
    preds[[m]] <- pred
  }

  # Process and save results
  colnames(accu.mat) <- colnames(auc.mat) <- names(preds) <- paste(nFeatures, sep = "")
  auc.vec <- colMeans(auc.mat)
  auc.cis <- GetCIs(auc.mat)
  best.model.inx <- which.max(auc.vec)
  
  rdtSet$analSet$multiROC <- list(
    type = rdtSet$analSet$type,
    train.mat = trainRuns,
    test.feats = nFeatures,
    pred.cv = perf.outp,
    true.cv = actualCls,
    imp.cv = feat.outp,
    pred.list = preds,
    accu.mat = accu.mat,
    auc.vec = auc.vec,
    auc.ci = auc.cis,
    best.model.inx = best.model.inx,
    target_group = target_group,
   other_group=other_group
  )
  
  return(.set.rdt.set(rdtSet))
}

PlotImpBiomarkers <- function(imgName, format = "png", dpi = 72, mdl.inx, measure = "freq", feat.num = 15) {
  rdtSet <- .get.rdt.set()
  imgName <- paste(imgName, "dpi", dpi, ".", format, sep = "")
  w <- 8; h <- 8
  rdtSet$imgSet$roc.imp.plot <- imgName
  rdtSet$imgSet$roc.imp.name <- mdl.inx
  
  Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, width = w, height = h, type = format, bg = "white")
  
  # Use norm_explore if available, otherwise default to norm data
  data <- if (is.null(rdtSet$dataSet$norm_explore)) t(rdtSet$dataSet$roc.norm) else t(rdtSet$dataSet$norm_explore)
  cls <- rdtSet$dataSet$roc.cls
  
  target_group <- rdtSet$analSet$multiROC$target_group
  other_group <- rdtSet$analSet$multiROC$other_group
  # Adjust for one-against-all if necessary
 
  if (nrow(data) != length(cls)) {
    stop("Mismatch between number of samples in data and cls. Ensure they have the same length.")
  }
  
  if (mdl.inx == -1) {
    mdl.inx <- rdtSet$analSet$multiROC$best.model.inx
  }
  
  imp.mat <- GetImpFeatureMat(rdtSet, rdtSet$analSet$multiROC$imp.cv, rdtSet$analSet$multiROC$test.feats[mdl.inx])
  imp.nms <- rownames(imp.mat)
  hit.nms <- imp.nms[imp.nms %in% colnames(data)]
  data <- data[, hit.nms, drop = FALSE]
  
  # Check if `data` has columns after filtering
  if (ncol(data) == 0) {
    stop("No matching features in data after filtering with importance matrix names. Check `imp.mat` and column names of `data`.")
  }
  
  # Compute median for each group
  mds <- apply(data, 2, function(x) {
    tapply(x, cls, median)
  })
  
  # Assign High or Low based on median values
  lowhigh <- t(apply(mds, 2, function(x) {
    ifelse(rank(x) == 1, "Low", "High")
  }))
  
  temp.dat <- data.frame(imp.mat, lowhigh)
  colnames(temp.dat) <- c(colnames(imp.mat), levels(cls))
  imp.fileNm <- "imp_features_cv.csv"
  fast.write.csv(temp.dat, file = imp.fileNm)
  
  rdtSet$analSet$imp.mat <- imp.mat
  rdtSet$analSet$lowhigh <- lowhigh
  rdtSet$analSet$roc.sig.nm <- imp.fileNm
  
  if (measure == "freq") {
    imp.vec <- sort(imp.mat[, 1], decreasing = TRUE)
    xlbl <- "Selected Frequency (%)"
  } else {
    imp.vec <- sort(imp.mat[, 2], decreasing = TRUE)
    xlbl <- "Average Importance"
  }
  
  feat.num <- min(feat.num, length(imp.vec))
  imp.vec <- rev(imp.vec[1:feat.num])
  
  nms.orig <- names(imp.vec)
  vip.nms <- substr(nms.orig, 1, 20)
  names(imp.vec) <- NULL
  
  # Save original graphical parameters and adjust plot margins
  op <- par(mar = c(6, 10, 3, 7))
  
  xlim.ext <- GetExtendRange(imp.vec, 12)
  dotchart(imp.vec, bg = "blue", xlab = xlbl, xlim = xlim.ext)
  mtext(side = 2, at = 1:feat.num, vip.nms, las = 2, line = 1)
  names(imp.vec) <- nms.orig
  axis.lims <- par("usr")
  shift <- 2 * par("cxy")[1]
  lgd.x <- axis.lims[2] + shift
  x <- rep(lgd.x, feat.num)
  y <- 1:feat.num
  par(xpd = TRUE)
  
  # Synchronize lowhigh with imp.vec and set colors
  lowhigh <- lowhigh[nms.orig, , drop = FALSE]
  col <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(2)
  bg <- ifelse(lowhigh == "High", col[1], col[2])
  
  # Plot with updated labels and colors
  points(x, y, pch = 22, bg = bg, cex = 3)
  text(x[1], axis.lims[4], target_group, srt = 45, adj = c(0.2, 0.5))
  
  # Continuous gradient legend for "High" and "Low" as a gradient rectangle
  gradient_col <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(100)
  rect_x <- x[1] + shift
  rect_y_bottom <- axis.lims[3] + 0.2 * diff(axis.lims[3:4])
  rect_y_top <- axis.lims[3] + 0.8 * diff(axis.lims[3:4])
  
  # Draw multiple small rectangles to create the gradient effect
  gradient_steps <- length(gradient_col)
  rect_height <- (rect_y_top - rect_y_bottom) / gradient_steps
  for (i in 1:gradient_steps) {
    rect(
      xleft = rect_x,
      ybottom = rect_y_bottom + (i - 1) * rect_height,
      xright = rect_x + shift / 4,
      ytop = rect_y_bottom + i * rect_height,
      col = gradient_col[i],
      border = NA
    )
  }
  
  # Add labels for "High" and "Low"
  text(rect_x + shift / 2, rect_y_top, "High", pos = 3)
  text(rect_x + shift / 2, rect_y_bottom, "Low", pos = 1)
  
  # Reset graphical parameters
  par(op)
  dev.off()
  
  return(.set.rdt.set(rdtSet))
}

PlotAccuracy<-function(rdtSet=NA, imgName, format="png", dpi=72){
  
rdtSet <- .get.rdt.set();;
  anal.mode <- rdtSet$analSet$mode;
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 9; h <- 7;
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  if(is.null(rdtSet$analSet$multiROC$accu.mat)){
    accu.mat <- rdtSet$analSet$ROCtest$accu.mat;
  }else{
    accu.mat <- rdtSet$analSet$multiROC$accu.mat;
  }
  
  mn.accu <- apply(accu.mat, 2, mean);
  ylabl <- 'Predictive Accuracy';
  ylim <- c(0,1);
  title <- 'Predictive accuracies with different features';
  txt.lbls <- paste(100*round(mn.accu,3),'%');
  
  matplot(t(accu.mat),type='l', lty=2, col="grey", xlab='Number of features',ylab=ylabl, ylim=ylim,
          axes=F,main=title);
  
  lines(1:ncol(accu.mat), mn.accu, lwd=2);
  points(mn.accu, pch=19, col=ifelse(1:length(mn.accu)==which.max(mn.accu),"red","blue"));
  text(mn.accu,labels=txt.lbls, adj=c(-0.3, -0.5), srt=45, xpd=T)
  axis(2);
  
  lbls <- colnames(accu.mat);
  axis(1, 1:length(mn.accu), labels=lbls);
  
  rdtSet$imgSet$roc.pred <- imgName;
  
  dev.off();
  return(.set.rdt.set(rdtSet));
}

GetFeatureNumbers <- function(feat.len){
  if(feat.len > 100){
    nFeatures <- c(5, 10, 15, 25, 50, 100);
  }else if(feat.len > 50){
    nFeatures <- c(3, 5, 10, 20, round(feat.len/2), feat.len);
  }else if(feat.len > 20){
    nFeatures <- c(2, 3, 5, 10, 20, feat.len);
  }else if(feat.len > 10){
    nFeatures <- c(2, 3, 5, 7, 10, feat.len);
  }else if(feat.len > 1){
    nFeatures <- sort(sample(2:feat.len));
  }else{
    print("Feature number is less than 2!")
    return();
  }
  nFeatures;
}

MergeDuplicates <- function(data, dim=2){
  
  if(is.null(dim(data))){ # a vector
    if(is.null(names(data))){
      print("Cannot detect duplicate data without names!!!");
      return();
    }
    nm.cls <- as.factor(names(data));
    uniq.len <- length(levels(nm.cls));
    if(uniq.len == length(data)){
      return(data);
    }
    new.data <- vector (mode="numeric",length=uniq.len);
    for(i in 1:uniq.len){
      dup.inx <- nm.cls == levels(nm.cls)[i];
      new.data[i] <- mean(data[dup.inx]);
    }
    names(new.data) <- levels(nm.cls);
    rem.len <- length(data) - length(new.data);
  }else{
    if(dim == 1){
      data <- t(data);
    }
    if(is.null(colnames(data))){
      print("Cannot detect duplicate data without var names!!!");
      return();
    }
    
    nm.cls <- as.factor(colnames(data));
    uniq.len <- length(levels(nm.cls));
    
    if(uniq.len == ncol(data)){
      if(dim == 1){
        data <- t(data);
      }
      return(data);
    }
    
    new.data <- matrix (nrow=nrow(data), ncol=uniq.len);
    for(i in 1:uniq.len){
      dup.inx <- which(nm.cls == levels(nm.cls)[i]);
      new.data[,i] <- apply(data[,dup.inx, drop=F], 1, mean);
    }
    rownames(new.data) <- rownames(data);
    colnames(new.data) <- levels(nm.cls);
    
    rem.len <- ncol(data) - ncol(new.data);
    if(dim == 1){
      new.data <- t(new.data);
    }
  }
  print(paste(rem.len, "duplicates are merged to their average"));
  new.data;
}


GetImpFeatureMat <- function(rdtSet, feat.outp, bestFeatNum){
  
  anal.mode <- rdtSet$analSet$mode; 
  
  # first order each run by cmpd names so that can be combined to a data frame
  feat.outp <- lapply(feat.outp, function(x) x[order(names(x))]);
  
  ####################################################################
  # First rank by frequencies of being selected in the given model ###
  ####################################################################
  # obtain their ranks
  freqRank <- lapply(feat.outp, function(x) rank(-x));
  runRanksMat <- do.call("cbind", freqRank);
  
  # order by their median rank across runs
  ordRunRanksMat <- as.data.frame(runRanksMat[order(apply(runRanksMat, 1, median)),]);
  
  # Then rank by mean importance measures
  impsMat <- as.data.frame(do.call("cbind", feat.outp));
  impsVec <- apply(impsMat, 1, mean);
  
  # if(anal.mode == "explore"){
    # now count the number being selected in the bestFeatNum
    selectedMat <- apply(ordRunRanksMat, 2, function(x) x <= bestFeatNum);
  # }else{
  #   selectedMat <- ordRunRanksMat;
  # }
  
  # calculate percentage of being selected in the best subsets
  percentVec <- apply(selectedMat, 1, sum)/ncol(ordRunRanksMat);
  
  # remove ones never being selected
  percentVec <- percentVec[percentVec > 0];
  
  # reorder the imps to percentVec
  impsVec <- impsVec[names(percentVec)];
  
  ###################################
  # combine and return the result
  ####################################
  imp.mat <- cbind(percentVec, impsVec);
  ord.inx <- order(imp.mat[,1], imp.mat[,2], decreasing=T);
  imp.mat <- imp.mat[ord.inx,];
  colnames(imp.mat) <- c("Rank Freq.", "Importance");
  
  return(imp.mat);
}

GetModelNames <- function(){
  rdtSet <- .get.rdt.set();
  test.nums <- rdtSet$analSet$multiROC$test.feats;
  paste("Model", 1:length(test.nums), "(", test.nums, "features)");
}

GetGroups <- function(group){
  rdtSet <- .get.rdt.set();
  levels(rdtSet$dataSet$meta.info[[group]])
}

GetGrpSampleNames <- function(grp){
  rdtSet <- .get.rdt.set();
  as.character(rownames(rdtSet$dataSet$meta.info)[rdtSet$dataSet$meta.info==grp])
}

 
GetTrainTestSplitMat <- function(y, propTraining = 2/3, nRuns = 30){
  
  nTotalSample <- length(y);
  
  smallestClass <- names(sort(table(y)))[1];
  nSmallest <- sum(y == smallestClass);
  
  nSmallestTrain <- round(propTraining * nSmallest);
  nBiggestTrain <- nSmallestTrain;
  
  nSmallestTest <- nSmallest - nSmallestTrain;
  nBiggestTest <- nTotalSample - (nSmallestTest + nSmallestTrain + nBiggestTrain);
  
  # sanity check for very large number of samples
  # for each run max 600 - 400 train, 200 test 
  big.samples <- FALSE;
  if(nSmallestTrain > 400){
    big.samples <- TRUE;
    nSmallestTrain <- nBiggestTrain <- 400;
    nSmallestTest <- nBiggestTest <- 200;
  }
  
  # split up in smallest class indices and biggest class indices
  smallestIndices <- which(y == smallestClass)
  biggestIndices <- seq(along = y)[-smallestIndices]
  
  nTrainingSample <- nSmallestTrain + nBiggestTrain;
  nTestSample <- nSmallestTest + nBiggestTest;
  
  trainingSampleAllRuns <- matrix(0, nrow = nRuns, ncol = nTrainingSample)
  testSampleAllRuns  <- matrix(0, nrow = nRuns, ncol = nTestSample);
  
  for (irun in 1:nRuns) {
    sampleSmallestTrain <- sample(smallestIndices, nSmallestTrain);
    sampleBiggestTrain <- sample(biggestIndices, nBiggestTrain);
    trainingSampleRun <- c(sampleSmallestTrain, sampleBiggestTrain);
    indicesTrainingSample <- rep(FALSE, length = nTotalSample);
    indicesTrainingSample[trainingSampleRun] <- TRUE;
    trainingSampleAllRuns[irun, ] <- which(indicesTrainingSample);
    if(big.samples){
      testSampleAllRuns[irun, ] <- sample(which(!indicesTrainingSample), 200);
    }else{
      testSampleAllRuns[irun, ] <- which(!indicesTrainingSample);
    }
  }
  
  # return the results
  list(
    training.mat = trainingSampleAllRuns,
    testing.mat = testSampleAllRuns
  );
}

RankFeatures <- function(x.in, y.in, method, lvNum){
  if(method == "rf"){ # use randomforest mean decrease accuracy
    rf <- randomForest::randomForest(x = x.in,y = y.in,importance=TRUE, keep.forest=F);
    return(randomForest::importance(rf)[ ,"MeanDecreaseAccuracy"]);
  }else if (method == "pls"){
    ncls <- as.numeric(y.in)-1;
    datmat <- as.matrix(x.in);
    pls <- pls::plsr(ncls~datmat,method='oscorespls', ncomp=lvNum);
    return(Get.VIP(pls, lvNum));
  }else if(method == "svm"){
    svmres <- e1071::svm(x.in, y.in, type = 'C', kernel="linear");
    imp.vec <- (t(svmres$coefs) %*% svmres$SV)[1,]^2;
    names(imp.vec) <- colnames(x.in);
    return(imp.vec);
  }else if(method == "auroc"){ # univariate based ou area under ROC
    imp.vec <- caTools::colAUC(x.in, y.in, plotROC=F)[1,];
    return(imp.vec);
  }else if(method == "tt"){ # univariate based ou area under ROC
    imp.vec <- Get.Fstat(x.in, as.factor(y.in)); # all f-stats
    names(imp.vec) <- colnames(x.in);
    return(imp.vec);
  }else if(method == "fisher"){ # univariate based ou area under ROC
    imp.vec <- Get.Fisher(x.in, as.factor(y.in));
    names(imp.vec) <- colnames(x.in);
    return(imp.vec);
  }else{
    print("Not supported!");
    return(NULL);
  }
}


Predict.class <- function(x.train, y.train, x.test, clsMethod="pls", lvNum, imp.out=F){
  
  # first convert class label to 0/1 so convert prob-> class will be easier
  y.train <- as.factor(as.numeric(y.train)-1);
  
  # note, we only need prob for class 1, pred can be inferred
  if (clsMethod == "rf"){
    model <- randomForest::randomForest(x.train, y.train, ntree=300, importance=T);
    prob.out <- predict(model, x.test, type="prob")[,"1"];
    if(imp.out){
      imp.vec <- randomForest::importance(model)[ ,"MeanDecreaseAccuracy"];
      return(list(imp.vec = imp.vec, prob.out = prob.out));
    }
    return(prob.out);
  }else if(clsMethod == "pls"){ # plsda
    model <- caret::plsda(x.train, y.train, method='oscorespls', ncomp=ifelse(ncol(x.train)>lvNum, lvNum, 2));
    prob.out <- predict(model, x.test, type="prob")[,"1",1];
    if(imp.out){
      imp.vec <- Get.VIP(model, lvNum);
      return(list(imp.vec = imp.vec, prob.out = prob.out));
    }
    return(prob.out);
  }else if(clsMethod == "lr"){ # logistic regression with selected variables (only in Test)
    x <- x.train;
    y <- y.train;
    xx.test <- x.test;
    
    names.x.origin <- names(x);
    names(x) <- paste0("V", 1:(ncol(x)));
    names(xx.test) <- names(x);
    
    data.xy <- data.frame(y, x);
    model <- logisticReg(data.xy);
    prob.out <- predict(model, xx.test, type="response");
    if(imp.out){
      imp.vec <- names(model$coefficients)[-1]
      return(list(imp.vec = imp.vec, prob.out = prob.out));
    }
    return(prob.out);
  }else{ # svm
    model <- e1071::svm(x.train, y.train, type = 'C', kernel="linear", probability=TRUE);
    prob.out <- attr(predict(model, x.test,  probability=TRUE), "probabilities")[,"1"];
    if(imp.out){
      imp.vec <- (t(model$coefs) %*% model$SV)[1,]^2;
      names(imp.vec) <- colnames(x.train);
      return(list(imp.vec = imp.vec, prob.out = prob.out));
    }
    return(prob.out);
  }
}

logisticReg <- function(d.xy) {
  fmla <- as.formula(paste("y ~", paste(names(d.xy)[-1], collapse="+")));
  model <- glm(fmla, data=d.xy, family=binomial(link="logit"), na.action="na.omit")
  return(model);
}

genLREquation <- function(coef.mdl){
  coef.mdl <- round(coef.mdl, 3);
  
  eq <- coef.mdl[[1]];
  for(i in 2:length(coef.mdl)) {
    eq <- paste(eq, ifelse(sign(coef.mdl[[i]])==1,"+","-"), abs(round(coef.mdl[[i]],3)), names(coef.mdl)[i]);
  }
  return(eq);
}


Get.VIP <- function(pls.obj, comp=2){
  # only use the top two comps
  b <- c(pls.obj$Yloadings)[1:comp];
  T <- pls.obj$scores[,1:comp, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- pls.obj$loading.weights[,1:comp, drop = FALSE]
  Wnorm2 <- colSums(W^2);
  SSW <- sweep(W^2, 2, SS / Wnorm2, "*")
  vips <- sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS));
  if(is.null(dim(vips))){
    vip.mns<-vips;
  }else{
    vip.mns<-apply(vips, 2, mean);      
  }
  vip.mns;
}

Get.Accuracy <- function(cm) {
    sum(diag(cm)) / sum(cm);
}

GetCIs <- function(data, param=F){
  means <- colMeans(data, na.rm=T);
  if(param){
    sds <- apply(data, 2, sd, na.rm=T);
    nn <- nrow(data);
    std.err <- sds/sqrt(nn);
    LCL <- round(means-1.96*std.err,3);
    LCL[LCL<0] <- 0;
    UCL <- round(means+1.96*std.err, 3);
    UCL[UCL>1] <- 1;
    res <- paste(LCL, "-", UCL, sep="");
  }else{
    cis <- apply(data, 2, quantile, probs=c(0.025, 0.975));
    cis <- round(cis,3);
    cis[cis<0] <- 0;
    cis[cis>1] <- 1;
    res <- paste(cis[1,], "-", cis[2,], sep="");
  }
  res;
}

ComputeAverageCurve<-function(perf, avg.method){
  # now get the average curve
  perf.avg = perf;
  if(avg.method == "vertical"){
    x.values <- seq(min(unlist(perf@x.values)), max(unlist(perf@x.values)),
                    length=max( sapply(perf@x.values, length)))
    for (i in 1:length(perf@y.values)) {
      perf.avg@y.values[[i]] <-
        approxfun(perf@x.values[[i]], perf@y.values[[i]],
                  ties=mean, rule=2)(x.values)
    }
    perf.avg@y.values <- list(rowMeans(data.frame(perf.avg@y.values )))
    perf.avg@x.values <- list(x.values)
  }else if(avg.method == "horizontal"){
    y.values <- seq(min(unlist(perf@y.values)), max(unlist(perf@y.values)),
                    length=max(sapply(perf@y.values, length)))
    for (i in 1:length(perf@x.values)) {
      perf.avg@x.values[[i]] <- approxfun(perf@y.values[[i]],
                                          perf@x.values[[i]],
                                          ties=mean, rule=2)(y.values)
    }
    perf.avg@x.values <- list(rowMeans( data.frame( perf.avg@x.values )));
    perf.avg@y.values <- list(y.values);
  }else{ # threshold
    all.alphas <- unlist(perf@alpha.values);
    min.alpha <- min(all.alphas);
    if(min.alpha == -Inf){
      min.alpha <- 0;
    }
    max.alpha <- max(all.alphas);
    if(max.alpha == Inf){
      max.alpha <- 1.0;
    }
    
    alpha.values <- rev(seq(min.alpha, max.alpha,length=max(sapply(perf@alpha.values, length))));
    perf.sampled <- perf;
    for (i in 1:length(perf.sampled@y.values)) {
      perf.sampled@x.values[[i]] <-
        approxfun(perf@alpha.values[[i]],perf@x.values[[i]],
                  rule=2, ties=mean)(alpha.values)
      perf.sampled@y.values[[i]] <-
        approxfun(perf@alpha.values[[i]], perf@y.values[[i]],
                  rule=2, ties=mean)(alpha.values)
    }
    ## compute average curve
    perf.avg <- perf.sampled
    perf.avg@x.values <- list(rowMeans(data.frame(perf.avg@x.values)))
    perf.avg@y.values <- list(rowMeans(data.frame(perf.avg@y.values)))
  }
  return(perf.avg);
}

ComputeHighLow <- function(perf){
  all.alphas <- unlist(perf@alpha.values);
  min.alpha <- min(all.alphas);
  if(min.alpha == -Inf){
    min.alpha <- 0;
  }
  max.alpha <- max(all.alphas);
  if(max.alpha == Inf){
    max.alpha <- 1.0;
  }
  
  alpha.values <- rev(seq(min.alpha, max.alpha,length=max(sapply(perf@alpha.values, length))));
  perf.sampled <- perf;
  for (i in 1:length(perf.sampled@y.values)) {
    perf.sampled@x.values[[i]] <-
      approxfun(perf@alpha.values[[i]],perf@x.values[[i]],
                rule=2, ties=mean)(alpha.values)
    perf.sampled@y.values[[i]] <-
      approxfun(perf@alpha.values[[i]], perf@y.values[[i]],
                rule=2, ties=mean)(alpha.values)
  }
  ## compute average curve
  y.data <- data.frame(perf.sampled@y.values)
  con.low <- apply(y.data, 1, quantile, 0.05);
  con.high <- apply(y.data, 1, quantile, 0.95);
  res <- list( 
    con.low = con.low,
    con.high = con.high
  );
  return (res);
}

Get.Fstat <-  function(x, fac, var.equal=TRUE) {
  
  x = t(x);
  sqr = function(x) x*x;
  stopifnot(length(fac)==ncol(x), is.factor(fac), is.matrix(x))
  x   <- x[,!is.na(fac), drop=FALSE]
  fac <- fac[!is.na(fac)]
  
  ## Number of levels (groups)
  k <- nlevels(fac)
  
  ## xm: a nrow(x) x nlevels(fac) matrix with the means of each factor level
  xm <- matrix(
    sapply(levels(fac), function(fl) rowMeans(x[,which(fac==fl), drop=FALSE])),
    nrow = nrow(x),
    ncol = nlevels(fac))
  
  ## x1: a matrix of group means, with as many rows as x, columns correspond to groups
  x1 <- xm[,fac, drop=FALSE]
  
  ## degree of freedom 1
  dff    <- k - 1
  
  ## x0: a matrix of same size as x with overall means
  x0 <- matrix(rowMeans(x), ncol=ncol(x), nrow=nrow(x))
  
  ## degree of freedom 2
  dfr    <- ncol(x) - dff - 1
  
  ## mean sum of squares
  mssf   <- rowSums(sqr(x1 - x0)) / dff
  mssr   <- rowSums(sqr( x - x1)) / dfr
  
  ## F statistic
  fstat  <- mssf/mssr
  return(fstat)
}

Get.Fisher <- function(x, fac, var.equal=TRUE) {
  inx1 <- which(y==levels(y)[1]);
  inx2 <- which(y==levels(y)[2]);
  p.value <- apply(as.matrix(x), 2,
                   function(x) {
                     tmp <- try(fisher.test(x[inx1], x[inx2]));
                     if(class(tmp) == "try-error") {
                       return(NA);
                     }else{
                       return(tmp$p.value);
                     }
                   });
  -log10(p.value);
}

PlotProbView <- function(imgName, format = "png", dpi = 72, mdl.inx, show, showPred) {
#save.image("probview.RData");
  rdtSet <- .get.rdt.set()
  smpl.nms <- rownames(t(rdtSet$dataSet$roc.norm))
  prob.vec <- rep(0.5, length(smpl.nms))
  names(prob.vec) <- smpl.nms
  target_group <- as.character(rdtSet$analSet$multiROC$target_group)
  other_group <- as.character(rdtSet$analSet$multiROC$other_group)
  if (mdl.inx == -1) {
    mdl.inx <- rdtSet$analSet$multiROC$best.model.inx
  }
  probs <- MergeDuplicates(unlist(rdtSet$analSet$multiROC$pred.cv[[mdl.inx]]))
  prob.vec[names(probs)] <- probs
  
  nms <- names(prob.vec)
  ord.inx <- order(nms)
  prob.vec <- prob.vec[ord.inx]
  # Convert `cls` to character vector to prevent factor-related issues
  cls <- rdtSet$dataSet$roc.cls[ord.inx]

  # Apply one-against-all transformation based on the target_group
 

  # Check the transformed `cls` vector
  print(cls)  # Confirm it shows only "IGT" and "other"
  
 
  # Proceed with plotting
  pred.out <- factor(ifelse(prob.vec > 0.5, target_group, other_group), levels = c(target_group, other_group))
  act.cls <- factor(cls, levels = c(target_group, other_group))
  
  prob.res <- data.frame(Probability = prob.vec, Predicted = pred.out, Actual = act.cls)
  write.table(prob.res, file = "roc_pred_prob.csv", sep = ",", col.names = TRUE)
  
  conf.res <- table(act.cls,pred.out);
  rdtSet$analSet$conf.table <- xtable::xtable(conf.res, caption = "Confusion Matrix (Cross-Validation)")
  rdtSet$analSet$conf.mat <- print(rdtSet$analSet$conf.table, type = "html", print.results = FALSE, caption.placement = "top", html.table.attributes = "border=1 width=150")
  
  imgName <- paste(imgName, "dpi", dpi, ".", format, sep = "")
  w <- 9; h <- 8
  rdtSet$imgSet$roc.prob.plot <- imgName
  rdtSet$imgSet$roc.prob.name <- mdl.inx
  
  Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, width = w, height = h, type = format, bg = "white")
  
  set.seed(123)
  y <- rnorm(length(prob.vec))
  max.y <- max(abs(y))
  ylim <- max.y * c(-1.05, 1.05)
  xlim <- c(0, 1.0)
  
  op <- par(mar = c(4, 4, 3, 6))
  pchs <- ifelse(cls == target_group, 19, 1)
  
  colors <- ifelse(show == 1, "darkgrey", "black")
  plot(prob.vec, y, pch = pchs, col = colors, xlim = xlim, ylim = ylim, xlab = "Predicted Class Probabilities", ylab = "Samples")
  abline(h = 0, lty = 2, col = "grey")
  abline(v = 0.5, lty = 2, col = "grey")
  
  par(xpd = TRUE)
  legend("right", inset = c(-0.11, 0), legend = c(target_group, other_group), pch = c(19, 1))
  
  if (showPred && !is.null(rdtSet$analSet$multiROC$test.res)) {
    test.y <- rnorm(length(rdtSet$analSet$multiROC$test.res))
    test.x <- rdtSet$analSet$multiROC$test.res
    pchs <- ifelse(rdtSet$dataSet$test.cls == target_group, 19, 1)
    points(test.x, test.y, pch = pchs, cex = 1.5, col = "red")
  }
  
  if (show == 1) {
    pred.vec <- ifelse(prob.vec > 0.5, target_group, other_group)
    misclassified <- (pred.vec != as.character(cls))
    text(prob.vec[misclassified], y[misclassified], nms[misclassified], pos = ifelse(pred.vec[misclassified] == target_group, 4, 2))
    
    if (showPred && !is.null(rdtSet$dataSet$test.cls)) {
      nms <- rownames(rdtSet$dataSet$test.data)
      pred.vec <- ifelse(test.x > 0.5, target_group, other_group)
      misclassified_test <- (pred.vec != as.character(rdtSet$dataSet$test.cls))
      text(test.x[misclassified_test], test.y[misclassified_test], nms[misclassified_test], pos = ifelse(pred.vec[misclassified_test] == target_group, 4, 2), cex = 0.9)
    }
  }
  par(op)
  dev.off()
  
  return(.set.rdt.set(rdtSet))
}


 PlotROC <- function(imgName, format="png", dpi=72, mdl.inx, avg.method, show.conf, show.holdout, focus="fpr", cutoff = 1.0){
  
  rdtSet <- .get.rdt.set();
  anal.mode <- rdtSet$analSet$mode;
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 8; h <- 8;

  rdtSet$imgSet$roc.multi.plot <- imgName;
  rdtSet$imgSet$roc.multi.model <- mdl.inx;

  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  op <- par(mar=c(5,4,3,3));
  
  if(mdl.inx == 0){ # compare all models
    preds <- rdtSet$analSet$multiROC$pred.list;
    auroc <- rdtSet$analSet$multiROC$auc.vec;
    perf <- ROCR::performance(preds[[1]], "tpr", "fpr");
    perf.avg <- ComputeAverageCurve(perf, avg.method);
    
    cols <- (1:length(preds))+1;
    ROCR::plot(perf.avg@x.values[[1]], perf.avg@y.values[[1]], type="n", axes=F,
         xlim=c(0,1), ylim=c(0,1),
         xlab="1-Specificity (False positive rate)",
         ylab="Sensitivity (True positive rate)"
    );
    
    box()
    axis(side=2)
    lab.at <- seq(0, 1, 0.2);
    grid.at <- seq(0, 1, 0.1);
    lab.labels <- lab.at
    axis(side=1, at=lab.at, labels=as.graphicsAnnot(sprintf( "%.1f", lab.labels)));
    abline(v=grid.at, lty=3, col="lightgrey");
    abline(h=grid.at, lty=3, col="lightgrey");
    lines(perf.avg@x.values[[1]], perf.avg@y.values[[1]], col=cols[1]);
    for(i in 2:length(preds)){
      perf <- ROCR::performance(preds[[i]], "tpr", "fpr");
      avg <- ComputeAverageCurve(perf, avg.method);
      lines(avg@x.values[[1]], avg@y.values[[1]], col=cols[i]);
    }
    
    best.inx <- which.max(auroc);
    
    # now add and format legends to the bottom right corner
    feats <- c("Var.", names(preds));
    feats <- substr(feats, 1, 8);
    feats <- sprintf("%-5s", feats);
    
    vals <- c("AUC", round(auroc, 3));
    
    vals <- sprintf("%-8s", vals);
    
    cis <- rdtSet$analSet$multiROC$auc.ci;
    cis <- c("CI", cis);
    legends <- paste(feats, vals, cis, sep="");
    
    pch <- c(NA, rep(15, length(preds)));
    cols <- c(NA, cols);
    
    legend("bottomright", legend = legends, pch=15, col=cols);
    
  }else if(mdl.inx > 0){ 
    
    preds <- ROCR::prediction(rdtSet$analSet$multiROC$pred.cv[[mdl.inx]], rdtSet$analSet$multiROC$true.cv);
    auroc <- round(rdtSet$analSet$multiROC$auc.vec[mdl.inx],3);
    auc.ci <- rdtSet$analSet$multiROC$auc.ci[mdl.inx];
    perf <- ROCR::performance(preds, "tpr", "fpr");
    perf.avg <- ComputeAverageCurve(perf, avg.method);
    y.all <- perf.avg@y.values[[1]];
    x.all <- perf.avg@x.values[[1]];
    lgd <- paste("Area under the curve (AUC) =", auroc, "\n",
                 "95% CI:", auc.ci);
    
    ROCR::plot(x.all, y.all, type="n", axes=F,
         xlim=c(0,1), ylim=c(0,1),
         xlab="1-Specificity (False positive rate)",
         ylab="Sensitivity (True positive rate)"
    );
    
    box()
    axis(side=2)
    lab.at <- seq(0, 1, 0.2);
    grid.at <- seq(0, 1, 0.1);
    lab.labels <- lab.at
    axis(side=1, at=lab.at, labels=as.graphicsAnnot(sprintf( "%.1f", lab.labels)));
    abline(v=grid.at, lty=3, col="lightgrey");
    abline(h=grid.at, lty=3, col="lightgrey");
    lines(x.all, y.all, type="l", lwd=2, col="blue");
    
    if(show.conf){
      res <- ComputeHighLow(perf);
      suppressWarnings(polygon(c(x.all, rev(x.all)), c(res$con.low, rev(res$con.high)), col="#0000ff22"))
    }
    
    legend("center", legend = lgd, bty="n");
    
  }
  dev.off();
  return(.set.rdt.set(rdtSet));
}

GetCurrentConfMat <- function(){
  rdtSet <- .get.rdt.set();
  return(rdtSet$analSet$conf.mat);
}

GetCurrentConfMatTest <- function(){
  rdtSet <- .get.rdt.set();
  return(rdtSet$analSet$conf.mat.test);
}

GetBestModelIndex <- function(){
  rdtSet <- .get.rdt.set();
  rdtSet$analSet$multiROC$best.model.inx;
}

#'Set custom data
#'@description The "selected.cmpds" should be for extraction
#'@param rdtSet Input the name of the created mSetObj (see InitDataObjects)
#'@param selected.cmpds Input the vector containing the compounds
#'@param selected.smpls Input the vector containing the samples
#'@export
#'
SetCustomData <- function(rdtSet=NA, selected.cmpds, selected.smpls){
 
  rdtSet <- .get.rdt.set();
  
  data.norm.orig <- rdtSet$dataSet$roc.norm;
  cls <- rdtSet$dataSet$roc.cls;
  
  if(length(selected.cmpds) > 0){
    data.norm <- data.norm.orig[rownames(data.norm.orig)%in% selected.cmpds, , drop=F];
    if(!is.null(rdtSet$dataSet$new.data)){
      rdtSet$dataSet$new.data <- rdtSet$dataSet$new.data[rownames(data.norm.orig)%in% selected.cmpds, , drop=F];
    }
  }
  
  if(length(selected.smpls) > 0){
    hit.inx <- colnames(data.norm) %in% selected.smpls;
    rdtSet$dataSet$test.data <- data.norm[, hit.inx,drop=F];
    rdtSet$dataSet$test.cls <- cls[hit.inx];
    data.norm <- data.norm[,!hit.inx ,drop=F];
    cls <- cls[!hit.inx];
  }else{
    rdtSet$dataSet$test.data <- NULL;
    rdtSet$dataSet$test.cls <- NULL;
  }
  
  rdtSet$dataSet$roc.norm.orig <- data.norm.orig;
  rdtSet$dataSet$roc.norm <- data.norm;
  rdtSet$dataSet$selected.cmpds <- paste(selected.cmpds, collapse="; ");
  rdtSet$dataSet$roc.cls <- cls;
  return(.set.rdt.set(rdtSet));
}




#'Perform MCCV for manually selected features
#'@description MCCV for manually selected features (no additional feature selection)
#'@usage PerformCV.test(mSetObj, method, lvNum, propTraining=2/3, nRuns=100)
#'@param rdtSet Input the name of the created mSetObj (see InitDataObjects)
#'@param method Select the classification method, "rf" for random forest classification, "pls" for PLS-DA,
#'and "svm" for support vector machine 
#'@param lvNum Input the number of latent variables to include in the analyis, only for PLS-DA classification  
#'@param propTraining Input the proportion of samples to use for training, by default it is 2/3 
#'@param nRuns Input the number of MCCV runs, by default it is 100 
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PerformCV.test <- function(rdtSet=NA, method='svm', lvNum=2, propTraining=2/3, nRuns=100){
  
  rdtSet <-  .get.rdt.set();
  rdtSet$analSet$tester.method <- method;
  rdtSet$analSet$tester.lvNum <- lvNum;
  data <- t(rdtSet$dataSet$roc.norm);
  cls <- rdtSet$dataSet$roc.cls;    
  target_group <- levels(cls)[1]
  other_group <- levels(cls)[2]
  if(method == "lr") {
    # check cls (response variable) whether it is number 0/1 or not
    if( length(levels(rdtSet$dataSet$roc.cls)) == 2 ) {
      rdtSet$dataSet$roc.cls.lbl <- levels(rdtSet$dataSet$roc.cls);  # already sorted
      cls <- as.numeric(rdtSet$dataSet$roc.cls)-1;
      rdtSet$dataSet$roc.cls.lbl.new <- sort(unique(cls));
      rdtSet$dataSet$cls01 <- cls; # integer values for response of logistic regression
    } else {
      AddErrMsg("level of response variable y/class has more than 2!");
      return(0);
    }
  } else {
    if(method == "pls"){# make sure the feature # > comp #
      if(lvNum > ncol(data)){
        AddErrMsg("PLS components cannot be more than total features selected!");
        return(0);
      }
    }
    #cls <- rdtSet$dataSet$roc.cls;
  }
  
  splitMat <- GetTrainTestSplitMat(cls, propTraining, nRuns);
  trainRuns <- splitMat$training.mat;
  testRuns <- splitMat$testing.mat;
  
  feat.outp <- actualCls <- perf.outp <- vector(length = nRuns, mode = "list");
  auc.vec <- accu.vec <- vector(mode="numeric", length=nRuns);
  
  for (irun in 1:nRuns){
    trainingSampleRun <- trainRuns[irun, ]
    testSampleRun <- testRuns[irun, ];
    x.train <- data[trainingSampleRun, ,drop=F];
    x.test <- data[testSampleRun, ,drop=F];
    y.train <- cls[trainingSampleRun];
    actualCls[[irun]] <- y.test <-  cls[testSampleRun];
    y.test.numeric <- as.numeric(y.test)-1
    y.train.numeric <- factor(as.numeric(y.train)-1)

    res <- Predict.class(x.train, y.train, x.test, method, lvNum, imp.out=T);
    feat.outp[[irun]] <- res$imp.vec;
    prob.out <- res$prob.out;
    
    # calculate AUC for each
    pred <- ROCR::prediction(prob.out, y.test);
    auc.vec[irun] <- slot(ROCR::performance(pred, "auc"), "y.values")[[1]];
   perf.outp[[irun]] <- prob.out;
     if (!all(levels(y.test) ==levels(factor(pred@labels[[1]]))) ) {
  auc.vec[irun] <- 1-auc.vec[irun] 
   perf.outp[[irun]] <- 1-prob.out;
}
    pred.out <- factor(ifelse(prob.out > 0.5, target_group, other_group ),levels=c(target_group, other_group));

  #print(c(levels(pred.out),"pred.out"))
  accu.vec[irun] <- Get.Accuracy(table(pred.out, y.test));
  }
 
  #############################################################################
  ## prepare results for default plotting
  ## 1) get best model based on AUROC for prob.view and imp.feature calculation
  ## 2) calculate accuracy and roc for all models under comparison
  #############################################################################
  
  prob.vec <- unlist(perf.outp);
  act.vec <- unlist(actualCls);
  preds <- ROCR::prediction(prob.vec, act.vec);
  auc <- mean(auc.vec);
  auc.ci <- GetCIs(as.matrix(auc.vec));
  
  #########################################
  # if there is holdout sample, do prediction
  if(!is.null(rdtSet$dataSet$test.data)){
    test.res <- Predict.class(data, cls, rdtSet$dataSet$test.data, method, lvNum);
  }else{
    test.res <- NULL;
  }
  #######################
  #########################################
  # if there is new samples, do prediction
  if(!is.null(rdtSet$dataSet$new.data)){
    new.res <<- Predict.class(data, cls, rdtSet$dataSet$new.data, method, lvNum);
  }else{
    new.res <- NULL;
  }
  
  #########################################
  # method = Logistic Regression
  # generate report tables for model with 10-fold Cross-Validation
  if( method == "lr") {
    x.train.all <- data;
    x.test.all  <- data;
    y.train.all <- as.character(cls);
    y.test.all  <- as.character(cls);
    
    # generate LR model and AUROC statistics. Then assign to the analSet list
    # return (list(model = mdlSummary, LR.auc=LR.auc));
    
    tbl.cls <- table(cls);
    if( (tbl.cls[[1]] < 10) & (tbl.cls[[2]] < 10) ) {
      AddErrMsg("The sample size of each group should be more than 10!");
      return (0);
    } else {
      doLogisticRegMdl(x.train.all, y.train.all, x.test.all, y.test.all);
    }
  }
  #######################
  if(!is.null(rdtSet$analSet$ROCtest$perm.res)){
    rdtSet$analSet$ROCtest <-  list(
      type = rdtSet$analSet$type,  # note, do not forget to carry the type "roc" when overwrite analSet
      train.mat = trainRuns,
      pred.cv = perf.outp,
      true.cv = actualCls,
      imp.cv = feat.outp,
      perm.res = rdtSet$analSet$ROCtest$perm.res,
      # record processed output
      pred.list = preds,
      accu.mat = t(as.matrix(accu.vec)),
      auc.ci = auc.ci,
      auc.vec = auc,
      test.res = test.res,
      new.res = new.res,
      target_group= target_group,
      other_group= other_group
    );
    
  } else {
    rdtSet$analSet$ROCtest <-  list(
      type = rdtSet$analSet$type,  # note, do not forget to carry the type "roc" when overwrite analSet
      train.mat = trainRuns,
      pred.cv = perf.outp,
      true.cv = actualCls,
      imp.cv = feat.outp,
      
      # record processed output
      pred.list = preds,
      accu.mat = t(as.matrix(accu.vec)),
      auc.ci = auc.ci,
      auc.vec = auc,
      test.res = test.res,
      new.res = new.res,
      target_group = target_group,
      other_group = other_group
    );
  }
  
  return(.set.rdt.set(rdtSet));
}


#'Plot ROC for the ROC Curve Based Model Creation and Evaluation module
#'@description Plot the ROC curve of the biomarker model created using a user-selected subset of features.
#'Pred and auroc are lists containing predictions and labels from different cross-validations. 
#'@usage PlotROCTest(mSetObj=NA, imgName, format="png", 
#'dpi=72, mdl.inx, avg.method, show.conf, show.holdout, focus="fpr", cutoff = 1.0)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param format Select the image format, "png", of "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", 
#'users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param mdl.inx Model index, 0 means to compare all models, 
#'input 1-6 to plot a ROC curve for one of the top six models  
#'@param avg.method Input the method to compute the average ROC curve, 
#'either "threshold", "vertical" or "horizontal"
#'@param show.conf Logical, if 1, show confidence interval, if 0 do not show
#'@param show.holdout Logical, if 1, show the ROC curve for hold-out validation, if 0 do not show 
#'@param focus "fpr" 
#'@param cutoff Input the threshold to limit the calculation of the ROC curve, the number must be between 0 and 1.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export 

PlotROCTest<-function(rdtSet=NA, imgName, format="png", dpi=72, mdl.inx, avg.method, show.conf, show.holdout, focus="fpr", cutoff = 1.0){
  
    rdtSet <-  .get.rdt.set();
  anal.mode <- rdtSet$analSet$mode;
  target_group <- rdtSet$analSet$ROCtest$target_group
  other_group <- rdtSet$analSet$ROCtest$other_group
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 8; h <- 8;
  rdtSet$imgSet$roc.testcurve.plot <- imgName;
  rdtSet$imgSet$roc.testcurve.name <- mdl.inx
  rdtSet$imgSet$roc.testcurve.method <- avg.method
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  op <- par(mar=c(5,4,3,3));
  
  if(anal.mode=="explore" && mdl.inx == 0){ # compare all models
    preds <- rdtSet$analSet$ROCtest$pred.list;
    auroc <- rdtSet$analSet$ROCtest$auc.vec;
    perf <- ROCR::performance(preds[[1]], "tpr", "fpr");
    perf.avg <- ComputeAverageCurve(perf, avg.method);
    
    cols <- (1:length(preds))+1;
    ROCR::plot(perf.avg@x.values[[1]], perf.avg@y.values[[1]], type="n", axes=F,
         xlim=c(0,1), ylim=c(0,1),
         xlab="1-Specificity (False positive rate)",
         ylab="Sensitivity (True positive rate)"
    );
    
    box()
    axis(side=2)
    lab.at <- seq(0, 1, 0.2);
    grid.at <- seq(0, 1, 0.1);
    lab.labels <- lab.at
    axis(side=1, at=lab.at, labels=as.graphicsAnnot(sprintf( "%.1f", lab.labels)));
    abline(v=grid.at, lty=3, col="lightgrey");
    abline(h=grid.at, lty=3, col="lightgrey");
    lines(perf.avg@x.values[[1]], perf.avg@y.values[[1]], col=cols[1]);
    for(i in 2:length(preds)){
      perf <- ROCR::performance(preds[[i]], "tpr", "fpr");
      avg <- ComputeAverageCurve(perf, avg.method);
      lines(avg@x.values[[1]], avg@y.values[[1]], col=cols[i]);
    }
    
    best.inx <- which.max(auroc);
    
    # now add and format legends to the bottom right corner
    feats <- c("Var.", names(preds));
    feats <- substr(feats, 1, 8);
    feats <- sprintf("%-5s", feats);
    
    vals <- c("AUC", round(auroc, 3));
    
    vals <- sprintf("%-8s", vals);
    
    cis <- rdtSet$analSet$multiROC$auc.ci;
    cis <- c("CI", cis);
    legends <- paste(feats, vals, cis, sep="");
    
    pch <- c(NA, rep(15, length(preds)));
    cols <- c(NA, cols);
    
    legend("bottomright", legend = legends, pch=15, col=cols);
    
  }else if(mdl.inx > 0 && anal.mode=="explore"){ 
    
    preds <- ROCR::prediction(rdtSet$analSet$ROCtest$pred.cv[[mdl.inx]], rdtSet$analSet$ROCtest$true.cv);
    auroc <- round(rdtSet$analSet$ROCtest$auc.vec[mdl.inx],3);
    auc.ci <- rdtSet$analSet$ROCtest$auc.ci[mdl.inx];
    perf <- ROCR::performance(preds, "tpr", "fpr");
    perf.avg <- ComputeAverageCurve(perf, avg.method);
    y.all <- perf.avg@y.values[[1]];
    x.all <- perf.avg@x.values[[1]];
    lgd <- paste("Area under the curve (AUC) =", auroc, "\n",
                 "95% CI:", auc.ci);
    
    ROCR::plot(x.all, y.all, type="n", axes=F,
         xlim=c(0,1), ylim=c(0,1),
         xlab="1-Specificity (False positive rate)",
         ylab="Sensitivity (True positive rate)"
    );
    
    box()
    axis(side=2)
    lab.at <- seq(0, 1, 0.2);
    grid.at <- seq(0, 1, 0.1);
    lab.labels <- lab.at
    axis(side=1, at=lab.at, labels=as.graphicsAnnot(sprintf( "%.1f", lab.labels)));
    abline(v=grid.at, lty=3, col="lightgrey");
    abline(h=grid.at, lty=3, col="lightgrey");
    lines(x.all, y.all, type="l", lwd=2, col="blue");
    
    if(show.conf){
      res <- ComputeHighLow(perf);
      suppressWarnings(polygon(c(x.all, rev(x.all)), c(res$con.low, rev(res$con.high)), col="#0000ff22"))
    }
    
    legend("center", legend = lgd, bty="n");
    
  }else{ # plot ROC of specific model and save the table for details
    
    preds <- ROCR::prediction(rdtSet$analSet$ROCtest$pred.cv, rdtSet$analSet$ROCtest$true.cv);
    auroc <- round(rdtSet$analSet$ROCtest$auc.vec[1],3)
    auc.ci <- rdtSet$analSet$ROCtest$auc.ci;
    
    perf <- ROCR::performance(preds, "tpr", "fpr");
    perf.avg <- ComputeAverageCurve(perf, avg.method);
    y.all <- perf.avg@y.values[[1]];
    x.all <- perf.avg@x.values[[1]];
    # to draw a roc curve line from (0,0)
    y.all <- c(0.0, y.all);
    x.all <- c(0.0, x.all);
    
    lgd <- paste("Area under the curve (AUC) =", auroc, "\n",
                 "95% CI:", auc.ci);
    
    ROCR::plot(x.all, y.all, type="n", axes=F,
         xlim=c(0,1), ylim=c(0,1),
         xlab="1-Specificity (False positive rate)",
         ylab="Sensitivity (True positive rate)"
    );
    
    box()
    axis(side=2)
    lab.at <- seq(0, 1, 0.2);
    grid.at <- seq(0, 1, 0.1);
    lab.labels <- lab.at
    axis(side=1, at=lab.at, labels=as.graphicsAnnot(sprintf( "%.1f", lab.labels)));
    abline(v=grid.at, lty=3, col="lightgrey");
    abline(h=grid.at, lty=3, col="lightgrey");
    lines(x.all, y.all, type="l", lwd=2, col="blue");
    
    if(show.conf){
      res <- ComputeHighLow(perf);            
      # to draw a roc curve line from (0,0)
      # suppressWarnings(polygon(c(x.all, rev(x.all)), c(res$con.low, rev(res$con.high)), col="#0000ff22"))
      suppressWarnings(polygon(c(x.all, rev(x.all)), c(c(0,res$con.low), c(rev(res$con.high),0)), col="#0000ff22"))
    }
    if(show.holdout){

      roc.obj <- pROC::roc(rdtSet$dataSet$test.cls, rdtSet$analSet$ROCtest$test.res, percent = F);
      test.x <- 1-roc.obj$spec;
      test.y <- roc.obj$sens;
      
      lbls <- c("Type", "CV","Holdout");
      lbls <- sprintf("%-10s",lbls);
      
      test.auc <- round(roc.obj$auc[[1]],3);
      aucs <- c("AUC", auroc, test.auc);
      
      lgd <- paste(lbls, aucs, sep="");
      lines(test.x, test.y, type="l", lwd=2, col="magenta");
      legend("bottomright", legend = lgd, pch=c(NA, 15, 15), col=c(NA, "blue", "magenta"));
    }else{
      legend("center", legend = lgd,  bty="n");
    }       
  }
  dev.off();
  return(.set.rdt.set(rdtSet));
}


#'Plot a summary view of the classification result of tester prediction
#'@description Plot of predicted class probabilities. On the x-axis is the proability, 
#'and the y-axis is the index of each predicted sample based on the probility. 
#'The samples are turned into separations at the x-axis.
#'This plot can be created for multivariate ROC curve analysis using SVM, PLS, and RandomForest.
#'Please note that sometimes, not all samples will be tested, instead they will be plotted
#'at the 0.5 neutral line. 
#'@usage PlotProbViewTest(mSetObj=NA, imgName, format="png", dpi=72, mdl.inx, show, showPred) 
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param format Select the image format, "png", of "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300. 
#'@param mdl.inx Model index, 0 means to compare all models, -1 means to use the best model, input 1-6 to plot a ROC curve for one of the top six models
#'@param show 1 or 0, if 1, label samples classified to the wrong groups 
#'@param showPred Show predicted samples
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotProbViewTest <- function(rdtSet=NA, imgName, format="png", dpi=72, mdl.inx, show, showPred) {
  
   rdtSet <-  .get.rdt.set(); 
  anal.mode <- rdtSet$analSet$mode;
   target_group <- rdtSet$analSet$ROCtest$target_group
  other_group <- rdtSet$analSet$ROCtest$other_group

  smpl.nms <- colnames(rdtSet$dataSet$roc.norm);
  prob.vec <- rep(0.5, length(smpl.nms));
  names(prob.vec) <- smpl.nms;
  
  probs <- MergeDuplicates(unlist(rdtSet$analSet$ROCtest$pred.cv));
 
  prob.vec[names(probs)] <- probs;
  
  nms <- names(prob.vec);
  ord.inx <- order(nms);
  prob.vec <- prob.vec[ord.inx];
    cls <-  rdtSet$dataSet$roc.cls[ord.inx];
  # remember to update the nms itself!
  nms <- names(prob.vec);
  
  # get html confusion matrix
  pred.out <- factor(ifelse(prob.vec > 0.5, target_group, other_group), levels = c(target_group, other_group));
    act.cls <- factor(cls, levels = c(target_group, other_group));
  
  prob.res <- data.frame(Probability=prob.vec, Predicted=pred.out, Actual=act.cls);
  
  write.table(prob.res, file="roc_pred_prob1.csv", sep=",", col.names = TRUE);
  
  conf.res <- table(act.cls,pred.out);
  rdtSet$analSet$conf.table <- xtable::xtable(conf.res, caption="Confusion Matrix (Cross-Validation)");
  rdtSet$analSet$conf.mat <- print(rdtSet$analSet$conf.table, type = "html", print.results=F, caption.placement="top", html.table.attributes="border=1 width=150" )     
  
  if(anal.mode == "test"){
    if(!is.null(rdtSet$dataSet$test.data)){
      
      test.pred <-  factor(ifelse(rdtSet$analSet$ROCtest$test.res > 0.5, target_group, other_group), levels = c(target_group, other_group));
      test.cls <-  factor(rdtSet$dataSet$test.cls, levels = c(target_group, other_group));
      
      test.df <- data.frame(Prob_HoldOut=rdtSet$analSet$ROCtest$test.res, Predicted_HoldOut=test.pred, Actual_HoldOut=test.cls);
      suppressMessages(write.table(test.df, file="roc_pred_prob1.csv", sep=",", append=TRUE, col.names = TRUE));
      
      test.res <- table(test.pred, test.cls);
      rdtSet$analSet$conf.mat.test <- print(xtable::xtable(test.res, 
                                                            caption="Confusion Matrix (Hold-out)"),
                                             type = "html", print.results=F, xtable.width=120, caption.placement="top",
                                             html.table.attributes="border=1 width=150" );
    }
  }
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 9; h <- 8;
  
  rdtSet$imgSet$roc.testprob.plot <- imgName;
  rdtSet$imgSet$roc.testprob.name <- mdl.inx

  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  set.seed(123);
  y <- rnorm(length(prob.vec));
  max.y <- max(abs(y));
  ylim <- max.y*c(-1.05, 1.05);
  
  xlim <- c(0, 1.0);
  
  op <- par(mar=c(4,4,3,6));
  pchs <- ifelse( cls  == other_group, 19, 1);
  
  colors <- ifelse(show==1, "darkgrey", "black");
  ROCR::plot(prob.vec, y, pch=pchs, col=colors, xlim=xlim, ylim= ylim, xlab = "Predicted Class Probabilities", ylab="Samples");
  abline(h = 0, lty = 2, col="grey");
  abline(v = 0.5, lty = 2, col="grey");
  
  par(xpd=T);
 
 legend("right", inset = c(-0.12, 0), legend = c(target_group, other_group), pch = c(1, 19))

  test.y <- test.x <- 0;
  if(showPred){
    test.y <- rnorm(length(rdtSet$analSet$ROCtest$test.res));
    test.x <- rdtSet$analSet$ROCtest$test.res;
   
    pchs <- ifelse(rdtSet$dataSet$test.cls==other_group, 19, 1);
    points(test.x, test.y, pch=pchs, cex=1.5, col="red");
  }
  
  if(show == 1){ 
    
        pred.vec <- ifelse(prob.vec > 0.5, target_group, other_group)
    misclassified <- (pred.vec != as.character(cls))
    text(prob.vec[misclassified], y[misclassified], nms[misclassified], pos = ifelse(pred.vec[misclassified] == other_group, 4, 2))
    
    if (showPred && !is.null(rdtSet$dataSet$test.cls)) {
      nms <- rownames(rdtSet$dataSet$test.data)
      pred.vec <- ifelse(test.x > 0.5, target_group, other_group)
      misclassified_test <- (pred.vec != as.character(rdtSet$dataSet$test.cls))
      text(test.x[misclassified_test], test.y[misclassified_test], nms[misclassified_test], pos = ifelse(pred.vec[misclassified_test] == other_group, 4, 2), cex = 0.9)
    }
  }
  par(op)
  dev.off();
  return(.set.rdt.set(rdtSet));
}


#'Plot classification performance using different features for Biomarker Tester
#'@description Plot of the accuracy of classification with an increasing number of features.
#'@usage PlotTestAccuracy(mSetObj=NA, imgName, format="png", dpi=72)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param format Select the image format, "png", of "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotTestAccuracy<-function(rdtSet=NA, imgName, format="png", dpi=72){
  
  rdtSet <-  .get.rdt.set();
  anal.mode <- rdtSet$analSet$mode;
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 9; h <- 7;
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  y.vec <- rdtSet$analSet$ROCtest$accu.mat[1,];
  ylim.ext <- GetExtendRange (y.vec, 12); # first increase ylim by 1/12
  boxplot(y.vec, col="#0000ff22", ylim=ylim.ext, outline=FALSE, boxwex=c(0.5, 0.5), ylab="Predictive Accuracy");
  stripchart(t(rdtSet$analSet$ROCtest$accu.mat), method = "jitter", vertical=T, add = T, pch=19);
  
  accu.info <- paste ("The average accuracy based on 100 cross validations is", round(mean(y.vec), 3));
  
  rdtSet$imgSet$roc.testpred <- imgName;
  
  if(!is.null(rdtSet$dataSet$test.cls)){
    test.pred <- ifelse(rdtSet$analSet$ROCtest$test.res > 0.5, target_group, other_group);
    test.cls <- as.numeric(rdtSet$dataSet$test.cls)-1;
    
    hit <- sum(test.pred == test.cls);
    percent <- round(hit/length(test.cls), 3);
    accu.info <- paste(accu.info, ". The accuracy for hold out data prediction is ",  percent,
                       "(", hit, "/",  length(test.cls), ").", sep="");
  }
  
  rdtSet$analSet$ROCtest$accu.info <- accu.info;
  
  dev.off();
  return(.set.rdt.set(rdtSet));
}


#'Export biomarker accuracy information
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@export
GetAccuracyInfo<-function(rdtSet=NA){
  rdtSet <-  .get.rdt.set();
  return(rdtSet$analSet$ROCtest$accu.info);
}


#'Perform permutation tests only for ROC Tester 
#'@description Perform permutation tests for the ROC Curve Based Model Creation and Evaluation module
#'@usage Perform.Permut(mSetObj=NA, perf.measure, perm.num, propTraining = 2/3)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param perf.measure Input the performance measure to rate the performance of the model, either the area
#'under the ROC curve ("auroc") or the predictive accuracy ("accu")
#'@param perm.num Input the number of permutations to perform
#'@param propTraining Numeric, input the fraction of samples to set aside for training. Default is set to 2/3.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
Perform.Permut<-function(rdtSet=NA, perf.measure, perm.num, propTraining = 2/3){
  
  rdtSet <-  .get.rdt.set();
  
  cvRuns = perm.num;
  propTraining = propTraining
  
  cls <- rdtSet$dataSet$roc.cls;
  datmat <- t(rdtSet$dataSet$roc.norm);

  if(rdtSet$analSet$mode == "test"){
    clsMethod <- rdtSet$analSet$tester.method;
  }else{
    clsMethod <- rdtSet$analSet$exp.method;
  }

  splitMat <- GetTrainTestSplitMat(cls, propTraining, cvRuns);
  trainInx <- splitMat$training.mat;
  testInx <- splitMat$testing.mat;
  
  # now, do the permutation
  perf.outp <- actualCls <- vector(length = cvRuns, mode = "list");
  # for (irun in 1:cvRuns){
  irun <- 1;
  while(irun < cvRuns){
    cls <- cls[order(runif(length(cls)))];
    trainingSampleRun <- trainInx[irun, ]
    testSampleRun <- testInx[irun, ];
    y.in <- cls[trainingSampleRun];
    # make sure the y.in contain only one group
    if(length(unique(as.numeric(y.in))) > 1){
      y.out <- cls[testSampleRun];
      actualCls[[irun]] <- as.numeric(y.out);
      perf.outp[[irun]] <- Get.pred(datmat[trainingSampleRun,], y.in,
                                    datmat[testSampleRun, ], y.out, clsMethod);
      irun <- irun + 1;
    }else{
      print("redo....");
    }
  }
  
  # get the AUROC for permuted data
  pred <- ROCR::prediction(perf.outp, actualCls);
  aucs <- try(unlist(slot(ROCR::performance(pred, "auc"), "y.values")));
  if (class(aucs)=="try-error"){
     AddErrMsg("Not enough distinct predictions to compute area under the ROC curve. Increase sample size or reduce permutation number.");
     return(0);
  }
  accs <- unlist(slot(ROCR::performance(pred, "acc"), "y.values"));
  perf.obj <- try(ROCR::performance(pred, "tpr", "fpr"));
  if (class(perf.obj)=="try-error"){
     AddErrMsg("Something wrong with the task. Increase sample size or reduce permutation number.");
     return(0);
  }
  # now, insert the average value of the original performance
  accs <- c(mean(rdtSet$analSet$ROCtest$accu.mat[1,]), accs);
  aucs <- c(mean(rdtSet$analSet$ROCtest$auc.vec), aucs);
  rdtSet$analSet$ROCtest$perm.res <- list(perf.measure=perf.measure, perf.obj=perf.obj, auc.vec=aucs, acc.vec=accs);
 return(.set.rdt.set(rdtSet));
}





.do.CVTest.LRmodel <- function(data.in, fmla.in, kfold=10, run.stepwise=FALSE){


  dw.case <- data.in[which(data.in$y == 1), ]; nrow(dw.case)
  dw.ctrl <- data.in[which(data.in$y == 0), ]; nrow(dw.ctrl)
  
  random.seed <- 10063927;
  grpIndice.case <- createCVset(nrow(dw.case), kfold, random.seed)
  grpIndice.ctrl <- createCVset(nrow(dw.ctrl), kfold, random.seed)
  
  all.trainning.y <- NULL
  all.trainning.fit <- NULL
  all.validation.y <- NULL
  all.validation.fit <- NULL
  
  for (i in 1:kfold) 
  {
    d.train.case <- dw.case[-grpIndice.case[[i]], ]; nrow(d.train.case)
    d.train.ctrl <- dw.ctrl[-grpIndice.ctrl[[i]], ]; nrow(d.train.ctrl)
    d.train <- rbind(d.train.case, d.train.ctrl); names(d.train)
    
    d.validation.case <- dw.case[grpIndice.case[[i]], ]; nrow(d.validation.case)
    d.validation.ctrl <- dw.ctrl[grpIndice.ctrl[[i]], ]; nrow(d.validation.ctrl)
    d.validation <- rbind(d.validation.case, d.validation.ctrl)  
    
    vnames <- names(d.train); 
    
    mdl <- glm(fmla.in, data=d.train, family=binomial(link="logit"))
    if(run.stepwise) { mdl <- step(mdl) }
    
    dval.pred <- predict(mdl, d.validation, type="response");  
    
    all.validation.y <- c(all.validation.y, d.validation$y)
    all.validation.fit <- c(all.validation.fit, dval.pred)
    all.trainning.y <- c(all.trainning.y, mdl$y)
    all.trainning.fit <- c(all.trainning.fit, mdl$fitted.values)
  }
  
  # 10-fold cross validation
  cv.r <- pROC::roc(all.validation.y ~ all.validation.fit, ci=T, ci.se=T, sp=seq(0,1,0.01)) # of: se (sensitivity), sp (specificity), thresholds, auc 
  # cv.threshold <- coords(cv.r, "best", ret=c("threshold"), best.method="youden", best.weights=c(5, nrow(dw.case)/nrow(data.in)) ); 
  cv.threshold <- as.numeric(pROC::coords(cv.r, "best", ret=c("threshold"), best.method="closest.topleft", transpose=TRUE)); 
  cv.rstat <- multi.stat(all.validation.fit > cv.threshold, all.validation.y);
  if(length(cv.rstat) == 1 && cv.rstat == 0){
    return(0);
  }
  cv.tbl <- table(all.validation.fit > cv.threshold, all.validation.y);
  
  # training performance
  train.r <- pROC::roc(all.trainning.y ~  all.trainning.fit, ci=T, ci.se=T, sp=seq(0,1,0.01)) # of: se (sensitivity), sp (specificity), thresholds, auc 
  train.threshold <- pROC::coords(train.r, "best", ret=c("threshold"), best.method="youden", transpose=TRUE); 
  train.rstat <- multi.stat(all.trainning.fit > cv.threshold, all.trainning.y) 
  if(length(train.rstat) == 1 && train.rstat == 0){
    return(0);
  } 
  return(list(
    train.r = train.r$ci,
    train.threshold = train.threshold,
    train.rstat = train.rstat,
    cv.robj = cv.r,
    cv.r = cv.r$ci,
    cv.threshold = cv.threshold,
    cv.rstat = cv.rstat,
    cv.tbl = cv.tbl,
    cv.y = all.validation.y,
    cv.pred = all.validation.fit,
    cv.threshold = cv.threshold
  ));
}


doLogisticRegMdl <- function(x.train, y.train, x.test, y.test){
  
  # use 10-fold CV as default; or LOOCV if sample size < 10
  x <- x.train;
  y <- y.train;
  x.cv.test <- x.test;
  y.cv.test <- y.test;
  
  #names.x.origin <- colnames(x);
  #names(x) <- paste0("V", 1:(ncol(x)));
  #names(x.cv.test) <- names(x);
  
  # glm can only take numeric + factor
  if(class(y) == "character"){
    y <- as.factor(y)
  }
  
  # LR Model generation
  data.xy <- data.frame(y, x);
  fmla <- as.formula(paste("y ~", paste(names(data.xy)[-1], collapse="+")));
  
  model <- glm(fmla, data=data.xy, family=binomial(link="logit"), na.action="na.omit");
  
  # if model with selected variables is failed, the use the stepwise variable selection
  if((!model$converged) | (model$deviance < 1)){
    model <- step(model, trace=0);
    fmla <- model$formula;
  }
  
  mdlSummary <- summary(model)$coefficients;
  #dimnames(mdlSummary)[[1]][-1] <- names.x.origin;
  # prob.out <- predict(model, x.cv.test, type="response");
  
  Odds <- round(exp(cbind(OR=coef(model), confint(model))), 2)[,1];
  # Odds <- round(exp(OR=coef(model)), 2);
  Odds[1] <- "-";
  LRmodel <- cbind(round(mdlSummary,3), Odds);
  LRmodel[,4] <- ifelse(LRmodel[,4] < 0.001, "< 0.001", LRmodel[,4]);
  
  # LR model's summary table as html format
  LRmodel.xtable <<- print(xtable::xtable(LRmodel, align="r|rrrcc"),
                           type = "html", print.results=F, xtable.width=600,
                           html.table.attributes="border=1 width=600" );
  
  coefTable <- coef(model);
  #names(coefTable)[-1] <- names.x.origin;
  if (model$converged & model$deviance > 1) {
    LReq <<- genLREquation(coefTable);
    LRConverged <<- "TRUE";
  } else {
    # (Error: Model was not converged)
    LReq <<- "(Error: Model was not converged)";
    LRConverged <<- "FALSE"; 
  }
  
  # ROC analysis with 10-fold CV test set   
  rStat <- .do.CVTest.LRmodel(data.xy, fmla.in=fmla, kfold=10, run.stepwise=FALSE)     
  # r <- roc(y.cv.test ~ prob.out, ci=T, ci.se=T, sp=seq(0,1,0.01)) # of: se (sensitivity), sp (specificity), thresholds, auc
  
  # ROC plot with k-fold CV test set for ROC plot
  LR.r <<- rStat$cv.robj;
  LR.y.origin <<- rStat$cv.y;
  LR.y.pred <<- rStat$cv.pred;
  LR.threshold <<- rStat$cv.threshold;
  
  # training/discovery; 10-fold CV
  auc  <- c( sprintf("%.3f (%.3f ~ %.3f)", rStat$train.r[2], rStat$train.r[1], rStat$train.r[3]),
             sprintf("%.3f (%.3f ~ %.3f)", round(rStat$cv.r[2],3), round(rStat$cv.r[1],3), round(rStat$cv.r[3],3)) );
  sens <- c( sprintf("%.3f (%.3f ~ %.3f)", rStat$train.rstat$sens, rStat$train.rstat$sensCI[1], rStat$train.rstat$sensCI[2]),
             sprintf("%.3f (%.3f ~ %.3f)", rStat$cv.rstat$sens, rStat$cv.rstat[1], rStat$cv.rstat$sensCI[2]) );
  spec <- c( sprintf("%.3f (%.3f ~ %.3f)", rStat$train.rstat$spec, rStat$train.rstat$specCI[1], rStat$train.rstat$specCI[2]),
             sprintf("%.3f (%.3f ~ %.3f)", rStat$cv.rstat$spec, rStat$cv.rstat$specCI[1], rStat$cv.rstat$specCI[2]) );
  
  LRperf <- cbind(auc, sens, spec);
  colnames(LRperf) <- c("AUC","Sensitivity","Specificity");
  rownames(LRperf) <- c("Training/Discovery","10-fold Cross-Validation");
  
  LRperf.xtable <<- print(xtable::xtable(LRperf, align="r|ccc"), include.rownames=TRUE, floating=FALSE,
                          type = "html", print.results=F, xtable.width=600,
                          html.table.attributes="border=1 width=600");
}





PlotROC.LRmodel <- function(rdtSet=NA, imgName, format="png", dpi=72, show.conf=FALSE, sp.bin=0.01) {

  rdtSet <- .get.rdt.set();
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- h <- 8;
  rdtSet$imgSet$roc.lr <- imgName;
  
  roc.object <- LR.r;
  y.origin <- LR.y.origin;
  y.pred <- LR.y.pred;
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  op <- par(mar=c(5,4,3,3));
  
  # auc.ci <- ci.auc(roc.object, method="bootstrap", boot.n=500, progress="none");
  auc.lbl <- paste("Area under the curve (AUC) = ", round(roc.object$ci[2],3), sep="");
  ci.lbl  <- paste("95% CI: ", round(roc.object$ci[1],3), " - ", round(roc.object$ci[3],3),sep="");
  perform.lbl <- paste(auc.lbl,"\n", ci.lbl, sep="");
  
  # plot ROC of specific model and save the table for details
  if(show.conf){
    # of: se (sensitivity), sp (specificity), thresholds, auc
    r.new <- roc(y.origin ~ y.pred, ci=T, of="se", sp=seq(0,1,sp.bin));
    pROC::plot.roc(r.new, add=F, ci.type="bars", ci.col="green", col="darkolivegreen", lty="dotted", legacy.axes=T,
             xlab="1 - Specificity (False positive rate)", ylab="Sensitivity (True positive rate)",
             cex.axis=1.0, cex.lab=1.0);
    pROC::plot.roc(pROC::smooth(r.new, method="density", n=512), add=T, col="blue", lwd=1);
    legend("bottomright", legend=c("Empirical","Smoothed", "95% CI"), cex=1.0, col=c(par("fg"), "blue", "green"), lwd=3, lty=c(3,1,1) );
    
    text(0.7, 0.4, perform.lbl, adj=c(0,1), col="black", cex=1.0);
  } else {
    pROC::plot.roc(roc.object, add=F, ci.type="no", col="blue", legacy.axes=T,
             xlab="1 - Specificity (False positive rate)", ylab="Sensitivity (True positive rate)",
             cex.axis=1.0, cex.lab=1.0, lwd=2, lty=1 );
    text(0.7, 0.4, perform.lbl, adj=c(0,1), col="black", cex=1.0);
  }
  
  dev.off();
return(.set.rdt.set(rdtSet));
}

Get.pred <- function(x.train, y.train, x.test, y.test, clsMethod="pls"){
  
  # first convert class label to 0/1 so convert prob-> class will be easier
  y.train <- as.factor(as.numeric(y.train)-1);
  y.test <- as.factor(as.numeric(y.test)-1);
  
  # note, we only need prob for class 1, pred can be inferred
  if (clsMethod == "rf"){
    model <- randomForest::randomForest(x.train, y.train, ntree=100, importance=F);
    prob.out <- predict(model, x.test, type="prob")[,"1"];
    return(prob.out);
  }else if(clsMethod == "pls"){ # plsda
    model <- caret::plsda(x.train, y.train, method='oscorespls');
    prob.out <- predict(model, x.test, type="prob")[,"1",1];
    return(prob.out);
  }else{ # svm
    model <- e1071::svm(x.train, y.train, type = 'C', kernel="linear", probability=TRUE);
    prob.out <- attr(predict(model, x.test,  probability=TRUE), "probabilities")[,"1"];
    return(prob.out);
  }
}


Plot.Permutation<-function(rdtSet=NA, imgName, format="png", dpi=72){
  
  rdtSet <- .get.rdt.set();
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 8; h <- 8;
  rdtSet$imgSet$roc.perm.plot <- imgName;
  
  if(rdtSet$analSet$ROCtest$perm.res$perf.measure == "auroc"){
    rdtSet$imgSet$roc.perm.method <- "auroc"
  }else{
    rdtSet$imgSet$roc.perm.method <- "predictive accuracy"
  }
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  if(rdtSet$analSet$ROCtest$perm.res$perf.measure == "auroc"){
    ROCR::plot(rdtSet$analSet$ROCtest$perm.res$perf.obj, col="grey", lwd=1, lty=2,
         xlab="1-Specificity (False Positive rate)",
         ylab="Sensitivity (True Positive rate)");
    
    # now add the original ROC
    
    preds <- ROCR::prediction(rdtSet$analSet$ROCtest$pred.cv, rdtSet$analSet$ROCtest$true.cv);
    auroc <- round(rdtSet$analSet$ROCtest$auc.vec[1],3)
    perf <- ROCR::performance(preds, "tpr", "fpr");
    # need to replace Inf with 1
    alpha.vals <- perf@alpha.values;
    perf@alpha.values <- lapply(alpha.vals, function(x){
      x[x==Inf] <- 1;
      x[x==-Inf] <- 0;
      x
    });
    
    
    ROCR::plot(perf,lwd=2,avg="threshold", col="blue", add=T);
    
    # calculate p value
    perm.vec <- rdtSet$analSet$ROCtest$perm.res$auc.vec;
    better.hits <- sum(perm.vec[-1]>=perm.vec[1]);
    num <- length(perm.vec);
    if(better.hits == 0) {
      p <- paste("p <", 1/num);
    }else{
      p <- better.hits/num;
      p <- paste("p =", signif(p, digits=5));
    }
    legend("center", legend = paste('Empirical p-value: ', p), bty="n", cex=1.5);
  }else{ # accuracies
    perms <- PreparePermResult(rdtSet$analSet$ROCtest$perm.res$acc.vec);
    perm.vec <- perms$permut;
    perm.p <- perms$permut.p;
    
    op<-par(mar=c(5,5,2,4));
    
    xlim.ext <- GetExtendRange (perm.vec, 10);
    hst <- hist(perm.vec[-1], breaks = "FD", freq=T, ylab="Frequency", xlim=xlim.ext, xaxt="n", xlab= 'Permutation test statistics', col="lightblue", main="");
    axis(1);
    
    # add the indicator using original label
    h <- max(hst$counts)
    arrows(perm.vec[1], h/5, perm.vec[1], 0, col="red", lwd=2);
    text(perm.vec[1], h/3.5, paste('Observed \n statistic \n', perm.p), xpd=T);
    par(op);
  }
  dev.off();
return(.set.rdt.set(rdtSet));
}