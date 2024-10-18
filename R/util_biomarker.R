#'Prepare data for ROC analysis
#'@description Prepare data for ROC analysis
#'@param rdtSet Input the name of the created rdtSet (see InitDataObjects)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PrepareROCData <- function(type="NA"){
  save.image("prepare.RData");
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
  save.image("calculate.RData");
  
  rdtSet <- .get.rdt.set()
  LRConverged <<- "FALSE"; 
  
  x <- t(rdtSet$dataSet$norm);
  y <- rdtSet$dataSet$cls;
  
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
    colMeans(data[which(rdtSet$dataSet$cls == cls), , drop=FALSE])
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
  if(nrow(rdtSet$dataSet$norm) < 1000){
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
  
  save.image("univroc.RData");
  rdtSet <- .get.rdt.set();
  
  imgName <- rdtSet$dataSet$url.var.nms[feat.nm];
  imgName = paste("roc_univ_", imgName, "_", version, "_dpi", dpi, ".", format, sep="");
  
  data_ori_norm <- rdtSet$dataSet$norm
  if(!is.null(rdtSet$dataSet$norm.orig)){
    data_ori_norm <- rdtSet$dataSet$norm.orig
  }
  
  x <- unname(unlist(data_ori_norm[feat.nm, ]));
  y <- rdtSet$dataSet$cls;
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
  
  data_ori_norm <- rdtSet$dataSet$norm
  if(!is.null(rdtSet$dataSet$norm.orig)){
    data_ori_norm <- rdtSet$dataSet$norm.orig
  }
  
  imgName <- rdtSet$dataSet$url.var.nms[feat.nm];
  imgName = paste("roc_boxplot_", imgName, "_", version, "_dpi", dpi, ".", format, sep="");
  
  x <- unname(unlist(data_ori_norm[feat.nm,]));
  y <- rdtSet$dataSet$cls;
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
  save.image("coords.RData");
  rdtSet <- .get.rdt.set();
  roc.obj <- rdtSet$analSet$roc.obj
  
  if (length(unique(rdtSet$dataSet$cls)) > 2) {
    # Multiclass case
    if (is.null(classLabel) || is.na(classLabel)) {
      stop("Please provide a class label to get ROC coordinates for a specific class.")
    }
    
    # Find the index of the class label
    class_labels <- unique(rdtSet$dataSet$cls)
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
  
  x <- unname(unlist(rdtSet$dataSet$norm[feat.nm, ]));
  y <- rdtSet$dataSet$cls;
  
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

PerformCV.explore <- function(cls.method, rank.method="auroc", lvNum=2, propTraining=2/3){
  rdtSet <- .get.rdt.set();
  
  rdtSet$analSet$exp.method <- cls.method;
  rdtSet$analSet$rank.method <- rank.method;
  rdtSet$analSet$exp.lvNum <- lvNum;

  data <- t(rdtSet$dataSet$norm);
  cls <- rdtSet$dataSet$cls;
  
  # number of subsampling to produce smooth curve
  if(nrow(data) > 500){
    nRuns <- 10;
  }else if(nrow(data) > 200){
    nRuns <- 20;
  }else if(nrow(data) > 100){
    nRuns <- 30;
  }else{
    nRuns <- 50;
  }
  
  nFeatures <- GetFeatureNumbers(ncol(data));
  
  feat.outp <- actualCls <- vector(length = nRuns, mode = "list");
  
  perf.outp <- vector(length = length(nFeatures), mode = "list");
  perf.outp <- lapply(perf.outp, function(x){x <- vector(length = nRuns, mode = "list"); return(x)});
  auc.mat <- accu.mat <- matrix(nrow=nRuns, ncol=length(nFeatures));
  
  splitMat <- GetTrainTestSplitMat(cls, propTraining, nRuns);
  trainRuns <- splitMat$training.mat;
  testRuns <- splitMat$testing.mat;
  
  for (irun in 1:nRuns){
    trainingSampleRun <- trainRuns[irun, ]
    testSampleRun <- testRuns[irun, ];
    x.in <- data[trainingSampleRun, ];
    y.train <- cls[trainingSampleRun];
    actualCls[[irun]] <- y.test <- cls[testSampleRun];
    
    # Feature ranking using only training data, 
    imp.vec <- RankFeatures(x.in, y.train, rank.method, lvNum);
    
    feat.outp[[irun]] <- imp.vec;
    ord.inx <- order(imp.vec, decreasing = TRUE);
    ordData <- data[, ord.inx];
    
    # buliding classifiers for each number of selected features and test on the test data
    for (inum in seq(along = nFeatures)){
      x.train <- ordData[trainingSampleRun, 1:nFeatures[inum]];
      x.test <- ordData[testSampleRun, 1:nFeatures[inum]];
      prob.out <- Predict.class(x.train, y.train, x.test, cls.method, lvNum);
      
      # calculate AUC for each
      pred <- ROCR::prediction(prob.out, y.test);
      auc.mat[irun, inum] <- slot(ROCR::performance(pred, "auc"), "y.values")[[1]];
      
      perf.outp[[inum]][[irun]] <- prob.out;
      pred.out <- as.factor(ifelse(prob.out > 0.5, 1, 0));
      accu.mat[irun, inum] <- Get.Accuracy(table(pred.out, y.test));
    }
  }
  
  #############################################################################
  ## prepare results for default plotting 
  ## 1) get best model based on AUROC for prob.view and imp.feature calculation
  ## 2) calculate accuracy and roc for all models under comparison
  #############################################################################
  
  preds <- vector(length = length(nFeatures), mode = "list");
  act.vec <- unlist(actualCls); # same for all subsets
  for(m in 1:length(nFeatures)){
    prob.vec <- unlist(perf.outp[[m]]);
    pred <- ROCR::prediction(prob.vec, act.vec);
    preds[[m]] <- pred; # prediction obj
    #auc.vec[m] <- slot(performance(pred, "auc"), "y.values")[[1]];
  }
  
  # accu.mat and preds are for all models
  colnames(accu.mat) <- colnames(auc.mat) <- names(preds) <- paste(nFeatures, sep = "");
  auc.vec <- colMeans(auc.mat);
  
  auc.cis <- GetCIs(auc.mat);
  # get the best based on AUROC
  best.model.inx <- which.max(auc.vec);
  
  rdtSet$analSet$multiROC <- list(
    type = rdtSet$analSet$type, # note, do not forget to carry the type "roc"
    train.mat = trainRuns,
    
    # record raw output
    test.feats = nFeatures,
    pred.cv = perf.outp, 
    true.cv = actualCls, 
    imp.cv = feat.outp,
    
    # record processed output
    pred.list = preds,
    accu.mat = accu.mat,
    auc.vec = auc.vec,
    auc.ci = auc.cis,
    best.model.inx = best.model.inx
  )
  return(.set.rdt.set(rdtSet));
}

PlotProbView <- function(imgName, format="png", dpi=72, mdl.inx, show, showPred) {
  
  rdtSet <- .get.rdt.set();
  anal.mode <- rdtSet$analSet$mode;
  
  smpl.nms <- rownames(rdtSet$dataSet$norm);
  prob.vec <- rep(0.5, length(smpl.nms));
  names(prob.vec) <- smpl.nms;
  
  if(mdl.inx == -1){
    mdl.inx <- rdtSet$analSet$multiROC$best.model.inx;
  }
  probs <- MergeDuplicates(unlist(rdtSet$analSet$multiROC$pred.cv[[mdl.inx]]));

  prob.vec[names(probs)] <- probs;
  
  nms <- names(prob.vec);
  ord.inx <- order(nms);
  prob.vec <- prob.vec[ord.inx];
  cls <- rdtSet$dataSet$cls[ord.inx];
  # remember to update the nms itself!
  nms <- names(prob.vec);
  
  # get html confusion matrix
  pred.out <- as.factor(ifelse(prob.vec > 0.5, 1, 0));
  act.cls <- as.numeric(cls)-1;
  
  prob.res <- data.frame(Probability=prob.vec, Predicted=pred.out, Actual=act.cls);
  
  write.table(prob.res, file="roc_pred_prob.csv", sep=",", col.names = TRUE);
  
  conf.res <- table(pred.out, act.cls);
  rdtSet$analSet$conf.table <- xtable::xtable(conf.res, caption="Confusion Matrix (Cross-Validation)");
  rdtSet$analSet$conf.mat <- print(rdtSet$analSet$conf.table, type = "html", print.results=F, caption.placement="top", html.table.attributes="border=1 width=150" )     
  

  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 9; h <- 8;

  rdtSet$imgSet$roc.prob.plot <- imgName;
  rdtSet$imgSet$roc.prob.name <- mdl.inx


  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  set.seed(123);
  y <- rnorm(length(prob.vec));
  max.y <- max(abs(y));
  ylim <- max.y*c(-1.05, 1.05);
  
  xlim <- c(0, 1.0);
  
  op <- par(mar=c(4,4,3,6));
  pchs <- ifelse(as.numeric(cls) == 1, 1, 19);
  
  colors <- ifelse(show==1, "darkgrey", "black");
  ROCR::plot(prob.vec, y, pch=pchs, col=colors, xlim=xlim, ylim= ylim, xlab = "Predicted Class Probabilities", ylab="Samples");
  abline(h = 0, lty = 2, col="grey");
  abline(v = 0.5, lty = 2, col="grey");
  
  par(xpd=T);
  legend("right",inset=c(-0.11,0), legend = unique(as.character(cls)), pch=unique(pchs));
  
  test.y <- test.x <- 0;
  if(showPred){
    test.y <- rnorm(length(rdtSet$analSet$multiROC$test.res));
    test.x <- rdtSet$analSet$multiROC$test.res;

    pchs <- ifelse(as.numeric(rdtSet$dataSet$test.cls) == 1, 1, 19);
    points(test.x, test.y, pch=pchs, cex=1.5, col="red");
  }
  
  if(show == 1){ 
    
    # add labels for sample classified wrong
    # the idea is to add labels to the left side for those with prob < 0.5
    # and add labels to the right side of the point for those with prob > 0.5
    # leave 0.5 out 
    
    # first wrong pred as 1 (right side)
    act.ones <- as.numeric(cls)-1 == 1;
    pred.vec <- ifelse(prob.vec > 0.5, 1, 0);
    
    wrong.inx <- (pred.vec != as.numeric(cls)-1) & pred.vec == 1;
    if(sum(wrong.inx) > 0){
      text(prob.vec[wrong.inx], y[wrong.inx], nms[wrong.inx], pos=4);
    }
    
    # first wrong pred as 0 (left side)
    act.zeros <- as.numeric(cls)-1 == 0;
    pred.vec <- ifelse(prob.vec < 0.5, 0, 0.5);
    wrong.inx <- pred.vec != as.numeric(cls)-1 & pred.vec == 0;
    if(sum(wrong.inx) > 0){
      text(prob.vec[wrong.inx], y[wrong.inx], nms[wrong.inx], pos=2);
    }
    
    if(showPred){
      nms <- rownames(rdtSet$dataSet$test.data);
      
      act.ones <- as.numeric(rdtSet$dataSet$test.cls)-1 == 1;
      act.zeros <- as.numeric(rdtSet$dataSet$test.cls)-1 == 0;
      
      # leave 0.5 there
      pred.vec <- ifelse(test.x > 0.5, 1, 0.5);
      wrong.inx <- (pred.vec != as.numeric(rdtSet$dataSet$test.cls)-1) & act.ones;
      if(sum(wrong.inx) > 0){
        text(test.x[wrong.inx], test.y[wrong.inx], nms[wrong.inx], pos=4, cex=0.9);
      }
      
      pred.vec <- ifelse(test.x < 0.5, 0, 0.5);
      wrong.inx <- pred.vec != as.numeric(rdtSet$dataSet$test.cls)-1 & act.zeros;
      if(sum(wrong.inx) > 0){
        text(test.x[wrong.inx], test.y[wrong.inx], nms[wrong.inx], pos=2, cex=0.9);
      }
    }
  }
  par(op)
  dev.off();
  return(.set.rdt.set(rdtSet));
}


PlotImpBiomarkers <- function(imgName, format="png", dpi=72, 
                        mdl.inx, measure = "freq", feat.num = 15) {
  
  rdtSet <- .get.rdt.set();

  anal.mode <- rdtSet$analSet$mode;
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 8; h <- 8;
  rdtSet$imgSet$roc.imp.plot <- imgName;
  rdtSet$imgSet$roc.imp.name <- mdl.inx
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  if(is.null(rdtSet$dataSet$norm_explore)){
    data <- rdtSet$dataSet$norm;
    rdtSet$dataSet$norm_explore <- rdtSet$dataSet$norm;
  } else {
    data <- rdtSet$dataSet$norm_explore
  }
 
  cls <- rdtSet$dataSet$cls;
  op <- par(mar=c(6,10,3,7));
  
  # if(anal.mode == "explore"){
    if(mdl.inx == -1){
      mdl.inx <- rdtSet$analSet$multiROC$best.model.inx;
    }
    imp.mat <- GetImpFeatureMat(rdtSet, rdtSet$analSet$multiROC$imp.cv, rdtSet$analSet$multiROC$test.feats[mdl.inx]);
  # }else{
  #   imp.mat <- GetImpFeatureMat(rdtSet, rdtSet$analSet$multiROC$imp.cv, null);
  # }
  
  imp.nms <- rownames(imp.mat);
  hit.nms <- imp.nms[imp.nms %in% colnames(data)];
  data <- data[, hit.nms];
  
  # note, tapply can only be applied to a vector, we need
  # to combine with apply in order to used on a data frame
  mds <- apply(data, 2,
               function(x){
                 tapply(x, cls, median);
               });

  lowhigh <- apply(mds, 2,
                   function(x){
                     ifelse(rank(x)==1, "Low", "High")
                   });
  lowhigh <- t(lowhigh);

  temp.dat <- data.frame(imp.mat, lowhigh);
  colnames(temp.dat) <- c(colnames(imp.mat), levels(cls));
  imp.fileNm <- "imp_features_cv.csv";
  fast.write.csv(temp.dat, file=imp.fileNm);
  temp.dat <- NULL;
  
  # record the imp.mat for table show
  rdtSet$analSet$imp.mat <- imp.mat;
  rdtSet$analSet$lowhigh <- lowhigh;
  rdtSet$analSet$roc.sig.nm <- imp.fileNm;

  if(measure=="freq"){
    imp.vec <- sort(imp.mat[,1], decreasing=T);
    xlbl <- "Selected Frequency (%)";
  }else{ # default sort by freq, need to reorder
    imp.vec <- sort(imp.mat[,2], decreasing=T);
    xlbl <- "Average Importance";
  }
  var.len <- length(imp.vec);
  
  if(feat.num <= 0){
    feat.num = 15;
  }else if(feat.num > var.len){
    feat.num <- var.len;
  }
  
  imp.vec<-rev(imp.vec[1:feat.num]);
  
  nms.orig <- names(imp.vec);
  vip.nms <-substr(nms.orig, 1, 20);
  names(imp.vec) <- NULL;
  
  xlim.ext <- GetExtendRange(imp.vec, 12);
  dotchart(imp.vec, bg="blue", xlab= xlbl, xlim=xlim.ext);
  
  mtext(side=2, at=1:feat.num, vip.nms, las=2, line=1)
  names(imp.vec) <- nms.orig;
  
  axis.lims <- par("usr"); # x1, x2, y1 ,y2
  
  # get character width
  shift <- 2*par("cxy")[1];
  lgd.x <- axis.lims[2] + shift;
  
  x <- rep(lgd.x, feat.num);
  y <- 1:feat.num;
  
  par(xpd=T);

  # now synchronize lowhigh with imp.vec
  lowhigh <- lowhigh[nms.orig, ,drop=FALSE];

  nc <- ncol(lowhigh);
  col <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(nc);
  
  # calculate background
  bg <- matrix("", nrow(lowhigh), nc);
  for (m in 1:nrow(lowhigh)){
    bg[m,] <- ifelse(lowhigh[m,]=="High", col[1], col[2]);
  } 

  cls.lbl <- levels(cls);

  for (n in 1:ncol(lowhigh)){
    points(x,y, bty="n", pch=22, bg=bg[,n], cex=3);
    # now add label
    text(x[1], axis.lims[4], cls.lbl[n], srt=45, adj=c(0.2,0.5));
    # shift x, note, this is good for current size
    x <- x + shift/1.25;
  }
  
  # now add color key, padding with more intermediate colors for contiuous band
  col <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(20)
  nc <- length(col);
  x <- rep(x[1] + shift, nc);
  
  shifty <- (axis.lims[4]-axis.lims[3])/3;
  starty <- axis.lims[3] + shifty;
  endy <- axis.lims[3] + 2*shifty;
  y <- seq(from = starty, to = endy, length = nc);
  
  points(x,y, bty="n", pch=15, col=rev(col), cex=2);
  
  text(x[1], endy+shifty/8, "High");
  text(x[1], starty-shifty/8, "Low");
  par(op);
  dev.off();
  
  return(.set.rdt.set(rdtSet));
  
}


PlotAccuracy<-function(rdtSet=NA, imgName, format="png", dpi=72){
  
  rdtSet <- .get.rdt.set();
  anal.mode <- rdtSet$analSet$mode;
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 9; h <- 7;
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  if(is.null(rdtSet$analSet$multiROC$accu.mat)){
    accu.mat <- rdt.set$analSet$ROCtest$accu.mat;
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