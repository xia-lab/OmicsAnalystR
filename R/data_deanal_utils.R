##################################################
## R script for OmicsAnalyst
## Description: functions for id annotation
##
## Authors: 
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

DoStatComparison <- function(filenm, alg="ttest", meta, selected, meta.vec, normOpt, p.lvl=0.05, fc.lvl=0.05, nonpar=FALSE){
  if(meta == "null"){
    meta = 1;
  }
  
  dataSet <- readRDS(filenm);
  
  if(normOpt != "none"){
    dataSet$data.comparison <- NormalizingDataOmics(dataSet$data.filtered, normOpt, "NA", "NA");
    
  }else{
    dataSet$data.comparison <- dataSet$data.filtered;
  }
  colnames(dataSet$data.comparison) <- rownames(dataSet$meta)
 
  if(exists("m2m",dataSet)){
    dataSet$data.comparison.taxa <- dataSet$data.filt.taxa
  }
  
  if(alg == "deseq2" || alg == "edger"){
    if(dataSet$isValueNormalized == "false"){
      dataSet$data.comparison <- round(dataSet$data.comparison)
    }
  }
  
  
  if(dataSet$de.method == alg && dataSet$de.norm == normOpt){
    return(UpdateDE(filenm, p.lvl, fc.lvl));
  }
  
  if(selected == "NA"){ # process page
    if(meta == ""){
      meta <- 1;
    }
    metavec <- dataSet$meta[,meta]
    sel <- unique(metavec)
  }else{
    metavec <- dataSet$meta[,meta]
    sel <- strsplit(selected, "; ")[[1]];
  }
  
  dataSet$meta$newcolumn = metavec
  metadf = dataSet$meta
  
  sel_meta1 = metadf[which(metadf[,"newcolumn"] %in% sel[1]),]
  sel_meta2 = metadf[which(metadf[,"newcolumn"] %in% sel[2]),] 
  nms1 <- rownames(sel_meta1)
  nms2 <- rownames(sel_meta2)
  sel_meta_more_than_2 = metadf[which(metadf[,"newcolumn"] %in% sel),]
  nms <- rownames(sel_meta_more_than_2)
  
  if(alg=="ttest"){
    if(length(unique(sel))>2){
      res <- GetFtestRes(dataSet, nms,"newcolumn", F, TRUE, F);
    }else{
      res <- GetTtestRes(dataSet, nms1,nms2,"newcolumn", F, TRUE, F);
    }
  }else if(alg =="kruskal"){
    if(length(unique(sel))>2){
      res <- GetFtestRes(dataSet, nms,"newcolumn", F, TRUE, T);
    }else{
      res <- GetTtestRes(dataSet, nms1,nms2,"newcolumn", F, TRUE, T);
    }
  }else if(alg =="limma"){
    res <- performLimma(dataSet, nms, "newcolumn")
  }else if(alg=="edger"){
    res <- performEdgeR(dataSet, nms, "newcolumn")
  }else if(alg =="deseq2"){
    performDeseq2(dataSet, nms, "newcolumn")
    .perform.computing();
    res <- .save.deseq.res(dataSet)
  }
  
  res <- res[,c(1,2)]
  rownames(res) <- rownames(dataSet$data.comparison)
  colnames(res) <- c("stat", "p_value")
  
  res <- na.omit(res)
  res <- res[order(res[,2], decreasing=FALSE),]
  pvals <- p.adjust(res[,"p_value"],method="BH");
  res <- cbind(res, pvals)
  res <- cbind(res, rownames(res))
  
  res[res == "NaN"] = 1
  pv <- as.numeric(res[,"p_value"])
  pv_no_zero <- pv[pv != 0]
  minval <- min(pv_no_zero)
  pv[pv == 0] = minval/2
  pvals <- -log10(pv);
  colorb<- ComputeColorGradient(pvals, "black", F, F);
  sizes <- as.numeric(rescale2NewRange(-log10(pv), 15, 35));
  res <- cbind(res, colorb);
  res <- cbind(res, sizes);
  ids <- names(dataSet$enrich_ids[rownames(res)])
  res <- cbind(res, ids);
  colnames(res) = c("stat", "p_value", "p_adj", "ids", "color", "size", "name");
  res <- as.matrix(res)
  de <- res
  
  
  
  if(exists("m2m",dataSet)){
    
    res2 <- GetTaxaCompRes(dataSet,nms1,nms2, nms, sel,"newcolumn", paired=FALSE,
                           equal.var=TRUE, nonpar=F,alg)
    res2 <- lapply(res2, function(x) x[,c(1,2)])
    for(i in 1:length(res2)){
      rownames(res2[[i]]) <- rownames(dataSet$data.comparison.taxa[[i]])
      
      colnames(res2[[i]]) <- c("stat", "p_value")
      res2[[i]] <- na.omit(res2[[i]])
      res2[[i]] <- res2[[i]][order(res2[[i]][,2], decreasing=FALSE),]
      pvals <- p.adjust(res2[[i]][,"p_value"],method="BH");
      res2[[i]] <- cbind(res2[[i]], pvals)
      res2[[i]] <- cbind(res2[[i]], rownames(res2[[i]]))
      
      
      res2[[i]][res2[[i]] == "NaN"] = 1
      
      pv <- as.numeric(res2[[i]][,"p_value"])
      pv_no_zero <- pv[pv != 0]
      minval <- min(pv_no_zero)
      pv[pv == 0] = minval/2
      pvals <- -log10(pv);
      colorb<- ComputeColorGradient(pvals, "black", F, F);
      sizes <- as.numeric(rescale2NewRange(-log10(pv), 15, 35));
      res2[[i]] <- cbind(  res2[[i]], colorb);
      res2[[i]] <- cbind(  res2[[i]], sizes);
      ids <- rownames(res2[[i]])
      res2[[i]] <- cbind(  res2[[i]], ids);
      colnames(  res2[[i]]) = c("stat", "p_value", "p_adj", "ids", "color", "size", "name");
      res2[[i]] <- as.matrix(res2[[i]])
    }
    dataSet$comp.res.taxa = res2;
  }
  dataSet$de.norm = normOpt;
  dataSet$de.method = alg;
  dataSet$comp.res = de;
  
  RegisterData(dataSet);
  
  return(UpdateDE(filenm, p.lvl, fc.lvl))
}


UpdateDE<-function(dataName, p.lvl = 0.05, fc.lvl = 1){
  dataSet <- readRDS(dataName);
  res <- dataSet$comp.res
  hit.inx <- as.numeric(res[, "p_value"]) <= p.lvl #pval
  
  if(sum(hit.inx) == 0){
    return (c(1, 0, nrow(res)));
  }
  # note, hit.inx can contain NA, not T/F
  hit.inx <- which(hit.inx);
  
  res.sig<-res[hit.inx, , drop=F];
  
  hit.inx <- abs(as.numeric(res.sig[, "stat"])) > fc.lvl #foldchange
  
  if(sum(hit.inx) == 0){
    return (c(1, 0, nrow(res)));
  }
  # note, hit.inx can contain NA, not T/F
  hit.inx <- which(hit.inx);
  
  res.sig<-res.sig[hit.inx, , drop=F];
  
  sig.count <- nrow(res.sig);
  de.genes <- rownames(res.sig);
  non.sig.count <- nrow(res)-sig.count;
  
  dataSet$sig.mat = res.sig;
  
  
  if(exists("m2m",dataSet)){
    
    res2 <- dataSet$comp.res.taxa
    
    hit.inx <- lapply(res2,function(x) as.numeric(x[, "p_value"]) <= p.lvl) #pval
    
    # note, hit.inx can contain NA, not T/F
    hit.inx <- lapply(hit.inx, function(x) which(x) );
    res.sig2 <- list()
    for(i in 1:length(res2)){
      if(length(hit.inx[i])==0){
        res.sig2[[i]] <- c(1, 0, nrow(res2[i]))
        
      }else{
        
        res.sig2[[i]] <- res2[[i]][hit.inx[[i]], , drop=F]
      }
      
    }
    
    
    hit.inx <- lapply(res.sig2,function(x) abs(as.numeric(x[, "stat"])) > fc.lvl) #foldchange
    hit.inx <- lapply(hit.inx, function(x) which(x) );
    for(i in 1:length(res.sig2)){
      if(length(hit.inx[i])==0){
        res.sig2[[i]] <- c(1, 0, nrow(res2[i]));
      }else{
        res.sig2[[i]] <- res.sig2[[i]][hit.inx[[i]], , drop=F];
      }
    }
    
    dataSet$sig.count.tax <- lapply(res.sig2, function(x)nrow(x) );
    # dataSet$de.genes <- lapply(res.sig2, function(x) rownames(x) );
    non.sig.count2 <- list()
    for(i in 1:length(res2)){
      non.sig.count2[[i]] <- nrow(res2[[i]])-dataSet$sig.count.tax[[i]]; 
    }
    dataSet$non.sig.count.tax <- non.sig.count2;
    names(res.sig2)<- colnames(dataSet$taxa_table);
    dataSet$sig.mat.tax = res.sig2;
  }
  
  RegisterData(dataSet);
  
  return(c(1, sig.count, non.sig.count));
}