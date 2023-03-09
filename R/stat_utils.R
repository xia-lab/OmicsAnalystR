##################################################
## R scripts for OmicsAnalyst
## Statistical comparison
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################


DoStatComparison <- function(filenm, alg="limma", meta, selected, meta.vec, normOpt, p.lvl=0.05, fc.lvl=0, nonpar=FALSE){
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
  
  # NOTE JE: rethink when these options are exposed: showing correctly for example data (ie. certain methods for counts vs. metabolomics)
  # but not for new, protocol data

  trimmed.data <- as.matrix(dataSet$data.comparison[,which(colnames(dataSet$data.comparison) %in% nms)])
  trimmed.meta.df <- dataSet$meta[which(rownames(dataSet$meta) %in% nms), ]
  trimmed.meta <- dataSet$meta[,"newcolumn"][which(rownames(dataSet$meta) %in% nms)]
  trimmed.meta <- make.names(trimmed.meta);
  if(min(trimmed.data) < 0){
    trimmed.data = trimmed.data + abs(min(trimmed.data))
  }

  res <- performLimma(trimmed.data, trimmed.meta);
  
  # NOTE JE: we need FC values here
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

# NOTE JE: unclear if FC values exist, since we didn't select them in initial DEA function
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
      res.sig2[[i]] <- c(1, 0, nrow(res2[i]))
      
    }else{
      
      res.sig2[[i]] <- res.sig2[[i]][hit.inx[[i]], , drop=F]
    }
    
  }
  
  
  dataSet$sig.count.tax <- lapply(res.sig2, function(x)nrow(x) );
 # dataSet$de.genes <- lapply(res.sig2, function(x) rownames(x) );
  non.sig.count2 <- list()
  for(i in 1:length(res2)){
    non.sig.count2[[i]] <- nrow(res2[[i]])-dataSet$sig.count.tax[[i]]; 
  }
  dataSet$non.sig.count.tax <- non.sig.count2
  
  names(res.sig2)<- colnames(dataSet$taxa_table)
  dataSet$sig.mat.tax = res.sig2;
}

  RegisterData(dataSet);

  return(c(1, sig.count, non.sig.count))
}

DoStatComparisonVis <- function(filenm, alg, meta, selected, meta.vec, omicstype, nonpar=FALSE){

  if(meta == "null"){
    meta = 1;
  }
  if(meta.vec == "NA"){ # process page
    #if(dataSet$name != filenm){
    dataSet <- readRDS(filenm);
    #}
  }else{
    if(omicstype != "NA"){
      sel.nms <- names(mdata.all)
      for(i in 1:length(sel.nms)){
        dat = readRDS(sel.nms[i])
        if(dat$type == omicstype){
          dataSet = dat;
        }
      }
    }else{
      dataSet <- .get.rdt.set();
    }
  }
  
  if(meta.vec == "NA"){ # process page
    if(meta == ""){
      meta=1;
    }
    metavec = dataSet$meta[,meta]
    sel = unique(metavec)
  }else{
    metavec <- strsplit(meta.vec, "; ")[[1]];
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
  }else{ 

    trimmed.data = as.matrix(dataSet$data.comparison[,which(colnames(dataSet$data.comparison) %in% nms)])
    trimmed.meta.df <- dataSet$meta[which(rownames(dataSet$meta) %in% nms), ]
    trimmed.meta = dataSet$meta[,"newcolumn"][which(rownames(dataSet$meta) %in% nms)]
    trimmed.meta <- make.names(trimmed.meta);
    if(min(trimmed.data) < 0){
      trimmed.data = trimmed.data + abs(min(trimmed.data))
    }

    res <- performLimma(trimmed.data, trimmed.meta);

}
  res = res[,c(1,2)]
  rownames(res) = rownames(dataSet$data.comparison)
  colnames(res) = c("stat", "p_value")
  
  res = na.omit(res)
  res = res[order(res[,2], decreasing=FALSE),]
  pvals <- p.adjust(res[,"p_value"],method="BH");
  res = cbind(res, pvals)
  res = cbind(res, rownames(res))
  de = res
  de[de == "NaN"] = 1
  pv = as.numeric(de[,"p_value"])
  pv_no_zero = pv[pv != 0]
  minval = min(pv_no_zero)
  pv[pv == 0] = minval/2
  pvals <- -log10(pv);
  colorb<- ComputeColorGradient(pvals, "black", F, F);
  sizes <- as.numeric(rescale2NewRange(-log10(pv), 15, 35));
  res = cbind(res, colorb);
  res = cbind(res, sizes);

  #ids <- names(dataSet$enrich_ids[order(match(dataSet$enrich_ids,rownames(res)))])
  enrich_ids <- dataSet$enrich_ids[dataSet$enrich_ids %in% rownames(res)]
  ids <- names(enrich_ids[order(match(enrich_ids,rownames(res)))])
  
  res = cbind(res, ids);
  colnames(res) = c("stat", "p_value", "p_adj", "ids", "color", "size", "name");
  res= as.matrix(res)
  library(RJSONIO)
  sink(filenm);
  cat(toJSON(res));
  sink();
  
  if(meta.vec == "NA"){
    filenm = "OK";
    dataSet$comp.res = de;
    dataSet$sel.meta = meta
    RegisterData(dataSet);
  }
  
  return(filenm)
}

GetTtestRes <- function(dataSet, nms1,nms2, meta, paired=FALSE, equal.var=TRUE, nonpar=F){
  trimmed.data = as.matrix(dataSet$data.comparison[,which(colnames(dataSet$data.comparison) %in% c(nms1,nms2))])
  dataSet$data.comparison = as.matrix(dataSet$data.comparison)
  inx1 = which(colnames(dataSet$data.comparison) %in% nms1)
  inx2 = which(colnames(dataSet$data.comparison) %in% nms2)
  if(nonpar){
    res <- apply(trimmed.data, 1, function(x) {
      tmp <- try(wilcox.test(x[inx1], x[inx2], paired = F));
      if(class(tmp) == "try-error") {
        return(c(NA, NA));
      }else{
        return(c(tmp$statistic, tmp$p.value));
      }
    })
    res = t(res)
  }else{#repalce with PerformFastUniv
    metad = dataSet$meta
    metad =metad[which(rownames(metad) %in% c(nms1,nms2)),]
    res <- PerformFastUnivTests(trimmed.data, factor(metad[,"newcolumn"]),T, F)

  }

  return(res);
}


GetFtestRes <- function(dat,nms, meta,paired=FALSE, equal.var=TRUE, nonpar=F){
  trimmed.data = as.matrix(dat$data.comparison[,which(colnames(dat$data.comparison) %in% nms)])
  dat$data.comparison = as.matrix(dat$data.comparison)
  metad = dat$meta
  
  metad =metad[which(rownames(metad) %in% nms),]
  if(nonpar){
    kwtest <- function(x, cls) {kruskal.test(x ~ cls);}
    aov.res <- apply(trimmed.data, 1, kwtest, cls=factor(metad[,"newcolumn"]));
    #extract all p values
    res <- unlist(lapply(aov.res, function(x) {c(x$statistic, x$p.value)}));
    res <- data.frame(matrix(res, nrow=length(aov.res), byrow=T), stringsAsFactors=FALSE);
  }else{
    res <- PerformFastUnivTests(trimmed.data, factor(metad[,"newcolumn"]),TRUE, F);
  }
  return(res);
}

#####to get result for each taxonomy level

GetTaxaCompRes <- function(dataSet,nms1,nms2, nms,sel, meta="newcolumn", paired=FALSE,
                           equal.var=TRUE, nonpar=F,alg){
  
  if(length(sel)==2){

    res <- lapply(dataSet$data.comparison.taxa, function(x){
      trimmed.data = as.matrix(x[,which(colnames(x) %in% c(nms1,nms2))])
      x = as.matrix(x)
      inx1 = which(colnames(x) %in% nms1)
      inx2 = which(colnames(x) %in% nms2)
      trimmed.meta.df = dataSet$meta;
      trimmed.meta.df =trimmed.meta.df[which(rownames(trimmed.meta.df) %in% c(nms1,nms2)),];
      trimmed.meta <- factor(trimmed.meta.df[,"newcolumn"]);
      res <- performLimma(trimmed.data, trimmed.meta);
    });
    return(res)
  }
}

#'ANOVA
#'@description Perform anova and only return p values and MSres (for Fisher's LSD)
#'@param x Input the data to perform ANOVA
#'@param cls Input class labels
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
aof <- function(x, cls) {
  aov(x ~ cls);
}

gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

performLimma <-function(trimmed.data, trimmed.meta){
  library(edgeR);
  cls <- as.factor(trimmed.meta); 
  inx = 0;
  myargs <- list();
  grp.nms <- levels(cls);
  
  for(m in 1:(length(grp.nms)-1)){
    for(n in (m+1):length(grp.nms)){
      inx <- inx + 1;
      myargs[[inx]] <- paste(grp.nms[m], "-", grp.nms[n], sep="")
    }
  }
  
  design <- model.matrix(~ 0 + cls) # no intercept
  colnames(design) <- levels(cls);
  myargs[["levels"]] <- design
  contrast.matrix <- do.call(makeContrasts, myargs)
  dataSet$design <- design;
  y=DGEList(counts = trimmed.data, group = cls)
  y <- calcNormFactors(y)
  y <- estimateCommonDisp(y, verbose = FALSE)
  y <- estimateTagwiseDisp(y)
  v <- voom(y, design, plot = F)
  vfit <- lmFit(v, design)
  vfit <- contrasts.fit(vfit, contrasts=contrast.matrix)
  efit <- eBayes(vfit)
  topFeatures <- topTable(efit, number = Inf, adjust.method = "fdr");
  if(length(unique(cls)) == 2){
    res = data.frame(stat=topFeatures[,"t"], P.Value=topFeatures[,"P.Value"], adj.P.Value=topFeatures[,"adj.P.Val"])
  }else{
    res = data.frame(stat=topFeatures[,"F"], P.Value=topFeatures[,"P.Value"], adj.P.Value=topFeatures[,"adj.P.Val"])
  }
  return(res)
}

PerformClusteringScatter <- function(filenm, type, nclust){
  nclust = as.numeric(nclust)
  Sys.setenv(RGL_USE_NULL = TRUE)
  library(rgl)
  library(cluster)
  reductionSet <- .get.rdt.set();
  pos.xyz = reductionSet$pos.xyz;
  if(type == "density"){
    library(ADPclust)
    if(init == "true"){
      ans = adpclust(pos.xyz)
    }else{
      clusterNum = as.numeric(clusterNum)
      ans = adpclust(pos.xyz, nclust=clusterNum)
    }
    cluster = ans$clusters
    dataSet$clusterObject = ans
  }else if (type == "kmeans"){
    ans = kmeans(pos.xyz, nclust)
    cluster = ans$cluster
  }else{
    distMat = dist(pos.xyz)
    hc = hclust(distMat, "complete")
    cluster <- cutree(hc, k = nclust)
  }
  dataSet$meta$cluster=unname(cluster)
  RegisterData(dataSet)
  
  netData <- list(cluster=unname(cluster), objects=a$objects, ellipse=meshes, meta=dataSet$meta);
  library(RJSONIO)
  sink(filenm);
  cat(toJSON(netData));
  sink();
  return(filenm)
}

PerformClusteringMeta <-function(filenm, meta, type, opt){ 
  Sys.setenv(RGL_USE_NULL = TRUE)
  library(rgl)
  library(cluster)
  reductionSet <- .get.rdt.set();
  pos.xyz = reductionSet$pos.xyz
  if(reductionOptGlobal  %in% c(loadingOpts, "procrustes")){
    metadf = reductionSet$newmeta
  }else{
    metadf = reductionSet$meta
  }
  selMetas = unique(metadf[,meta])
  clustersHolder = list()
  library(ADPclust)
  library(factoextra)
  if(opt == "meta"){
    for(i in 1:length(unique(selMetas))){
      inx = metadf[,meta] == selMetas[i]
      clusList = list()
      kmax=10
      r=fviz_nbclust(pos.xyz[inx,], kmeans, method = "silhouette", k.max = kmax)
      maxScore = max(r$data$y)
      inx2 = r$data$y == maxScore
      optimalNumber = as.numeric(r$data$clusters[inx2][1])
      if(nrow(pos.xyz[inx,])<optimalNumber){
        optimalNumber = nrow(pos.xyz[inx,])
      }
      if(type == "kmeans"){
        ans = kmeans(pos.xyz[inx,], optimalNumber)
        clusList$cluster = unname(ans$cluster)
      }else if(type == "peak"){
        ans = adpclust(pos.xyz[inx,], nclust=optimalNumber)
        clusList$cluster = ans$clusters
      }else if(type == "meanshift"){
        library(ks)
        ans = kms(pos.xyz[inx,])
        clusList$cluster = ans$label
      }
      clusList$inxs = which(metadf[,meta] == selMetas[i])
      clustersHolder[[i]]=clusList
    }
  }else{
    clusList = list()
    kmax=10
    
    if(nrow(pos.xyz)<10){
      kmax = nrow(pos.xyz)
    }     
    r=fviz_nbclust(pos.xyz, kmeans, method = "silhouette", k.max = kmax)
    maxScore = max(r$data$y);
    inx2 = r$data$y == maxScore;
    optimalNumber = as.numeric(r$data$clusters[inx2][1])
    if(nrow(pos.xyz)<optimalNumber){
      optimalNumber = nrow(pos.xyz);
    }
    if(type == "kmeans"){
      ans = kmeans(pos.xyz, optimalNumber);
      clusList$cluster = unname(ans$cluster);
    }else if(type == "peak"){
      ans = adpclust(pos.xyz, nclust=optimalNumber);
      clusList$cluster = ans$clusters;
    }else if(type == "meanshift"){
      library(ks)
      ans = kms(pos.xyz, min.clust.size=10);
      clusList$cluster = ans$label;
    }
    clusList$inxs = seq.int(1, length(clusList$cluster));
    clustersHolder[[1]]=clusList;
    cluster <- clusList$cluster;
  }
  
  if(opt == "meta"){
    netData <- list(cluster=clustersHolder, metaNm=selMetas, objects="NA", ellipse="NA", meta=metadf);
  }else{
    netData <- list(cluster=clustersHolder, metaNm=paste0("Cluster ", unique(cluster)), objects="NA", ellipse="NA", meta=metadf);  
  }
  
  library(RJSONIO);
  sink(filenm);
  cat(toJSON(netData));
  sink();
  return(filenm);
}


PerformCustomClustering <- function(filenm, type, ids){
  Sys.setenv(RGL_USE_NULL = TRUE)
  library(rgl)
  library(cluster)
  library(ADPclust)
  library(factoextra)
  
  idsvec <- strsplit(ids, "; ")[[1]];
  reductionSet <- .get.rdt.set();
  pos.xyz = reductionSet$pos.xyz
  if(reductionOptGlobal  %in% loadingOpts){
    metadf = reductionSet$newmeta
  }else{
    metadf = reductionSet$meta
  }
  
  clustersHolder = list()
  
  inx = which(rownames(pos.xyz) %in% idsvec)
  clusList = list()
  kmax=10
  r=fviz_nbclust(pos.xyz[inx,], kmeans, method = "silhouette", k.max = kmax)
  maxScore = max(r$data$y)
  inx2 = r$data$y == maxScore
  optimalNumber = as.numeric(r$data$clusters[inx2][1])
  if(nrow(pos.xyz[inx,])<optimalNumber){
    optimalNumber = nrow(pos.xyz[inx,])
  }
  if(type == "kmeans"){
    ans = kmeans(pos.xyz[inx,], optimalNumber)
    clusList$cluster = unname(ans$cluster)
  }else if(type == "peak"){
    ans = adpclust(pos.xyz[inx,], nclust=optimalNumber)
    clusList$cluster = ans$clusters
  }else if(type == "meanshift"){
    library(ks)
    ans = kms(pos.xyz[inx,])
    clusList$cluster = ans$label
  }
  clusList$inxs = inx
  
  netData <- list(cluster=clusList);  
  library(RJSONIO)
  sink(filenm);
  cat(toJSON(netData));
  sink();
  return(filenm)
}
