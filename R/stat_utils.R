##################################################
## R scripts for OmicsAnalyst
## Statistical comparison
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

UpdateDE <- function(dataName, fc.lvl, p.lvl = 0.05){
  save.image("updateDE.RData")
  dataSet <- readDataset(dataName)
  res     <- dataSet$comp.res

  if (is.null(res)){
    current.msg <<- "There is no result generated. Have you performed analysis?"
    return(c(0, 0, 0))
  }

  ## ---------- figure out per-gene fold-change -------------------------
  if ("coefficient" %in% colnames(res)){              # old layout
    fc.vec <- abs(as.numeric(res$coefficient))
  } else {                                            # new layout
    ave.idx <- match("AveExpr", colnames(res))
    if (is.na(ave.idx) || ave.idx <= 1){
      stop("Cannot locate comparison columns (before AveExpr).")
    }
    comp.mat <- as.matrix(res[ , 1:(ave.idx-1), drop = FALSE])
    mode(comp.mat) <- "numeric"
    fc.vec  <- apply(comp.mat, 1, function(x) max(abs(x), na.rm = TRUE))
  }

  ## ---------- hit index (p-value & FC) --------------------------------
  hit.inx <- (as.numeric(res[ , "P.Value"]) <= p.lvl) & (fc.vec >= fc.lvl)

  if (sum(hit.inx, na.rm = TRUE) == 0){
    return(c(1, 0, nrow(res)))
  }
  hit.rows   <- which(hit.inx)
  res.sig    <- res[hit.rows, , drop = FALSE]

  sig.count      <- nrow(res.sig)
  non.sig.count  <- nrow(res) - sig.count
  dataSet$sig.mat <- res.sig

  ## ---------- taxa branch (unchanged except kept original p_adj filter) ----
  if (exists("m2m", dataSet)){
    res2      <- dataSet$comp.res.taxa
    hit.inx2  <- lapply(res2, function(x) as.numeric(x[ , "p_adj"]) <= p.lvl)
    hit.inx2  <- lapply(hit.inx2, function(x) which(x))
    res.sig2  <- lapply(seq_along(res2), function(i){
                        if (length(hit.inx2[[i]]) == 0) return(res2[[i]][0 , , drop = FALSE])
                        res2[[i]][hit.inx2[[i]], , drop = FALSE]
                      })
    names(res.sig2)           <- colnames(dataSet$taxa_table)
    dataSet$sig.mat.tax       <- res.sig2
    dataSet$sig.count.tax     <- lapply(res.sig2, nrow)
    dataSet$non.sig.count.tax <- Map(function(r,s) nrow(r)-s, res2, dataSet$sig.count.tax)
  }

  RegisterData(dataSet)
  return(c(1, sig.count, non.sig.count))
}


#####to get result for each taxonomy level

GetTaxaCompRes <- function(dataSet,nms1,nms2, nms,sel, meta="newcolumn", paired=FALSE,
                           equal.var=TRUE, nonpar=F,alg){
  
  if(length(sel)==2){

    res <- lapply(dataSet$data.proc.taxa, function(x){
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

performLimma <- function(trimmed.data, trimmed.meta){
  
  library(limma);
  
  # process class labels
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
  
  # create design matrix
  design <- model.matrix(~ 0 + cls) # no intercept
  colnames(design) <- levels(cls);
  myargs[["levels"]] <- design
  contrast.matrix <- do.call(makeContrasts, myargs)
  dataSet$design <- design;
  
  # perform differential analysis
  fit <- lmFit(trimmed.data, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  topFeatures <- topTable(fit2, number = Inf, adjust.method = "fdr");
  
  # process results
  if(length(unique(cls)) == 2){
    res = data.frame(stat=topFeatures[,"t"], P.Value=topFeatures[,"P.Value"], adj.P.Value=topFeatures[,"adj.P.Val"])
  }else{
    res = data.frame(stat=topFeatures[,"F"], P.Value=topFeatures[,"P.Value"], adj.P.Value=topFeatures[,"adj.P.Val"])
  }
  
  rownames(res) <- rownames(topFeatures)
  
  return(res)
}

PerformClusteringScatter <- function(filenm, type, nclust){
  nclust = as.numeric(nclust)
  Sys.setenv(RGL_USE_NULL = TRUE)
  library(rgl)
  library(cluster)
  reductionSet <- .get.rdt.set();
  pos.xyz = reductionSet[[reductionSet$reductionOpt ]]$pos.xyz;
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
  pos.xyz = reductionSet[[reductionSet$reductionOpt ]]$pos.xyz;
    metadf = reductionSet$newmeta

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
  pos.xyz = reductionSet[[reductionSet$reductionOpt]]$pos.xyz;
  metadf = reductionSet$newmeta
  
  clustersHolder = list()
  
  inx = which(rownames(pos.xyz) %in% idsvec)
  clusList = list()
  kmax = 10
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
