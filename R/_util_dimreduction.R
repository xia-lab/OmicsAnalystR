##################################################
## R script for OmicsAnalyst
## Description: Compute dimension reduction
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################

my.reduce.dimension <- function(omicsType, reductionOpt, method="globalscore", dimn){
  if(method == ""){
    method="globalscore"
  }
  
  dimn = as.numeric(dimn);
  if(dimn == 2){
    dimn = 3;
  }
  sel.nms <- names(mdata.all)[mdata.all==1];
  if(length(sel.nms) != 2 && !reductionOpt %in% c("mcia", "diablo", "mbpca")){
    current.msg <<- "Exactly two datasets are required for the selected method";
    return(0)
  }else{
    data.list = list()
    d.list = list()
    omics.type = vector();
    featureNms <- vector();
    for(i in 1:length(sel.nms)){
      dataSet = readRDS(sel.nms[i])
      omics.type <- c(omics.type, dataSet$type)

      if(dataSet$type %in% c("prot", "rna_b", "mirna")){
      #ids <- dataSet$comp.res[, "ids"];
      #name <- dataSet$comp.res[,"name"];
      #dataSet$comp.res[, "ids"] <- unname(name);
      #dataSet$comp.res[, "name"] <- unname(ids);
      }

      d.list[[dataSet$type]] = list()
      data.list[[dataSet$type]] <- dataSet$data.proc
      d.list[[dataSet$type]][["data.proc"]] = dataSet$data.proc#dat1
      d.list[[dataSet$type]][["meta"]] = dataSet$meta#meta1
      d.list[[dataSet$type]][["comp.res"]] = dataSet$comp.res#comp.res
      d.list[[dataSet$type]][["enrich.nms"]] = dataSet$enrich_ids
      if(i == 1){
        
        newmeta <- dataSet$meta
        comp.res1 = dataSet$comp.res
        enrich.nms1 = dataSet$enrich_ids
        comp.res.inx1 = rep(1, nrow(comp.res1));
        featureNms <- rownames(dataSet$data.proc);
      }else{
        newmeta <- rbind(newmeta, dataSet$meta)
        comp.res1 = rbind(comp.res1, dataSet$comp.res)
        enrich.nms1 = c(enrich.nms1, dataSet$enrich_ids);
        comp.res.inx1 = c(comp.res.inx1, rep(i, nrow(dataSet$comp.res)));
        featureNms <- c(featureNms, rownames(dataSet$data.proc));
      }
      
    }
    dataSet$comp.res = comp.res1
    dataSet$enrich_ids = enrich.nms1
    dataSet$comp.res.inx = comp.res.inx1
  }

  if(reductionOpt == "mbpca"){
    library(mogsa)
    moa <- mbpca(data.list, ncomp = dimn, k = "all", method = "blockScore", option = "inertia", 
                 center=TRUE, scale=TRUE, 
                 unit.p = T, unit.obs=T);
    pos.xyz <- moaScore(moa)[,c(1:3)];
    loading.pos.xyz <- moa@loading
    for(i in 1:length(names(data.list))){
      rownames(loading.pos.xyz) <- gsub(paste0("_",names(data.list)[i]), "", rownames(loading.pos.xyz))
    }
    loadingNames=rownames(loading.pos.xyz)
    names = rownames(pos.xyz)
    newmeta = dataSet$meta
    rownames(dataSet$comp.res) = loadingNames   
    dataSet$misc$pct = signif(moa@eig,4)*100
    dataSet$dim.res <- moa
  }else if(reductionOpt == "diablo"){
    library(mixOmics)
    Y <- dataSet$meta[,1]
    dats = lapply(data.list, function(x){
      x <- data.frame(x, stringsAsFactors = FALSE);
      x <- t(x);
    })
    design = matrix(1, ncol = length(dats), nrow = length(dats), 
                    dimnames = list(names(dats), names(dats)))
    
    res = mixOmics:::block.splsda(X = dats, Y = Y, ncomp =dimn, design = design)
    pos.xyz <- res$variates[[1]]
    pos.xyz2 <- res$variates[[2]]
    
    for(i in 1:length(res$loadings)){
      pos = as.data.frame(res$loadings[[i]])
      rn <- rownames(res$loadings[[i]])
      pos <- unitAutoScale(pos);
      res$loadings[[i]] <- pos
      rownames(res$loadings[[i]]) <- rn
    }
    
    loading.pos.xyz <- rbind(res$loadings[[1]], res$loadings[[2]])
    rownames(loading.pos.xyz) = c(rownames(data.list[[1]]), rownames(data.list[[2]]))
    loadingNames=rownames(loading.pos.xyz)
    
    names = rownames(pos.xyz)
    newmeta = dataSet$meta
    rownames(dataSet$comp.res) = dataSet$comp.res[, ncol(dataSet$comp.res)]
    dataSet$dim.res <- res
    
    for(i in 1:length(names(data.list))){
    if("prop_expl_var" %in% names(res)){
      var.vec <- res$prop_expl_var
    }else if("explained_variance" %in% names(res)){
      var.vec <- res$explained_variance
    }else{
      var.vec <- 0;
    }
    dataSet$misc$pct2[[names(data.list)[i]]] = unname(signif(as.numeric(var.vec[[i]],4)))*100
    if(i==1){
    loading.pos.xyz <- res$loadings[[i]]
    loadingNames=rownames(data.list[[i]])
    }else{
    loading.pos.xyz <- rbind(loading.pos.xyz, res$loadings[[i]])
    loadingNames=c(loadingNames,rownames(data.list[[i]]))
    }
    }
    rownames(loading.pos.xyz) <- loadingNames;

    if(! "diablo" %in% dim.res.methods){
      dim.res.methods <<- c(dim.res.methods, "diablo")
    }
  }else if(reductionOpt == "rcca"){
    library(mixOmics)
    dats = lapply(data.list, function(x){
      x <- data.frame(x, stringsAsFactors = FALSE);
      x <- t(x);
    })
    if(nrow(dats[[1]])>nrow(dats[[2]])){
      X=dats[[2]];
      Y=dats[[1]]
    }else{
      X=dats[[1]];
      Y=dats[[2]]
    }
    if(method == "cv"){
      grid1 <- seq(0, 0.2, length = 5) 
      grid2 <- seq(0.0001, 0.2, length = 5)
      cv <- tune.rcc(X, Y, grid1 = grid1, grid2 = grid2, validation = "loo")
      res <- mixOmics:::rcc(X, Y, ncomp = dimn, lambda1 = cv$opt.lambda1, lambda2 = cv$opt.lambda2)
    }else{
      res <- mixOmics:::rcc(X, Y, method="shrinkage", ncomp = dimn)
    }
    pos.xyz <- res$variates
    loading.pos.xyz <- rbind(res$loadings[[1]], res$loadings[[2]])
    rownames(loading.pos.xyz) = c(rownames(data.list[[1]]), rownames(data.list[[2]]))
    loadingNames=rownames(loading.pos.xyz)
    
    names = rownames(pos.xyz)
    newmeta = dataSet$meta
    rownames(dataSet$comp.res) = dataSet$comp.res[, ncol(dataSet$comp.res)]
    dataSet$dim.res <- res
    if(! "rcca" %in% dim.res.methods){
      dim.res.methods <<- c(dim.res.methods, "rcca")
    }
  }else if(reductionOpt == "spls"){
    library(mixOmics)
    dats = lapply(data.list, function(x){
      x <- data.frame(x, stringsAsFactors = FALSE);
      x <- t(x);
    })
    res <- mixOmics:::spls(dats[[1]], dats[[2]], ncomp = 10, mode = "regression")
    
    for(i in 1:length(res$variates)){
      pos = as.data.frame(res$variates[[i]])
      rn <- rownames(res$variates[[i]])
      pos <- unitAutoScale(pos);
      res$variates[[i]] <- pos
      rownames(res$variates[[i]]) <- paste0(rn, "_", names(data.list)[i])
    }
    
    pos.xyz <- res$variates[[1]]
    pos.xyz2 <- res$variates[[2]]
    
    for(i in 1:length(res$loadings)){
      pos = as.data.frame(res$loadings[[i]])
      rn <- rownames(res$loadings[[i]])
      pos <- unitAutoScale(pos);
      res$loadings[[i]] <- pos
      rownames(res$loadings[[i]]) <- rn
    }
    
    loading.pos.xyz <- rbind(res$loadings[[1]], res$loadings[[2]])
    rownames(loading.pos.xyz) = c(rownames(data.list[[1]]), rownames(data.list[[2]]))
    loadingNames=rownames(loading.pos.xyz)
    names = rownames(pos.xyz)
    newmeta = dataSet$meta
    
    rownames(dataSet$comp.res) = dataSet$comp.res[, ncol(dataSet$comp.res)]
    dataSet$dim.res <- res

    if("prop_expl_var" %in% names(res)){
      var.vec <- res$prop_expl_var
    }else if("explained_variance" %in% names(res)){
      var.vec <- res$explained_variance
    }else{
      var.vec <- 0;
    }

    dataSet$misc$pct2[[names(data.list)[1]]] = unname(signif(as.numeric(var.vec$X,4)))*100
    dataSet$misc$pct2[[names(data.list)[2]]] = unname(signif(as.numeric(var.vec$Y,4)))*100
    if(! "spls" %in% dim.res.methods){
      dim.res.methods <<- c(dim.res.methods, "rcca")
    }
  }else if(reductionOpt == "mcia"){
    library(omicade4)
    # check dims sapply(dfs, dim)
    
    mcoin <- mcia(data.list, cia.nf=3)
    pos.xyz = mcoin$mcoa$Tl1
    loading.pos.xyz = mcoin$mcoa$axis #not really loading, score plot of second omics
    rownames(loading.pos.xyz) = featureNms
    loadingNames = rownames(loading.pos.xyz)
    names = rownames(pos.xyz)
    ndata <- length(mcoin$coa)
    syn <- mcoin$mcoa$SynVar
    sync <- c()
    for (i in 1:(ndata)) {
      sync <- rbind(sync, syn)
    }
    sync <- sync[,c(1:3)]
    colnames(sync) = c("Axis1", "Axis2", "Axis3");
    
    seg.names = rownames(sync)
    sync2 = rbind(sync, pos.xyz)
    maxVal <- max(abs(sync2))
    sync2 <- sync2/maxVal
    sync= sync2[c(1:length(seg.names)),]
    rownames(sync)=seg.names;
    dataSet$mcia.seg.points = sync;
    dataSet$misc$Correlation = mcoin$mcoa$RV[1,2]
    dataSet$misc$pct = signif(mcoin$mcoa$pseudoeig,4)*100
    pos.xyz = sync2[c((length(seg.names)+1):nrow(sync2)),]
    rownames(pos.xyz) = names
    
    type.vec<-vector()
    for(i in 1:length(d.list)){
      if(i==1){
        type.vec <- rep(names(d.list)[i], ncol(d.list[[i]][["data.proc"]]))
      }else{
        type.vec <- c(type.vec, rep(names(d.list)[i], ncol(d.list[[i]][["data.proc"]]))  )      
      }
    }
    newmeta$omics <- type.vec;
    dataSet$dim.res <- mcoin
  }else if(reductionOpt == "procrustes"){
    library(vegan)
    ndat1 <- decostand(t(data.list[[1]]), method = "standardize")
    pca.dat1 <- rda(ndat1)
    ndat2 <- decostand(t(data.list[[2]]), method = "standardize")
    pca.dat2 <- rda(ndat2)
    
    if(dimn == 2){
      choicesVec = c(1,2)
    }else{
      choicesVec = c(1,2,3)
    }
    res= procrustes(pca.dat1, pca.dat2, choices=choicesVec, symmetric = T, scale = TRUE)
    res2=protest(X = pca.dat1, Y = pca.dat2, scores = "sites", permutations = 999) 
    dataSet$misc["Sum of Squares"] = res2$ss
    dataSet$misc$Significance = res2$signif
    dataSet$misc$Correlation = res2$scale
    pos.xyz = rbind(res$X, res$Yrot)
    names = make.unique(as.character(rownames(pos.xyz)))
    newmeta$omics[c(1:(length(names)/2))] = names(d.list)[1]
    newmeta$omics[c( ((length(names)/2)+1) : (length(names)))] = names(d.list)[2]
    dataSet$dim.res <- list(res, res2)
  }
  if(reductionOpt != "mcia"){ # had to rescale with segment points
    pos.xyz <- as.data.frame(pos.xyz)
    pos.xyz <- unitAutoScale(pos.xyz);
    rownames(pos.xyz) <- names;
  }
  if(reductionOpt %in% c("diablo", "spls")){
    names2 <- rownames(pos.xyz2)
    pos.xyz2 <- as.data.frame(pos.xyz2)
    pos.xyz2 <- unitAutoScale(pos.xyz2);
    rownames(pos.xyz2) <- names2;
    dataSet$pos.xyz2 <- pos.xyz2
  }
  dataSet$pos.xyz = pos.xyz
  dataSet$newmeta = newmeta
  if(reductionOpt != "procrustes"){
    hit.inx <- match(loadingNames, unname(enrich.nms1));
    loadingSymbols <- names(enrich.nms1[hit.inx]);
  }
  if(reductionOpt %in% loadingOpts){
    loading.pos.xyz = as.data.frame(loading.pos.xyz)
    loading.pos.xyz <- unitAutoScale(loading.pos.xyz);
    
    rownames(loading.pos.xyz) = loadingNames
    dataSet$loading.names = loadingNames
    for(i in 1:length(names(data.list))){
      loadingNames <- gsub(paste0("_", names(data.list)[i]), "", loadingNames)
    }
    dataSet$loading.enrich = loadingSymbols
    dataSet$loading.pos.xyz = loading.pos.xyz
    if(reductionOpt == "o2pls"){ #second loading plot
      loading.pos.xyz2 <- as.data.frame(loading.pos.xyz2)
      loading.pos.xyz2 <- unitAutoScale(loading.pos.xyz2);
      rownames(loading.pos.xyz2) = loadingNames2
      dataSet$loading.names2 = loadingNames2
      dataSet$loading.pos.xyz2 = loading.pos.xyz2
    }
  }
  reductionOptGlobal <<- reductionOpt
  dataSet$omicstype <-names(data.list)
  .set.rdt.set(dataSet);
  
  if(dimn == 2){
    netData <-list()
    if(reductionOpt == "mcia"){
      pos.xyz = rbind(pos.xyz, dataSet$mcia.seg.points);
    }
    netData$score=pos.xyz
    if(reductionOpt %in% loadingOpts){
      netData$loading = list(pos=loading.pos.xyz, name=loadingSymbols);
    }else{
      netData$loading = "NA";
    }
    sink(TwoDJsonName);
    cat(toJSON(netData));
    sink();
  }
  return(1)
}

