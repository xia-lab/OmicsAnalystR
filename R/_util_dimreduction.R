##################################################
## R script for OmicsAnalyst
## Description: Compute dimension reduction
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################

my.reduce.dimension <- function(omicsType, reductionOpt, dimn){\  
  dimn = as.numeric(dimn);
  if(dimn == 2){
    dimn = 3;
  }

  sel.nms <- names(mdata.all)[mdata.all==1];
    data.list = list()
    d.list = list()
    omics.type = vector();
    featureNms <- vector();

    for(i in 1:length(sel.nms)){

      dataSet = readRDS(sel.nms[i])
      omics.type <- c(omics.type, dataSet$type)

      d.list[[dataSet$type]] = list()
      data.list[[dataSet$type]] <- dataSet$data.proc
      d.list[[dataSet$type]][["data.proc"]] = dataSet$data.proc #dat1
      d.list[[dataSet$type]][["meta"]] = dataSet$meta #meta1
      d.list[[dataSet$type]][["comp.res"]] = dataSet$comp.res #comp.res
      d.list[[dataSet$type]][["enrich.nms"]] = dataSet$enrich_ids

      if(i == 1){
        
        newmeta <- dataSet$meta
        comp.res1 = dataSet$comp.res
        enrich.nms1 = dataSet$enrich_ids
        comp.res.inx1 = rep(1, nrow(comp.res1));
        featureNms <- rownames(dataSet$data.proc);

      } else {
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

  if(reductionOpt == "diablo"){
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
  } else if(reductionOpt == "mcia") {
    library(omicade4)
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

  }

  if(reductionOpt != "mcia"){ # had to rescale with segment points
    pos.xyz <- as.data.frame(pos.xyz)
    pos.xyz <- unitAutoScale(pos.xyz);
    rownames(pos.xyz) <- names;
  }

  if(reductionOpt == "diablo"){
    names2 <- rownames(pos.xyz2)
    pos.xyz2 <- as.data.frame(pos.xyz2)
    pos.xyz2 <- unitAutoScale(pos.xyz2);
    rownames(pos.xyz2) <- names2;
    dataSet$pos.xyz2 <- pos.xyz2
  }

  dataSet$pos.xyz = pos.xyz
  dataSet$newmeta = newmeta

  hit.inx <- match(loadingNames, unname(enrich.nms1));
  loadingSymbols <- names(enrich.nms1[hit.inx]);

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

