##################################################
## R script for OmicsAnalyst
## Description: Compute dimension reduction
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################

unsupervised.reduce.dimension <- function(reductionOpt){  
  ncomps = 10;
  
  sel.nms <- names(mdata.all)[mdata.all==1];
  data.list = list()
  omics.type = vector();
  featureNms <- vector();
  
  for(i in 1:length(sel.nms)){
    
    dataSet = readRDS(sel.nms[i])
    omics.type <- c(omics.type, dataSet$type)
    data.list[[dataSet$type]] <- dataSet$data.proc
    
    if(i == 1){       
      comp.res1 = dataSet$comp.res
      enrich.nms1 = dataSet$enrich_ids
      comp.res.inx1 = rep(1, nrow(comp.res1));
      featureNms <- rownames(dataSet$data.proc);
    } else {
      comp.res1 = rbind(comp.res1, dataSet$comp.res)
      enrich.nms1 = c(enrich.nms1, dataSet$enrich_ids);
      comp.res.inx1 = c(comp.res.inx1, rep(i, nrow(dataSet$comp.res)));
      featureNms <- c(featureNms, rownames(dataSet$data.proc));
    }
  }
  
  reductionSet <- list()
  reductionSet$comp.res = comp.res1
  reductionSet$enrich_ids = enrich.nms1
  reductionSet$comp.res.inx = comp.res.inx1
  reductionSet$meta = dataSet$meta
  
  if(reductionOpt == "mcia") {
    
    library(omicade4)
    mcoin <- mcia(data.list, cia.nf=ncomps)
    
    pos.xyz = mcoin$mcoa$SynVar;
    loading.pos.xyz = mcoin$mcoa$Tco;
    rownames(loading.pos.xyz) = featureNms
    
    # get sample and weight names
    loadingNames = featureNms
    names = rownames(pos.xyz)
    
    reductionSet$misc$pct = signif(mcoin$mcoa$pseudoeig,4)*100;
    reductionSet$pos.xyz = pos.xyz;
    reductionSet$loading.pos.xyz = loading.pos.xyz;
  }
  
  hit.inx <- match(loadingNames, unname(enrich.nms1));
  loadingSymbols <- names(enrich.nms1[hit.inx]);
  reductionSet$loading.enrich = loadingSymbols
  reductionSet$loading.names = loadingNames
  reductionSet$omicstype <-names(data.list)
  
  reductionOptGlobal <<- reductionOpt
  .set.rdt.set(reductionSet);
  
  return(1)
}


supervised.reduce.dimension <- function(reductionOpt){ 
  ncomps = 10;
  
  sel.nms <- names(mdata.all)[mdata.all==1];
  data.list = list()
  omics.type = vector();
  featureNms <- vector();
  
  for(i in 1:length(sel.nms)){
    
    dataSet = readRDS(sel.nms[i])
    omics.type <- c(omics.type, dataSet$type)
    data.list[[dataSet$type]] <- dataSet$data.proc
    
    if(i == 1){       
      comp.res1 = dataSet$comp.res
      enrich.nms1 = dataSet$enrich_ids
      comp.res.inx1 = rep(1, nrow(comp.res1));
      featureNms <- rownames(dataSet$data.proc);
    } else {
      comp.res1 = rbind(comp.res1, dataSet$comp.res)
      enrich.nms1 = c(enrich.nms1, dataSet$enrich_ids);
      comp.res.inx1 = c(comp.res.inx1, rep(i, nrow(dataSet$comp.res)));
      featureNms <- c(featureNms, rownames(dataSet$data.proc));
    }
  }
  
  reductionSet <- list()
  reductionSet$comp.res = comp.res1
  reductionSet$enrich_ids = enrich.nms1
  reductionSet$comp.res.inx = comp.res.inx1
  reductionSet$meta = dataSet$meta
  
  if(reductionOpt == "diablo"){
    library(mixOmics)
    Y <- reductionSet$meta[,1]
    
    dats = lapply(data.list, function(x){
      x <- data.frame(x, stringsAsFactors = FALSE);
      x <- t(x);
    });
    
    design = matrix(0.2, ncol = length(dats), nrow = length(dats), 
                    dimnames = list(names(dats), names(dats)));
    diag(design) = 0;
    
    res = mixOmics:::block.splsda(X = dats, Y = Y, ncomp = 5, design = design)
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
    
    for(i in 1:length(names(data.list))){
      if("prop_expl_var" %in% names(res)){
        var.vec <- res$prop_expl_var
      }else if("explained_variance" %in% names(res)){
        var.vec <- res$explained_variance
      }else{
        var.vec <- 0;
      }
      reductionSet$misc$pct[[names(data.list)[i]]] = unname(signif(as.numeric(var.vec[[i]],4)))*100
      if(i==1){
        loading.pos.xyz <- res$loadings[[i]]
        loadingNames=rownames(data.list[[i]])
      }else{
        loading.pos.xyz <- rbind(loading.pos.xyz, res$loadings[[i]])
        loadingNames=c(loadingNames,rownames(data.list[[i]]))
      }
    }
    rownames(loading.pos.xyz) <- loadingNames;
    reductionSet$pos.xyz = pos.xyz;
    reductionSet$loading.pos.xyz = loading.pos.xyz;
    
    if(! "diablo" %in% dim.res.methods){
      dim.res.methods <<- c(dim.res.methods, "diablo")
    }
  }
  
  hit.inx <- match(loadingNames, unname(enrich.nms1));
  loadingSymbols <- names(enrich.nms1[hit.inx]);
  reductionSet$loading.enrich = loadingSymbols
  reductionSet$loading.names = loadingNames
  reductionSet$omicstype <-names(data.list)
  
  reductionOptGlobal <<- reductionOpt
  .set.rdt.set(reductionSet);
  
    return(1);
}
