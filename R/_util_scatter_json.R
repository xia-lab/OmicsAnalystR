##################################################
## R script for OmicsAnalyst
## Description: Compute 3D scatter plot from dimension reduction results
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################

my.json.scatter <- function(filenm){
  omicstype.vec <- list();
  sig.mats <- list();
  seeds <- vector()
  sel.nms <- names(mdata.all)[mdata.all==1];

  for(i in 1:length(sel.nms)){
    dataSet = readRDS(sel.nms[i])
    sig.mats[[dataSet$type]] <- dataSet$sig.mat
    omicstype.vec <- c(omicstype.vec, dataSet$type)
    if(i == 1){
      seeds <- rownames(dataSet$sig.mat) 
    }else{
      seeds <- c(seeds, rownames(dataSet$sig.mat))
    }
    meta <- dataSet$meta;
    sel.meta <-dataSet$sel.meta;
  }
  
  reductionSet <- .get.rdt.set();
  Sys.setenv(RGL_USE_NULL = TRUE)
  library(rgl)
  library(igraph)
  pos.xyz =  reductionSet$pos.xyz
  nodes <- vector(mode="list");
  names <- c(rownames(pos.xyz))
  metadf = reductionSet$meta

  
  a=list();
  a$objects = "NA";
  meshes="NA"
  
  col = vector();
  if("omics" %in% colnames(metadf)){
    sel.meta = "omics"
  }
  
  # can be selected meta as well if = reductionSet$sel.meta
  meta.vec = as.vector(metadf[,1])
  meta.vec.num = as.integer(as.factor(metadf[,1]))
  col.s <- gg_color_hue(length(unique(meta.vec)), "green")
  for(i in 1:length(meta.vec.num)){
    col[i] = col.s[meta.vec.num[i]];
  }
  color = col
  nodeSize = 18;

  if(length(names)>200){
    nodeSize = 16;
  }
  
  for(i in 1:length(names)){
    nodes[[i]] <- list(
      id=names[i],
      label=names[i],
      size=nodeSize,
      meta=meta.vec[i],
      cluster=meta.vec.num[i],
      fx = unname(pos.xyz[i,1])*1000,
      fy = unname(pos.xyz[i,2])*1000,
      fz = unname(pos.xyz[i,3])*1000,
      colorb=color[i],
      colorw=color[i],
      topocolb=color[i],
      topocolw=color[i],
      expcolb=color[i],
      expcolw=color[i],
      attributes=list(
        expr = 1,
        degree=1,
        between=1
      )
    );
  }

  modules = "NA"

  # save node table
  ellipse ="NA"
  
  library(RJSONIO)
    
    loading.data = reductionSet$loading.pos.xyz
    cluster = reductionSet$loadingCluster
    aLoading=list();
    aLoading$objects = "NA";

    names = reductionSet$loading.enrich
    ids = reductionSet$loading.names
    rownames(loading.data) = names
    de = reductionSet$comp.res[which(reductionSet$comp.res[,"ids"] %in% ids),]
    de[de == "NaN"] = 1
    pv = as.numeric(de[,"p_value"])
    pv_no_zero = pv[pv != 0]
    minval = min(pv_no_zero)
    pv[pv == 0] = minval/2
    pvals <<- -log10(pv);
    type.vec <- pvals;
    if(reductionSet$comp.res.inx[1] != "NA"){
      for(i in 1:length(unique(reductionSet$comp.res.inx))){
        inx = reductionSet$comp.res.inx == i
        type.vec[inx] <- omicstype.vec[[i]]
      }
    }
    colors<- ComputeColorGradient(pvals, "black", F, F);
    colorb <- colors;
    sizes <- as.numeric(rescale2NewRange(-log10(pv), 15, 25));
    loading.nodes <- vector(mode="list"); #is node2 loading?


    seed.inx <- names %in% unique(seeds);
    seed_arr <- rep("notSeed",length(names));
    seed_arr[seed.inx] <- "seed";

    for(i in 1:length(pvals)){
      loading.nodes[[i]] <- list(
        id=ids[i],
        label=names[i],
        size=sizes[i],
        cluster=1,
        omicstype=type.vec[i],
        fx = unname(loading.data[i,1])*1000,
        fy = unname(loading.data[i,2])*1000,
        fz = unname(loading.data[i,3])*1000,
        seedArr = seed_arr[i],
        colorb=colorb[i],
        colorw=colorb[i],
        topocolb="#ffa500",
        topocolw="#ffa500",
        expcolb="#ffa500",
        expcolw="#ffa500",
        attributes=list(
          expr = 1,
          degree=1,
          between=1
        )
      );
    }

    netData <- list(omicstype=omicstype.vec, nodes=nodes, edges="", modules=modules, objects=a$objects, ellipse=meshes, meta=metadf, loading=loading.nodes, reductionOpt=reductionOptGlobal , objectsLoading=aLoading$objects, sigMat=sig.mats);
  
  # user can compare dimred results to single-omics PCA
  pca.scatter <- qs::qread("pca.scatter.qs");
    for(i in 1:length(omicstype.vec)){
      pos<-pca.scatter[[paste0("pca_", omicstype.vec[i])]]$score
      
      pca_nodes <- nodes[c(1:nrow(pos))]
      for(j in 1:nrow(pos)){
        pca_nodes[[j]][["id"]] <-rownames(pos)[j]
        pca_nodes[[j]][["label"]] <-rownames(pos)[j]
        pca_nodes[[j]][["fx"]] <-pos[j,1]
        pca_nodes[[j]][["fy"]] <-pos[j,2]
        pca_nodes[[j]][["fz"]] <-pos[j,3]
      }
      nm <- paste0("pca_", omicstype.vec[[i]])
      netData[[nm]] <- pca_nodes;
      
      loading.data<-pca.scatter[[paste0("pca_", omicstype.vec[i])]]$loading
      loadingNames <- rownames(loading.data);
      loading.enrich = reductionSet$enrich_ids[order(match(reductionSet$enrich_ids, loadingNames))]
      ids = loadingNames
      names = names(loading.enrich)
      
      
      pca_loading <- loading.nodes[c(1:nrow(loading.data))];
      count = 1
      for(k in 1:length(loading.nodes)){
        if(loading.nodes[[k]][["id"]] %in% rownames(loading.data) || loading.nodes[[k]][["label"]] %in% rownames(loading.data)){
          pca_loading[[count]] <- loading.nodes[[k]]
          pca_loading[[count]][["id"]]=ids[k]
          pca_loading[[count]][["label"]]=names[k]
          pca_loading[[count]][["omicstype"]]=type.vec[k]
          inx = which(rownames(loading.data) == loading.nodes[[k]][["id"]]);
          if(length(inx) == 0){
            inx = which(rownames(loading.data) == loading.nodes[[k]][["label"]]);
          }
          pca_loading[[count]][["fx"]] = loading.data[inx,1]
          pca_loading[[count]][["fx"]] = loading.data[inx,2]
          pca_loading[[count]][["fx"]] = loading.data[inx,3]
          count = count +1;

        }
      }
      nm <- paste0("pca_", omicstype.vec[[i]], "_loading")
      netData[[nm]] <- pca_loading;
    }
  
  reductionSet$misc$pct2 <- c(reductionSet$misc$pct2, pca.scatter$pct2);
  
  netData[["misc"]] <- reductionSet$misc
  jsonNms$scatter <<- filenm;
  
  sink(filenm);
  cat(toJSON(netData));
  sink();
  
  .set.rdt.set(reductionSet);
  return(1);
}
