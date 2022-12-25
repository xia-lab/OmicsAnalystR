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
  
  dataSet <- .get.rdt.set();
  reductionSet <- dataSet
  Sys.setenv(RGL_USE_NULL = TRUE)
  library(rgl)
  library(igraph)
  pos.xyz =  dataSet$pos.xyz
  nodes <- vector(mode="list");
  names <- c(rownames(pos.xyz))
  if(reductionOptGlobal %in% c(loadingOpts,"procrustes")){
    metadf = dataSet$newmeta
  }else{
    metadf = meta
  }
  
  a=list();
  a$objects = "NA";
  meshes="NA"
  
  col = vector();
  if("omics" %in% colnames(metadf)){
    sel.meta = "omics"
  }
  
  # can be selected meta as well if = dataSet$sel.meta
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
  
  if(reductionOptGlobal == "mcia"){
    names2 = rownames(dataSet$mcia.seg.points)
    nodes2 = list()
    for(i in 1:length(names2)){
      nodes2[[i]] <- list(
        id=names2[i],
        label=names2[i],
        size=0,
        meta="mcia.seg",
        cluster=1,
        fx = unname(dataSet$mcia.seg.points[i,1])*1000,
        fy = unname(dataSet$mcia.seg.points[i,2])*1000,
        fz = unname(dataSet$mcia.seg.points[i,3])*1000,
        colorb="#D3D3D3",
        colorw="#D3D3D3",
        topocolb="#D3D3D3",
        topocolw="#D3D3D3",
        expcolb="#D3D3D3",
        expcolw="#D3D3D3",
        attributes=list(
          expr = 1,
          degree=1,
          between=1
        )
      );
    }
    #nodes = c(nodes, nodes2); # do not use linking points now
  }else if(reductionOptGlobal %in% c("diablo", "spls")){
    names <- rownames(dataSet$pos.xyz2)
    nodes_samples2 <- vector(mode="list");
    pos.xyz2 =  dataSet$pos.xyz2;
    for(i in 1:length(names)){
      nodes_samples2[[i]] <- list(
        id=names[i],
        label=names[i],
        size=nodeSize,
        meta=meta.vec[i],
        cluster=1,
        fx = unname(pos.xyz2[i,1])*1000,
        fy = unname(pos.xyz2[i,2])*1000,
        fz = unname(pos.xyz2[i,3])*1000,
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
  }
  
  if(reductionOptGlobal == "procrustes"){
    pos.xyz.length = nrow(pos.xyz)
    edge.mat <- cbind(id=c(1:(pos.xyz.length)), source=names[c(1:(pos.xyz.length/2))], target=names[c(((pos.xyz.length/2)+1):pos.xyz.length) ], opacity = 0);
  #mcia hides linking edges
  #}else if(reductionOptGlobal == "mcia"){
   # pos.xyz.length = nrow(pos.xyz)
   # edge.mat <- cbind(id=c(1:(pos.xyz.length)), source=names, target=c(rownames(dataSet$mcia.seg.points), rownames(dataSet$mcia.seg.points)), opacity = 0);
  }else{
    edge.mat = "";
  }
  modules = "NA"
  # save node table
  ellipse ="NA"
  dataSet$pos.xyz = pos.xyz
  
  library(RJSONIO)
  
  if(!(reductionOptGlobal %in% loadingOpts) ){
    netData <- list(nodes=nodes, edges=edge.mat, modules=modules, objects=a$objects, ellipse=meshes, meta=metadf, loading="NA", reductionOpt=reductionOptGlobal, misc=dataSet$misc, omicstype=dataSet$omicstype);
  }else{
    if(reductionOptGlobal == "pca" ){
      metadf = dataSet$meta
    }else{
      metadf = dataSet$newmeta
    }
    
    loading.data = dataSet$loading.pos.xyz
    cluster = dataSet$loadingCluster
    aLoading=list();
    aLoading$objects = "NA";

    names = dataSet$loading.enrich
    ids = dataSet$loading.names
    rownames(loading.data) = names
    de = dataSet$comp.res[which(dataSet$comp.res[,"ids"] %in% ids),]
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
    nodes2 <- vector(mode="list");

    #loading.data = loading.data[which(rownames(loading.data) %in% as.character(ids)),];


    seed.inx <- names %in% unique(seeds);
    seed_arr <- rep("notSeed",length(names));
    seed_arr[seed.inx] <- "seed";

    for(i in 1:length(pvals)){
      nodes2[[i]] <- list(
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
    netData <- list(omicstype=omicstype.vec, nodes=nodes, edges=edge.mat, modules=modules, objects=a$objects, ellipse=meshes, meta=metadf, loading=nodes2, reductionOpt=reductionOptGlobal , objectsLoading=aLoading$objects, sigMat=sig.mats);
  }
  
  if(reductionOptGlobal %in% c("diablo", "spls")){
    type <- omicstype.vec[[2]]
    netData[[ type]] <- nodes_samples2;
    if(length(omicstype.vec)>2){
      sel.nms <- names(mdata.all)[mdata.all==1];
      
      for(j in 3:length(sel.nms)){
        names <- rownames(pos.xyz)
        
        nodes_samples3 <- vector(mode="list");
        pos <- reductionSet$dim.res$variates[[j]]
        pos<- unitAutoScale(pos[,c(1:3)]);
        for(i in 1:length(names)){
          nodes_samples3[[i]] <- list(
            id=names[i],
            label=names[i],
            size=nodeSize,
            meta=meta.vec[i],
            cluster=1,
            fx = unname(pos[i,1])*1000,
            fy = unname(pos[i,2])*1000,
            fz = unname(pos[i,3])*1000,
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
        type <- omicstype.vec[[j]]
        netData[[type]] <- nodes_samples3;
        
      }
    }
  }
  
  pca.scatter <- qs::qread("pca.scatter.qs");
  if(reductionOptGlobal != "procrustes"){
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
      #print(head(loading.data))
      loadingNames <- rownames(loading.data);
      loading.enrich = dataSet$enrich_ids[order(match(dataSet$enrich_ids, loadingNames))]
      ids = loadingNames
      names = names(loading.enrich)
      
      
      pca_loading <- nodes2[c(1:nrow(loading.data))];
      count = 1
      for(k in 1:length(nodes2)){
        if(nodes2[[k]][["id"]] %in% rownames(loading.data) || nodes2[[k]][["label"]] %in% rownames(loading.data)){
          pca_loading[[count]] <- nodes2[[k]]
          pca_loading[[count]][["id"]]=ids[k]
          pca_loading[[count]][["label"]]=names[k]
          pca_loading[[count]][["omicstype"]]=type.vec[k]
          inx = which(rownames(loading.data) == nodes2[[k]][["id"]]);
          if(length(inx) == 0){
            inx = which(rownames(loading.data) == nodes2[[k]][["label"]]);
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
  }
  
  dataSet$misc$pct2 <- c(dataSet$misc$pct2, pca.scatter$pct2);
  
  netData[["misc"]] <- dataSet$misc
  partialToBeSaved <<- c(partialToBeSaved, c(filenm));
  jsonNms$scatter <<- filenm;
  
  sink(filenm);
  cat(toJSON(netData));
  sink();
  
  .set.rdt.set(dataSet);
  return(1);
}
