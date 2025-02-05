my.enrich.net<-function( netNm="abc", type="list", overlapType="mixed"){
  
  infoSet <- readSet(infoSet, "infoSet");
 
  enr.mat <- infoSet$imgSet$enrTables[[type]]$table

  hits <-  enr.mat[,"Hits"];
  pvals <- enr.mat[,"P.Value"];
  
  pvalue <- pvals;
  names(pvalue)<-enr.mat$Pathway
  id <- enr.mat$Pathway;
  
  
  if(is.null(enr.mat)){
    return(0);
  }
  
  require(igraph);
  require(reshape);
  
  current.geneset <- infoSet$imgSet$enrTables[[type]]$current.geneset;
  hits.query <- infoSet$imgSet$enrTables[[type]]$hits.query;
  hits.query <- hits.query[enr.mat$Pathway];
  geneSets <- hits.query;
  
  n <- nrow(enr.mat);
  w <- matrix(NA, nrow=n, ncol=n);
  colnames(w) <- rownames(w) <- id;
  for (i in 1:n) {
    for (j in i:n) {
      w[i,j] <- overlap_ratio(geneSets[id[i]], geneSets[id[j]], overlapType)
    }
  }
  wd <- reshape::melt(w);
  wd <- wd[wd[,1] != wd[,2],];
  wd <- wd[!is.na(wd[,3]),];
  
  g <- graph_from_data_frame(wd[,-3], directed=F);
  if(type == "list"){
    g <- delete_edges(g, E(g)[wd[,3] < 0.3]);
  }else{
    g <- delete_edges(g, E(g)[wd[,3] < 0.3]);
  }
  idx <- unlist(sapply(V(g)$name, function(x) which(x == id)));
  
  # define local function
  my.normalize <- function(x){
    return((x- min(x)) /(max(x)-min(x)))
  }
  my.rescale <- function(x, from, to){
    (x - min(x)) / max(x - min(x)) * (to - from) + from
  }
  
  
  # Compute normalized p-values for color gradient
  normalized_pvalues <- -log(my.normalize(pvalue) + min(pvalue/2))
  
  # Ensure that you compute colors only for existing vertices in the graph
  existing_vertices <- V(g)$name
  vertex_colors <- ComputeColorGradient(normalized_pvalues[names(normalized_pvalues) %in% existing_vertices], "black", F, F)
  vertex_colorsw <- ComputeColorGradient(normalized_pvalues[names(normalized_pvalues) %in% existing_vertices], "white", F, F)
  
  # Assign colors only to existing vertices
  V(g)$color <- vertex_colors
  V(g)$colorw <- vertex_colorsw
  
  cnt <- hits;
  names(cnt) <- id;
  cnt2 <- cnt[V(g)$name];
  
  if (all(cnt2 == cnt2[1])){
    V(g)$size <- rep(16, length(cnt2))
  }else{
    V(g)$size <- my.rescale(log(cnt2+1, base=10), 8, 32);
  }
  
  # layout
  pos.xy <- layout_nicely(g);
  
  # now create the json object
  nodes <- vector(mode="list");
  node.nms <- V(g)$name;
  node.sizes <- V(g)$size;
  node.cols <- V(g)$color;
  node.colsw <- V(g)$colorw;
  
  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id = node.nms[i],
      label=node.nms[i],
      size = node.sizes[i],
      true_size=node.sizes[i], 
      molType="set",
      colorb=node.cols[i],
      colorw=node.colsw[i],
      posx = pos.xy[i,1],
      posy = pos.xy[i,2]
    );
  }
  
  edge.mat <- as_edgelist(g);
  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2]);
  
  # covert to json
  bedges <- stack(hits.query);
  b.mat <- matrix(NA, nrow=nrow(bedges), ncol=2);
  b.mat[,1] <- bedges[,"values"];
  b.mat[,2] <- as.character(bedges[,"ind"]);
  b.mat <- b.mat[complete.cases(b.mat),]
  colnames(b.mat) <- c("source", "target");
  bg <- graph.data.frame(b.mat, directed=F);
  idx <- unlist(sapply(V(bg)$name, function(x) which(x == id)));
  cols <- color_scale("red", "#E5C494");
  V(bg)$label<- V(bg)$name
  
  V(bg)$color[V(bg)$name %in% enr.mat$Pathway] <- ComputeColorGradient(-log(pvalue), "black", F, F);
  V(bg)$colorw[V(bg)$name %in% enr.mat$Pathway] <- ComputeColorGradient(-log(pvalue), "white", F, F);
  node.nms <- V(bg)$name;
  if(type == "limma"){
    tbl <- infoSet$imgSet$enrTables[[type]]$sig.mat
    V(bg)$label[!V(bg)$name %in% enr.mat$Pathway] <- tbl$label[match(V(bg)$name[!V(bg)$name %in% enr.mat$Pathway], tbl$ids)]
 
    tbl <- tbl[which(rownames(tbl) %in% V(bg)$name),]
    expr.val <- tbl[,1];
    expvals <- expr.val;
    names(expvals) <- tbl$label
    expvals <- expvals[V(bg)$label]
    V(bg)$color[!V(bg)$name %in% enr.mat$Pathway] <- ComputeColorGradient(unname(expvals), "black", T,T);
    V(bg)$colorw[!V(bg)$name %in% enr.mat$Pathway] <- ComputeColorGradient(unname(expvals), "black", T, T);
  }else if(anal.type == "genelist" && sum(as.numeric(paramSet$all.prot.mat[,1])) != 0){
    tbl <- paramSet$all.prot.mat
    gene.nms <- V(bg)$name[which(!V(bg)$name %in% rownames(enr.mat))]
    tbl <- tbl[which(tbl[,2] %in% gene.nms),]
 
    expr.val <- tbl[,1];
    names(expr.val) <- tbl
    expvals <- expr.val
    expvals <- expvals[node.nms]
    expvals <- expvals[!is.na(expvals)]
    V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T, T);
    V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T, T);
    
  }else{
    expvals <- rep(0,length(V(bg)$color)); 
    V(bg)$color[!V(bg)$name %in% enr.mat$Pathway] <- "#00FFFF";
    V(bg)$colorw[!V(bg)$name %in% enr.mat$Pathway] <- "#668B8B"
  }
  node.dgr2 <- as.numeric(degree(bg));
  V(bg)$size <- my.rescale(log(node.dgr2, base=10), 8, 24); 
  
  # layout
  pos.xy <- layout_nicely(bg);
  
  # now create the json object
  bnodes <- vector(mode="list");
  node.sizes <- V(bg)$size;
  node.cols <- V(bg)$color;
  node.colsw <- V(bg)$colorw;
  
  shapes <- rep("circle", length(node.nms));
  hit.inx <- node.nms %in% b.mat[,"source"];
  shapes[hit.inx] <- "gene";
  node.lbls <- V(bg)$label
  for(i in 1:length(node.sizes)){
    bnodes[[i]] <- list(
      id = node.nms[i],
      label=node.lbls[i], 
      size=node.sizes[i], 
      colorb=node.cols[i],
      colorw=node.colsw[i],
      true_size=node.sizes[i], 
      molType=shapes[i],
      exp= unname(expvals[node.lbls[i]]),
      posx = pos.xy[i,1],
      posy = pos.xy[i,2]
    );
  }
  
  ppi.comps <- vector(mode="list");
  infoSet$paramSet$current.net.nm <- netNm
  ppi.comps[[netNm]] <- bg;
  infoSet$analSet$ppi.comps <- ppi.comps
  
  bedge.mat <- get.edgelist(bg);
  bedge.mat <- cbind(id=paste0("b", 1:nrow(bedge.mat)), source=bedge.mat[,1], target=bedge.mat[,2]);
  initsbls <-  setNames(tbl$label,tbl$ids)
  
  #for rjson generation
  edge.mat <- apply(edge.mat, 1, as.list)
  bedge.mat <- apply(bedge.mat, 1, as.list)
  enr.mat <- apply(enr.mat, 1, as.list)
   
  
  netData <- list(nodes=nodes, 
                  edges=edge.mat, 
                  bnodes=bnodes, 
                  bedges=bedge.mat, 
                  enr=unname(enr.mat), 
                  id=names(enr.mat), 
                  #sizes=analSet$listSizes, 
                  hits=hits.query, 
                  genelist=initsbls, 
                  analType=type, 
                  #org=paramSet$data.org, 
                  backgroundColor=list("#514F6A", "#222222"),
                  dat.opt = infoSet$imgSet$enrTables[[type]]$selDataNm,
                  naviString = "Enrichment Network");
 
  netName <- paste0(netNm, ".json");
  print(netName)
  partialToBeSaved <- c( infoSet$paramSet$partialToBeSaved, c(netName));
  infoSet$paramSet$jsonNms <- netName;
  infoSet$enrichNet <- netData;
  saveSet(infoSet);
  sink(netName);
  cat(rjson::toJSON(netData));
  sink();
  
  return(1);
}

overlap_ratio <- function(x, y, type) {
  x <- unlist(x)
  y <- unlist(y)
  if(type == "mixed"){
    res <- 0.5 * length(intersect(x, y))/length(unique(y)) + 0.5 * length(intersect(x, y))/length(unique(x));
  }else if(type == "overlap"){
    if(length(x)>length(y)){
      res=length(intersect(x, y))/length(unique(y))
    }else{
      res=length(intersect(x, y))/length(unique(x))
    }
  }else{
    res=length(intersect(x, y))/length(unique(c(x,y)))
  }
  return(res)
}
