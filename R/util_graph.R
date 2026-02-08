
my.convert.igraph <- function(net.nm, fileNm, idType="NA"){
  library(igraph);
  reductionSet <- .get.rdt.set();
  g <- ppi.comps[[net.nm]];
  expr.vec <- rep(0, length(V(g)$name))
  names(expr.vec) = V(g)$name
  sel.nms <- names(mdata.all)[mdata.all==1];
  seeds <- vector();
  if(length(sel.nms) != 2){
    return(0)
  }else{
    data.list = list()
    if(!(is.null(reductionSet$taxlvl)) & reductionSet$taxlvl != "Feature"){
      
      micidx <-reductionSet$micidx 
      dataSet = readDataset(sel.nms[micidx])
      dat1 = dataSet$data.proc.taxa[[reductionSet$taxlvl]]
      meta1 = dataSet$meta
      comp.res1 = dataSet$comp.res.taxa[[reductionSet$taxlvl]]
      enrich.nms1 = unique(dataSet$taxa_table[,reductionSet$taxlvl])
      enrich.nms1 <- setNames(enrich.nms1,enrich.nms1)
      rownames(meta1) = paste0(rownames(meta1), ".", dataSet$type)
      comp.res.inx1 = rep(micidx, nrow(comp.res1))
      data.list[[dataSet$type]] = dat1
      comp.inx = which(names(expr.vec) %in% rownames(comp.res1))
      expr.vec[comp.inx] = as.numeric(comp.res1[which(rownames(comp.res1) %in% names(expr.vec)) ,"coefficient"])
      seeds <- rownames(dataSet$sig.mat.tax[[reductionSet$taxlvl]]) 
      
      residx <-reductionSet$residx
      dataSet2 = readDataset(sel.nms[residx])
      dat2 = dataSet2$data.proc
      meta2 = dataSet2$meta
      # OPTIMIZED: Direct rbind for 2 datasets instead of growing
      comp.res1 = rbind(comp.res1, dataSet2$comp.res)
      enrich.nms2 = dataSet2$enrich_ids
      rownames(meta2) = paste0(rownames(meta2), ".", dataSet2$type)
      newmeta = rbind(meta1, meta2);
      comp.res.inx1 = c(comp.res.inx1, rep(residx, nrow(dataSet2$comp.res)))
      enrich.nms1 = c(enrich.nms1,  enrich.nms2);
      data.list[[dataSet2$type]] = dat2
      comp.inx = which(names(expr.vec) %in% rownames(comp.res1))
      expr.vec[comp.inx] = as.numeric(comp.res1[which(rownames(comp.res1) %in% names(expr.vec)) , "coefficient"])
      seeds <- c(seeds, rownames(dataSet2$sig.mat)) 
      
    }else{
      
      for(i in 1:length(sel.nms)){
        if(i == 1){
          dataSet = readDataset(sel.nms[i]);
          dat1 = dataSet$data.proc;
          meta1 = dataSet$meta;
          comp.res1 = dataSet$comp.res;
          enrich.nms1 = dataSet$enrich_ids;
          rownames(meta1) = paste0(rownames(meta1), ".", dataSet$type);
          comp.res.inx1 = rep(1, nrow(comp.res1));
          data.list[[dataSet$type]] = dat1;
          comp.inx = which(names(expr.vec) %in% rownames(comp.res1));
          expr.vec[comp.inx] = as.numeric(comp.res1[which(rownames(comp.res1) %in% names(expr.vec)), "coefficient"]);
          seeds <- rownames(dataSet$sig.mat);
          idTypes <- dataSet$idType;
          fileNms <- dataSet$name;
        }else{
          dataSet2 = readDataset(sel.nms[i]);
          dat2 = dataSet2$data.proc;
          meta2 = dataSet2$meta;
          comp.res1 = rbind(comp.res1, dataSet2$comp.res);
          enrich.nms2 = dataSet2$enrich_ids;
          rownames(meta2) = paste0(rownames(meta2), ".", dataSet2$type);
          newmeta = rbind(meta1, meta2);
          comp.res.inx1 = c(comp.res.inx1, rep(i, nrow(dataSet2$comp.res)));
          enrich.nms1 = c(enrich.nms1,  enrich.nms2);
          data.list[[dataSet2$type]] = dat2;
          comp.inx = which(names(expr.vec) %in% rownames(comp.res1));
          expr.vec[comp.inx] = as.numeric(comp.res1[which(rownames(comp.res1) %in% names(expr.vec)), "coefficient"]);
          seeds <- c(seeds, rownames(dataSet2$sig.mat));
          idTypes <- c(idTypes, dataSet2$idType);
          fileNms <- c(fileNms,dataSet2$name);
        }
      }
    }
    
  }
  
  node.exp <- unname(expr.vec);
  
  current.net.nm <<- net.nm;
  # annotation
  nms <- V(g)$name;
  hit.inx <- match(nms, enrich.nms1);
  lbls <- V(g)$label;
  
  # setup shape (gene circle, other squares)
  shapes <- rep("circle", length(nms));
  
  if(!is.null(E(g)$weight)){
    E(g)$correlation <- E(g)$weight;
    E(g)$weight <- abs(E(g)$weight);
    ppi.comps[[net.nm]] <- g;
  }
  
  # get edge data
  edge.mat <- as_edgelist(g);
  if(!is.null(E(g)$correlation)){
    edge.sizes <- as.numeric(rescale2NewRange(abs(E(g)$correlation), 0.5, 3));
    edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], weight=signif(E(g)$correlation,4), size=edge.sizes, init_size=edge.sizes);
  }else if(!is.null(E(g)$weight)){
    edge.sizes <- as.numeric(rescale2NewRange(abs(E(g)$correlation), 0.5, 3));
    edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], weight=signif(E(g)$correlation,4), size=edge.sizes, init_size=edge.sizes);
  }else{
    edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2]);
  }
  
  # now get coords
  pos.xy <- PerformLayOut(net.nm, "backbone");
  
  # get the note data
  node.btw <- as.numeric(betweenness(g));
  node.dgr <- as.numeric(degree(g));
  # node size to degree values
  if(vcount(g)>500){
    min.size <- 2;
  }else if(vcount(g)>200){
    min.size <- 3;
  }else{
    min.size <- 4;
  }
  node.sizes <- as.numeric(rescale2NewRange((log(node.btw + 1))^2, min.size, 12));
  centered <- T;
  notcentered <- F;
  # update node color based on betweenness
  require("RColorBrewer");
  topo.val <- log(node.btw+1);
  topo.colsb <- ComputeColorGradient(topo.val, "black", notcentered, FALSE);
  topo.colsw <- ComputeColorGradient(topo.val, "white", notcentered, FALSE);
  topo.colsc <- ComputeColorGradient(topo.val, "colorblind", notcentered, TRUE);
  
  # color based on expression
  bad.inx <- is.na(node.exp) | node.exp==0;
  if(!all(bad.inx)){
    exp.val <- node.exp;
    if(names(table(expr.vec > 0))[1] == "TRUE"){
      centerBool = F;
    }else{
      centerBool = T;
    }
    node.colsb.exp <- ComputeColorGradient(exp.val, "black", centerBool, FALSE); 
    node.colsw.exp <- ComputeColorGradient(exp.val, "white", centerBool, FALSE);
    node.colsc.exp <- ComputeColorGradient(exp.val, "colorblind", centerBool, TRUE);
    
    node.colsb.exp[bad.inx] <- "#d3d3d3"; 
    node.colsw.exp[bad.inx] <- "#c6c6c6"; 
    node.colsc.exp[bad.inx] <- "#c6c6c6"; 
    # node.colsw.exp[bad.inx] <- "#b3b3b3";
  }else{
    node.colsb.exp <- rep("#d3d3d3",length(node.exp)); 
    node.colsw.exp <- rep("#c6c6c6",length(node.exp)); 
    node.colsc.exp <- rep("#b3b3b3",length(node.exp)); 
  }
  
  # now update for bipartite network
  mol.types <- rep("NA",length(node.exp)); 
  #mol.types <- rep(names(data.list)[1],length(node.exp)); 
  mol.types <- V(g)$type;
  mir.inx <- V(g)$type %in% unique(V(g)$type)[2]
  shapes[mir.inx] <- "square";
  mol.types[mir.inx] <- names(data.list)[2]
  
  freq = table(mol.types);
  
  duplicated.types=mol.types;
  for(i in 1:length(unique(mol.types))){
    duplicated.types[duplicated.types == names(freq[i])]=order(freq)[i];
  }

  
  node.cols <- topo.colsb;
  ntype <- names(freq);
  color.vec <- generate_colors(length(ntype), "colorblind");
  for(i in 1:length(ntype)){
    if(ntype[i] %in% c("Protein", "Seed", "Gene") ){
      color.vec[i] = "#BD0313";
    }else{
      node.cols[which(mol.types ==ntype[i])]=color.vec[i];
    }
  }
  
  
  colVec <- color.vec;
  
  V(g)$moltype <- V(g)$type;
  V(g)$layers <- as.numeric(as.factor(V(g)$type));
  numOfTypes <- vector();
  
  
  if(!(is.null(reductionSet$taxlvl)) & reductionSet$taxlvl != "Feature"){
    
    
    dataSet <- readDataset(sel.nms[[reductionSet$micidx]]);   
    dat.nms <- unique(dataSet$taxa_table[,reductionSet$taxlvl]);
    inx <- which(V(g)$name %in% dat.nms);
    inxNum <- V(g)$name %in% dat.nms;
    
    numOfTypes <- c(numOfTypes, sum(inxNum));
    names(numOfTypes)[reductionSet$micidx] <- sel.nms[[reductionSet$micidx]];
    V(g)$layers[inx] <- sel.nms[[reductionSet$micidx]];
    
    dataSet <- readDataset(sel.nms[[reductionSet$residx]]);   
    dat.nms <- unique(unname(dataSet$enrich_ids));
    inx <- which(V(g)$name %in% dat.nms);
    inxNum <- V(g)$name %in% dat.nms;
    
    numOfTypes <- c(numOfTypes, sum(inxNum));
    names(numOfTypes)[reductionSet$residx] <- sel.nms[[reductionSet$residx]];
    V(g)$layers[inx] <- sel.nms[[reductionSet$residx]];
    hl.node <- unique(c(reductionSet[["toHighlight"]]$from,reductionSet[["toHighlight"]]$to));
    hl.idx <- match(hl.node,names(expr.vec));
    color.vec <- c(color.vec,"#FFD700");
    node.cols[hl.idx] <- "#FFD700";
    
    
    edge.modify <- merge(reductionSet[["toHighlight"]][,-4],as.data.frame(edge.mat),all.x = TRUE,,by.x = c("from","to"),by.y = c("source","target"));
    
    edge.modify$weight <-  as.numeric(edge.modify$weight)+  as.numeric(edge.modify$potential);
    edge.mat<- cbind(edge.mat,highlight=0);
    edge.mat[match(edge.modify$id,edge.mat[,"id"]),"weight"] <- edge.modify$weight;
    edge.mat[match(edge.modify$id,edge.mat[,"id"]),"highlight"] <- 1;
    
    
  }else{
    
    for( i in 1:length(sel.nms)){
      
      dataSet <- readDataset(sel.nms[[i]]);   
      dat.nms <- unique(unname(dataSet$enrich_ids));
      inx <- which(V(g)$name %in% dat.nms)
      inxNum <- V(g)$name %in% dat.nms
      
      numOfTypes <- c(numOfTypes, sum(inxNum))
      names(numOfTypes)[i] <- sel.nms[[i]]
      V(g)$layers[inx] <- sel.nms[[i]];
    }
  }
  
  
  if(length(ntype)==1){
    # seed.inx <- nms %in% unique(seed.proteins);
    # mol.types[seed.inx] <- "Seed"
    # mol.types[!seed.inx] <- "Protein"
    
  }else{
    topo.colsw <- node.cols; # dark blue
    topo.colsb <- node.cols;
    topo.colsc <- node.cols;
  }
  numOfTypes <- sort(numOfTypes)
  for( i in 1:length(names(numOfTypes))){
    V(g)$layers[ which(V(g)$layers == names(numOfTypes)[i])] = as.numeric(i);
  }
  V(g)$layers <- as.numeric(V(g)$layers)
  
  ppi.comps[[net.nm]] <<- g;
  
  seed.inx <- nms %in% unique(seeds);
  seed_arr <- rep("notSeed",length(node.dgr));
  seed_arr[seed.inx] <- "seed";
  
  # now create the json object
  nodes <- vector(mode="list");

  ids_and_omicstype = paste0(nms, "_", mol.types);
  for(i in 1:length(node.sizes)){
    
    nodes[[i]] <- list(
      id=nms[i],
      featureId=V(g)$featureId[i],
      idnb = i, 
      label=lbls[i],
      x = pos.xy[i,1],
      y = pos.xy[i,2],
      molType = mol.types[i],
      size=node.sizes[i], 
      seedArr = seed_arr[i],
      type=shapes[i],
      colorb=topo.colsb[i],
      colorw=topo.colsw[i],
      colorc=topo.colsc[i],
      expcolb=node.colsb.exp[i],
      expcolw=node.colsw.exp[i],
      expcolc=node.colsc.exp[i],
      highlight = 0,
      attributes=list(
        expr = node.exp[i],
        degree=node.dgr[i], 
        between=node.btw[i])
    );
  }
  
  # save node table
  nd.tbl <- data.frame(Id=nms, Label=lbls, Degree=node.dgr, Betweenness=round(node.btw,2), Expression=node.exp);
  # order 
  ord.inx <- order(nd.tbl[,3], nd.tbl[,4], decreasing = TRUE)
  nd.tbl <- nd.tbl[ord.inx, ];
  fast.write(nd.tbl, file="node_table.csv", row.names=FALSE);
  reductionSet$imgSet$node_table <- nd.tbl;

  # Export to Arrow for Java DataTable zero-copy access
  tryCatch({
    ExportNodeTableArrow()
  }, error = function(e) {
    warning(paste("NodeTable Arrow export failed:", e$message))
  })

  if(length(V(g)$name)>100 && ppi.net$db.type != "uploaded"){
    modules <- FindCommunities("walktrap", FALSE);
  }else{
    modules <- "NA"
  }
  node.res <- ppi.net$node.res
  if("domain" %in% colnames(node.res) || ppi.net$db.type == "signal"){
    node.info <- edge.infoU;
    node.types <- node.types
  }else{
    node.info <- "";
    node.types <- "";
  }
  
  # covert to json
  require(rjson);
  #formattin json file because of rjson
  edges.list <- apply(edge.mat, 1, as.list)
  netData <- list(nodes=nodes, edges=edges.list, idTypes=idTypes,fileNms = fileNms, org=data.org, analType=anal.type, naviString = "network", modules=modules, tblNm="", nodeTypes= unique(mol.types),omicstype = unique(mol.types),nodeColors = unique(color.vec) ,idType="entrez");
  
  if(!is.null(E(g)$correlation)){
    netData[["maxCorrelation"]] <- max(E(g)$correlation)
    netData[["minCorrelation"]] <- min(abs(E(g)$correlation))
  }

  network_img <- saveNetworkImage(g, net.nm, pos.xy, node.cols, node.sizes, shapes, E(g)$correlation);
  if(!is.null(network_img)){
    netData[["network.image"]] <- network_img;
  }

  infoSet <- readSet(infoSet, "infoSet");
  infoSet$paramSet$jsonNms$network <- fileNm
  saveSet(infoSet);
  .set.rdt.set(reductionSet)

  sink(fileNm);
  cat(rjson::toJSON(netData));
  sink();
  return(1);
}

saveNetworkImage <- function(graph, netName, layoutCoords, vertexColors, vertexSizes, vertexShapes, edgeCorr){
  if(is.null(graph) || is.null(netName) || is.null(layoutCoords)){
    return(NULL);
  }
  load_cairo();
  load_pheatmap(); # ensures Cairo available? (no, just width)
  safeName <- gsub("[^A-Za-z0-9_]+", "_", netName);
  if(safeName == ""){
    safeName <- "network";
  }
  imgName <- paste0("network_top50_", safeName, ".png");
  width <- max(6, min(14, nrow(layoutCoords)/10 + 6));
  height <- max(6, min(10, nrow(layoutCoords)/15 + 5));
  Cairo::Cairo(file = imgName, unit="in", dpi=150, width=width, height=height, type="png", bg="#0f0f0f");
  edgeCols <- if(!is.null(edgeCorr)){
    ifelse(edgeCorr >= 0, "#46a6ff", "#ff7f45")
  } else {
    rep("#2f4f4f", length(E(graph)))
  }
  edgeWidths <- if(!is.null(edgeCorr)){
    pmax(0.6, abs(edgeCorr) * 2)
  } else {
    rep(1, length(E(graph)))
  }
  vertex.shape <- ifelse(vertexShapes == "square", "square", "circle");
  plot(graph,
       layout=layoutCoords,
       vertex.color=vertexColors,
       vertex.size=vertexSizes,
       vertex.shape=vertex.shape,
       vertex.label=NA,
       edge.color=edgeCols,
       edge.width=edgeWidths,
       edge.curved=0.2,
       vertex.frame.color="black",
       margin=c(0,0,0,0));
  dev.off();
  return(imgName);
}

my.perform.layout <- function(net.nm, algo, focus=""){
  library(igraph);
  g <- ppi.comps[[net.nm]];
  vc <- vcount(g);
  if(algo == "Default"){
    if(vc > 5000) {
      pos.xy <- layout_with_lgl(g);
    }else if(vc < 100){
      pos.xy <- layout_with_kk(g);
    }else{
      pos.xy <- layout_with_fr(g);
    }
  }else if(algo == "FrR"){
    pos.xy <- layout_with_fr(g, area=34*vc^2);
  }else if(algo == "random"){
    pos.xy <- layout_randomly (g);
  }else if(algo == "lgl"){
    pos.xy <- layout_with_lgl(g);
  }else if(algo == "gopt"){
    pos.xy <- layout_with_graphopt(g)
  }else if(algo == "fr"){
    pos.xy <- layout_with_fr(g, dim=3, niter=500)
  }else if(algo == "kk"){
    pos.xy <- layout_with_kk(g, dim=3, maxiter=500)
  }else if(algo == "tree"){
    l <- layout_with_sugiyama(g, vgap=vc/4)
    pos.xy <- -l$layout
  }else if(algo == "circular_tripartite"){
    library(ggforce)
    l <- layout_with_sugiyama(g, layers = as.numeric(V(g)$layers)*(vc/3) +30)
    layout <- l$layout
    
    radial <- radial_trans(
      r.range = rev(range(layout[,2])),
      a.range = range(layout[,1]),
      offset = 0
    )
    coords <- radial$transform(layout[,2], layout[,1])
    layout[,1] <- coords$x
    layout[,2] <- coords$y
    pos.xy <-layout
  }else if(algo == "tripartite"){
    l <- layout_with_sugiyama(g, layers = as.numeric(V(g)$layers)*(vc/4))
    pos.xy <- -l$layout[,2:1] 
  }else if(algo == "concentric"){
    library(graphlayouts)
    # the fist element in the list for concentric is the central node.
    if(focus==""){
      inx=1;
    }else{
      inx = which(V(g)$name == focus)
    }
    coords <- layout_with_focus(g,inx)
    pos.xy <- coords$xy
  }else if(algo == "backbone"){
    library(graphlayouts)
    if(length(V(g)$name)<2000){
      coords = layout_with_stress(g)
      pos.xy = coords
    }else{
      coords = layout_with_sparse_stress(g,pivots=100)
      pos.xy = coords
    }
  }
  pos.xy;
}
