
my.convert.igraph <- function(net.nm, filenm, idType="NA"){
  reductionSet <- .get.rdt.set();
  g <- ppi.comps[[net.nm]];
  
  V(g)$type <- V(g)$name;
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
      dataSet = qs::qread(sel.nms[micidx])
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
      dataSet2 = qs::qread(sel.nms[residx])
      dat2 = dataSet2$data.proc
      meta2 = dataSet2$meta
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
          dataSet = qs::qread(sel.nms[i]);
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
        }else{
          dataSet2 = qs::qread(sel.nms[i]);
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
        }
      }
    }
    
  }
  
  node.exp <- unname(expr.vec);
  
  current.net.nm <<- net.nm;
  # annotation
  nms <- V(g)$name;
  hit.inx <- match(nms, enrich.nms1);
  lbls <- names(enrich.nms1[hit.inx]);
  
  # setup shape (gene circle, other squares)
  shapes <- rep("circle", length(nms));
  
  if(!is.null(E(g)$weight)){
    E(g)$correlation <- E(g)$weight;
    E(g)$weight <- abs(E(g)$weight);
    ppi.comps[[net.nm]] <- g;
  }
  
  # get edge data
  edge.mat <- get.edgelist(g);
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
  if(is.null(V(g)$moltype)){
    mol.types <- rep("NA",length(node.exp)); 
  }else{
    mol.types <- V(g)$moltype; 
  }
  mol.types <- rep(names(data.list)[1],length(node.exp)); 
  mir.inx <- nms %in% enrich.nms2
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
  
  V(g)$moltype <- mol.types;
  V(g)$layers <- as.numeric(as.factor(mol.types));
  numOfTypes <- vector();
  
  
  if(!(is.null(reductionSet$taxlvl)) & reductionSet$taxlvl != "Feature"){
    
    
    dataSet <- qs::qread(sel.nms[[reductionSet$micidx]]);   
    dat.nms <- unique(dataSet$taxa_table[,reductionSet$taxlvl]);
    inx <- which(V(g)$name %in% dat.nms);
    inxNum <- V(g)$name %in% dat.nms;
    
    numOfTypes <- c(numOfTypes, sum(inxNum));
    names(numOfTypes)[reductionSet$micidx] <- sel.nms[[reductionSet$micidx]];
    V(g)$layers[inx] <- sel.nms[[reductionSet$micidx]];
    
    dataSet <- qs::qread(sel.nms[[reductionSet$residx]]);   
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
      
      dataSet <- qs::qread(sel.nms[[i]]);   
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
      featureId=nms[i],
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
  netData <- list(nodes=nodes, edges=edges.list, idType=idType, org=data.org, analType=anal.type, naviString = "network", modules=modules, tblNm="", nodeTypes= unique(mol.types), nodeColors = unique(color.vec) ,idType="entrez");
  
  if(!is.null(E(g)$correlation)){
    netData[["maxCorrelation"]] <- max(E(g)$correlation)
    netData[["minCorrelation"]] <- min(abs(E(g)$correlation))
  }
  
  jsonNms$network <<- filenm
  
  sink(filenm);
  cat(rjson::toJSON(netData));
  sink();
  return(1);
}

my.perform.layout <- function(net.nm, algo, focus=""){
  g <- ppi.comps[[net.nm]];
  vc <- vcount(g);
  if(algo == "Default"){
    if(vc > 3000) {
      pos.xy <- layout.lgl(g, maxiter = 100);
    }else if(vc > 2000) {
      pos.xy <- layout.lgl(g, maxiter = 150);
    }else if(vc > 1000) {
      pos.xy <- layout.lgl(g, maxiter = 200);
    }else if(vc < 150){
      pos.xy <- layout.kamada.kawai(g);
    }else{
      pos.xy <- layout.fruchterman.reingold(g);
    }
  }else if(algo == "FrR"){
    pos.xy <- layout.fruchterman.reingold(g);
  }else if(algo == "random"){
    pos.xy <- layout.random(g);
  }else if(algo == "lgl"){
    if(vc > 3000) {
      pos.xy <- layout.lgl(g, maxiter = 100);
    }else if(vc > 2000) {
      pos.xy <- layout.lgl(g, maxiter = 150);
    }else {
      pos.xy <- layout.lgl(g, maxiter = 200);
    }
  }else if(algo == "gopt"){
    # this is a slow one
    if(vc > 3000) {
      maxiter <- 50;
    }else if(vc > 2000) {
      maxiter <- 100;
    }else if(vc > 1000) {
      maxiter <- 200;
    }else{
      maxiter <- 500;
    }
    pos.xy <- layout.graphopt(g, niter=maxiter);
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
