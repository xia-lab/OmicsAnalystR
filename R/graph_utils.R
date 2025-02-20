##################################################
## R script for OmicsAnalyst
## Description: General graph manipulation functions 
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################

GetColorSchema <- function(my.grps){
  # test if total group number is over 9
  my.grps <- as.factor(my.grps);
  grp.num <- length(levels(my.grps));
  if(grp.num > 9){
    pal12 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
               "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
               "#FFFF99", "#B15928");
    dist.cols <- colorRampPalette(pal12)(grp.num);
    lvs <- levels(my.grps);
    colors <- vector(mode="character", length=length(my.grps));
    for(i in 1:length(lvs)){
      colors[my.grps == lvs[i]] <- dist.cols[i];
    }
  }else{
    colors <- as.numeric(my.grps)+1;
  }
  return (colors);
}


DecomposeGraph <- function(gObj, minNodeNum = 2, jsonBool = F){

  # now decompose to individual connected subnetworks
    if(jsonBool == "netjson"){
        comps <-list(gObj)
    }else{
        comps <-decompose.graph(gObj, min.vertices=minNodeNum);
    }
  if(length(comps) == 0){
    current.msg <<- paste("No subnetwork was identified with at least", minNodeNum, "nodes!");
    return(NULL);
  }
  
  # first compute subnet stats
  net.stats <- ComputeSubnetStats(comps);
  ord.inx <- order(as.numeric(net.stats[,2]), decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  comps <- comps[ord.inx];
  names(comps) <- rownames(net.stats) <- paste("subnetwork", 1:length(comps), sep="");
  
  # note, we report stats for all nets (at least 3 nodes);
  hit.inx <- net.stats$Node >= minNodeNum;
  comps <- comps[hit.inx];
  #overall <- list();
  #overall[["overall"]] <- g
  #ppi.comps <- append(overall, ppi.comps, after=1);
 
  # now record
  ppi.comps <<- comps;
  net.stats <<- net.stats;
  
  sub.stats <- unlist(lapply(comps, vcount)); 
  return(sub.stats);
}


ComputeSubnetStats <- function(comps){
  net.stats <- as.data.frame(matrix(0, ncol = 3, nrow = length(comps)));
  colnames(net.stats) <- c("Node", "Edge", "Query");
  reductionSet <- .get.rdt.set();
  if(exists("corr.mat.inter.taxa",reductionSet) & reductionSet$taxlvl !="Feature"){

      for(j in 1:length(comps)){
        g <- comps[[j]];
        nd.queries <- V(g)$name;
        sel.nms <- names(mdata.all)[mdata.all==1];
        nd.res <- "";
        dataSet <- readDataset(sel.nms[reductionSet$micidx]);
        lbl = reductionSet$taxlvl
     
        nms <- unique(dataSet$taxa_table[,reductionSet$taxlvl]);
          if(sum(nms %in% nd.queries)>0 && !grepl(lbl, nd.res)){
            nd.res <- paste0(lbl,": ", sum(nms %in% nd.queries), "; ", nd.res)
          }
           dataSet <- readDataset(sel.nms[reductionSet$residx]);
        lbl = "Metabolite"
        
        nms <- unique(unname(dataSet$enrich_ids));
        if(sum(nms %in% nd.queries)>0 && !grepl(lbl, nd.res)){
          nd.res <- paste0(lbl,": ", sum(nms %in% nd.queries), "; ", nd.res)
        }
        
          
        }
        net.stats[j,] <- c(nd.res ,ecount(g),0);
      

  }else{
      
    for(j in 1:length(comps)){
      g <- comps[[j]];
      nd.queries <- V(g)$name;
      sel.nms <- names(mdata.all)[mdata.all==1];
      nd.res <- "";
      for( i in 1:length(sel.nms)){
           dataSet <- readDataset(sel.nms[[i]]);
          lbl = dataSet$readableType;
        
        if(sum(V(g)$type == dataSet$type) && !grepl(lbl, nd.res)){
          nd.res <- paste0(lbl,": ", sum(V(g)$type == dataSet$type), "; ", nd.res)
        }
        
      }
      net.stats[j,] <- c(nd.res ,ecount(g),0);
    }
  
    }

  return(net.stats);
}


ComputeColorGradient <- function(nd.vec, background="black", centered, colorblind){
  require("RColorBrewer");
  
  minval <- min(nd.vec, na.rm=TRUE);
  maxval <- max(nd.vec, na.rm=TRUE);
  res <- maxval-minval;
  
  if(res == 0){
    return(rep("#FF0000", length(nd.vec)));
  }
  color <- GetColorGradient(background, centered, colorblind);
  breaks <- generate_breaks(nd.vec, length(color), center = centered);
  return(scale_vec_colours(nd.vec, col = color, breaks = breaks));
}

generate_breaks <- function(x, n, center = F){
  if(center){
    m <- max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    res <- seq(-m, m, length.out = n + 1)
  }
  else{
    res <- seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
  }
  return(res)
}

scale_vec_colours <- function(x, col = rainbow(10), breaks = NA){
  breaks <- sort(unique(breaks));
  return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
}

GetColorGradient <- function(background, center, colorblind=F) {
  if (background == "black") {
    if (center) {
      if (colorblind) {
        return(c(colorRampPalette(c("#3182bd", "#6baed6", "#bdd7e7", "#eff3ff"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      } else {
        return(c(colorRampPalette(c("#31A231", "#5BC85B", "#90EE90", "#C1FFC1"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      }
    } else {
      if (colorblind) {
        return(c(colorRampPalette(c("#3182bd", "#bbfdff"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      } else {
        return(colorRampPalette(rev(heat.colors(9)))(100))
      }
    }
  } else {
    if (center) {
      if (colorblind) {
        return(c(colorRampPalette(c("#3182bd", "#bbfdff"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      } else {
        return(c(colorRampPalette(c("#137B13", "#31A231", "#5BC85B", "#90EE90"))(50), colorRampPalette(c("#FF7783", "#E32636", "#BD0313", "#96000D"))(50)))
      }
    } else {
      return(colorRampPalette(hsv(h = seq(0.72, 1, 0.035), s = 0.72, v = 1))(100))
    }
  }
}

# also save to GraphML
ExportNetwork <- function(fileName){
  current.net <- ppi.comps[[current.net.nm]];
  write.graph(current.net, file=fileName, format="graphml");
}


ExtractModule<- function(nodeids){
  set.seed(8574);
  nodes <- strsplit(nodeids, "; ")[[1]];
  g <- ppi.comps[[current.net.nm]];
  # try to see if the nodes themselves are already connected
  hit.inx <- V(g)$name %in% nodes; 
  gObj <- induced.subgraph(g, V(g)$name[hit.inx]);
  
  # now find connected components
  comps <-decompose.graph(gObj, min.vertices=1);
  
  if(length(comps) == 1){ # nodes are all connected
    g <- comps[[1]];
  }else{
    # extract modules
    paths.list <-list();
    sd.len <- length(nodes);
    for(pos in 1:sd.len){
      paths.list[[pos]] <- get.shortest.paths(g, nodes[pos], nodes[-(1:pos)])$vpath;
    }
    nds.inxs <- unique(unlist(paths.list));
    nodes2rm <- V(g)$name[-nds.inxs];
    g <- simplify(delete.vertices(g, nodes2rm), edge.attr.comb="first");
  }
  nodeList <- get.data.frame(g, "vertices");
  if(nrow(nodeList) < 3){
    return ("NA");
  }

  if(ncol(nodeList) == 1){
    nodeList <- data.frame(Id=nodeList[,1], Label=nodeList[,1]);
  }
  
  module.count <- module.count + 1;
  module.nm <- paste("module", module.count, sep="");
  colnames(nodeList) <- c("Id", "Label");
  ndFileNm <- paste(module.nm, "_node_list.csv", sep="");
  fast.write(nodeList, file=ndFileNm, row.names=F);
  
  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  edgFileNm <- paste(module.nm, "_edge_list.csv", sep="");
  fast.write(edgeList, file=edgFileNm, row.names=F);
  
  filenm <- paste(module.nm, ".json", sep="");
  
  # record the module 
  ppi.comps[[module.nm]] <<- g;
  UpdateSubnetStats();
  
  module.count <<- module.count
  
  convertIgraph2JSON(module.nm, filenm);
  return (filenm);
}

PerformLayOut <- function(net.nm, algo, focus=""){
    if(!exists("my.perform.layout")){ # public web on same user dir
        compiler::loadcmp("../../rscripts/OmicsAnalystR/R/util_graph.Rc");    
    }
    return(my.perform.layout(net.nm, algo, focus));
}

UpdateNetworkLayout <- function(algo, filenm, focus=""){
  current.net <- ppi.comps[[current.net.nm]];
  pos.xy <- PerformLayOut(current.net.nm, algo, focus);
  nms <- V(current.net)$name;
  nodes <- vector(mode="list");
  if(algo %in% c("fr", "kk")){
    for(i in 1:length(nms)){
      nodes[[i]] <- list(
        id=nms[i],
        x=pos.xy[i,1],
        y=pos.xy[i,2],
        z=pos.xy[i,3]
      );
    }
  }else{
    for(i in 1:length(nms)){
      nodes[[i]] <- list(
        id=nms[i], 
        x=pos.xy[i,1], 
        y=pos.xy[i,2]
      );
    }
  }
  # now only save the node pos to json
  require(rjson);
  netData <- list(nodes=nodes);
  sink(filenm);
  cat(rjson::toJSON(netData));
  sink();
  return(filenm);
}

GetNetsNameString <- function(){
  nms <- paste(rownames(net.stats), collapse="||");
  return(nms);
}

# support walktrap, infomap and lab propagation
FindCommunities <- function(method="walktrap", use.weight=FALSE){
  # make sure this is the connected
  current.net <- ppi.comps[[current.net.nm]];
  g <- current.net;
  if(!is.connected(g)){
    g <- decompose.graph(current.net, min.vertices=2)[[1]];
  }
  total.size <- length(V(g));
  
  if(use.weight){ # this is only tested for walktrap, should work for other method
    # now need to compute weights for edges
    egs <- get.edges(g, E(g)); #node inx
    nodes <- V(g)$name;
    # conver to node id
    negs <- cbind(nodes[egs[,1]],nodes[egs[,2]]);
    
    # get min FC change
    base.wt <- min(abs(seed.expr))/10;
    
    # check if user only give a gene list without logFC or all same fake value
    if(length(unique(seed.expr)) == 1){
      seed.expr <- rep(1, nrow(negs));
      base.wt <- 0.1; # weight cannot be 0 in walktrap
    }
    
    wts <- matrix(base.wt, ncol=2, nrow = nrow(negs));
    for(i in 1:ncol(negs)){
      nd.ids <- negs[,i];
      hit.inx <- match(names(seed.expr), nd.ids);
      pos.inx <- hit.inx[!is.na(hit.inx)];
      wts[pos.inx,i]<- seed.expr[!is.na(hit.inx)]+0.1;
    }
    nwt <- apply(abs(wts), 1, function(x){mean(x)^2})    
  }
  
  if(method == "walktrap"){
    fc <- walktrap.community(g);
  }else if(method == "infomap"){
    fc <- infomap.community(g);
  }else if(method == "labelprop"){
    fc <- label.propagation.community(g);
  }else{
    print(paste("Unknown method:", method));
    return ("NA||Unknown method!");
  }
  
  if(length(fc) == 0 || modularity(fc) == 0){
    return ("NA||No communities were detected!");
  }
  
  # only get communities
  communities <- communities(fc);
  community.vec <- vector(mode="character", length=length(communities));
  gene.community <- NULL;
  qnum.vec <- NULL;
  pval.vec <- NULL;
  rowcount <- 0;
  nms <- V(g)$name;
  hit.inx <- match(nms, ppi.net$node.data[,1]);
  sybls <- ppi.net$node.data[hit.inx,2];
  names(sybls) <- V(g)$name;
  for(i in 1:length(communities)){
    # update for igraph 1.0.1 
    path.ids <- communities[[i]];
    psize <- length(path.ids);
    if(psize < 5){
      next; # ignore very small community
    }

    hits <- seed.proteins %in% path.ids;
    qnums <- sum(hits);
    
    if(qnums == 0){
      next; # ignor community containing no queries
    }
    
    rowcount <- rowcount + 1;
    pids <- paste(path.ids, collapse="->");
    #path.sybls <- V(g)$Label[path.inx];
    path.sybls <- sybls[path.ids];
    com.mat <- cbind(path.ids, path.sybls, rep(i, length(path.ids)));
    gene.community <- rbind(gene.community, com.mat);
    qnum.vec <- c(qnum.vec, qnums);
    
    # calculate p values (comparing in- out- degrees)
    #subgraph <- induced.subgraph(g, path.inx);
    subgraph <- induced.subgraph(g, path.ids);
    in.degrees <- degree(subgraph);
    #out.degrees <- degree(g, path.inx) - in.degrees;
    out.degrees <- degree(g, path.ids) - in.degrees;
    ppval <- wilcox.test(in.degrees, out.degrees)$p.value;
    ppval <- signif(ppval, 3);
    pval.vec <- c(pval.vec, ppval);
    
    # calculate community score
    community.vec[rowcount] <- paste(c(psize, qnums, ppval, pids), collapse=";");
  }
  pvall <<- pval.vec
  if(length(pval.vec)>1){
    ord.inx <- order(pval.vec, decreasing=F);
    community.vec <- community.vec[ord.inx];
    qnum.vec <- qnum.vec[ord.inx];
    ord.inx <- order(qnum.vec, decreasing=T);
    community.vec <- community.vec[ord.inx];
  }
  
  all.communities <- paste(community.vec, collapse="||");
  if(!is.null(gene.community)){
    colnames(gene.community) <- c("Id", "Label", "Module");
    fast.write(gene.community, file="module_table.csv", row.names=F);
    return(all.communities);
  }else{
    return("NA");
  }
  
}

community.significance.test <- function(graph, vs, ...) {
  subgraph <- induced.subgraph(graph, vs)
  in.degrees <- degree(subgraph)
  out.degrees <- degree(graph, vs) - in.degrees
  wilcox.test(in.degrees, out.degrees, ...)
}


convertIgraph2JSON <- function(net.nm, filenm, idType="NA"){

    if(!exists("my.convert.igraph")){ # public web on same user dir
        compiler::loadcmp("../../rscripts/OmicsAnalystR/R/util_graph.Rc");    
    }
   
    return(my.convert.igraph(net.nm, filenm, idType));
}


# from to should be valid nodeIDs
GetShortestPaths <- function(from, to){
  current.net <- ppi.comps[[current.net.nm]];
  paths <- get.all.shortest.paths(current.net, from, to)$res;
  if(length(paths) == 0){
    return (paste("No connection between the two nodes!"));
  }
  
  path.vec <- vector(mode="character", length=length(paths));
  for(i in 1:length(paths)){
    path.inx <- paths[[i]]; 
    path.ids <- V(current.net)$name[path.inx];
    path.sybls <- path.ids;
    pids <- paste(path.ids, collapse="->");
    psbls <- paste(path.sybls, collapse="->");
    path.vec[i] <- paste(c(pids, psbls), collapse=";")
  }
  
  if(length(path.vec) > 50){
    path.vec <- path.vec[1:50];
  }
  
  all.paths <- paste(path.vec, collapse="||");
  return(all.paths);
}

GetNetsName <- function(){
  rownames(net.stats);
}

GetNetsEdgeNum <- function(){
  as.numeric(net.stats$Edge);
}

GetNetsNodeNum <- function(){
  as.character(net.stats$Node);
}

GetNetsQueryNum <- function(){
  as.numeric(net.stats$Query);
}

ProcessGraphFile <- function(graph=new_g, labels, typeList=type.list, generateJson = F){  
  overall.graph <<- graph
  nms <- V(graph)$name;
  if(length(nms)<1){
    nms <- V(graph)$id;
    graph = set_vertex_attr(graph, "name", value=nms)
  }
  lblsNm <- names(labels)
  names(lblsNm) <- unname(labels)
  lbls <- unname(lblsNm[V(graph)$featureId]);
  node.data = data.frame(nms, lbls);
  graph = set_vertex_attr(graph, "label", value=lbls)
  seed.proteins <<- nms;
  
  if(!is.null(typeList) && is.null(V(graph)$type)){
    typeVec <- rep("NA", length(nms))
    inx.list <- list();
    if(!is.null(typeList)){
      for( i in 1:length(typeList)){
        inx <- nms %in% typeList[[i]]
        typeVec[inx] <- names(typeList)[i]
      }
    }
    graph = set_vertex_attr(graph, "type", value=typeVec)
  }
  
  
  seed.genes <<- seed.proteins;
  e=get.edgelist(graph)
  edge.data= data.frame(Source=e[,1], Target=e[,2])
  
  seed.expr <<- rep(0, length(node.data));
  substats <- DecomposeGraph(graph);
  
  if(is.null(substats)){
    msg.vec <<- "No subnetworks containing at least 3 edges are identified"
    return(0);
  }
  
  net.nm <- names(ppi.comps)[1];
  net.nmu <<- net.nm;
  current.net.nm <<- net.nm;
  #ppi.comps[["overall"]] <- graph
  ppi.comps <<- ppi.comps
  g <- ppi.comps[[net.nm]];
  ppi.net <<- list(db.type="abc",
                   db.type="ppi", 
                   order=1, 
                   seeds=nms, 
                   table.nm=" ", 
                   node.data = node.data,
                   edge.data = edge.data
  );
  data.idType <<- "NA"; 
  if(generateJson){
    convertIgraph2JSON(current.net.nm , "omicsanalyst_net_0.json");
  }
  return(1);
}
