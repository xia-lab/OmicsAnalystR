##################################################
## R scripts for NetworkAnalyst
## Description: biological network analysis methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################


FilterBipartiNet <- function(nd.type="all", min.dgr, min.btw){

    all.nms <- V(overall.graph)$name;
    edge.mat <- as_edgelist(overall.graph);
    dgrs <- degree(overall.graph);
    nodes2rm.dgr <- nodes2rm.btw <- NULL;

    if(nd.type == "gene"){
        hit.inx <- all.nms %in% edge.mat[,1];
    }else if(nd.type=="other"){
        hit.inx <- all.nms %in% edge.mat[,2];
    }else{ # all
        hit.inx <- rep(TRUE, length(all.nms));
    }

    if(min.dgr > 0){
        rm.inx <- dgrs <= min.dgr & hit.inx;
        nodes2rm.dgr <- V(overall.graph)$name[rm.inx];
    }
    if(min.btw > 0){
        btws <- betweenness(overall.graph);
        rm.inx <- btws <= min.btw & hit.inx;
        nodes2rm.btw <- V(overall.graph)$name[rm.inx];
    }

    nodes2rm <- unique(c(nodes2rm.dgr, nodes2rm.btw));
    overall.graph <- simplify(delete.vertices(overall.graph, nodes2rm));
    current.msg <<- paste("A total of", length(nodes2rm) , "was reduced.");
    substats <- DecomposeGraph(overall.graph);
    if(!is.null(substats)){
        overall.graph <<- overall.graph;
        return(c(length(seed.genes),length(seed.proteins), vcount(overall.graph), ecount(overall.graph), length(ppi.comps), substats));
    }else{
        return(0);
    }
}

PrepareNetwork <- function(net.nm, json.nm){
   my.ppi <- ppi.comps[[net.nm]];
   nd.nms <- V(my.ppi)$name;


   convertIgraph2JSON(net.nm, json.nm);
   current.net.nm <<- net.nm;
   return(1);
}

GetNodeIDs <- function(){
    V(overall.graph)$name;
}

GetNodeNames <- function(){
    V(overall.graph)$Label;
}

GetNodeDegrees <- function(){
    degree(overall.graph);
}

GetNodeBetweenness <- function(){
    round(betweenness(overall.graph, directed=F, normalized=F), 2);
}

# for a given graph, obtain the smallest subgraphs that contain
# all the seed nodes. This is acheived by iteratively remove
# the marginal nodes (degree = 1) that are not in the seeds
GetMinConnectedGraphs <- function(max.len = 200){
    set.seed(8574);
    # first get shortest paths for all pair-wise seeds
    my.seeds <- seed.proteins;
    sd.len <- length(my.seeds);
    paths.list <-list();

    # first trim overall.graph to remove no-seed nodes of degree 1
    dgrs <- degree(overall.graph);
    keep.inx <- dgrs > 1 | (names(dgrs) %in% my.seeds);
    nodes2rm <- V(overall.graph)$name[!keep.inx];
    overall.graph <-  simplify(delete.vertices(overall.graph, nodes2rm));

    # need to restrict the operation b/c get.shortest.paths is very time consuming
    # for top max.len highest degrees
    if(sd.len > max.len){
        hit.inx <- names(dgrs) %in% my.seeds;
        sd.dgrs <- dgrs[hit.inx];
        sd.dgrs <- rev(sort(sd.dgrs));
        # need to synchronize all (seed.proteins) and top seeds (my.seeds)
        seed.proteins <- names(sd.dgrs);
        if(max.len>table(hit.inx)[["TRUE"]]){
            sd.len <-  table(hit.inx)[["TRUE"]];
        }else{
            sd.len <-  max.len;
        }
        my.seeds <- seed.proteins[1:sd.len];
        current.msg <<- paste("The minimum connected network was computed using the top", sd.len, "seed proteins in the network based on their degrees.");
    }else{
        current.msg <<- paste("The minimum connected network was computed using all seed proteins in the network.");
    }
    # now calculate the shortest paths for
    # each seed vs. all other seeds (note, to remove pairs already calculated previously)
    for(pos in 1:sd.len){
            paths.list[[pos]] <- get.shortest.paths(overall.graph, my.seeds[pos], seed.proteins[-(1:pos)])$vpath;
    }
    nds.inxs <- unique(unlist(paths.list));
    nodes2rm <- V(overall.graph)$name[-nds.inxs];
    g <- simplify(delete.vertices(overall.graph, nodes2rm));

    nodeList <- get.data.frame(g, "vertices");
    colnames(nodeList) <- c("Id", "Label");
    fast.write(nodeList, file="orig_node_list.csv", row.names=F);

    edgeList <- get.data.frame(g, "edges");
    edgeList <- cbind(rownames(edgeList), edgeList);
    colnames(edgeList) <- c("Id", "Source", "Target");
    fast.write(edgeList, file="orig_edge_list.csv", row.names=F);

    path.list <- NULL;
    substats <- DecomposeGraph(g);
    net.stats<<-net.stats
    if(!is.null(substats)){
        overall.graph <<- g;
        return(c(length(seed.genes),length(seed.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), substats));
    }else{
        return(0);
    }
}

UpdateSubnetStats <- function(){
    old.nms <- names(ppi.comps);
    net.stats <- ComputeSubnetStats(ppi.comps);
    ord.inx <- order(as.numeric(net.stats[,2]), decreasing=TRUE);
    net.stats <- net.stats[ord.inx,];
    rownames(net.stats) <- old.nms[ord.inx];
    net.stats <<- net.stats;
}

# exclude nodes in current.net (networkview)
ExcludeNodes <- function(nodeids, filenm){

    nodes2rm <- strsplit(nodeids, ";")[[1]];
    current.net <- ppi.comps[[current.net.nm]];
    current.net <- delete.vertices(current.net, nodes2rm);

    # need to remove all orphan nodes
    bad.vs<-V(current.net)$name[degree(current.net) == 0];
    current.net <- delete.vertices(current.net, bad.vs);

    # return all those nodes that are removed
    nds2rm <- paste(c(bad.vs, nodes2rm), collapse="||");

    # update topo measures
    node.btw <- as.numeric(betweenness(current.net));
    node.dgr <- as.numeric(degree(current.net));
    node.exp <- as.numeric(get.vertex.attribute(current.net, name="abundance", index = V(current.net)));
    nms <- V(current.net)$name;
    hit.inx <- match(nms, ppi.net$node.data[,1]);
    lbls <- ppi.net$node.data[hit.inx,2];

    nodes <- vector(mode="list");
    for(i in 1:length(nms)){
        nodes[[i]] <- list(
                  id=nms[i],
                  label=lbls[i],
                  degree=node.dgr[i],
                  between=node.btw[i],
                  expr = node.exp[i]
                );
    }
    # now only save the node pos to json
    require(RJSONIO);
    netData <- list(deletes=nds2rm,nodes=nodes);
    sink(filenm);
    cat(toJSON(netData));
    sink();

    ppi.comps[[current.net.nm]] <<- current.net;
    UpdateSubnetStats();

    return(filenm);
}

# exclude nodes in overall net (network builder)
ExcludeNodesOverall <- function(nodeids, id.type, vismode){
    # all convert to uniprot ID
    lines <- strsplit(nodeids, "\r|\n|\r\n")[[1]];
    lines<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", lines, perl=TRUE);
    if(vismode != "network"){
        prot.anots <- convertIdToEntrez(lines, id.type);
        nodes2rm <- unique(prot.anots$accession);
    }else{
        prot.anots <- lines
        nodes2rm <- unique(lines);
    }

    # now find the overlap
    nodes2rm <- nodes2rm[nodes2rm %in% V(overall.graph)$name];
    g <- delete.vertices(overall.graph, nodes2rm);

    nodeList <- get.data.frame(g, "vertices");
    colnames(nodeList) <- c("Id", "Label");

    fast.write(nodeList, file="orig_node_list.csv", row.names=F);

    edgeList <- get.data.frame(g, "edges");
    edgeList <- cbind(rownames(edgeList), edgeList);
    colnames(edgeList) <- c("Id", "Source", "Target");
    fast.write(edgeList, file="orig_edge_list.csv", row.names=F);

    substats <- DecomposeGraph(g);
    if(!is.null(substats)){
        overall.graph <<- g;
        return(c(length(seed.genes),length(seed.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), substats));
    }else{
        return(0);
    }
}

PrepareSubnetDownloads <- function(nm){
  g <- ppi.comps[[nm]];
  # need to update graph so that id is compound names rather than ID
  V(g)$name <- as.character(doID2LabelMapping(V(g)$name));
  saveNetworkInSIF(g, nm);
}

# adapted from BioNet
saveNetworkInSIF <- function(network, name){
    edges <- .graph.sif(network=network, file=name);
    sif.nm <- paste(name, ".sif", sep="");
    if(length(list.edge.attributes(network))!=0){
	edge.nms <- .graph.eda(network=network, file=name, edgelist.names=edges);
        sif.nm <- c(sif.nm, edge.nms);

    }
    if(length(list.vertex.attributes(network))!=0){
	node.nms <- .graph.noa(network=network, file=name);
        sif.nm <- c(sif.nm, node.nms);
    }
    # need to save all sif and associated attribute files into a zip file for download
    zip(paste(name,"_sif",".zip", sep=""), sif.nm);
}

# internal function to write cytoscape .sif file
.graph.sif <- function(network, file){
    edgelist.names <- igraph::as_edgelist(network, names=TRUE)
    edgelist.names <- cbind(edgelist.names[,1], rep("pp", length(E(network))), edgelist.names[,2]);
    write.table(edgelist.names, row.names=FALSE, col.names=FALSE, file=paste(file, ".sif", sep=""), sep="\t", quote=FALSE)
    return(edgelist.names)
}

doID2LabelMapping <- function(entrez.vec){
    if(exists("nodeListu")){
        hit.inx <- match(entrez.vec, nodeListu[, "Id"]);
        symbols <- nodeListu[hit.inx, "Label"];

        # if not gene symbol, use id by itself
        na.inx <- is.na(symbols);
        symbols[na.inx] <- entrez.vec[na.inx];
        return(symbols);
    }else{ # network upload
        return(entrez.vec);
    }
}

# internal method to write cytoscape node attribute files
.graph.noa <- function(network, file){
  all.nms <- c();
  attrib <- list.vertex.attributes(network)
  for(i in 1:length(attrib)){
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "character")
    {
      type <- "String"
    }
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "integer")
    {
      type <- "Integer"
    }
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "numeric")
    {
      type <- "Double"
    }
    noa <- cbind(V(network)$name, rep("=", length(V(network))), get.vertex.attribute(network, attrib[i]))
    first.line <- paste(attrib[i], " (class=java.lang.", type, ")", sep="")
    file.nm <- paste(file, "_", attrib[i], ".NA", sep="");
    write(first.line, file=file.nm, ncolumns = 1, append=FALSE, sep=" ")
    write.table(noa, row.names = FALSE, col.names = FALSE, file=file.nm, sep=" ", append=TRUE, quote=FALSE);
    all.nms <- c(all.nms, file.nm);
  }
  return(all.nms);
}

# internal method to write cytoscape edge attribute files
.graph.eda <- function(network, file, edgelist.names){
  all.nms <- c();
  attrib <- list.edge.attributes(network)
  for(i in 1:length(attrib)){
    if(is(get.edge.attribute(network, attrib[i]))[1] == "character")
    {
      type <- "String"
    }
    if(is(get.edge.attribute(network, attrib[i]))[1] == "integer")
    {
      type <- "Integer"
    }
    if(is(get.edge.attribute(network, attrib[i]))[1] == "numeric")
    {
      type <- "Double"
    }
    eda <- cbind(cbind(edgelist.names[,1], rep("(pp)", length(E(network))), edgelist.names[,3]), rep("=", length(E(network))), get.edge.attribute(network, attrib[i]))
    first.line <- paste(attrib[i], " (class=java.lang.", type, ")", sep="");
    file.nm <- paste(file, "_", attrib[i], ".EA", sep="");
    write(first.line, file=file.nm, ncolumns=1, append=FALSE, sep =" ")
    write.table(eda, row.names = FALSE, col.names = FALSE, file=file.nm, sep=" ", append=TRUE, quote=FALSE);
    all.nms <- c(all.nms, file.nm);
  }
  return(all.nms);
}

SetCellCoexNumber <-function(num){
    cellCoexNumber <<- num;
}

SetDiffNetName <-function(nm){
    diffNetName <<- nm;
}

SetDiffFilter <-function(pct){
    diffPct <<- pct/10;
}

ResetNetworkSelection <- function(){
    dataSet$nets <- list();
    RegisterData(dataSet)
}

GetSelectedNetworks <- function(){
    return(names(dataSet$nets));
}

ClearNetSelection <- function(type){
    if(type %in% names(dataSet$nets)){
        inx = which( names(dataSet$nets) == type)
        #print(inx);
         dataSet$nets[[inx]] = NULL;
    }
    RegisterData(dataSet)
    return(1)
}

# === libraries
doUniprot2EntrezMapping<-function(uniprot.vec){
  db.map <-  queryGeneDB("entrez_uniprot", data.org);
  hit.inx <- match(uniprot.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  na.inx <- is.na(entrezs);
  entrezs[na.inx] <- uniprot.vec[na.inx];
  return(entrezs);
}

doEmblGene2EntrezMapping <- function(emblgene.vec){
  # db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
  # db.map <-  readRDS(db.path);
  db.map <-  queryGeneDB("entrez_embl_gene", data.org);
  hit.inx <- match(emblgene.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doEmblProtein2EntrezMapping <- function(emblprotein.vec){
  # db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
  # db.map <-  readRDS(db.path);
  db.map <-  queryGeneDB("entrez_embl_protein", data.org);
  hit.inx <- match(emblprotein.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}
