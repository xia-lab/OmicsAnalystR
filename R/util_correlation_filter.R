my.correlation.filter <- function(corSign="both", crossOmicsOnly="false", networkInfer="NA", threshold.inter=0.5, 
                                threshold.intra=0.9, numToKeep=2000, updateRes="false", taxlvl="genus", datagem="agora"){

  reductionSet <- .get.rdt.set();
  reductionSet$threshold.inter <- threshold.inter
  reductionSet$threshold.intra <- threshold.intra
  reductionSet$crossOmicsOnly <- crossOmicsOnly;
   load_igraph();
  if(!(exists("selDatsCorr.taxa",reductionSet))){
    sel.nms <- names(mdata.all)[mdata.all == 1];
    dataSetList <- lapply(sel.nms, readDataset);
    labels <- unlist(lapply(dataSetList, function(x) x$enrich_ids))
    types <- unlist(lapply(dataSetList, function(x) rep(x$type, length(x$enrich_ids))))
    type_df <- data.frame(name=c(labels),
                          type=c(types))
    
    # Create a lookup table
    type_lookup <- setNames(type_df$type, type_df$name)
   
    corr.mat <- qs::qread(reductionSet$corr.mat.path);
    corr.p.mat<- qs::qread("corr.p.mat.qs");
    corr.p.mat <- reshape2::melt(corr.p.mat) 
    sel.dats <- reductionSet$selDatsCorr;
    rowlen <- nrow(corr.mat);

    g <- igraph::graph_from_adjacency_matrix(corr.mat,mode = "undirected", diag = FALSE, weighted = 'correlation')
 
   # Assign types to nodes using the lookup table
    toMatch <- unlist(lapply(dataSetList, function(x) x$type));
    pattern <- paste0("^(.*)((", paste(toMatch, collapse = "|"), "))$")
    V(g)$type <- gsub(pattern, "\\2", V(g)$name);
    V(g)$label <- gsub(pattern, "", V(g)$name);
    
    edge_list <- igraph::as_data_frame(simplify(g, remove.loops = TRUE,edge.attr.comb="max"), "edges")
 
    v1 <- V(g)$name[V(g)$type == unique(types)[1]]
    v2 <- V(g)$name[V(g)$type == unique(types)[2]]
    # Get the edges that connect nodes of different types
    inter_inx <- V(g)[ends(g, E(g))[, 1]]$type != V(g)[ends(g, E(g))[, 2]]$type;
    intra_inx <- V(g)[ends(g, E(g))[, 1]]$type == V(g)[ends(g, E(g))[, 2]]$type;
    inter_g <- delete_edges(g, E(g)[intra_inx])
    intra_g <- delete_edges(g, E(g)[inter_inx])
   
    # for histogram only
    reductionSet$corr.graph.path <- "corr.graph.qs";
    qs::qsave(list(corr.graph.inter=inter_g, corr.graph.intra=intra_g), "corr.graph.qs")
    
    cor.list <- list(all = NULL, inter = NULL, intra = NULL)
    
    if (corSign == "both") {
      toRm.inter <- E(inter_g)[!abs(correlation) >  threshold.inter]
      toRm.intra <- E(intra_g)[!abs(correlation) >  threshold.intra]
    } else if (corSign == "positive") {
      toRm.inter <- E(inter_g)[!correlation >  threshold.inter]
      toRm.intra <- E(intra_g)[!correlation >  threshold.intra]
    } else {
      toRm.inter <- E(inter_g)[!correlation <  -threshold.inter]
      toRm.intra <- E(intra_g)[!correlation <  -threshold.intra]
    }
    
    inter_g_sub <- delete_edges(inter_g, E(inter_g)[toRm.inter])
    intra_g_sub <- delete_edges(intra_g, E(intra_g)[toRm.intra])
    
    edge_list_inter <- get.edgelist(inter_g_sub)
    edge_list_intra <- get.edgelist(intra_g_sub)
    
    # add correlation weights
    cor.list$inter <- data.frame(edge_list_inter,  as.numeric(E(inter_g_sub)$correlation))
    cor.list$intra <- data.frame(edge_list_intra,  as.numeric(E(intra_g_sub)$correlation))
    colnames(cor.list$inter) <- c("source", "target","correlation");
    colnames(cor.list$intra) <- c("source", "target","correlation");
    cor.list$all  <- rbind(cor.list$inter, cor.list$intra)

    reductionSet$cor.list.path <- "cor.list.qs";
    qs::qsave(cor.list, file="cor.list.qs");
    if (crossOmicsOnly == "true") {
      cor_edge_list <- cor.list$inter
    } else {
      cor_edge_list <- dplyr::bind_rows(cor.list$inter, cor.list$intra)
    }
    
    # filter interomics only
    numToKeep <- as.numeric(numToKeep)
    if (numToKeep > length(unique(cor_edge_list$correlation))) {
      numToKeep <- length(unique(cor_edge_list$correlation))
    }

    if (nrow(cor_edge_list) >= 3) {
      top.edge <- sort(abs(unique(cor_edge_list$correlation)))[c(1:numToKeep)];
      top.inx <- match(abs(cor_edge_list$correlation), top.edge);
      cor_edge_list <- cor_edge_list[!is.na(top.inx), ,drop=F];
      new_g <- igraph::graph_from_data_frame(cor_edge_list, directed = FALSE)
      new_g <- igraph::simplify(new_g, edge.attr.comb = "mean")
      
      pattern <- paste0("^(.*)((", paste(toMatch, collapse = "|"), "))$")
      V(new_g)$type <- gsub(pattern, "\\2", V(new_g)$name);
      toMatch2 <- unlist(lapply(dataSetList, function(x) paste0("_", x$type)));
      pattern2 <- paste0("(", paste(toMatch2, collapse = "|"), ")$")
      V(new_g)$featureId <- gsub(pattern2, "", V(new_g)$name);
      type.list <- list();
      for(i in 1:length(sel.nms)){
        type.list[[sel.nms[[i]]]] <- unique(cor_edge_list[,i]);
      }

      cor_edge_list <- cor_edge_list %>%
          dplyr::left_join(corr.p.mat, by = c("source" = "Var1", "target" = "Var2")) %>%
           dplyr::mutate(pval = value)
      cor_edge_list$value =NULL 

      cor_edge_list$label1 <- gsub(paste(paste0("_",toMatch,"$"), collapse = "|"), "", cor_edge_list$source )
      cor_edge_list$label2 <- gsub(paste(paste0("_",toMatch,"$"), collapse = "|"), "", cor_edge_list$target )
      cor_edge_list$label1 <- names(reductionSet$labels)[match( cor_edge_list$label1,reductionSet$labels)]
      cor_edge_list$label2 <- names(reductionSet$labels)[match( cor_edge_list$label2,reductionSet$labels)]
       reductionSet$corNet <- cor_edge_list;
      write.csv(corr.mat,"corNet.csv",row.names=F)

      reductionSet$taxlvl <-"Feature"
      .set.rdt.set(reductionSet);
      intres <- ProcessGraphFile(new_g, labels, type.list);
      return(intres);
    } else {
      msg.vec <<- "Less than 3 correlations have been identified using current parameters. Failed to create correlation network";
      .set.rdt.set(reductionSet);
      return(0)
    }
    
  } else  {
    
    taxlvl <- gsub("(^[[:alpha:]])", "\\U\\1", taxlvl, perl=TRUE)
    
    if(exists("corr.mat.inter.taxa",reductionSet) & taxlvl==reductionSet$taxlvl){
      cor_edge_list_inter <- reductionSet$corr.mat.inter.taxa[[taxlvl]] 
      cor_edge_list_intra <- reductionSet$corr.mat.intra.taxa[[taxlvl]] 
      
    }else{
 
      sel.inx <- mdata.all==1;
      sel.nms <- names(mdata.all)[sel.inx];
      micidx <- reductionSet$micidx
      residx <- reductionSet$residx
      
      dataSet <- readDataset(sel.nms[micidx]);
      
      labels <- unique(dataSet$taxa_table[,taxlvl])
      labels <- setNames(labels, labels)
      dataSet <- readDataset(sel.nms[residx]);
      labels <- c(labels, dataSet$enrich_ids);
      
      corr.mat <- reductionSet$corr.mat.taxa[[taxlvl]]
      sel.dats <- list()
      sel.dats[[micidx]] <- reductionSet$selDatsCorr.taxa[[taxlvl]];
      sel.dats[[residx]] <- reductionSet$selDatsCorr[[residx]];
      if(networkInfer != "NA"){
        library(parmigene);
        if(networkInfer == "aracne"){
          corr.mat = aracne.m(corr.mat);
        }else if (networkInfer == "clr"){
          corr.mat = clr(corr.mat);
        }else if (networkInfer == "mrnet"){
          corr.mat = mrnet(corr.mat);
        }
      }
      
      rowlen <- nrow(corr.mat);
      corr.mat.intra <- corr.mat[-c(which(rownames(corr.mat) %in% rownames(sel.dats[[residx]]))),-c(which(colnames(corr.mat) %in% rownames(sel.dats[[micidx]])))]; #between-omics
      corr.mat.inter1 <- corr.mat[c(which(rownames(corr.mat) %in% rownames(sel.dats[[1]]))),c(which(colnames(corr.mat) %in% rownames(sel.dats[[1]])))]; #within-omics
      corr.mat.inter2 <- corr.mat[c(which(rownames(corr.mat) %in% rownames(sel.dats[[2]]))),c(which(colnames(corr.mat) %in% rownames(sel.dats[[2]])))]; #within-omics
      
      
      cor_g_inter1 <- graph_from_incidence_matrix(corr.mat.inter1 , directed="false", weighted = 'correlation');
      cor_g_inter2 <- graph_from_incidence_matrix(corr.mat.inter2 , directed="false", weighted = 'correlation');
      cor_g_intra <- graph_from_incidence_matrix(corr.mat.intra , directed="false", weighted = 'correlation');  
      
      cor_edge_list_inter1 <- igraph:::as_data_frame(cor_g_inter1, 'edges');
      cor_edge_list_inter1 <- cor_edge_list_inter1[!is.na(cor_edge_list_inter1$correlation), ]  
      cor_edge_list_inter2 <- igraph:::as_data_frame(cor_g_inter2, 'edges');
      cor_edge_list_inter2 <- cor_edge_list_inter2[!is.na(cor_edge_list_inter2$correlation), ]     
      cor_edge_list_inter <- rbind(cor_edge_list_inter1, cor_edge_list_inter2)
      
      cor_edge_list_intra <- igraph:::as_data_frame(cor_g_intra, 'edges');
      cor_edge_list_intra <- cor_edge_list_intra[!is.na(cor_edge_list_intra$correlation), ]  
      
      # filter self loops
      cor_edge_list_inter <- cor_edge_list_inter[cor_edge_list_inter$correlation != 1, ]
      cor_edge_list_intra <- cor_edge_list_intra[cor_edge_list_intra$correlation != 1, ]
      cor.list <- list();
      cor.list[["all"]] <- rbind(cor_edge_list_inter, cor_edge_list_intra)
      cor.list[["inter"]] <- cor_edge_list_inter
      cor.list[["intra"]] <- cor_edge_list_intra
      qs::qsave(cor.list, file= paste0("cor.list",taxlvl,".qs"));
      
      if(!(exists("corr.mat.inter.taxa",reductionSet))){
        reductionSet$corr.mat.inter.taxa <- list()
        reductionSet$corr.mat.intra.taxa <- list()
        
      }
      
      reductionSet$corr.mat.inter.taxa[[taxlvl]] <- cor_edge_list_inter[cor_edge_list_inter$correlation != 1, ]
      reductionSet$corr.mat.intra.taxa[[taxlvl]] <- cor_edge_list_intra[cor_edge_list_intra$correlation != 1, ]
      
    }
    
    
    if(corSign == "both"){
      cor.inx.inter <- abs(cor_edge_list_inter$correlation) > threshold.inter
      cor.inx.intra <- abs(cor_edge_list_intra$correlation) > threshold.intra
    }else if(corSign == "positive"){
      cor.inx.inter <- cor_edge_list_inter$correlation > threshold.inter
      cor.inx.intra <- cor_edge_list_intra$correlation > threshold.intra
    }else{
      cor.inx.inter <- cor_edge_list_inter$correlation < -threshold.inter
      cor.inx.intra <- cor_edge_list_intra$correlation < -threshold.intra
    }
    
    
    cor_edge_list_inter <- cor_edge_list_inter[cor.inx.inter, ];
    cor_edge_list_intra <- cor_edge_list_intra[cor.inx.intra, ];
    
    cor_edge_list <- cor_edge_list_inter
    
    
    if (crossOmicsOnly == "true") {
      cor_edge_list <- cor_edge_list_inter
    } else {
      cor_edge_list <- dplyr::bind_rows(cor_edge_list_inter, cor_edge_list_intra)
    }
    #filter interomics only
    numToKeep <- as.numeric(numToKeep)
    if(numToKeep > length(unique(cor_edge_list$correlation)) ){
      numToKeep <-length(unique(cor_edge_list$correlation))
    }
    
    top.edge <- sort(abs(unique(cor_edge_list$correlation)))[c(1:numToKeep)]; #default only show top 20% significant edges when #edges<1000
    top.inx <- match(abs(cor_edge_list$correlation), top.edge);
    cor_edge_list <- cor_edge_list[!is.na(top.inx), ,drop=F];
    
    
    new_g <- graph_from_data_frame(cor_edge_list, F);
    new_g <-simplify(new_g, edge.attr.comb="mean")
    
    if(nrow(cor_edge_list)<3){
      msg.vec <<- paste0("Less than 3 correlations have been identified using an inter-omics correlation threshold of ", threshold.inter, 
                         " and intra-omics correlation threshold of ", threshold.intra);
    }
    
    new_g <- graph_from_data_frame(cor_edge_list, F);
    new_g <-simplify(new_g, edge.attr.comb="mean")
    type.list <- list();
    for(i in 1:length(sel.nms)){
      type.list[[sel.nms[[i]]]] <- unique(cor_edge_list[,i]);
    }
    
    reductionSet$taxlvl <- taxlvl
    reductionSet$datagem <- datagem
    
    qvec <- unique(rownames(reductionSet$selDatsCorr.taxa[["Genus"]]))
    mvec <- unique(rownames(reductionSet[["selDatsCorr"]][[residx]]))
    m2mscore<-M2Mscore(qvec,mvec,taxlvl,dataGem=datagem)
    toHighlight <- merge(m2mscore,cor_edge_list_intra,all.x = TRUE)
    toHighlight <- unique(toHighlight[!(is.na(toHighlight$correlation)),])
    toHighlight <- toHighlight[which(toHighlight$correlation>0),]
    reductionSet$m2mscore<-m2mscore
    reductionSet$toHighlight <- toHighlight
    .set.rdt.set(reductionSet);

    
    intres <- ProcessGraphFile(new_g, labels, type.list);
    return(intres)
    
  }  
}
