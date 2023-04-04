##################################################
## R scripts for OmicsAnalyst
## Description: Related to correlation analysis
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#default feature selection based on sig genes
DoFeatSelectionForCorr <- function(type="default", retainedNumber=20, retainedComp=3){

  sel.dats <- list();
  labels <- vector();
  reductionSet <- .get.rdt.set()
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  
  
  if(type %in% c("default","custom")){
    for(i in 1:length(sel.nms)){
      nm = sel.nms[i]
      
      dataSet <- qs::qread(nm);
      
      if(i==1){
        all.mat <- dataSet$data.proc
      }else{
        all.mat <- rbind(all.mat, dataSet$data.proc)
      }
      
      if(type == "default"){
        sig.mat <- dataSet$sig.mat
      }else{
        sig.mat <- dataSet$custom.sig.mat
      }
      
      inxAll = which(rownames(all.mat) %in% rownames(sig.mat));
      inx = which(rownames(dataSet$data.proc) %in% rownames(sig.mat));
      
      all.mat <- all.mat[inxAll, ]
      dat <- dataSet$data.proc[inx, ]
      sel.dats[[i]] <- dat
      labels <- c(labels, rownames(dat));
      
      
      if(exists("m2m",dataSet)){ 
        
        all.mat.taxa <- dataSet$data.proc.taxa
        
        if(type == "default"){
          sig.mat.taxa <- dataSet$sig.mat.tax
        }else{
          sig.mat <- dataSet$custom.sig.mat
        }
        inxAll <- list()
        inx <- list()
        sel.dats.taxa <- list()
        labels.taxa <- list()
        for(j in 1:length(all.mat.taxa)){
          inxAll[[j]] = which(rownames(all.mat.taxa[[j]]) %in% rownames(sig.mat.taxa[[j]]));
          inx[[j]] = which(rownames(dataSet$data.proc.taxa[[j]]) %in% rownames(sig.mat.taxa[[j]]));
          all.mat.taxa[[j]] <- all.mat.taxa[[i]][inxAll[[j]], ]
          sel.dats.taxa[[j]] <- dataSet$data.proc.taxa[[j]][inx[[j]], ]
          
        }
        names(sel.dats.taxa)<- colnames(dataSet$taxa_table)
        reductionSet$selDatsCorr.taxa <- sel.dats.taxa
        reductionSet$micidx<-i
        reductionSet$residx<-3-i
      }
      
    }
  }else{
    sel.dats <- list();
    reductionSet$corr.axis.nms <- list();
    for(i in 1:length(sel.nms)){
      nm = sel.nms[i]
      dataSet <- qs::qread(nm);
      
      inx = which(reductionSet$loading.pos.xyz$ids %in% rownames(dataSet$data.proc));
      loading.df <- reductionSet$loading.pos.xyz[inx, ]
      
      if(retainedNumber > nrow(loading.df)){
        numToKeep <- nrow(loading.df);
      }else{
        numToKeep <- retainedNumber
      }
      
      #scores <- GetDist3D(loading.df, c(0,0,0))
      #scores <- scores * dataSet$misc$pct;
      for(j in 1:retainedComp){
        if(j == 1){
          loading <- loading.df[,1]
          names(loading) <- rownames(loading.df)
          loading <- loading[order(-abs(loading))]
          reductionSet$corr.axis.nms[[j]] <-names(loading)[c(1:numToKeep)]
          toKeep <- names(loading)[c(1:numToKeep)]
        }else{
          loading <- loading.df[,j]
          names(loading) <- rownames(loading.df)
          loading <- loading[order(-abs(loading))]
          reductionSet$corr.axis.nms[[j]] <-names(loading)[c(1:numToKeep)]
          toKeep <- c(toKeep, names(loading)[c(1:numToKeep)])
        }
      }
      
      dat <- dataSet$data.proc
      dat <- dat[rownames(dat) %in% toKeep, ]
      
      sel.dats[[i]] <- dat
    }
  }
  reductionSet$selDatsCorr <- sel.dats
  .set.rdt.set(reductionSet);
  return(1)
}

DoCorrelationFilter <- function(sign="both", crossOmicsOnly="false",networkInfer="NA", threshold.inter=0.9, 
                                threshold.intra=0.5, numToKeep=2000, update="false"){
                 DoCorrelationFilterTaxa(sign, crossOmicsOnly, networkInfer,threshold.inter,threshold.intra, numToKeep,"genus","agora","false");                 
}

DoCorrelationFilterTaxa <- function(sign="both", crossOmicsOnly="false",networkInfer="NA", threshold.inter=0.9, 
                                threshold.intra=0.5, numToKeep=2000,taxlvl="genus",datagem="agora",update="false"){

  reductionSet <- .get.rdt.set();
  
  if(update=="false" | !(exists("selDatsCorr.taxa",reductionSet))){
    
    sel.inx <- mdata.all==1;
    sel.nms <- names(mdata.all)[sel.inx];
    for(i in 1:length(sel.nms)){
      nm = sel.nms[i]
      dataSet <- qs::qread(nm);
      labels <- c(labels, dataSet$enrich_ids);
    } 
    
    corr.mat <- reductionSet$corr.mat
    sel.dats <- reductionSet$selDatsCorr;
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
    corr.mat.intra <- corr.mat[-c(which(rownames(corr.mat) %in% rownames(sel.dats[[2]]))),-c(which(colnames(corr.mat) %in% rownames(sel.dats[[1]])))]; #between-omics
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
    qs::qsave(cor.list, file="cor.list.qs");
    
    reductionSet$corr.mat.inter <- cor_edge_list_inter[cor_edge_list_inter$correlation != 1, ]
    reductionSet$corr.mat.intra <- cor_edge_list_intra[cor_edge_list_intra$correlation != 1, ]
    
    if(sign == "both"){
      cor.inx.inter <- abs(cor_edge_list_inter$correlation) > threshold.inter
      cor.inx.intra <- abs(cor_edge_list_intra$correlation) > threshold.intra
    }else if(sign == "positive"){
      cor.inx.inter <- cor_edge_list_inter$correlation > threshold.inter
      cor.inx.intra <- cor_edge_list_intra$correlation > threshold.intra
    }else{
      cor.inx.inter <- cor_edge_list_inter$correlation < -threshold.inter
      cor.inx.intra <- cor_edge_list_intra$correlation < -threshold.intra
    }
    
    cor_edge_list_inter <- cor_edge_list_inter[cor.inx.inter, ];
    cor_edge_list_intra <- cor_edge_list_intra[cor.inx.intra, ];
    
    cor_edge_list <- cor_edge_list_inter
    
    if(crossOmicsOnly == "true"){
      cor_edge_list <- cor_edge_list_intra
    }else{
      cor_edge_list <- rbind(cor_edge_list, cor_edge_list_intra)
    }
    
    #filter interomics only
    numToKeep <- as.numeric(numToKeep)
    if(numToKeep > length(unique(cor_edge_list$correlation)) ){
      numToKeep <-length(unique(cor_edge_list$correlation))
    }
    
    #if(sort(abs(unique(cor_edge_list$correlation)))[1:numToKeep]) > 1000){
    #top.edge <- sort(abs(unique(cor_edge_list$correlation)))[1:1000];
    #}else{                                                    #default only show top 100 significant edges when #edges>1000   
    top.edge <- sort(abs(unique(cor_edge_list$correlation)))[c(1:numToKeep)]; #default only show top 20% significant edges when #edges<1000
    #}
    top.inx <- match(abs(cor_edge_list$correlation), top.edge);
    
    cor_edge_list <- cor_edge_list[!is.na(top.inx), ,drop=F];
    
    new_g <- graph_from_data_frame(cor_edge_list, F);
    new_g <-simplify(new_g, edge.attr.comb="mean")
    
    if(nrow(cor_edge_list)<3){
      msg.vec <<- paste0("Less than 3 correlations have been identified using an inter-omics correlation threshold of ", threshold.inter, " and intra-omics correlation threshold of ", threshold.intra);
    }
    
    new_g <- graph_from_data_frame(cor_edge_list, F);
    new_g <-simplify(new_g, edge.attr.comb="mean")
    type.list <- list();
    for(i in 1:length(sel.nms)){
      type.list[[sel.nms[[i]]]] <- unique(cor_edge_list[,i]);
    }
    reductionSet$taxlvl <-"Feature"
    .set.rdt.set(reductionSet);
    
    intres <- ProcessGraphFile(new_g, labels, type.list);
  }else  {
    
    taxlvl <- gsub("(^[[:alpha:]])", "\\U\\1", taxlvl, perl=TRUE)
    
    if(exists("corr.mat.inter.taxa",reductionSet) & taxlvl==reductionSet$taxlvl){
      cor_edge_list_inter <- reductionSet$corr.mat.inter.taxa[[taxlvl]] 
      cor_edge_list_intra <- reductionSet$corr.mat.intra.taxa[[taxlvl]] 
      
    }else{
      
      
      sel.inx <- mdata.all==1;
      sel.nms <- names(mdata.all)[sel.inx];
      micidx <- reductionSet$micidx
      residx <- reductionSet$residx
      
      dataSet <- qs::qread(sel.nms[micidx]);
      
      labels <- unique(dataSet$taxa_table[,taxlvl])
      labels <- setNames(labels, labels)
      dataSet <- qs::qread(sel.nms[residx]);
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
    
    
    if(sign == "both"){
      cor.inx.inter <- abs(cor_edge_list_inter$correlation) > threshold.inter
      cor.inx.intra <- abs(cor_edge_list_intra$correlation) > threshold.intra
    }else if(sign == "positive"){
      cor.inx.inter <- cor_edge_list_inter$correlation > threshold.inter
      cor.inx.intra <- cor_edge_list_intra$correlation > threshold.intra
    }else{
      cor.inx.inter <- cor_edge_list_inter$correlation < -threshold.inter
      cor.inx.intra <- cor_edge_list_intra$correlation < -threshold.intra
    }
    
    
    cor_edge_list_inter <- cor_edge_list_inter[cor.inx.inter, ];
    cor_edge_list_intra <- cor_edge_list_intra[cor.inx.intra, ];
    
    cor_edge_list <- cor_edge_list_inter
    
    
    if(crossOmicsOnly == "true"){
      cor_edge_list <- cor_edge_list_intra
    }else{
      cor_edge_list <- rbind(cor_edge_list, cor_edge_list_intra)
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
      msg.vec <<- paste0("Less than 3 correlations have been identified using an inter-omics correlation threshold of ", threshold.inter, " and intra-omics correlation threshold of ", threshold.intra);
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
    
  }
  
  
}

ExportOmicsPairs <- function(fileName, type){
  
  cor.list <- qs::qread("cor.list.qs");
  cor.obj <- cor.list[[type]];
  colnames(cor.obj) = c("Id1", "Id2", "Correlation");
  write.table(cor.obj,row.name=F, file=fileName);
  
}

DoOmicsCorrelation <- function(cor.method="univariate",cor.stat="pearson"){
  labels <- vector();
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  m2midx<-0
  for(i in 1:length(sel.nms)){
    nm = sel.nms[i]
    dataSet <- qs::qread(nm);
    labels <- c(labels, dataSet$enrich_ids);
    
    if(exists("m2m",dataSet)){
      
      labels.taxa <- lapply(dataSet$data.taxa, function(x) rownames(x))
      labels.taxa <-  lapply(labels.taxa, function(x) setNames(x,x))
    }else if(exists("labels.taxa")){
      labels.taxa <- lapply(labels.taxa, function(x) c(x,dataSet$enrich_ids))
    }
    
    
  }
  
  reductionSet <- .get.rdt.set();
  
  sel.dats <- reductionSet$selDatsCorr;
  
  require("igraph");
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  
  if(length(sel.nms) != 2){
    return(0)
  }
  residx <- reductionSet$residx
  
  if(cor.method == "univariate"){
    corr.mat <- cor(cbind(t(sel.dats[[1]]), t(sel.dats[[2]])), method=cor.stat);
    if(exists("selDatsCorr.taxa",reductionSet)){
      
      corr.mat.taxa <- lapply(reductionSet$selDatsCorr.taxa, function(x){
        cor(cbind(t(x), t(sel.dats[[residx]])), method=cor.stat)
      })
      
      reductionSet$corr.mat.taxa <- corr.mat.taxa
    }
    
  }else if(cor.method == "MI"){
    library(parmigene)
    res = knnmi.all(rbind(sel.dats[[1]], sel.dats[[2]]), k=5)
    scale = 1/max(res)
    corr.mat = res * scale
    
    if(exists("selDatsCorr.taxa",reductionSet)){
      corr.mat.taxa <- list()
      for(j in 1:length(reductionSet$selDatsCorr.taxa)){
        res = knnmi.all(rbind(reductionSet$selDatsCorr.taxa[[j]], sel.dats[[residx]]), k=5)
        scale = 1/max(res)
        corr.mat.taxa[[j]] = res * scale
        
      }
      reductionSet$corr.mat.taxa <- corr.mat.taxa
    }
    
    
  }else{
    library(ppcor);
    sel.res <- cbind(t(sel.dats[[1]]), t(sel.dats[[2]]))
    res <- pcor(sel.res, method=cor.stat);
    corr.mat <- res$estimate;
    rownames(corr.mat) <- colnames(sel.res)
    colnames(corr.mat) <- colnames(sel.res)
    
    if(exists("selDatsCorr.taxa",reductionSet)){
      
      sel.res <- lapply(reductionSet$selDatsCorr.taxa, function(x){
        cbind(t(x), t(sel.dats[[residx]]))
      })
      res <- lapply(sel.res, function(x){pcor(x, method=cor.stat)})
      corr.mat.taxa <- lapply(res, function(x) x$estimate);
      for(j in 1:length(corr.mat.taxa)){
        
        rownames(corr.mat.taxa[[j]]) <- colnames(sel.res[[j]])
        colnames(corr.mat.taxa[[j]]) <- colnames(sel.res[[j]])
      }
      reductionSet$corr.mat.taxa <- corr.mat.taxa
    }
  }
  reductionSet$corr.mat <- corr.mat
  .set.rdt.set(reductionSet);
  return(1);
}


PlotCorrHistogram <- function(imgNm, dpi=72, format="png"){
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  library(ggplot2)
  reductionSet <- .get.rdt.set();
  
  fig.list <- list();
  for( i in 1:2){
    if(i == 1){
      cors <- reductionSet$corr.mat.inter
      titleText <- "Between-omics correlation"
    }else{
      cors <- reductionSet$corr.mat.intra
      titleText <- "Intra-omics correlation"
    }
    cor.data <- cors[upper.tri(cors, diag = FALSE)]  # we're only interested in one of the off-diagonals, otherwise there'd be duplicates
    cor.data <- as.data.frame(cor.data)  # that's how ggplot likes it
    summary(cor.data)
    colnames(cor.data) <- "coefficient"
    coefficient.p <- function(r, n) {
      pofr <- ((1-r^2)^((n-4)/2))/beta(a = 1/2, b = (n-2)/2)
      return(pofr)
    }
    
    fig.list[[i]] <- ggplot(cors, aes(x=correlation)) + geom_histogram() +
      xlim(-1,1) +
      theme_bw()
  }
  library(Cairo)
  library(ggpubr)
  Cairo(file=imgNm, width=10, height=8, unit="in", type="png", bg="white", dpi=dpi);
  p1 <- ggarrange(plotlist=fig.list, ncol = 1, labels=c("Between-omics correlation", "Intra-omics correlation"))
  print(p1)
  dev.off();
  
}


PlotDegreeHistogram <- function(imgNm, netNm = "NA", dpi=72, format="png"){
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  Cairo(file=imgNm, width=400, height=400, type="png", bg="white");
  library(ggplot2)
  G.degrees <- degree(overall.graph)
  
  G.degree.histogram <- as.data.frame(table(G.degrees))
  G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])
  
  p <- ggplot(G.degree.histogram, aes(x = G.degrees, y = Freq)) +
    geom_point() +
    scale_x_continuous("Degree\n(nodes containing that amount of connections)",
                       breaks = c(1, 3, 10, 30, 100, 300),
                       trans = "log10") +
    scale_y_continuous("Frequency\n(number of nodes)",
                       breaks = c(1, 3, 10, 30, 100, 300, 1000),
                       trans = "log10") +
    ggtitle("Degree Distribution (log-log)") +
    theme_bw()  +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off();
}

PlotBetweennessHistogram <- function(imgNm, netNm = "NA", dpi=72, format="png"){
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  Cairo(file=imgNm, width=400, height=400, type="png", bg="white");
  library(ggplot2)
  if(netNm != "NA"){
    overall.graph <- ppi.comps[[netNm]];
  }
  G.degrees <- betweenness(overall.graph)
  
  G.degree.histogram <- as.data.frame(table(G.degrees))
  G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])
  
  p <- ggplot(G.degree.histogram, aes(x = G.degrees, y = Freq)) +
    geom_point() +
    scale_x_continuous("Betweenness\n(nodes with that amount of betweenness)",
                       breaks = c(1, 3, 10, 30, 100, 300,1000,3000,10000,30000),
                       trans = "log10") +
    scale_y_continuous("Frequency\n(number of nodes)",
                       breaks = c(1, 3, 10, 30, 100, 300, 1000),
                       trans = "log10") +
    ggtitle("Betweenness Distribution (log-log)") +
    theme_bw()  +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off();
}

GetNetworkTopology <- function(netnm){
  g <- ppi.comps[[netnm]];
  globalProperties <-list();
  globalProperties[["Diameter"]] <-diameter(g);
  globalProperties[["Radius"]] <-radius(g);
  globalProperties[["Average path length"]] <-signif(mean_distance(g), 3);
  globalProperties[["Clustering coefficient"]] <- transitivity(g, type="global");
  propertiesVector <- c(globalProperties[[1]], globalProperties[[2]], globalProperties[[3]], globalProperties[[4]]);
  return(propertiesVector)
}


