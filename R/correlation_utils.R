##################################################
## R scripts for OmicsAnalyst
## Description: Related to correlation analysis
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#default feature selection based on sig genes

DoFeatSelectionForCorr <- function(type="default", retainedNumber=20, retainedComp=3) {

if(!exists("mem.featSelectionForCorr")){
    require("memoise");
    mem.featSelectionForCorr <<- memoise(.do.feat.selectionForCorr);
  }
  return(mem.featSelectionForCorr(type,retainedNumber,retainedComp));
 
}


.do.feat.selectionForCorr <- function(type="default", retainedNumber=20, retainedComp=3){
  print(c(type,"DoFeatSelectionForCorr"))
  sel.dats <- list();
  labels <- vector();
  reductionSet <- .get.rdt.set()
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];

  if(type %in% c("default","custom")){
    for(i in 1:length(sel.nms)){
      dataName = sel.nms[i];
      dataSet <- readDataset(dataName);
      
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
      
      #if more than 1000 sig features, then limit to top 1000 only
      if(nrow(sig.mat) > 1000){
        sig.mat <- sig.mat[c(1:1000),];
      }
      
      inxAll = which(rownames(all.mat) %in% rownames(sig.mat));
      inx = which(rownames(dataSet$data.proc) %in% rownames(sig.mat));
      
      all.mat <- all.mat[inxAll, ];
      dat <- dataSet$data.proc[inx, ];
      rownames(dat) <- paste0(rownames(dat), "_", dataSet$type);
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
      dataName = sel.nms[i]
      dataSet <- readDataset(dataName);

      inx = which(reductionSet[[type]]$loading.pos.xyz$ids %in% rownames(dataSet$data.proc));
      loading.df <- reductionSet[[type]]$loading.pos.xyz[inx, ]
      
      if(retainedNumber > nrow(loading.df)){
        numToKeep <- nrow(loading.df);
      }else{
        numToKeep <- retainedNumber
      }
      
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
      
      library(stringr)
      
      # Check if all elements start with "X"
      all_start_with_x <- all(str_detect(toKeep, "^X"))
      if(all_start_with_x){
        toKeep <- substring(toKeep, 2, nchar(toKeep))
      }
      
      dat <- dataSet$data.proc
      dat <- dat[rownames(dat) %in% toKeep, ]
      rownames(dat) <- paste0(rownames(dat), "_", dataSet$type);
      sel.dats[[i]] <- dat
    }
   
  }
 
  names(sel.dats) = sel.nms
  reductionSet$selDatsCorr <- sel.dats;
  reductionSet$feat.sel.type <- type;
  .set.rdt.set(reductionSet);
  return(1)
}


DoCorrelationFilter <- function(corSign="both", crossOmicsOnly="false", networkInfer="NA", threshold.inter=0.5, 
                                threshold.intra=0.9, numToKeep=2000, updateRes="false", taxlvl="genus", datagem="agora"){

    if(!exists("my.correlation.filter")){ # public web on same user dir
        compiler::loadcmp("../../rscripts/OmicsAnalystR/R/util_correlation_filter.Rc");    
    }
    return(my.correlation.filter(corSign, crossOmicsOnly, networkInfer, threshold.inter, threshold.intra, numToKeep, updateRes, taxlvl, datagem));
}

GenerateNetworkJson <- function(fileName="omicsanalyst_net_0.json"){
    intres <- convertIgraph2JSON(current.net.nm , fileName);
    return(intres);
}

ExportOmicsPairs <- function(fileName, type){
  
  cor.list <- qs::qread("cor.list.qs");
  cor.obj <- cor.list[[type]];
  colnames(cor.obj) = c("Id1", "Id2", "Correlation");
  write.table(cor.obj,row.name=F, file=fileName);
  
}


DoOmicsCorrelation <- function(cor.method="univariate",cor.stat="pearson",ifAll,metaSel,group){
  labels <- vector();
  print(c(ifAll,metaSel,group))
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
   
  m2midx<-0
  for(i in 1:length(sel.nms)){
    dataName = sel.nms[i]
    dataSet <- readDataset(dataName);
    labels <- c(labels, dataSet$enrich_ids);
    if(exists("m2m",dataSet)){
      labels.taxa <- lapply(dataSet$data.taxa, function(x) rownames(x))
      labels.taxa <-  lapply(labels.taxa, function(x) setNames(x,x))
    }else if(exists("labels.taxa")){
      labels.taxa <- lapply(labels.taxa, function(x) c(x,dataSet$enrich_ids))
    }
  }
  
  reductionSet <- .get.rdt.set();
  reductionSet$cor.stat <- cor.stat;
  reductionSet$cor.method <- cor.method;
  sel.dats <- reductionSet$selDatsCorr;
  load_igraph();


  if(ifAll != "true"){
    meta.info = reductionSet$dataSet$meta.info
    sampleInclude = row.names(meta.info[which(meta.info[[metaSel]]==group),])
  

    if(length(sampleInclude)<4){
         return(-1)
    }else{
     sel.dats <- lapply(sel.dats,function(x){
         return(x[,sampleInclude])
     })
   
   }
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

  library(Hmisc)
 
    res = rcorr(as.matrix(cbind(t(sel.dats[[1]]), t(sel.dats[[2]]))),type = cor.stat);
   corr.mat <- res$r
   corr.p.mat <- res$P
      
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
        corr.p.mat <- res$p.value;
    rownames(corr.mat) <-    colnames(corr.mat) <- colnames(sel.res)
    rownames(corr.p.mat) <-    colnames(corr.p.mat) <-colnames(sel.res)
    
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
  reductionSet$labels <- labels
  reductionSet$corr.mat.path <- "corr.mat.qs"
  qs::qsave(corr.mat, "corr.mat.qs");
  qs::qsave(corr.p.mat,"corr.p.mat.qs")
  .set.rdt.set(reductionSet);
  
  return(1);
}

PlotCorrViolin <- function(imgNm, dpi=72, format="png"){
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  load_ggplot();
  library(scales)
  
  reductionSet <- .get.rdt.set();
  graphs <- qs::qread(reductionSet$corr.graph.path);

  
  fig.list <- list();
  for( i in 1:2){
    if(i == 1){
      g <- graphs$corr.graph.inter 
      titleText <- "Between-omics correlation"
      threshold <- reductionSet$threshold.inter
    }else{
      g <- graphs$corr.graph.intra
      titleText <- "Intra-omics correlation"
      threshold <- reductionSet$threshold.intra
      
    }
    
    df_res <- data.frame(get.edgelist(g),  as.numeric(E(g)$correlation))
    df_res <- df_res[!duplicated(df_res), ]
    colnames(df_res) <- c("source", "target","correlation");
    df_res$type <- titleText;
    
    sig.pos <- oob_censor(df_res$correlation, c(threshold, 1))
    sig.neg <- oob_censor(df_res$correlation, c(-1, -threshold))
    sig.pos.num <- length(na.omit(sig.pos))
    sig.neg.num <- length(na.omit(sig.neg))
    
    fig.list[[i]] <- ggplot2::ggplot(df_res, aes(x = type, y = correlation, fill = type)) +
      geom_violin(trim = FALSE, fill = "#d3d3d3", show.legend = FALSE) + 
      labs(x = titleText) +
      geom_jitter(height = 0, width = 0.05, alpha=0,show.legend = FALSE) +
      theme(legend.position = "none") +  #xlab(dataSet$analysisVar) +
      #stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE) +
      scale_x_discrete(labels = NULL) +
      scale_y_continuous(limits = c(-1, 1)) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "red", size = 0.5) +
      annotate("text", x =0.5, y = threshold, label = sig.pos.num, vjust = -1) +
      geom_hline(yintercept = -threshold, linetype = "dashed", color = "red", size = 0.5) +
      annotate("text", x =0.5, y = -threshold, label = sig.neg.num, vjust = 1.5) +
      theme_bw()+
      theme(text=element_text(size=13), plot.title = element_text(size = 11, hjust=0.5), panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  }
  load_cairo();
  library(ggpubr)
  Cairo(file=imgNm, width=10, height=8, unit="in", type="png", bg="white", dpi=dpi);
  p1 <- ggarrange(plotlist=fig.list, ncol = 2)
  print(p1)
  dev.off();
  
  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$correlation_distribution <- imgNm;
  saveSet(infoSet);

}


PlotDegreeHistogram <- function(imgNm, netNm = "NA", dpi=72, format="png"){
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  Cairo(file=imgNm, width=400, height=400, type="png", bg="white");
  load_ggplot();
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

  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$degree_distribution <- imgNm;
  saveSet(infoSet);
}

PlotBetweennessHistogram <- function(imgNm, netNm = "NA", dpi=72, format="png"){
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  Cairo(file=imgNm, width=400, height=400, type="png", bg="white");
  load_ggplot();
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

  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$betweenness_distribution <- imgNm;
  saveSet(infoSet);
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

#expand non-square adjancy matrix to square
expand.matrix <- function(A){
  m <- nrow(A)
  n <- ncol(A)
  B <- matrix(0,nrow = m, ncol = m)
  C <- matrix(0,nrow = n, ncol = n)
  cbind(rbind(B,t(A)),rbind(A,C))
}



###########generate chordgram
GenerateChordGram <- function(thresh=0.5,maxN,pval,imgName = "chordgram", format = "png", dpi = 300){
  print(c(maxN,pval))
  plotjs <- paste0(imgName, ".json");
  reductionSet <- .get.rdt.set(); 
  corr.mat <- qs::qread("corr.mat.qs");
  corr.p.mat<- qs::qread("corr.p.mat.qs");
  corr.mat <- reshape2::melt(corr.mat)
  corr.p.mat <- reshape2::melt(corr.p.mat) 
   corr.mat$pval <- corr.p.mat$value
   sel.inx <- mdata.all==1; 
 sel.nms <- names(mdata.all)[sel.inx];
   
 sel.dats <- reductionSet$selDatsCorr[sel.nms]; 

  #corr.mat <-   corr.mat[corr.mat$Var1!=corr.mat$Var2,]
   corr.mat <- corr.mat[abs(corr.mat$value)>thresh & corr.p.mat$value < pval,]
  corr.mat <- corr.mat[corr.mat$Var1 %in% rownames(sel.dats[[1]]) & corr.mat$Var2 %in% rownames(sel.dats[[2]]),]
  corr.mat <- corr.mat[order(abs(corr.mat$value),decreasing = T),]
 
  if(nrow(corr.mat)>maxN){
    corr.mat <- corr.mat[1:maxN,]
   }
 
    dataSetList <- lapply(sel.nms, readDataset);

corr.mat$Var1 <- gsub(paste0("_", dataSetList[[1]]$type), "", corr.mat$Var1);
corr.mat$Var1 <- names(dataSetList[[1]]$enrich_ids)[match(corr.mat$Var1,dataSetList[[1]]$enrich_ids)]

corr.mat$Var2 <- gsub(paste0("_", dataSetList[[2]]$type), "", corr.mat$Var2);
corr.mat$Var2 <- names(dataSetList[[2]]$enrich_ids)[match(corr.mat$Var2,dataSetList[[2]]$enrich_ids)]


  reductionSet$chordGram <- corr.mat
  write.csv(corr.mat,"chord_diagram.csv",row.names=F)
  library(jsonlite)
  write_json(corr.mat,plotjs, pretty = TRUE)
  .set.rdt.set(reductionSet);
  return(nrow(corr.mat))
}


###############
DoOmicsDiffCorrelation <- function(cor.method="univariate",cor.stat="pearson",comp.meta="Diagnosis",selnm1,selnm2){
  labels <- vector();
  sel.nms <- c(selnm1,selnm2);
   
  m2midx<-0
  for(i in 1:length(sel.nms)){
    dataName = sel.nms[i]
    dataSet <- readDataset(dataName);
    labels <- c(labels, dataSet$enrich_ids);
    if(exists("m2m",dataSet)){
      labels.taxa <- lapply(dataSet$data.taxa, function(x) rownames(x))
      labels.taxa <-  lapply(labels.taxa, function(x) setNames(x,x))
    }else if(exists("labels.taxa")){
      labels.taxa <- lapply(labels.taxa, function(x) c(x,dataSet$enrich_ids))
    }
  }
  
  reductionSet <- .get.rdt.set();
  reductionSet$cor.stat <- cor.stat;
  reductionSet$cor.method <- cor.method;
  sel.dats <- reductionSet$selDatsCorr;
  sel.dats <- sel.dats[names(sel.dats) %in% sel.nms]
  sel.dats <- lapply(sel.dats, function(x){
    if(nrow(x)>500){
      return(x[1:500,])
    }else{
      return(x)
    }
  })
  load_igraph();
 
  residx <- reductionSet$residx
  
  meta.info =  reductionSet$dataSet$meta.info
  corr.ls <- split(rownames(meta.info),meta.info[[comp.meta]])  

  corr.ls <- corr.ls[unlist(lapply(corr.ls, length))>4]
  
  
 
  if(cor.method == "univariate"){
    library(Hmisc)
    corr.mat.ls <- lapply(corr.ls, function(samp){
       res = rcorr(as.matrix(cbind(t(sel.dats[[1]][,samp]), t(sel.dats[[2]][,samp]))),type = cor.stat);
     })

      if(exists("selDatsCorr.taxa",reductionSet)){
      corr.mat.taxa <- lapply(reductionSet$selDatsCorr.taxa, function(x){
        cor(cbind(t(x), t(sel.dats[[residx]])), method=cor.stat)
      })
      
      reductionSet$corr.mat.taxa <- corr.mat.taxa
    }
    
  }

  reductionSet$diffnet.mat.path <- "diffnet.mat.qs"
  qs::qsave(corr.mat.ls, "diffnet.mat.qs");
  .set.rdt.set(reductionSet);
  return(1);
}
 
 
GenerateDiffNet <- function(corr_thresh=0.7,p_thresh=0.05,imgName = "diffnet", format = "png", dpi = 300,dt1,dt2,layout="kk"){
    print(imgName)
  plotjs <- paste0(imgName, ".json");
  reductionSet <- .get.rdt.set();
  sel.dats <- reductionSet$selDatsCorr[c(dt1,dt2)];
  dataSetList <- lapply(c(dt1,dt2), readDataset);
  type1 <- dataSetList[[1]][["type"]]
  type2 <- dataSetList[[2]][["type"]]
  corr.mat.ls <- qs::qread("diffnet.mat.qs");

  corr.mat.ls <- lapply(corr.mat.ls, function(x){
    corr.mat<-reshape2::melt(x$r)
    p.mat<-reshape2::melt(x$P)
    corr.mat$pval <- p.mat$value
    corr.mat <- corr.mat[corr.mat$Var1 %in% rownames(sel.dats[[1]]) & corr.mat$Var2 %in% rownames(sel.dats[[2]]),]
    #corr.mat <- corr.mat[abs(corr.mat$value)>thresh,]
    return( corr.mat)
  })
  
  sig.idx <- unique(unlist(lapply(corr.mat.ls, function(x) return(which(abs(x$value)>corr_thresh&x$pval< p_thresh)) )))
  corr.mat.ls <- lapply(corr.mat.ls, function(x){
    x <- x[sig.idx,]
   # x$value[x$value<thresh] <-0
    names(x) <-c("from","to","corr","pval")
    x$weight <- abs(x$corr)
    x$from = as.character(x$from)
    x$to = as.character(x$to)
    return(x)
  })

corr.mat.ls <- lapply(corr.mat.ls,function(x){
   df=x
   df$type1 = type1
   df$type2 = type2
   df$source = gsub(paste0("_",type1,"$"),"",df$from)
   df$target = gsub(paste0("_",type2,"$"),"",df$to)
   df$source = names(dataSetList[[1]]$enrich_ids)[match(df$source,dataSetList[[1]]$enrich_ids)]
   df$target = names(dataSetList[[2]]$enrich_ids)[match(df$target,dataSetList[[2]]$enrich_ids)]
   df$corr[df$weight<corr_thresh & df$pval > p_thresh] <- 0
   df$weight[df$weight<corr_thresh & df$pval > p_thresh] <- 0
   return(df)
})
 
  reductionSet$diffList <- corr.mat.ls 
  .set.rdt.set(reductionSet);
  library(jsonlite)
  write_json(corr.mat.ls,plotjs, pretty = TRUE)
  return(1);


  library(igraph)
  library(ggplot2)
  library(ggraph)
  nodes <- data.frame(
    id=1:length( unique(c(corr.mat.ls[[1]]$from, corr.mat.ls[[1]]$to))),
   label = unique(c(corr.mat.ls[[1]]$from, corr.mat.ls[[1]]$to))
   )
  nodes$type <- ifelse(nodes$label %in% rownames(sel.dats[[1]]),type1,type2)
  
  
  edges.ls <-lapply(corr.mat.ls, function(x){
    x$from <- nodes$id[match(x$from,nodes$label)]
    x$to <- nodes$id[match(x$to,nodes$label)]
    return(x)
  })
  edges = edges.ls[[1]] 
  nodes$label <- gsub(type1,"",  nodes$label)
  nodes$label <- gsub(type2,"",  nodes$label)
  nodes$label <- gsub("_$","",  nodes$label)
  graphall <- graph_from_data_frame(edges[edges$weight>0,], vertices = nodes, directed = FALSE)

  min_edge_length = 1
   layout_fixed <- layout_with_kk(graphall)
  #layout_fixed <- layout.norm(layout_fixed, xmin = -min_edge_length, xmax = min_edge_length,  ymin = -min_edge_length, ymax = min_edge_length)
  unique_types <- unique(nodes$type)
  
  node_shapes <- setNames(c(21,22),unique_types)
  
   node_colors <-setNames(c("#00b300","#ffaa00"),unique_types)   # Node color by type
 
  layout_df <- data.frame(
    x = layout_fixed[, 1],
    y = layout_fixed[, 2],
    name = V(graphall)$name  # Node names for alignment
  )
  edges.ls <- lapply(edges.ls, function(x){
    x <- x[x$weight>corr_thresh & x$pval<p_thresh,]
    return(x)
  })
   
  edges.ls <- edges.ls[unlist(lapply(edges.ls,nrow))>0]
  subgraph_ls <- lapply(edges.ls, function(x){
    subg<-induced_subgraph(graphall, vids = V(graphall)[V(graphall) %in% x$from |V(graphall) %in% x$to])
    return(subg)
    })

   
  graph_res <-lapply(names(subgraph_ls), function(title){
    
    g<-plot_graph_with_fixed_layout(subgraph_ls[[title]],layout_df,title,node_colors,node_shapes)
    return(g)
  })

  combined_plot <-  patchwork::wrap_plots(graph_res, ncol = 2, guides = "collect") +
    patchwork::plot_layout(guides = "collect", heights = unit(c(1, 1), "null"), widths = unit(1, "null")) &
    theme(plot.margin = margin(15, 15, 15, 15))
 
  
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  
  if(length(graph_res)>2){
   h=8
  }else{

   h=6
  }
  Cairo::Cairo(file = imgName, type = format, dpi = dpi, width = 9, height = h, units = "in", bg = "white")
  print(combined_plot)
  dev.off()
  
  reductionSet$analSet$diffgraph.ls <- graph_res
  .set.rdt.set(reductionSet);
  return(1)
}

plot_graph_with_fixed_layout <- function(graph, layout_df, title,node_colors,node_shapes) {
  # Merge layout_df to preserve positions for all nodes
  layout_sub <- merge(data.frame(name = V(graph)$name), layout_df, by = "name", sort = FALSE)
  edge_colors <- ifelse(E(graph)$corr > 0, "red", "blue")  
  # Set 'name' explicitly as vertex labels
  E(graph)$color <- edge_colors  # Scale edge width
  V(graph)$color <- node_colors[V(graph)$type]  # Node fill color
  V(graph)$border.color <- V(graph)$color 
  V(graph)$degree <- degree(graph) 
  E(graph)$width <- E(graph)$weight 
 
  # Plot the graph
  ggraph(graph, layout = "manual", x = layout_sub$x, y = layout_sub$y) +
    # Edges: color and width based on correlation and weight
    geom_edge_link0(aes(edge_color = color, edge_width = width), alpha = 0.8) +
    
    # Nodes: Shape and fill dynamically mapped to 'type'
    geom_node_point(aes(shape = type, fill = type,size=4 ),  color = "black", stroke = 0.5) +
    
    # Node labels
     geom_node_text(aes(label = label), size = 2.5, color = "black", repel = TRUE) +
    
    # Scales for edge and node properties
    scale_edge_color_identity(guide = "none") +  # Use edge colors directly
    scale_edge_width(range = c(0.5, 1.2)) +  # Scale edge width
    scale_shape_manual(values = node_shapes) +  # Shapes by type
    scale_fill_manual(values = node_colors) +   # Fill colors by type
    scale_size_continuous(range = c(2, 6)) +
    ggtitle(title)+
    
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 13,margin = margin(t = 10, b = 10)),
          legend.position = "none" )
}


GetChordSymbols1 <- function() {
  rdtSet <- .get.rdt.set()
  symbols <- rdtSet$chordGram[,"Var1"]
   return(symbols)
}

GetChordSymbols2 <- function() {
  rdtSet <- .get.rdt.set()
  symbols <- rdtSet$chordGram[,"Var2"]
   return(symbols)
}

GetChordColNames <- function() {
  rdtSet <- .get.rdt.set()
   if (is.null(rdtSet$chordGram)) {
    stop("correlation reatult table not found.")
  }
 
 chord_colnames <- setdiff(colnames(rdtSet$chordGram),c("Var1","Var2")) # Exclude the symbol column
  
  return(chord_colnames)
}

 
GetChordFileName <- function() {
  rdtSet <- .get.rdt.set()
  
  if (is.null(rdtSet$chordGram)) {
    stop("correlation reatult table not found.")
  }

  return("chord_diagram.csv")
}

GetChordMat <- function() {
  rdtSet <- .get.rdt.set()
 
   if (is.null(rdtSet$chordGram)) {
    stop("correlation reatult table not found.")
  }
 
  chord_matrix <- as.matrix(subset(rdtSet$chordGram, select = -c(Var1,Var2))) # Removing the symbol column
  print(head( chord_matrix ))
  return(chord_matrix)
}

GetchordIds <- function() {
  rdtSet <- .get.rdt.set()
  
  if (is.null(rdtSet$chordGram)) {
    stop("correlation reatult table not found.")
  }
  ids <- rownames(rdtSet$chordGram)
  
  return(ids)
}

GetCorrNetSymbols1 <- function() {
  rdtSet <- .get.rdt.set()
  symbols <- rdtSet$corNet[,"source"]
  return(symbols)
}

GetCorrNetSymbols2 <- function() {
  rdtSet <- .get.rdt.set()
  symbols <- rdtSet$corNet[,"target"]
  return(symbols)
}

GetCorrNetColNames <- function() {
  rdtSet <- .get.rdt.set()
  if (is.null(rdtSet$corNet)) {
    stop("correlation reatult table not found.")
  }
  
  corNet_colnames <- setdiff(colnames(rdtSet$corNet),c("source","target")) # Exclude the symbol column
  
  return(corNet_colnames)
}


GetCorrNetFileName <- function() {
  rdtSet <- .get.rdt.set()
  
  if (is.null(rdtSet$corNet)) {
    stop("correlation reatult table not found.")
  }
  
  return("corNet.csv")
}

GetCorrNetMat <- function() {
  rdtSet <- .get.rdt.set()
  
  if (is.null(rdtSet$corNet)) {
    stop("correlation reatult table not found.")
  }
  
  corNet_matrix <- as.matrix(subset(rdtSet$corNet, select = -c(source,target))) # Removing the symbol column
  #print(head( corNet_matrix ))
  return(corNet_matrix)
}

GetCorrNetIds <- function() {
  rdtSet <- .get.rdt.set()
  
  if (is.null(rdtSet$corNet)) {
    stop("correlation reatult table not found.")
  }
  ids <- rownames(rdtSet$corNet)
  
  return(ids)
}

 
plotFeatCorr <- function(reductionSet=NA,imgName,feat1="PC.O.18.2.0.16.0.0",feat2="DNMT3A",dpi=72,format="png"){
  
  reductionSet <- .get.rdt.set();
  
  imgName <- paste(imgName, "dpi", dpi, ".", format, sep = "")
  library(ggplot2)
  sel.inx <- mdata.all==1; 
  sel.nms <- names(mdata.all)[sel.inx];
  dataSetList <- lapply(sel.nms, readDataset);
  id1 = paste0(dataSetList[[1]]$enrich_ids[feat1],"_",dataSetList[[1]]$type)
  id2 = paste0(dataSetList[[2]]$enrich_ids[feat2],"_",dataSetList[[2]]$type)
  selDatsCorr <- reductionSet$selDatsCorr[sel.nms]
  
  data<-data.frame(x= t(selDatsCorr[[1]])[,id1] ,y=t(selDatsCorr[[2]])[,id2])
    
    p<-ggplot(data, aes(x = x, y = y)) +
      geom_point(size = 3) + # Dots represent the samples
      geom_smooth(method = "lm", se = FALSE) +             
      labs(title = "",
           x = feat1,
           y = feat2 ) +
      theme_minimal()+
      theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), # Add panel border
            plot.margin = margin(10, 10, 10, 10),
            axis.title = element_text(size = 11, color = "black"),
            axis.text = element_text(size = 10, color = "black")) 
    
  w <- 7.5
  h<-6
  Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, width = w, height = h, type = format, bg = "white")
  print(p)
  dev.off()
  
}