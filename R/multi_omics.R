
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
      
      dataSet <- readRDS(nm);
      
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
      dataSet <- readRDS(nm);
      
      inx = which( rownames(reductionSet$loading.pos.xyz) %in% rownames(dataSet$data.proc));
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
          loading <- loading[order(-abs(loading))]
          reductionSet$corr.axis.nms[[j]] <-names(loading)[c(1:numToKeep)]
          toKeep <- names(loading)[c(1:numToKeep)]
        }else{
          loading <- loading.df[,j]
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
      dataSet <- readRDS(nm);
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
      
      dataSet <- readRDS(sel.nms[micidx]);
      
      labels <- unique(dataSet$taxa_table[,taxlvl])
      labels <- setNames(labels, labels)
      dataSet <- readRDS(sel.nms[residx]);
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
    dataSet <- readRDS(nm);
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

PlotDiagnostic <- function(alg, imgName, dpi=72, format="png"){
  dpi <- as.numeric(dpi);
  imgNm <- paste(imgName, "dpi", dpi, ".", format, sep="");
  require("Cairo");
  if(alg %in% c("snf", "spectrum") ){
    h=8
    fig.list <- list()
    library(ggpubr);
  }else{
    h=8
  }
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  reductionSet<-.get.rdt.set();
  
  Cairo(file=imgNm, width=10, height=h, type="png",unit="in", bg="white", dpi=dpi);
  if(alg == "procrustes"){
    res <- reductionSet$dim.res
    error = residuals(res[[1]])
    require("ggrepel")
    
    error.df = data.frame(Samples=names(error), Procrustes_residual=unname(error))
    
    rankres <- rank(-abs(error), ties.method="random");
    inx.x <- which(rankres < 6);
    inx.y <- error[inx.x];
    nms <- names(error)[inx.x];
    subsetdf <- error.df[which(error.df$Samples %in% nms),]
    
    p = ggplot(error.df, aes(x = Samples, y = Procrustes_residual)) + geom_col() + geom_hline(yintercept = summary(error)[c(2,3,5)])+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      geom_text_repel(
        data = subsetdf,
        aes(label = Samples),
        size = 5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      ) +
      theme_bw()
    print(p)
  }else if(alg == "rcca"){
    require(mixOmics)
    plot(reductionSet$dim.res, scree.type = "barplot")
    
  }else if(alg == "spls"){
    #save.image("diag.RData");
    require(mixOmics)
    tune.spls <- mixOmics:::perf(reductionSet$dim.res, validation = "Mfold", folds = 10, progressBar = FALSE, nrepeat = 1)
    if("Q2.total" %in% names(tune.spls)){
      plot(tune.spls$Q2.total)
    }else{
      r = data.frame(Comp=seq.int(nrow(tune.spls$measures$Q2.total$values)), Q2.total=tune.spls$measures$Q2.total$values$value)
      rownames(r) = paste("comp", seq.int(nrow(r)))
      plot(r)
    }
    abline(h = 0.0975)
    
  }else if(alg == "diablo"){
    require(mixOmics)
    res = reductionSet$dim.res    
    set.seed(123) # for reproducibility, only when the `cpus' argument is not used
    # this code takes a couple of min to run
    perf.res <- mixOmics:::perf(res, validation = 'Mfold', folds = 10, nrepeat = 1, dist="max.dist")
    diablo.comp <<- median(perf.res$choice.ncomp$WeightedVote)
    plot(perf.res) 
  }else if(alg == "mcia"){
    res = reductionSet$dim.res 
    p1<-plot.mcoin(type=3, res, phenovec=reductionSet$cls, sample.lab=FALSE, df.color=length(names(mdata.all)))   
    print(p1)
  }else if(alg == "mbpca"){
    res = reductionSet$dim.res 
    plotEig(3, res@eig)
    #mogsa::plot(res, value="eig", type=2, xlab="Component", ylab="Eigenvalue")    
  }else if(alg == "spectrum"){
    #res <- ComputeSilhouette("snf");
    #sp2 <- fviz_silhouette(res[[2]], titlenm=sel.nms[[2]])+
    #theme(axis.title.y=element_blank())
    #lims <- layer_scales(sp2)$y$get_limits()
    
    #sp1 <- fviz_silhouette(res[[1]], titlenm=sel.nms[[1]])+ ylim(lims)+
    #theme(legend.position="none")
    
    
    #fig.list.sub<-list();
    #fig.list.sub[[1]] <- sp1
    #fig.list.sub[[2]] <- sp2
    
    #p11 <- ggarrange(plotlist=fig.list.sub, ncol = 2, nrow = 1)
    #fig.list[[1]] <- p11
    #if(!is.null(reductionSet$clustRes$eigenvector_analysis)){
    #fig.list[[1]] <- function(){plotEig(length(unique(reductionSet$clustVec)), reductionSet$clustRes$eigenvector_analysis$evals)};
    #}else{
    #fig.list[[1]] <- function(){plotEig(length(unique(reductionSet$clustVec)), reductionSet$clustRes$eigensystem$values[1:10])};
    #}
    if(!is.null(reductionSet$clustRes$eigenvector_analysis)){
      plotEig(length(unique(reductionSet$clustVec)), reductionSet$clustRes$eigenvector_analysis[,2]);
    }else{
      plotEig(length(unique(reductionSet$clustVec)), reductionSet$clustRes$eigensystem$values[1:10]);
    }
  }else if(alg == "perturbation"){
    
    res <- reductionSet$clustRes
    library(ggpubr)
    library(ggplot2)
    xlabel="Number of clusters"
    ylabel="AUC"
    auc1 <- res$dataTypeResult[[1]]$Discrepancy$AUC[-1]
    auc.df1 <- data.frame(K=seq.int(length(auc1))+1, evals=auc1)
    auc2 <- res$dataTypeResult[[2]]$Discrepancy$AUC[-1]
    auc.df2 <- data.frame(K=seq.int(length(auc2))+1, evals=auc2)
    auc.df1$evals2 = auc.df2$evals
    colnames(auc.df1) = c("K", sel.nms[[1]], sel.nms[[2]])
    library("tidyverse")
    df <- auc.df1 %>%
      select(K,  sel.nms[[1]], sel.nms[[2]]) %>%
      gather(key = "variable", value = "value", -K)
    p1 <- ggplot(df, aes(x = K, y = value)) + 
      geom_point(aes(color = variable), size=2) + 
      geom_line(aes(color = variable)) + 
      geom_vline(xintercept = length(unique(res$cluster1)),linetype="dashed")+ xlab(xlabel) + 
      ylab(ylabel) +
      theme_bw()
    print(p1)
  }else if(alg == "snf"){
    #res <- ComputeSilhouette("snf");
    
    #sp2 <- fviz_silhouette(res[[2]], titlenm=sel.nms[[2]])+
    #theme(axis.title.y=element_blank())
    #lims <- layer_scales(sp2)$y$get_limits()
    
    #sp1 <- fviz_silhouette(res[[1]], titlenm=sel.nms[[1]])+ ylim(lims)+
    #theme(legend.position="none")
    
    #fig.list.sub<-list();
    #fig.list.sub[[1]] <- sp1
    #fig.list.sub[[2]] <- sp2
    
    #p11 <- ggarrange(plotlist=fig.list.sub, ncol = 2, nrow = 1)
    #fig.list[[1]] <- p11
    #fig.list[[1]] <-  function(){plotEig(length(unique(reductionSet$clustVec)), reductionSet$clustRes[[5]])}
    #p1 <- ggarrange(plotlist=fig.list, ncol = 1, nrow = 1)
    plotEig(length(unique(reductionSet$clustVec)), reductionSet$clustRes[[5]])
  }
  dev.off();
  
  return(1);
}

DoIntegrativeAnalysis <- function(method, sign="both", threshold=0.6, nComp){
  require("igraph");
  if(! method %in% dim.res.methods){
    intRes <- DoDimensionReductionIntegrative("omics", method, method="globalscore", 3);
  }
  
  threshold <- as.numeric(threshold)
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx]
  data <- list()
  labels <- vector();
  for(i in 1:length(sel.nms)){
    dataSet = readRDS(sel.nms[i])
    dat <- dataSet$data.proc
    df <- data.frame(dat, stringsAsFactors = FALSE)
    df <- t(df)
    data[[sel.nms[i]]] <- df
    labels <- c(labels, dataSet$enrich_ids)
  }
  reductionSet<-.get.rdt.set();
  res <- reductionSet$dim.res
  net.res <- mixOmics:::network(res, cutoff = threshold, save="jpeg")
  
  cor_edge_list <- igraph:::as_data_frame(net.res$gR, 'edges');
  if(sign == "both"){
    cor.inx <- abs(cor_edge_list$weight) > threshold
  }else if(sign == "positive"){
    cor.inx <- cor_edge_list$weight > threshold
  }else{
    cor.inx <- cor_edge_list$weigt < -threshold
  }
  
  only_sig <- cor_edge_list[cor.inx, ];
  new_g <- graph_from_data_frame(only_sig, F);
  
  type.list <- list()
  for(i in 1:length(sel.nms)){
    type.list[[sel.nms[[i]]]] <- unique(cor_edge_list[,i]);
  }
  
  ProcessGraphFile(new_g, labels, type.list);
}

NormalizeDataWrapper <-function (nm, opt, colNorm){
  if(nm == "NA"){
    sel.nms <- names(mdata.all)
    for(i in 1:length(sel.nms)){
      dataSet = readRDS(sel.nms[i])
      data <- NormalizingDataOmics(dataSet$data.filtered, opt, colNorm, "NA")
      dataSet$data.proc <- data;
      if(exists("m2m",dataSet)){
        data.norm.taxa <- lapply(dataSet$dataSet$data.filt.taxa, function(x) {
          NormalizingDataOmics(x, opt, colNorm, "NA")
        }
        )
        dataSet$data.proc.taxa <- data.norm.taxa
      }
      RegisterData(dataSet)
    }
    return(1)
  }else{
    
    dataSet <- readRDS(nm);
    data <- NormalizingDataOmics(dataSet$data.filtered, opt, colNorm, "NA")
    dataSet$data.proc <- data;

    if(exists("m2m", dataSet)){
      data.norm.taxa <- lapply(dataSet$data.filt.taxa, function(x) {
        NormalizingDataOmics(x, opt, colNorm, "NA")
      })
      dataSet$data.proc.taxa <- data.norm.taxa
    }

    RegisterData(dataSet)
    return(1)
  }
}

ScaleDataWrapper <-function (nm, scaleNorm){
  
  if(nm == "NA"){
    sel.nms <- names(mdata.all)
    for(i in 1:length(sel.nms)){
      dataSet = readRDS(sel.nms[i])
      data <- NormalizingDataOmics(dataSet$data.filtered, "NA", "NA", scaleNorm)
      dataSet$data.proc <- data;
      if(exists("m2m",dataSet)){
        data.norm.taxa <- lapply(dataSet$dataSet$data.filt.taxa, function(x) {
          NormalizingDataOmics(x, "NA", "NA", scaleNorm)
        }
        )
        dataSet$data.proc.taxa <- data.norm.taxa
      }
      RegisterData(dataSet)
    }
    return(1)
  }else{
    
    dataSet <- readRDS(nm);
    
    data <- NormalizingDataOmics(dataSet$data.filtered, "NA", "NA", scaleNorm)
    dataSet$data.proc <- data;
    if(exists("m2m",dataSet)){
      
      data.norm.taxa <- lapply(dataSet$data.filt.taxa, function(x) {
        NormalizingDataOmics(x, "NA", "NA", scaleNorm)
      }
      )
      
      
      dataSet$data.proc.taxa <- data.norm.taxa
      
      
    }
    
    RegisterData(dataSet)
    return(1)
  }
}

NormalizingDataOmics <-function (data, norm.opt="NA", colNorm="NA", scaleNorm="NA"){
  msg <- ""
  rnms <- rownames(data)
  cnms <- colnames(data)

  # column(sample)-wise normalization
  if(colNorm=="SumNorm"){
    data<-t(apply(data, 2, SumNorm));
    rownm<-"Normalization to constant sum";
  }else if(colNorm=="MedianNorm"){
    data<-t(apply(data, 2, MedianNorm));
    rownm<-"Normalization to sample median";
  }else{
    # nothing to do
    rownm<-"N/A";
  }
  
  if(norm.opt=="log"){
    min.val <- min(data[data>0], na.rm=T)/10;
    numberOfNeg = sum(data<=0, na.rm = TRUE) + 1; 
    totalNumber = length(data)
    data[data<=0] <- min.val;
    data <- log2(data);
    msg <- paste(msg, "Log2 transformation.", collapse=" ");
  }else if(norm.opt=="vsn"){
    require(limma);
    data <- normalizeVSN(data);
    msg <- paste(msg, "VSN normalization.", collapse=" ");
  }else if(norm.opt=="quantile"){
    require('preprocessCore');
    data <- normalize.quantiles(as.matrix(data), copy=TRUE);
    msg <- paste(msg, "Quantile normalization.", collapse=" ");
  }else if(norm.opt=="combined"){
    require(limma);
    data <- normalizeVSN(data);
    require('preprocessCore');
    data <- normalize.quantiles(as.matrix(data), copy=TRUE);
    msg <- paste(msg, "VSN followed by quantile normalization.", collapse=" ");
  }else if(norm.opt=="logcount"){ # for count data, do it in DE analysis, as it is dependent on design matrix
    require(edgeR);
    nf <- calcNormFactors(data);
    y <- voom(data,plot=F,lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    msg <- paste(msg, "Limma based on log2-counts per million transformation.", collapse=" ");
  } else if(norm.opt=="rle"){
    data <- edgeRnorm(data,method="RLE");
    msg <- c(msg, paste("Performed RLE Normalization"));
  }else if(norm.opt=="TMM"){
    data <- edgeRnorm(data,method="TMM");
    msg <- c(msg, paste("Performed TMM Normalization"));
  }else if(norm.opt=="clr"){
    data <- apply(data, 2, clr_transform);
    msg <- "Performed centered-log-ratio normalization.";
  }else if(norm.opt=='LogNorm'){
    min.val <- min(abs(data[data!=0]))/10;
    data<-apply(data, 2, LogNorm, min.val);
  }else if(norm.opt=='CrNorm'){
    norm.data <- abs(data)^(1/3);
    norm.data[data<0] <- - norm.data[data<0];
    data <- norm.data;
  }
  
  
  # scaling
  if(scaleNorm=='MeanCenter'){
    data<-apply(data, 1, MeanCenter);
    scalenm<-"Mean Centering";
  }else if(scaleNorm=='AutoNorm'){
    data<-apply(data, 1, AutoNorm);
    scalenm<-"Autoscaling";
  }else if(scaleNorm=='ParetoNorm'){
    data<-apply(data, 1, ParetoNorm);
    scalenm<-"Pareto Scaling";
  }else if(scaleNorm=='RangeNorm'){
    data<-apply(data, 1, RangeNorm);
    scalenm<-"Range Scaling";
  }else if(scaleNorm=="colsum"){
    data <- sweep(data, 2, colSums(data), FUN="/")
    data <- data*10000000;
    #msg <- c(msg, paste("Performed total sum normalization."));
  }else if(scaleNorm=="upperquartile"){
    suppressMessages(library(edgeR))
    otuUQ <- edgeRnorm(data,method="upperquartile");
    data <- as.matrix(otuUQ$counts);
    #msg <- c(msg, paste("Performed upper quartile normalization"));
  }else if(scaleNorm=="CSS"){
    suppressMessages(library(metagenomeSeq))
    #biom and mothur data also has to be in class(matrix only not in phyloseq:otu_table)
    data1 <- as(data,"matrix");
    dataMR <- newMRexperiment(data1);
    data <- cumNorm(dataMR,p=cumNormStat(dataMR));
    data <- MRcounts(data,norm = T);
    #msg <- c(msg, paste("Performed cumulative sum scaling normalization"));
  }else{
    scalenm<-"N/A";
  }
  if(scaleNorm %in% c('MeanCenter', 'AutoNorm', 'ParetoNorm', 'RangeNorm')){
    data <- t(data)
  }
  
  data <- as.data.frame(data)
  rownames(data) <- rnms;
  colnames(data) <- cnms;
  msg.vec <<- msg;
  return(data)
}

# Function to annotate metabolite data to internal database
AnnotateMetaboliteData <- function(dataName, idtype){
  #if(dataSet$name != dataName){
  dataSet <- readRDS(dataName);
  dataSet$name <- dataName
  #}
  
  data <- dataSet$data.raw;
  qvec <- rownames(data);
  
  # record all the data
  if(!exists("name.map", where = dataSet)){
    dataSet$name.map <- list();
  }
  
  # map to cpd db from metaboanalyst
  dataSet <- MetaboliteMappingExact(dataSet, qvec, idtype)
  
  # do some sanity check
  todo.inx <- which(is.na(dataSet$name.map$hit.inx));
  resint <- 1;
  print(length(todo.inx)/length(dataSet$name.map$hit.inx));
  if(length(todo.inx)/length(dataSet$name.map$hit.inx) > 0.5){
    msg <- c("Over half of the compound IDs could not be matched to our database. Please make 
             sure that correct compound IDs or common compound names are used.");
    resint <- 2;
  }else if (length(todo.inx) > 15){
    msg <- c("There are >15 compounds without matches. You can either proceed or if necessary, update these compound IDs and upload again.");   
    resint <- 1;     
  }else{
    msg <- paste0("Name matching OK!", " A total of ", nrow(data), " are mapped.");  
    resint <- 1; 
  }
  
  msg <- c(msg, paste0("A total of ", length(na.omit(dataSet$name.map$hit.inx)), " are mapped."))
  #dataSet$enrich_ids <- rownames(dataSet$data.proc)
  #dataSet$enrich_ids[dataSet$name.map$match.state == 1] <- dataSet$name.map$hit.values[dataSet$name.map$match.state == 1]
  #names(dataSet$enrich_ids) <- dataSet$enrich_ids
  data <- as.matrix(dataSet$data.raw);
  rownames(data) <- unname(dataSet$enrich_ids);
  data <- RemoveDuplicates(data, "mean", quiet=T); # remove duplicates
  data <- as.data.frame(data)
  dataSet$data.annotated <- data
  RegisterData(dataSet);
  
  msg.vec <<- msg
  
  return(resint)
}

ReadOmicsData <- function(fileName, omics.type=NA) {
  # need to handle reading .csv files too!
  
  data <- .readDataTable(fileName)
  dataSet <- list();
  
  meta.info <- list();
  cls.inx <- grep("^#CLASS", data[,1]);
  
  if(length(cls.inx) > 0){ 
    for(i in 1:length(cls.inx)){
      inx <- cls.inx[i];
      cls.nm <- substring(data[inx, 1],2); # discard the first char #
      if(nchar(cls.nm) > 6){
        cls.nm <- substring(cls.nm, 7); # remove class
      }
      cls.lbls <- data[inx, -1];
      # test NA
      na.inx <- is.na(cls.lbls);
      cls.lbls[na.inx] <- "NA";
      cls.lbls <- ClearFactorStrings(cls.nm, cls.lbls);
      
      meta.info[[cls.nm]] <- cls.lbls;
    }
    meta.info <- data.frame(meta.info);
    dataSet$meta <- meta.info
    data = data[-cls.inx,];
    rownames(dataSet$meta) <- colnames(data)[-1];
  }
  
  if(class(data) == "try-error" || ncol(data) == 1){
    AddErrMsg("Data format error. Failed to read in the data!");
    AddErrMsg("Make sure the data table is saved as comma separated values (.csv) format!");
    AddErrMsg("Please also check the followings: ");
    AddErrMsg("Either sample or feature names must in UTF-8 encoding; Latin, Greek letters are not allowed.");
    AddErrMsg("We recommend to use a combination of English letters, underscore, and numbers for naming purpose.");
    AddErrMsg("Make sure sample names and feature (peak, compound) names are unique.");
    AddErrMsg("Missing values should be blank or NA without quote.");
    AddErrMsg("Make sure the file delimeters are commas.");
    return(0);
  }
  
  var.nms <- data[,1];
  data[,1] <- NULL;
  smpl.nms <- colnames(data);
  data <- as.matrix(data);
  rownames(data) <- var.nms;
  
  data <- RemoveDuplicates(data, "mean", quiet=T); # remove duplicates
  data <- as.data.frame(data)
  var.nms <- rownames(data)
  
  msg <- paste("A total of ", ncol(data), " samples and ", nrow(data), " features were found")
  
  # Basic checks - no duplicate samples names
  # Check for uniqueness of sample name
  if(length(unique(smpl.nms))!=length(smpl.nms)){
    dup.nm <- paste(smpl.nms[duplicated(smpl.nms)], collapse=" ");
    AddErrMsg("Duplicate sample names are not allowed!");
    AddErrMsg(dup.nm);
    return(0);
  }
  
  # checking variable names - no duplicate variables for metabolites and microbiome?
  #if(length(unique(var.nms))!=length(var.nms)){
  #  dup.nm <- paste(var.nms[duplicated(var.nms)], collapse=" ");
  #  AddErrMsg("Duplicate feature names are not allowed!");
  #  AddErrMsg(dup.nm);
  #  return(0);
  #}
  
  # now check for special characters in the data labels
  if(sum(is.na(iconv(smpl.nms)))>0){
    na.inx <- is.na(iconv(smpl.nms));
    nms <- paste(smpl.nms[na.inx], collapse="; ");
    AddErrMsg(paste("No special letters (i.e. Latin, Greek) are allowed in sample names!", nms, collapse=" "));
    return(0);
  }
  
  if(sum(is.na(iconv(var.nms)))>0){
    na.inx <- is.na(iconv(var.nms));
    nms <- paste(var.nms[na.inx], collapse="; ");
    AddErrMsg(paste("No special letters (i.e. Latin, Greek) are allowed in feature names!", nms, collapse=" "));
    return(0);
  }
  
  # only keep alphabets, numbers, ",", "." "_", "-" "/"
  smpl.nms <- CleanNames(smpl.nms, "sample_name");
  
  # keep a copy of original names for saving tables 
  orig.var.nms <- var.nms;
  var.nms <- CleanNames(var.nms, "var_name"); # allow space, comma and period
  names(orig.var.nms) <- var.nms;
  
  current.msg <<- msg;
  # now create the dataSet
  dataSet$orig.var.nms <- orig.var.nms;
  data <- data.frame(apply(data, 2, function(x) as.numeric(as.character(x))))
  # now reassgin the dimension names
  colnames(data) <- smpl.nms;
  rownames(data) <- var.nms;
  dataSet$data.proc <- data
  dataSet$data.raw <- data
  dataSet$data.annotated <- ""
  dataSet$data.missed <- ""
  dataSet$data.filtered <- ""
  dataSet$name <- fileName;
  dataSet$de.method <- "NA"
  dataSet$type <- omics.type;
  if(omics.type == "rna_b"){
    readableType <- "Transcriptomics";
  }else if (omics.type == "met_t" || omics.type == "met_u"){
    readableType <- "Metabolomics";
  }else if (omics.type == "mic_m"){
    readableType <- "Microbiome";
  }else if (omics.type == "prot"){
    readableType <- "Proteomics";
  }else if (omics.type == "mirna"){
    readableType <- "miRNA";
  }else{
    readableType <-  omics.type;
  }
  dataSet$readableType <- readableType;
  dataSet$enrich_ids = rownames(dataSet$data.proc)
  names(dataSet$enrich_ids) = rownames(dataSet$data.proc)
  # update current dataset
  RegisterData(dataSet);
  partialToBeSaved <<- c(partialToBeSaved, c(fileName))
  return(1)
}

SetOmicsType <- function(fileName, omics.type=NA) {
  dataSet <- readRDS(fileName)
  dataSet$type <- omics.type
  if(omics.type == "rna_b"){
    readableType <- "Transcriptomics";
  }else if (omics.type == "met_t" || omics.type == "met_u"){
    readableType <- "Metabolomics";
    # dataSet$isMet <- TRUE
  }else if (omics.type == "mic_m"){
    readableType <- "Microbiome";
    #   dataSet$isMic <- TRUE
  }else if (omics.type == "prot"){
    readableType <- "Proteomics";
  }else if (omics.type == "mirna"){
    readableType <- "miRNA";
  }else{
    readableType <-  omics.type;
  }
  dataSet$readableType <- readableType;
  
  RegisterData(dataSet);
}

#'Plot PCA plot for multi-omics samples
#'@description 
#'@param imgNm name of the image to output
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'

PlotMultiPCA <- function(imgNm, dpi=72, format="png",factor="1"){
  
  require("Cairo");
  library(ggplot2)
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  sel.nms <- names(mdata.all)
  fig.list <- list()
  
  pca.list<- list()
  pct <- list();
  
  for(i in 1:length(sel.nms)){
    dataSet = readRDS(sel.nms[i])
    x <- dataSet$data.proc
    pca <- prcomp(t(na.omit(x)), center=T, scale=T);
    imp.pca<-summary(pca)$importance;
    xlabel <- paste0("PC1"," (", 100*round(imp.pca[2,][1], 3), "%)")
    ylabel <- paste0("PC2"," (", 100*round(imp.pca[2,][2], 3), "%)")
    names <- colnames(x);
    pca.res <- as.data.frame(pca$x);
    pca.res <- pca.res[,c(1,2)]
    
    xlim <- GetExtendRange(pca.res$PC1);
    ylim <- GetExtendRange(pca.res$PC2);
    
    Factor <- dataSet$meta[,1];
    pca.rest <- pca.res
    pca.rest$Conditions <- Factor
    pca.rest$names <- rownames(pca.res)
    
    pcafig <- ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions)) +
      geom_point(size=3, alpha=0.5) + 
      xlim(xlim) + 
      ylim(ylim) + 
      xlab(xlabel) + 
      ylab(ylabel) +
      theme_bw()
    fig.list[[i]] <- pcafig
    
    pos.xyz = data.frame(x=pca$x[,1], y=pca$x[,2], z=pca$x[,3]);
    loading.pos.xyz = data.frame(pca$rotation[,c(1:3)]);
    loadingNames = rownames(pca$rotation);
    
    pos.xyz <- as.data.frame(pos.xyz);
    pos.xyz <- unitAutoScale(pos.xyz);
    rownames(pos.xyz) = rownames(dataSet$meta);
    dataSet$pos.xyz = pos.xyz;
    
    loadingNames <- rownames(loading.pos.xyz);
    loading.pos.xyz <- as.data.frame(loading.pos.xyz);
    loading <- unitAutoScale(loading.pos.xyz);
    rownames(loading) <- loadingNames;
    nm <- paste0("pca_", dataSet$type);
    pca.list[[nm]][["score"]] <- pos.xyz * 1000;
    pca.list[[nm]][["loading"]] <- loading* 1000;    
    
    pct2Nm <- paste0("pca_", dataSet$type)
    pct[[pct2Nm]] <- unname(round(imp.pca[2,],3))[c(1:3)]*100;
  }
  pca.list$pct2 <- pct;
  qs::qsave(pca.list, file="pca.scatter.qs");
  
  h<-6*round(length(fig.list)/2)
  Cairo(file=imgNm, width=14, height=h, type=format, bg="white", unit="in", dpi=dpi);
  library("ggpubr")
  p1 <- ggarrange(plotlist=fig.list, ncol = 2, nrow = round(length(fig.list)/2), labels=sel.nms)
  print(p1)
  dev.off();
  
}


PlotMultiDensity <- function(imgNm, dpi=72, format="png",factor="1"){
  require("Cairo");
  library(ggplot2)
  dpi <- as.numeric(dpi)
  factor <- as.numeric(factor)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  sel.nms <- names(mdata.all)
  fig.list <- list()
  merged.df <- data.frame
  df.list <- list()
  for(i in 1:length(sel.nms)){
    dataSet <- readRDS(sel.nms[i])
    dat <- dataSet$data.proc
    if(i==1){
      st <- stack(as.data.frame(dat))
      st$Dataset <- rep(sel.nms[i], nrow(dat))
      merged.df <- st
    }else{
      st <- stack(as.data.frame(dat))
      st$Dataset <- rep(sel.nms[i], nrow(dat)) 
      merged.df <-rbind(merged.df, st)
    } 
    
    
  }
  
  type<-merged.df$type
  merged.df$ind <- paste0(merged.df$ind, "_", merged.df$type)
  
  Cairo(file=imgNm, width=10, height=6, type=format, bg="white", dpi=dpi, unit="in");
  g =ggplot(merged.df, aes(x=values)) + 
    geom_line(aes(color=Dataset, group=ind), stat="density", alpha=0.1) + 
    geom_line(aes(color=Dataset), stat="density", alpha=0.7, size=3) +
    theme_bw()
  print(g)
  dev.off();
} 

GetMultiSummary<- function(){
  sel.nms <- names(mdata.all);
  featureNum <- "";
  dat.nms <- "";
  for(i in 1:length(sel.nms)){
    dataSet = readRDS(sel.nms[i])
    dat <- dataSet$data.proc
    if(i == 1){
      cls.lbls <- dataSet$meta[,1]
      featureNum <- nrow(dat)
      sampleNum <- ncol(dat)
      dat.nms <- sel.nms[i]
    }else{
      featureNum <- c(featureNum, nrow(dat))
      dat.nms <- c(dat.nms, sel.nms[i])
    }
  }
  featureNum <- paste(featureNum, collapse="; ")
  dat.nms <- paste(dat.nms, collapse="; ")
  cls.lbls <- unique(cls.lbls)
  cls.lbls <- paste(cls.lbls, collapse="; ")
  return(c(sampleNum, featureNum, dat.nms, cls.lbls) )
}

GetOmicsDataDims <- function(dataName){
  #if(dataSet$name != dataName){
  dataSet <- readRDS(dataName);
  #}
  dm <- dim(dataSet$data.proc);
  naNum <- sum(is.na(dataSet$data.proc));
  return(c(dm, naNum));
} 

doIdMapping <- function(q.vec, type){
  if(type %in% c("mir_acc", "mir_id", "mirnet")){
    require('RSQLite');
    path <- paste0(sqlite.path, "mir2gene.sqlite")
    mir.db <- dbConnect(SQLite(), path);
    if(type == "mir_id"){
      q.vec <- tolower(q.vec)
    }
    query <- paste (shQuote(q.vec),collapse=",");
    table.n <- data.org
    statement <- paste("SELECT * FROM ", data.org, " WHERE ((",type," IN (",query,")) OR (mir_id IN (", query, ")))", sep="");
    
    mir.dic <- .query.sqlite(mir.db, statement);       
    
    entrezs <- mir.dic[c("mir_acc", type)];
    if(nrow(entrezs) == 0){
      return(0)
    }
    entrezs <- entrezs[!duplicated(entrezs[,"mir_acc"]),]
    rownames(entrezs) = seq.int(nrow(entrezs));
    entrezs <- data.frame(lapply(entrezs, as.character), stringsAsFactors=FALSE)
    colnames(entrezs) = c("gene_id", "accession");
    
    return(entrezs);
  }
}

ComputeHeatmap <- function(fileNm, type){
  sel.nms <- names(mdata.all);
  reductionSet <- .get.rdt.set();
  if(type == "NA"){
    reductionSet$clustVec <- "NA";
  }
  .set.rdt.set(reductionSet);
  res.list <- list()
  for(i in 1:length(sel.nms)){
    dataSet <- readRDS(sel.nms[i])
    res <- ComputePathHeatmapTable(dataSet);
    res.list[[i]] <- res;
  }
  require(RJSONIO);
  json.mat <- toJSON(res.list, .na='null');
  sink(fileNm);
  cat(json.mat);
  sink();
  current.msg <<- "Data is now ready for heatmap visualization!";
  jsonNms$heatmap <<- fileNm
  partialToBeSaved <<- c(partialToBeSaved, c(fileNm))
  return(1)
}

ComputePathHeatmapTable <- function(dataSet){
  data <- dataSet$data.proc
  sig.ids <- rownames(dataSet$data.proc);
  enrich.nms1 <- dataSet$enrich_ids
  
  metadf <- dataSet$meta
  meta.nms <- colnames(metadf)[-ncol(metadf)]
  metadf <- as.data.frame(metadf[,-which(colnames(metadf) == "newcolumn")])
  colnames(metadf) = meta.nms
  
  res <- dataSet$comp.res[rownames(data),c(1:2)]
  colnames(res) <- c("statistic", "p.value")
  stat.pvals <- unname(as.vector(res[,2]));
  
  # scale each gene 
  dat <- t(scale(t(data)));
  
  rankPval = order(as.vector(stat.pvals))
  stat.pvals = stat.pvals[rankPval]
  dat = dat[rankPval,]
  reductionSet <- .get.rdt.set();
  if(reductionSet$clustVec != "NA"){
    vec <- reductionSet$clustVec
    names(vec) <- colnames(dat);
    veco <- vec[order(vec)];
    smpl.int.rk <- match(names(veco), colnames(dat));
    dat <- dat[, smpl.int.rk];
    
  }
  
  # now pearson and euclidean will be the same after scaleing
  dat.dist <- dist(dat); 
  
  orig.smpl.nms <- colnames(dat);
  orig.gene.nms <- rownames(dat);
  
  rest <- t(apply(dat, 1, function(x){as.numeric(cut(x, breaks=30))}));
  # do clustering and save cluster info
  # convert order to rank (score that can used to sort)
  if(nrow(dat)> 1){
    dat.dist <- dist(dat);
    gene.ward.ord <- hclust(dat.dist, "ward.D")$order;
    gene.ward.rk <- match(orig.gene.nms, orig.gene.nms[gene.ward.ord]);
    gene.ave.ord <- hclust(dat.dist, "ave")$order;
    gene.ave.rk <- match(orig.gene.nms, orig.gene.nms[gene.ave.ord]);
    gene.single.ord <- hclust(dat.dist, "single")$order;
    gene.single.rk <- match(orig.gene.nms, orig.gene.nms[gene.single.ord]);
    gene.complete.ord <- hclust(dat.dist, "complete")$order;
    gene.complete.rk <- match(orig.gene.nms, orig.gene.nms[gene.complete.ord]);
    
    dat.dist <- dist(t(dat));
    smpl.ward.ord <- hclust(dat.dist, "ward.D")$order;
    smpl.ward.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ward.ord])
    smpl.ave.ord <- hclust(dat.dist, "ave")$order;
    smpl.ave.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ave.ord])
    smpl.single.ord <- hclust(dat.dist, "single")$order;
    smpl.single.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.single.ord])
    smpl.complete.ord <- hclust(dat.dist, "complete")$order;
    smpl.complete.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.complete.ord])
  }else{
    # force not to be single element vector which will be scaler
    #stat.pvals <- matrix(stat.pvals);
    gene.ward.rk <- gene.ave.rk <- gene.single.rk <- gene.complete.rk <- matrix(1);
    smpl.ward.rk <- smpl.ave.rk <- smpl.single.rk <- smpl.complete.rk <- 1:ncol(dat);
  }
  
  gene.cluster <- list(
    ward = gene.ward.rk,
    average = gene.ave.rk,
    single = gene.single.rk,
    complete = gene.complete.rk,
    pval = seq.int(c(1:length(gene.ward.rk)))
  );
  
  sample.cluster <- list(
    ward = smpl.ward.rk,
    average = smpl.ave.rk,
    single = smpl.single.rk,
    complete = smpl.complete.rk
  );
  
  
  # prepare meta info    
  # 1) convert meta.data info numbers
  # 2) match number to string (factor level)
  if(reductionSet$clustVec != "NA"){
    metadf$Cluster <- reductionSet$clustVec;
    sample.cluster[[reductionSet$clustType]] <- smpl.int.rk;
    metadf <- metadf[smpl.int.rk,]
  }
  
  meta <- metadf;
  grps <- colnames(metadf)
  nmeta <- meta.vec <- NULL;
  uniq.num <- 0;
  for (i in 1:ncol(meta)){
    cls <- meta[,i];
    grp.nm <- grps[i];
    meta.vec <- c(meta.vec, as.character(cls))
    # make sure each label are unqiue across multiple meta data
    ncls <- paste(grp.nm, as.numeric(cls)); # note, here to retain ordered factor
    nmeta <- c(nmeta, ncls);
    sample.cluster[[grps[i]]] <- order(cls)
  }
  
  # convert back to numeric
  nmeta <- as.numeric(as.factor(nmeta))+99;
  unik.inx <- !duplicated(nmeta)   
  
  # get corresponding names
  meta_anot <- meta.vec[unik.inx]; 
  names(meta_anot) <- nmeta[unik.inx]; # name annotatation by their numbers
  
  nmeta <- matrix(nmeta, ncol=ncol(meta), byrow=F);
  colnames(nmeta) <- grps;
  
  # for each gene/row, first normalize and then tranform real values to 30 breaks 
  
  # note, use {} will lose order; use [[],[]] to retain the order
  
  gene.id = orig.gene.nms; if(length(gene.id) ==1) { gene.id <- matrix(gene.id) };
  hit.inx <- match(gene.id , unname(enrich.nms1));
  lbls <- names(enrich.nms1[hit.inx]);

  # re-order by p-value
  gene.id <- gene.id[rankPval];
  lbls <- lbls[rankPval];
  
  json.res <- list(
    data.name = dataSet$name,
    data.type = dataSet$type,
    gene.id = gene.id,
    gene.name = lbls,
    gene.cluster = gene.cluster,
    sample.cluster = sample.cluster,
    sample.names = orig.smpl.nms,
    meta = data.frame(nmeta),
    meta.anot = meta_anot,
    data = rest,
    org = data.org,
    pval = stat.pvals
  );
  
  return(json.res);
}

PlotDiagnosticPca <- function(imgNm, dpi=72, format="png",type="diablo"){
  
  #save.image("TestDI.RData");
  require("Cairo");
  library(ggplot2);
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  fig.list <- list()
  reductionSet <- .get.rdt.set();
  if(type == "diablo"){ 
    library(grid)
    library(gridExtra)
    library(gridGraphics);
    library(cowplot);
    Cairo(file=imgNm, width=10, height=24, type=format, bg="white", unit="in", dpi=dpi);
    fig.list[[1]] <- as_grob(function(){
      plotDiablo(reductionSet$dim.res, ncomp = 1)
    })
    fig.list[[2]] <- as_grob(function(){
      plotDiablo(reductionSet$dim.res, ncomp = 2)
    })
    fig.list[[3]] <- as_grob(
      plotDiablo(reductionSet$dim.res, ncomp = 3)
    )
    grid.arrange(grobs =fig.list, nrow=length(fig.list))
    dev.off();  
  }else if(type == "rcca" || type == "spls"){
    res = reductionSet$dim.res 
    Factor <- as.factor(reductionSet$meta$newcolumn)
    library(ggplot2)
    library("ggpubr")
    scrs <- list()
    for(i in 1:length(res$variates)){
      
      pca.rest <- as.data.frame(res$variates[[i]][,c(1:3)])
      colnames(pca.rest) <- c("PC1", "PC2", "PC3")
      pca.rest$Conditions <- Factor
      pca.rest$names <- rownames(res$variates[[i]])
      
      xlim <- GetExtendRange(pca.rest[,1]);
      ylim <- GetExtendRange(pca.rest[,2]);
      if("prop_expl_var" %in% names(res)){
        var.vec <- res$prop_expl_var[[i]] 
      }else{
        var.vec <- res$explained_variance[[i]] 
      }
      # proe <- signif(as.numeric(var.vec), 4)
      proe <- signif(var.vec/sum(var.vec), 4)
      
      xlabel <- paste0("Variate 1 (", proe[1]*100,"%)")
      ylabel <- paste0("Variate 2 (", proe[2]*100,"%)")
      
      if(i == 1){
        fig <- ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions)) +
          geom_point(size=3, alpha=0.5) + 
          xlim(xlim) + 
          ylim(ylim) + 
          xlab(xlabel) + 
          ylab(ylabel) + 
          theme_bw() +
          theme(legend.position = "none") + 
          ggtitle(reductionSet$omicstype[[1]])
        fig.list[[i]] <- fig
      } else {
        fig <- ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions)) +
          geom_point(size=3, alpha=0.5) + 
          xlim(xlim) + 
          ylim(ylim) + 
          xlab(xlabel) + 
          ylab(ylabel) + 
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()) + 
          ggtitle(reductionSet$omicstype[[2]]) +
          theme_bw()
        fig.list[[i]] <- fig
      }
    }
    
    h<-8
    Cairo(file=imgNm, width=14, height=h, type=format, bg="white", unit="in", dpi=dpi);
    p1 <- ggarrange(plotlist=fig.list, ncol = 2, nrow = 1, widths=c(7,8));
    print(p1)
    dev.off();
  }else if(type == "mcia"){
    library(omicade4)
    mcoin <- reductionSet$dim.res
    h<-8
    Cairo(file=imgNm, width=10, height=h, type=format, bg="white", unit="in", dpi=dpi);
    plot.mcoin(type=1, mcoin, phenovec=reductionSet$meta$newcolumn, sample.lab=FALSE, df.color=length(names(mdata.all)))
    dev.off();
  } else if(type == "mbpca"){
    res = reductionSet$dim.res 
    scrs <- moaScore(res);
    scr <- as.data.frame(scrs[,c(1:3)])
    Factor <- as.factor(reductionSet$meta$newcolumn)
    pca.rest <- scr
    pca.rest$Conditions <- Factor
    pca.rest$names <- rownames(scr)
    
    xlim <- GetExtendRange(pca.rest$PC1);
    ylim <- GetExtendRange(pca.rest$PC2);
    
    xlabel <- paste0("PC1")
    ylabel <- paste0("PC2")
    
    library(ggplot2)
    pcafig <- ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions)) +
      geom_point(size=3, alpha=0.5) + xlim(xlim) + ylim(ylim)+ xlab(xlabel) + ylab(ylabel) +
      theme_bw()
    
    Cairo(file=imgNm, width=10, height = 8, type=format, bg="white", unit="in", dpi=dpi);
    print(pcafig)
    dev.off();
  } else if(type == "procrustes"){
    library(ggplot2)
    library(grid)
    pro.test <- reductionSet$dim.res[[1]]
    pct <- pro.test$svd$d
    ctest <- data.frame(rda1=pro.test$Yrot[,1], rda2=pro.test$Yrot[,2], xrda1=pro.test$X[,1],
                        xrda2=pro.test$X[,2],Type=reductionSet$newmeta[,"omics"], Conditions = reductionSet$newmeta[,1])
    xlabel <- paste0("Component 1 ", "(" , signif(pct[1],4), ")")
    ylabel <- paste0("Component 2 ", "(" , signif(pct[2],4), ")")
    
    p <- ggplot(ctest) +
      geom_point(aes(x=rda1, y=rda2, colour=Conditions, shape=Type)) +
      geom_point(aes(x=xrda1, y=xrda2, colour=Conditions, shape=Type)) +
      geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=Conditions), alpha=0.4,arrow=arrow(length=unit(0.1,"cm"))) + 
      xlab(xlabel) + ylab(ylabel) +
      theme_bw()
    Cairo(file=imgNm, width=10, height=10, type=format, bg="white", unit="in", dpi=dpi);
    print(p)
    dev.off();
  }
}

PlotDiagnosticLoading <- function(imgNm, dpi=72, format="png",type="diablo"){
  require("Cairo");
  library(ggplot2)
  reductionSet <- .get.rdt.set();
  dpi<-as.numeric(dpi);
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  if(type == "diablo"){
    library(grid)
    library(gridExtra)
    library(cowplot)
    fig.list <- list()
    fig.list[[1]] <- as_grob(function(){
      plotLoadings(reductionSet$dim.res, ndisplay=10, comp = 1, contrib="max", method="median", size.name=1.1, legend=T)
    })
    fig.list[[2]] <- as_grob(function(){
      plotLoadings(reductionSet$dim.res, ndisplay=10, comp = 2, contrib="max", method="median", size.name=1.1, legend=T)
    })
    fig.list[[3]] <-as_grob(function(){
      plotLoadings(reductionSet$dim.res, ndisplay=10, comp = 3, contrib="max", method="median", size.name=1.1, legend=T)
    })
    h <- 8*round(length(fig.list))
    Cairo(file=imgNm, width=13, height=h, type=format, bg="white", unit="in", dpi=dpi);
    grid.arrange(grobs =fig.list, nrow=length(fig.list))
    dev.off();
  }else if(type == "rcca" || type == "spls"){
    Cairo(file=imgNm, width=12, height=10, type=format, bg="white", unit="in", dpi=dpi);
    plotVar(reductionSet$dim.res, comp = 1:2, cutoff = 0.5, var.names = c(TRUE, TRUE),
            cex = c(4, 4), title = 'rCCA comp 1 - 2')
    dev.off();
  }else if(type == "mcia"){
    library(omicade4)
    mcoin <- reductionSet$dim.res
    
    Cairo(file=imgNm, width=10, height=10, type=format, bg="white", unit="in", dpi=dpi);
    plot.mcoin(type=2, mcoin, phenovec=reductionSet$cls, sample.lab=FALSE, df.color=1:length(names(mdata.all)))
    dev.off();
  }else if(type == "mbpca"){
    library(ggplot2)
    moa <- reductionSet$dim.res 
    loading <- moa@loading[,c(1:3)]
    loading <- as.data.frame(loading)
    colnames(loading) = c("PC1", "PC2", "PC3")
    d.types <- rownames(moa@RV)
    loading$Type <- gsub(".*_", "", rownames(loading))
    for(i in 1:length(d.types)){
      rownames(loading) = gsub(paste0("_", d.types[i]), "",rownames(loading))
    }
    loading$Type <- as.factor(loading$Type)
    
    xlim <- GetExtendRange(loading[,1]);
    ylim <- GetExtendRange(loading[,2]);
    
    xlabel <- paste0("PC1")
    ylabel <- paste0("PC2")
    
    pcafig <- ggplot(loading, aes(x=PC1, y=PC2,  color=Type)) +
      geom_point(size=3, alpha=0.5) + 
      xlim(xlim) + ylim(ylim) + xlab(xlabel) + ylab(ylabel) +
      theme_bw()
    
    Cairo(file=imgNm, width=10, height=10, type=format, bg="white", unit="in", dpi=dpi);
    print(pcafig)
    dev.off();
  }
}


PlotHeatmapDiagnosticPca <- function(imgNm, dpi=72, format="png",type="spectrum"){
  require("Cairo");
  library(ggplot2)
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  sel.nms <- names(mdata.all)
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readRDS(sel.nms[i])
    data.list[[i]] <- dat$data.proc
  }
  reductionSet <- .get.rdt.set();
  
  fig.list <- list();
  result <- reductionSet$clustRes
  clust <- reductionSet$clustVec
  
  fig.list <- list()
  for(i in 1:length(sel.nms)){
    dataSet = readRDS(sel.nms[i])
    x <- dataSet$data.proc
    pca <- prcomp(t(na.omit(x)));
    imp.pca<-summary(pca)$importance;
    xlabel <- paste0("PC1"," (", 100*round(imp.pca[2,][1], 3), "%)")
    ylabel <- paste0("PC2"," (", 100*round(imp.pca[2,][2], 3), "%)")
    names <- colnames(x);
    pca.res <- as.data.frame(pca$x);
    pca.res <- pca.res[,c(1,2)]
    
    xlim <- GetExtendRange(pca.res$PC1);
    ylim <- GetExtendRange(pca.res$PC2);
    
    Factor <- as.factor(clust)
    pca.rest <- pca.res
    pca.rest$Cluster <- Factor
    pca.rest$names <- rownames(pca.res)
    Conditions = as.factor(dataSet$cls)
    pca.rest$Conditions <- Conditions
    
    pcafig <- ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Cluster, shape=Conditions)) +
      geom_point(size=3, alpha=0.5) + 
      xlim(xlim) + ylim(ylim) + xlab(xlabel) + ylab(ylabel) +
      theme_bw()
    fig.list[[i]] <- pcafig
    
    pcafigcls <- ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions)) +
      geom_point(size=3, alpha=0.5) + 
      xlim(xlim) + ylim(ylim) + xlab(xlabel) + ylab(ylabel) +
      theme_bw()
    if(i == 1){
      
      #fig.list[[3]] <- pcafigcls
    }else if( i == 2){
      #fig.list[[4]] <- pcafigcls
    }
  }
  h <- 6*round(length(fig.list)/2)
  Cairo(file=imgNm, width=14, height=h, type=format, bg="white", unit="in", dpi=dpi);
  library("ggpubr")
  p1 <- ggarrange(plotlist=fig.list, ncol = 2, nrow = round(length(fig.list)/2), labels=sel.nms)
  print(p1)
  dev.off();
}

ComputeSpectrum <- function(method="1", clusterNum="-1"){
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readRDS(sel.nms[i])
    data.list[[i]] <- dat$data.proc
  }
  reductionSet <- .get.rdt.set();
  reductionSet$omicstype <- sel.nms
  
  clusterNum <- as.numeric(clusterNum)
  library("Spectrum")
  if(clusterNum == -1){
    clusterNum = NULL;
    if(method == "eigengap"){
      method=1;
    }else if (method == "multimodality"){
      method=2;
    }
  }else{
    method=3;
  }
  
  res <- Spectrum::Spectrum(data.list, method=method, fixk=clusterNum, maxk=10, missing=TRUE, fontsize=8, dotsize=2, showres=F, silent=T)
  clust <- res$assignments
  
  SNFNMI = .calNMI(clust, as.numeric(dat$cls));
  
  reductionSet$clustType <- "Spectrum"
  reductionSet$clustVec <- clust
  reductionSet$clustRes <- res
  reductionSet$clustDistMat <- res$similarity_matrix
  reductionSet$clustNmi <- SNFNMI
  .set.rdt.set(reductionSet);
  
  return(1)
}

ComputePins <- function(method="kmeans", clusterNum="auto"){
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readRDS(sel.nms[i])
    data.list[[i]] <- t(dat$data.proc)
  }
  reductionSet <- .get.rdt.set();
  reductionSet$omicstype <- sel.nms
  clusterNum <- as.numeric(clusterNum)
  if(clusterNum == -1){
    clusterNum = NULL;
  }
  
  library("PINSPlus")
  result <- SubtypingOmicsData(dataList = data.list, clusteringMethod = "hclust", k = clusterNum, verbose=F, agreementCutoff=0.8)
  clust <- result$cluster1
  
  SNFNMI = .calNMI(clust, as.numeric(dat$cls))
  
  reductionSet$clustType <- "Perturbation"
  reductionSet$clustVec <- clust
  reductionSet$clustRes <- result
  reductionSet$clustNmi <- SNFNMI
  
  .set.rdt.set(reductionSet);
  return(1)
}


plot.mcoin <- function(type=1, x, axes=1:2, sample.lab=TRUE, sample.legend=TRUE, sample.color=1, phenovec=NULL, df.color=1, df.pch=NA, gene.nlab=0, ...) {
  mcoin <- x
  if (!inherits(mcoin, "mcia"))
    stop("mcia object expected, please run mcia first")
  
  ndata <- length(mcoin$coa)
  eig <- mcoin$mcoa$pseudoeig
  cov2 <- mcoin$mcoa$cov2[, axes] #pseueig
  
  if (!length(sample.color) %in% c(1, nrow(mcoin$mcoa$SynVar)))
    stop("length of sample.color should be either 1 or # of samples")
  if (!length(df.color) %in% c(1, ndata))
    stop("length of df.color should be either 1 or # of samples")
  
  if (is.na(df.pch[1])) {
    pch <- c(16, 17, 15, 18, 1, 2, 0, 5)
    if (ndata > 8) {
      pch <- rep(pch, ceiling(ndata/8))[1:ndata]
      warning("more than 8 datasets in mcia, df.pch is recycled used")
    } else {
      pch <- pch[1:ndata]
    }
  } else {
    if (length(df.pch) != ndata)
      stop("the length of df.pch should be the same with datasets in mcia, recycled use df.pch is not allowed")
    pch <- df.pch
  }
  
  #layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
  if(type == 1){
    splot.sample.mcia(mcoin, axis1=axes[1], axis2=axes[2], 
                      col=sample.color, pch=pch,
                      sample.lab=sample.lab,
                      legend=sample.legend,
                      sub="sample space",
                      phenovec=phenovec)
  }else if (type == 2){
    
    splot.mol.mcia(mcoin, axis1=axes[1], axis2=axes[2], col=df.color, pch=pch, gene.nlab=gene.nlab)
  }else if (type == 3){
    # plot c(2, 1)
    # par(mar=c(3, 3, 3, 3))
    # barplot(eig, col="black",names.arg=paste("eig", 1:length(eig)))
    # legend(x="topright", pch=pch, col=df.color, legend=names(mcoin$coa), box.col="white")
    
    nkeig <- ncol(mcoin$mcoa$Tli)
    eig <- mcoin$mcoa$pseudoeig[c(1:10)]
    proe <- eig/sum(eig)
    
    par(mar=c(3, 4, 1, 4))
    neig <- length(eig)
    po <- barplot(eig, plot=FALSE)
    barplot(eig,names.arg= 1:neig, col=c(rep("cyan", nkeig), rep("gray", neig-nkeig)), xlab="", ylab="Eigen Value")
    
    par(new=T) 
    plot(cbind(po, proe), frame.plot=FALSE, pch=20, axes=FALSE, xlab="Eigen Vector", ylab="")
    points(x=po, y=proe, pch=20)
    lines(x=po, y=proe)
    axis(side=4)
    mtext(side=4, text="Percentage", line=2.5, cex=0.8)
    legend(x="topright", pch=pch, col=df.color, legend=names(mcoin$coa), box.col="white")
  }else{
    #plot c(2, 2)
    par(mar=c(4.5, 4.5, 0.5, 0.5))
    plot(cov2, pch=".", main="", axes=TRUE, col=NA, 
         xlab=paste("pseudoeig", axes[1]), ylab=paste("pseudoeig", axes[2]))
    scatterutil.grid(0)
    points(cov2, pch=pch, col=df.color)
  }
}

splot.sample.mcia <- function(x, axis1=1, axis2=2, 
                              col=1, pch=20, sample.lab=TRUE, 
                              legend=TRUE, sub="",phenovec=NULL) {
  
  # plot matched samples from mcia
  if (!inherits(x, "mcia"))
    stop("x should be an object of class mcia!")
  if (!is.null(phenovec))
    phenovec <- as.factor(phenovec)
  ndata <- length(x$coa)
  dfxy <- x$mcoa$Tl1
  syn <- x$mcoa$SynVar
  pl <- pch
  cl <- col
  
  if (!is.null(phenovec) && length(phenovec) != nrow(syn))
    stop("the length of phenovec should be the same with # of samples")
  if (!axis1 %in% 1:ncol(dfxy))
    stop("Uncharacterized axis selected by axis1")
  if (!axis2 %in% 1:ncol(dfxy))
    stop("Uncharacterized axis selected by axis1")
  if (!length(col) %in% c(1, nrow(syn), length(levels(phenovec))))
    stop("length of col should be either 1 or # of samples")
  if (!length(pch) %in% c(1, ndata))
    stop("length of pch should be either 1 or # of data frame")
  
  sync <- c()
  for (i in 1:(ndata)) {
    sync <- rbind(sync, syn)
  }
  c <- x$mcoa$TL$"T"
  
  if (length(col) == 1) {
    if (is.null(phenovec))
      col <- rep(col, length(c)) else
        col <- c(phenovec)
  } else if (length(col) == length(levels(phenovec))) {
    if (!is.null(phenovec))
      col <- col[c(phenovec)]
  } else
    col <- rep(col, ndata)
  
  if (length(pch) == 1)
    pch <- rep(pch, length(c)) else
      pch <- rep(pch, table(c))
  
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  coo <- scatterutil.base(dfxy = dfxy, xax = axis1, yax = axis2, sub = sub, 
                          xlim = NULL, ylim = NULL, grid = TRUE, 
                          addaxes = TRUE, cgrid = 1, include.origin = TRUE, 
                          origin = c(0, 0), csub = 1.25, possub = "bottomleft", 
                          pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE)
  
  points(dfxy[, c(axis1, axis2)], pch=pch, col=col)
  segments(sync[, axis1], sync[, axis2], dfxy[, axis1], dfxy[, axis2], col=col)
  
  if (sample.lab) {
    lab <- rownames(syn)
    text(x=syn[, axis1], y=syn[, axis2], lab)
  }
  
  if (legend && any(c(length(cl) != 1, length(pl) != 1, !is.null(phenovec)))) {
    
    ple <- c()
    if (length(pl) != 1) {
      pl <- pl
      ple <- names(x$coa)
    }
    
    cle <- c()
    if (length(cl) == nrow(syn)) {
      cl <- cl
      cle <- rownames(syn)
    } else if (length(cl) == 1 && !is.null(phenovec)) {
      cl <- sort(unique(c(phenovec)))
      cle <- levels(phenovec)
    } else if (length(cl) == length(levels(phenovec))) {
      cl <- cl
      cle <- levels(phenovec)
    }
    
    pch.i <- c(pl, rep(20, length(cl)))
    if (length(pl)==1)
      col.i <- cl else
        col.i <- c(rep(1, ndata), cl)
    le.i <- c(ple, cle)
    
    legend("topleft", fill=FALSE, col=col.i, pch=pch.i, legend=le.i, border=F, bty="n")
  }
  box()
}

splot.mol.mcia <- function(x, axis1, axis2, col, pch, gene.nlab) {
  ndata <- length(x$coa)
  co <- x$mcoa$Tco[, c(axis1, axis2)]
  c <- as.numeric(x$mcoa$TC$"T")
  if (length(col) == 1)
    col <- rep(col, ndata)
  if (length(pch) == 1)
    pch <- rep(pch, ndata)
  made4::plotgenes(co, colpoints=0, nlab=gene.nlab, pch=".", axis1=1, axis2=2, sub="variable space")
  for (i in 1:ndata) {
    points(co[c %in% i, ], col=col[i], pch=pch[i])
  }
}


GetDiagnosticSummary<- function(type){
  if(type %in% c("perturbation", "spectrum", "snf")){
    reductionSet <- .get.rdt.set();
    clustNum <- length(unique(reductionSet$clustVec))
    return(c(clustNum, signif(reductionSet$clustNmi)))
  }else if(type == "procrustes"){
    reductionSet <- .get.rdt.set();
    res <- reductionSet$dim.res[[2]];
    return(c(signif(res$ss,4), signif(res$scale,4)));
  }else{
    return(c("","") )
  }
}

plotEig <- function(nclust=3, eig.vec){
  nkeig <- as.numeric(nclust)
  eig <- eig.vec
  proe <- eig/sum(eig)
  par(mar=c(3, 4, 1, 4))
  neig <- length(eig)
  po <- barplot(eig, plot=FALSE)
  barplot(eig,names.arg= 1:neig, col=c(rep("cyan", nkeig), rep("gray", neig-nkeig)), xlab="", ylab="Eigen Value")
  
  par(new=T) 
  plot(cbind(po, proe), frame.plot=FALSE, pch=20, axes=FALSE, xlab="Eigen Vector", ylab="")
  points(x=po, y=proe, pch=20)
  lines(x=po, y=proe)
  axis(side=4)
  mtext(side=4, text="Percentage", line=2.5, cex=0.8)
}

widget_file_size <- function(p, filenm) {
  d <- tempdir()
  filenm <- paste0(filenm, ".html")
  withr::with_dir(d, htmlwidgets::saveWidget(p,filenm))
  f <- file.path(d, filenm)
}

estimateClusters <- function (W, NUMC = 2:10) 
{
  if (min(NUMC) == 1) {
    warning("Note that we always assume there are more than one cluster.")
    NUMC <- NUMC[NUMC > 1]
  }
  W <- (W + t(W))/2
  diag(W) <- 0
  if (length(NUMC) <= 0) {
    warning(paste("Invalid NUMC provided, must be an integer vector", 
                  "with atleast one other number than 1.", "Using default NUMC=c(2,3,4,5)", 
                  sep = ""))
    NUMC <- 2:5
  }
  degs <- rowSums(W)
  degs[degs == 0] <- .Machine$double.eps
  D <- diag(degs)
  L <- D - W
  Di <- diag(1/sqrt(degs))
  L <- Di %*% L %*% Di
  eigs <- eigen(L)
  eigs_order <- sort(eigs$values, index.return = T)$ix
  eigs$values <- eigs$values[eigs_order]
  eigs$vectors <- eigs$vectors[, eigs_order]
  eigengap <- abs(diff(eigs$values))
  quality <- list()
  for (c_index in 1:length(NUMC)) {
    ck <- NUMC[c_index]
    UU <- eigs$vectors[, 1:ck]
    EigenvectorsDiscrete <- .discretisation(UU)[[1]]
    EigenVectors <- EigenvectorsDiscrete^2
    temp1 <- EigenVectors[do.call(order, lapply(1:ncol(EigenVectors), 
                                                function(i) EigenVectors[, i])), ]
    temp1 <- t(apply(temp1, 1, sort, TRUE))
    quality[[c_index]] <- (1 - eigs$values[ck + 1])/(1 - 
                                                       eigs$values[ck]) * sum(sum(diag(1/(temp1[, 1] + .Machine$double.eps)) %*% 
                                                                                    temp1[, 1:max(2, ck - 1)]))
  }
  t1 <- sort(eigengap[NUMC], decreasing = TRUE, index.return = T)$ix
  K1 <- NUMC[t1[1]]
  K12 <- NUMC[t1[2]]
  t2 <- sort(unlist(quality), index.return = TRUE)$ix
  K2 <- NUMC[t2[1]]
  K22 <- NUMC[t2[2]]
  output <- list(`Eigen-gap best` = K1, `Eigen-gap 2nd best` = K12, 
                 `Rotation cost best` = K2, `Rotation cost 2nd best` = K22, "eigs" =  (1-eigs$values)[NUMC])
  return(output)
}


.discretisation <- function(eigenVectors) {
  
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  
  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
    
  }
  
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}


ComputeSNF <- function(method="1", clusterNum="auto"){
  K = 20;		# number of neighbors, usually (10~30)
  alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
  T = 10; 	# Number of Iterations, usually (10~20)
  
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readRDS(sel.nms[i])
    data.list[[i]] <- dat$data.proc
  }
  reductionSet <- .get.rdt.set();
  reductionSet$omicstype <- sel.nms
  
  library("SNFtool")
  truelabel = as.numeric(dat$cls)
  Data1 <- t(data.list[[1]])
  Data2 <- t(data.list[[2]])
  
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  
  W = SNF(list(W1,W2), K, T);
  res <- estimateClusters(W, 2:10);
  clusterNum <- as.numeric(clusterNum)
  if(clusterNum == -1){
    C = res[[1]]
  }else{
    C = clusterNum
  }
  group = spectralClustering(W, C); 
  SNFNMI = .calNMI(group, truelabel)
  
  reductionSet$clustType <- "Spectrum"
  reductionSet$clustVec <- group
  reductionSet$clustRes <- res
  reductionSet$clustDistMat <- W
  reductionSet$clustNmi <- SNFNMI
  
  .set.rdt.set(reductionSet);
  return(1)
}

.calNMI<-function (x, y) 
{
  x <- as.vector(x)
  y <- as.vector(y)
  mutual.info <- (.mutualInformation(x, y)/sqrt(.entropy(x) * 
                                                  .entropy(y)))
  return(max(0, mutual.info, na.rm = TRUE))
}

ComputeSilhouette <-function(type){
  
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readRDS(sel.nms[i])
    data.list[[i]] <- dat$data.proc
  }
  
  library(cluster)
  reductionSet <- .get.rdt.set();
  clusterVec <- reductionSet$clustVec
  # function to compute average silhouette for k clusters
  if(type == "spectrum"){
    clustMatrix <- reductionSet$clustRes$similarity_matrix
  }else if(type == "snf"){
    clustMatrix <- reductionSet$clustDistMat
  }
  
  ss1 <- cluster::silhouette(clusterVec, daisy(t(data.list[[1]])))
  ss2 <- cluster::silhouette(clusterVec, daisy(t(data.list[[2]])))
  
  return(list(ss1,ss2))
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

# normalize to zero mean and unit variance
AutoNorm<-function(x){
  (x - mean(x))/sd(x, na.rm=T);
}

CheckMetaIntegrity <- function(){
  sel.nms <- names(mdata.all)
  
  data.list <- list()
  cnms <- list()
  metas <- list();
  for(i in 1:length(sel.nms)){
    dat = readRDS(sel.nms[i])
    cnms[[i]] <- colnames(dat$data.proc);
    metas[[i]] <- dat$meta;
  }
  
  if(length(metas) == 0){
    msg.vec <<- paste0('Please make sure row(s) corresponding to meta-data start with "#CLASS" or to include a metadata file.' );
    return(0)
  }
  
  for(i in 1:length(sel.nms)){
    for(j in 1:length(sel.nms)){
      bool <- SameElements(cnms[[i]], cnms[[j]])
      if(!bool){
        msg.vec <<- paste0("Please make sure the sample names are consistent and reupload the data sets. ", sel.nms[i], ": ", unique(cnms[[i]]),"; ", sel.nms[j], ": ", unique(cnms[[j]]));
        return(0)
      }
      boolMeta <- identical(metas[[i]],metas[[j]])
      print(boolMeta);
      if(!boolMeta){
        msg.vec <<- "Please make sure the meta data is consistent across all uploaded data sets.";
        return(0)
      }
    }
  }
  return(1)
  
}

SameElements <- function(a, b) return(identical(sort(a), sort(b)));

# Function to annotate metabolite data to internal database
SkippingAnnotation <- function(dataName, idtype){
  
  #if(dataSet$name != dataName){
  dataSet <- readRDS(dataName);
  #}
  
  data <- dataSet$data.raw;
  qvec <- rownames(data);
  
  dataSet$enrich_ids <- rownames(dataSet$data.raw)
  names(dataSet$enrich_ids) <- rownames(dataSet$data.raw)
  dataSet$data.annotated <- dataSet$data.raw
  
  RegisterData(dataSet);
  
  
  return(1)
}


#'Plot t-sne plot for multi-omics samples
#'@description 
#'@param imgNm name of the image to output
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'

PlotMultiTsne <- function(imgNm, dpi=72, format="png",factor="1"){
  require("Cairo");
  library(ggplot2)
  library(Rtsne)
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  sel.nms <- names(mdata.all)
  fig.list <- list()
  
  pca.list<- list()
  pct <- list();
  
  for(i in 1:length(sel.nms)){
    dataSet = readRDS(sel.nms[i])
    
    x <- t(dataSet$data.proc)
    max.perx <- floor((nrow(x)-1)/3);
    if(max.perx > 30){
      max.perx <- 30;
    }
    tsne_out <- Rtsne(x,pca=FALSE,perplexity=max.perx,theta=0.0, check_duplicates=F) # Run TSNE
    xlabel <- "tsne1"
    ylabel <-"tsne2"
    names <- colnames(x);
    pca.res <- as.data.frame(tsne_out$Y);
    pca.res <- pca.res[,c(1,2)]
    colnames(pca.res) <- c("tsne1", "tsne2")
    xlim <- GetExtendRange(pca.res[,1]);
    ylim <- GetExtendRange(pca.res[,2]);
    
    Factor <- dataSet$meta[,1];
    pca.rest <- pca.res
    pca.rest$Conditions <- Factor
    pca.rest$names <- rownames(pca.res)
    
    pcafig <- ggplot(pca.rest, aes(x=tsne1, y=tsne2,  color=Conditions)) +
      geom_point(size=3, alpha=0.5) + 
      xlim(xlim) + 
      ylim(ylim) + 
      xlab(xlabel) + 
      ylab(ylabel) +
      theme_bw()
    fig.list[[i]] <- pcafig
  }
  
  h<-6*round(length(fig.list)/2)
  Cairo(file=imgNm, width=14, height=h, type=format, bg="white", unit="in", dpi=dpi);
  library("ggpubr")
  p1 <- ggarrange(plotlist=fig.list, ncol = 2, nrow = round(length(fig.list)/2), labels=sel.nms)
  print(p1)
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
  #print(propertiesVector);
  return(propertiesVector);
}


