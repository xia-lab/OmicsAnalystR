

DoIntegrativeAnalysis <- function(method, sign="both", threshold=0.6, nComp){
  load_igraph();
  intRes <- DoDimensionReductionIntegrative(method);
  threshold <- as.numeric(threshold)
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx]
  data <- list()
  labels <- vector();
  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i])
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
    cor.inx <- cor_edge_list$weigth < -threshold
  }
  
  only_sig <- cor_edge_list[cor.inx, ];
  new_g <- graph_from_data_frame(only_sig, F);
  
  type.list <- list()
  for(i in 1:length(sel.nms)){
    type.list[[sel.nms[[i]]]] <- unique(cor_edge_list[,i]);
  }
  ProcessGraphFile(new_g, labels, type.list, T);
}

NormalizeDataWrapper <-function (nm, opt, colNorm){
  if(nm == "NA"){
    sel.nms <- names(mdata.all)
  }else{
    sel.nms <- c(nm);
  }
  for(i in 1:length(sel.nms)){
    dataSet = readDataset(sel.nms[i])
    data <- NormalizingDataOmics(dataSet$data.filtered, opt, colNorm, "NA")
    dataSet$data.proc <- data;
    if(exists("m2m",dataSet)){
      data.norm.taxa <- lapply(dataSet$dataSet$data.filt.taxa, function(x) {
        NormalizingDataOmics(x, opt, colNorm, "NA")
      }
      )
      dataSet$data.proc.taxa <- data.norm.taxa
    }
    fast.write(data, paste0("table_", dataSet$name));
    RegisterData(dataSet)
  }
  return(1)
}

ScaleDataWrapper <-function (nm, scaleNorm){
  if(nm == "NA"){
    sel.nms <- names(mdata.all)
  }else{
    sel.nms <- c(nm);
  }
  for(i in 1:length(sel.nms)){
    dataSet = readDataset(sel.nms[i])
    data <- NormalizingDataOmics(dataSet$data.proc, "NA", "NA", scaleNorm)
    dataSet$data.proc <- data;
    if(exists("m2m",dataSet)){
      data.norm.taxa <- lapply(dataSet$dataSet$data.proc.taxa, function(x) {
        NormalizingDataOmics(x, "NA", "NA", scaleNorm)
      })
      dataSet$data.proc.taxa <- data.norm.taxa
    }
    RegisterData(dataSet)
  }
  return(1);
}

FilterDataMultiOmicsHarmonization <- function(dataName,filterMethod, filterPercent = 0){
  filterPercent <- as.numeric(filterPercent);
  if(dataName == "NA"){
    sel.nms <- names(mdata.all)
  } else {
    sel.nms <- c(dataName);
  }
  
  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i])
    int.mat <- qs::qread(dataSet$data.annotated.path);
    int.mat <- int.mat[,colnames(int.mat) %in% colnames(dataSet$data.proc)];

    if(filterMethod == "variance"){
      data <- FilterDataByVariance(int.mat, filterPercent);
    }else{
      featVar <- apply(data, 1, var);
      if(var(featVar) < 0.001){
        print("Detected autoscale");
        msg.vec <<- paste0(dataSet$name, " appears to be autoscaled. Filtering can not be performed on autoscaled dataset!");
        return(2);
      }
      res <- PerformFeatureFilter(t(int.mat), filterMethod, filterPercent, "", T);
      data <- t(res$data);
    }
    if(any(class(data) == "character")){
      msg.vec <<- paste0(dataSet$name, " appears to be autoscaled. Filtering can not be performed on autoscaled dataset!");
      print("Detected autoscale");
      return(2)
    }
    dataSet$data.proc <- data;
    if(exists("m2m",dataSet)){
      data.norm.taxa <- lapply(dataSet$dataSet$data.annotated.taxa, function(x) {
        if(filterMethod == "variance"){
          data <- FilterDataByVariance(x, filterPercent);
        }else{
          filterPercent <- filterPercent/100;
          res <- PerformFeatureFilter(t(x), filterMethod, filterPercent, "", T);
          data <- t(res$data);
        }
        return(data);
      })
      dataSet$data.proc.taxa <- data.norm.taxa
    }

    RegisterData(dataSet)
  }
  return(1)
}

FilterDataByVariance <- function(data, filterPercent){
  featVar <- apply(data, 1, var);
  if(var(featVar) < 0.001){
    return("Already autoscaled");
  }
  varThresh <- quantile(featVar, (filterPercent/100));
  featKeep <- which(featVar > varThresh);
  data <- data[featKeep, ];
  return(data);
}

#'Plot PCA plot for multi-omics samples
#'@description 
#'@param imgNm name of the image to output
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'

PlotMultiPCA <- function(imgNm, dpi=72, format="png",factor="1", interactive=F){
  load_cairo();
  load_ggplot();
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  sel.nms <- names(mdata.all)
  fig.list <- list()
  
  pca.list<- list()
  pct <- list();
  all_data <- list()

  for(i in 1:length(sel.nms)){
    print(i)
    dataSet = readDataset(sel.nms[i])
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
    pca.rest$dataset <- dataSet$name # This is the new column specifying the dataset/source
    all_data[[i]] <- pca.rest

    #pcafig <- ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions)) +
    #  geom_point(size=3, alpha=0.5) + 
    #  xlim(xlim) + 
    #  ylim(ylim) + 
    #  xlab(xlabel) + 
    #  ylab(ylabel) +
    #  ggtitle(dataSet$name)+
    #  theme_bw()+
    #  theme(text=element_text(size=13), plot.title = element_text(hjust = 0.5));
    
    #fig.list[[i]] <- pcafig
    
    #computing 3d pca
    #pos.xyz = data.frame(x=pca$x[,1], y=pca$x[,2], z=pca$x[,3]);
    #loading.pos.xyz = data.frame(pca$rotation[,c(1:3)]);
    #loadingNames = rownames(pca$rotation);
    
    #pos.xyz <- as.data.frame(pos.xyz);
    #pos.xyz <- unitAutoScale(pos.xyz);
    #rownames(pos.xyz) = rownames(dataSet$meta);
    #dataSet$pos.xyz = pos.xyz;
    
    #loadingNames <- rownames(loading.pos.xyz);
    #loading.pos.xyz <- as.data.frame(loading.pos.xyz);
    #loading <- unitAutoScale(loading.pos.xyz);
    #rownames(loading) <- loadingNames;
    #nm <- paste0("pca_", dataSet$type);
    #pca.list[[nm]][["score"]] <- pos.xyz * 1000;
    #pca.list[[nm]][["loading"]] <- loading* 1000;    
    
    #pct2Nm <- paste0("pca_", dataSet$type)
    #pct[[pct2Nm]] <- unname(round(imp.pca[2,],3))[c(1:3)]*100;
  }

  pca.list$pct2 <- pct;
  qs::qsave(pca.list, file="pca.scatter.qs");
  
  #h<-6*round(length(sel.nms)/2)
  h<-6
  #library("ggpubr")
  #p1 <- ggarrange(plotlist=fig.list, ncol = 2, nrow = round(length(fig.list)/2))
  combined_data <- do.call(rbind, all_data)
  
  p1 <- ggplot(combined_data, aes(x=PC1, y=PC2, color=Conditions)) +
    geom_point(size=3, alpha=0.5) + 
    facet_wrap(~ dataset, scales = "free") + # Use facet_wrap or facet_grid
    theme_bw() +
    theme(text=element_text(size=13))


  if(interactive){
    library(plotly);
        m <- list(
                l = 50,
                r = 50,
                b = 20,
                t = 20,
                pad = 0.5
            )

    ggp_build <- layout(ggplotly(p1), autosize = FALSE, width = 1000, height = 600, margin = m)
    return(ggp_build);
  }else{
  Cairo(file=imgNm, width=14, height=h, type=format, bg="white", unit="in", dpi=dpi);
  print(p1)
  dev.off();
  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$qc_multi_pca <- imgNm;
  saveSet(infoSet);
  }
}

PlotMultiDensity <- function(imgNm, dpi=72, format="png",factor="1", interactive=F){
  load_cairo();
  load_ggplot();
  dpi <- as.numeric(dpi)
  factor <- as.numeric(factor)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  sel.nms <- names(mdata.all)
  fig.list <- list()
  # Pre-allocate list to avoid sequential rbind (memory optimization)
  df.list <- vector("list", length(sel.nms))

  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i])
    dat <- dataSet$data.proc
    st <- stack(as.data.frame(dat))
    st$Dataset <- rep(sel.nms[i], nrow(dat))
    df.list[[i]] <- st  # Store in list instead of sequential rbind
  }

  # Combine all data frames at once (instead of sequential rbind)
  merged.df <- do.call(rbind, df.list)
  
  type<-merged.df$type
  merged.df$ind <- paste0(merged.df$ind, "_", merged.df$type)
    g =ggplot(merged.df, aes(x=values)) + 
    geom_line(aes(color=Dataset, group=ind), stat="density", alpha=0.1) + 
    geom_line(aes(color=Dataset), stat="density", alpha=0.7, size=3) +
    theme_bw() +
    theme(text=element_text(size=13))



  if(interactive){
    library(plotly);
        m <- list(
                l = 50,
                r = 50,
                b = 20,
                t = 20,
                pad = 0.5
            )

    ggp_build <- layout(ggplotly(g), autosize = FALSE, width = 1000, height = 600, margin = m)
    return(ggp_build);
  }else{
  Cairo(file=imgNm, width=10, height=6, type=format, bg="white", dpi=dpi, unit="in");
  print(g)
  dev.off();
  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$qc_multi_density <- imgNm;
  saveSet(infoSet);
  }
} 

CheckMetaIntegrity <- function(){
  sel.nms <- names(mdata.all)
  
  data.list <- list()
  cnms <- list()
  metas <- list();
  for(i in 1:length(sel.nms)){
    dat = readDataset(sel.nms[i])
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
      #print(boolMeta);
      if(!boolMeta){
        msg.vec <<- "Please make sure the meta data is consistent across all uploaded data sets.";
        return(0)
      }
    }
  }
  return(1)
  
}

SameElements <- function(a, b) return(identical(sort(a), sort(b)));

#'Plot t-sne plot for multi-omics samples
#'@description 
#'@param imgNm name of the image to output
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'

PlotMultiTsne <- function(imgNm, dpi=72, format="png",factor="1"){
  load_cairo();
  load_ggplot();
  library(Rtsne)
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  sel.nms <- names(mdata.all)
  fig.list <- list()
  
  pca.list<- list()
  pct <- list();
  
  for(i in 1:length(sel.nms)){
    dataSet = readDataset(sel.nms[i])
    
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
