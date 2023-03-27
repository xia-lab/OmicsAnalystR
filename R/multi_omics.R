

DoIntegrativeAnalysis <- function(method, sign="both", threshold=0.6, nComp){
  require("igraph");
  if(! method %in% dim.res.methods){
    intRes <- DoDimensionReductionIntegrative("omics", method, 3);
  }
  
  threshold <- as.numeric(threshold)
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx]
  data <- list()
  labels <- vector();
  for(i in 1:length(sel.nms)){
    dataSet = qs::qread(sel.nms[i])
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
      dataSet = qs::qread(sel.nms[i])
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
    
    dataSet <- qs::qread(nm);
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
      dataSet = qs::qread(sel.nms[i])
      data <- NormalizingDataOmics(dataSet$data.proc, "NA", "NA", scaleNorm)
      dataSet$data.proc <- data;
      if(exists("m2m",dataSet)){
        data.norm.taxa <- lapply(dataSet$dataSet$data.proc.taxa, function(x) {
          NormalizingDataOmics(x, "NA", "NA", scaleNorm)
        }
        )
        dataSet$data.proc.taxa <- data.norm.taxa
      }
      RegisterData(dataSet)
    }
    return(1)
  }else{
    
    dataSet <- qs::qread(nm);
    
    data <- NormalizingDataOmics(dataSet$data.proc, "NA", "NA", scaleNorm)
    dataSet$data.proc <- data;
    if(exists("m2m",dataSet)){
      
      data.norm.taxa <- lapply(dataSet$data.proc.taxa, function(x) {
        NormalizingDataOmics(x, "NA", "NA", scaleNorm)
      }
      )
      
      
      dataSet$data.proc.taxa <- data.norm.taxa
      
      
    }
    
    RegisterData(dataSet)
    return(1)
  }
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
      cls.lbls <- ClearFactorStrings(cls.lbls, cls.nm);
      
      meta.info[[cls.nm]] <- cls.lbls;
    }
    meta.info <- data.frame(meta.info);
    data = data[-cls.inx,];
    rownames(dataSet$meta) <- colnames(data)[-1];
  }

  dataSet$meta <- meta.info;
  
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
  return(1)
}

SetOmicsType <- function(fileName, omics.type=NA) {
  dataSet <- qs::qread(fileName)
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
    dataSet = qs::qread(sel.nms[i])
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
    dataSet <- qs::qread(sel.nms[i])
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
    dataSet = qs::qread(sel.nms[i])
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
  dataSet <- qs::qread(dataName);
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



PlotDimredVarexp <- function(imgNm, dpi=72, format="png"){

  require("Cairo");
  library(ggplot2)
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  reductionSet <- .get.rdt.set();

  df <- reductionSet$var.exp;
  df <- reshape2::melt(df)
  colnames(df) <- c("Factor", "Dataset", "value")
  min_r2 = 0
  max_r2 = max(df$value)
  
  p1 <- ggplot(df, aes_string(y="Factor", x="Dataset")) + 
    geom_tile(aes_string(fill="value"), color="black") +
    labs(x="", y="", title="") +
    scale_fill_gradientn(colors=c("gray97","darkblue"), guide="colorbar", limits=c(min_r2,max_r2)) +
    guides(fill=guide_colorbar("Var. (%)")) +
    theme(
      axis.text.x = element_text(size=rel(1.0), color="black"),
      axis.text.y = element_text(size=rel(1.1), color="black"),
      axis.line = element_blank(),
      axis.ticks =  element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size=rel(1.0))
    )
  
  Cairo(file=imgNm, width=10, height=8, type=format, bg="white", unit="in", dpi=dpi);
  print(p1)
  dev.off();
}

PlotDimredFactors <- function(meta, pc.num = 5, imgNm, dpi=72, format="png"){
  
  require("Cairo");
  library(ggplot2)
  library(GGally)
  library(grid)
  
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  sel.nms <- names(mdata.all)
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = qs::qread(sel.nms[i])
    data.list[[i]] <- dat$data.proc
  }
  reductionSet <- .get.rdt.set();
  
  pclabels <- paste0("Factor ", 1:pc.num);
  
  # draw plot
  Cairo::Cairo(file = imgNm, unit="in", dpi=dpi, width=10, height=10, type=format, bg="white");
  
  data <- as.data.frame(reductionSet$pos.xyz[,1:pc.num]);
  meta.info <- reductionSet$meta;
  meta.info <- meta.info[match(rownames(data), rownames(meta.info)),]
  

  inx <- which(colnames(meta.info) == meta)
  cls <- meta.info[, inx];
  #cls.type <- mSetObj$dataSet$meta.types[inx] ##### UPDATE THIS AFTER SUPPORT COMPLEX META
  cls.type <- "disc";
  
  if (cls.type == "disc"){ ## code to execute if metadata class is discrete
    
    uniq.cols <- GetColorSchema(unique(cls))
    
    # make plot
    p <- ggpairs(data, 
                 lower = list(continuous = wrap("points")), 
                 upper = list(continuous = wrap("density")),
                 diag = list(continuous = wrap("densityDiag", alpha = 0.5, color = NA)),
                 columnLabels = pclabels, mapping = aes(color = cls))
    
    auxplot <- ggplot(data.frame(cls = cls),aes(x=cls,y=cls,color=cls)) + 
      theme_bw() + geom_point(size = 6) + theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=11)) + 
      scale_color_manual(values = uniq.cols) + guides(col = guide_legend(nrow = 1))
    p <- p + theme_bw() + scale_color_manual(values = uniq.cols) + scale_fill_manual(values = uniq.cols) + 
      theme(plot.margin = unit(c(0.25, 0.25, 0.6, 0.25), "in"))
    mylegend <- grab_legend(auxplot)
    
  } else { ## code to excute if metadata class is continuous
    
    colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(20));
    num.cls <- as.numeric(as.character(cls));
    cols <- colors[as.numeric(cut(num.cls,breaks = 20))];
    
    # make plot
    p <- ggpairs(data, lower = list(continuous = wrap("points", color = cols)), 
                 upper = list(continuous = wrap("density", color = "#505050")),
                 diag = list(continuous = wrap("densityDiag", fill = "#505050", color = NA)),
                 columnLabels = pclabels)
    
    auxplot <- ggplot(data.frame(cls = num.cls), aes(x=cls, y=cls, color=cls)) + 
      theme_bw() + geom_point(size = 6) + 
      theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=11)) + 
      guides(col = guide_legend(nrow = 1))
    
    p <- p + theme_bw() + theme(plot.margin = unit(c(0.25, 0.25, 0.8, 0.25), "in"))
    mylegend <- grab_legend(auxplot)
    
  }
  
  grid.newpage()
  grid.draw(p)
  vp = viewport(x=5, y=0.3, width=.35, height=.3, default.units = "in") ## control legend position
  pushViewport(vp)
  grid.draw(mylegend)
  upViewport()
  dev.off()
}


CheckMetaIntegrity <- function(){
  sel.nms <- names(mdata.all)
  
  data.list <- list()
  cnms <- list()
  metas <- list();
  for(i in 1:length(sel.nms)){
    dat = qs::qread(sel.nms[i])
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
    dataSet = qs::qread(sel.nms[i])
    
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

