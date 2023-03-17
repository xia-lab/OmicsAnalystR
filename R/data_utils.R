##################################################
## R scripts for OmicsAnalyst
## Description: Data IO functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

.get.rdt.set <- function(){
  return(qs::qread("rdt.set.qs"));
}

.set.rdt.set <- function(my.set){
  qs::qsave(my.set, file="rdt.set.qs");
}

#'Initialize resources for analysis
#'@description call this function before performing any analysis
#'@param onWeb whether the script is running in local or on web
#'@author Guangyan Zhou \email{guangyan.zhou@mail.mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
Init.Data <- function(){ 
  # to control parallel computing for some packages
  Sys.setenv("OMP_NUM_THREADS" = 2); 
  Sys.setenv("OPENBLAS_NUM_THREADS" = 2);

  dim.res.methods <<- vector();
  jsonNms <<- list()
  reductionSet<- list()
  reductionSet$clustVec <- "NA";
  fileTypeu <<- "NA"
  partialToBeSaved <<- c("Rload.RData", "Rhistory.R")
  regids <<- vector()
  rcmd.vecu <<- vector()
  table.nmu <<- ""
  isKo <<- F
  integOpts <<- c("mcia")
  merged.reduction <<- F
  reductionOptGlobal <<- "pca"
  dataSet <- list(annotated=FALSE);
  dataSet$misc <- list();
  dataSet$de.method <- "NA";
  anal.type <<- "multiomics";
  dataSet <<- dataSet;
  mdata.all <<- list(); 
  msg.vec <<- vector(mode="character");
  current.msg <<- "";
  lib.path <<- "../../data/";
  data.org <<- NULL;
  module.count <<- 0;
  if(file.exists("/home/glassfish/sqlite/")){
    sqlite.path <<- "/home/glassfish/sqlite/";  #public server
  }else if(file.exists("/Users/xia/Dropbox/sqlite/")){
    sqlite.path <<- "/Users/xia/Dropbox/sqlite/";  #xia local
  }else if(file.exists("/Users/jeffxia/Dropbox/sqlite/")){
    sqlite.path <<- "/Users/jeffxia/Dropbox/sqlite/";  #xia local 2
  }else if(file.exists("/home/zgy/sqlite")){
    sqlite.path <<- "/home/zgy/sqlite/";  #zzggyy local
  }else if(file.exists("/media/zzggyy/disk/")){
    sqlite.path <<- "/media/zzggyy/disk/sqlite/";  #zzggyy local
  } else if(file.exists("/Users/lzy/sqlite")){
    sqlite.path <<- "/Users/lzy/sqlite/";  #ly local
  }else {
    sqlite.path <<- "/Users/jessicaewald/sqlite/sqlite/"; #jess local
  }

  cmdSet <- list(annotated=FALSE);
  saveSet(cmdSet, "cmdSet");

  .set.rdt.set(reductionSet);
}

# this should be fixed across data sets
SetOrganism <- function(org){
  data.org <<- org;
}

#'Sanity check individual dataset for meta-analysis 
#'@param fileName The filename of dataset in qs format
#'@author Guangyan Zhou \email{guangyan.zhou@mail.mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
SanityCheckData <- function(fileName){
  dataSet <- readRDS(fileName);
  
  # general sanity check then omics specific
  
  # use first column by default
  cls <- dataSet$meta[,1]
  
  # check class info
  cls.lbl <- as.factor(as.character(cls));

  min.grp.size <- min(table(cls.lbl));
  cls.num <- length(levels(cls.lbl));
  
  msg <- paste(cls.num, "groups were defined in samples based on the first metadata column.");
  dataSet$cls.num <- cls.num;
  dataSet$min.grp.size <- min.grp.size;
  
  # check numerical matrix  
  int.mat <- dataSet$data.proc;
  rowNms <- rownames(int.mat);
  colNms <- colnames(int.mat);
  naNms <- sum(is.na(int.mat));
  
  num.mat <- apply(int.mat, 2, as.numeric);
  if(sum(is.na(num.mat)) > naNms){
    # try to remove "," in thousand seperator if it is the cause
    num.mat <- apply(int.mat,2,function(x) as.numeric(gsub(",", "", x)));
    if(sum(is.na(num.mat)) > naNms){
      msg<-c(msg,"<font color=\"red\">Non-numeric values were found and replaced by NA.</font>");
    }else{
      msg<-c(msg,"All data values are numeric.");
    }
  }else{
    msg<-c(msg,"All data values are numeric.");
  }
  
  int.mat <- num.mat;
  rownames(int.mat) <- rowNms;
  colnames(int.mat)<- colNms;
  
  # check for feats with all constant (var =0)
  varCol <- apply(int.mat, 1, var, na.rm=T);
  constCol <- (varCol == 0 | is.na(varCol));
  constNum <- sum(constCol, na.rm=T);
  
  if(constNum > 0){
    msg<-c(msg, paste("<font color=\"red\">", constNum, "features with a constant or single value across samples were found and deleted.</font>"));
    int.mat <- int.mat[!constCol, , drop=FALSE];
  }
  
  # check zero, NA values
  totalCount <- nrow(int.mat)*ncol(int.mat);
  naCount <- sum(is.na(int.mat));
  naPercent <- round(100*naCount/totalCount,1)
  
  msg <- c(msg, paste("A total of ", naCount, " (", naPercent, "%) missing values were detected.", sep=""));
  
  # obtain original half of minimal positive value (threshold)
  minConc <- min(int.mat[int.mat>0], na.rm=T)/2;
  
  # remove smpls/exp with over half missing value
  good.inx<-apply(is.na(int.mat), 2, sum)/nrow(int.mat)<0.6;
  if(sum(!good.inx)>0){
    msg <- c(msg, paste(sum(!good.inx), "Low quality samples (>60% missing) removed."));
    int.mat <- int.mat[,good.inx, drop=FALSE];
    meta.info <- dataSet$meta;
    meta.info <- meta.info[good.inx, , drop=F];
    dataSet$meta <- meta.info;
    if(ncol(int.mat)/smpl.num < 0.5){
      AddErrMsg("Too many missing values - low quality data rejected!");
      return(0);
    }
  }
  
  # remove ffeatures/variables with over half missing value    
  gd.inx <- apply(is.na(int.mat), 1, sum)/ncol(int.mat)<0.75;
  gene.num <- nrow(int.mat);
  if(sum(!gd.inx) > 0){
    int.mat <- int.mat[gd.inx,];
    msg <- c(msg, paste(sum(!gd.inx), "low quality genes (>75% missing) removed"));
    if(nrow(int.mat)/gene.num < 0.25){
      AddErrMsg("Too many missing values - low quality data rejected.");
      return(0);
    }
  }
  
  current.msg <<- msg;
  dataSet$minConc <- minConc;
  dataSet$data.proc <- int.mat;
  dataSet$cls <- cls.lbl
  
  RegisterData(dataSet);
  return(1);
}


# When multiple genelists/datasets/results, record their name and save the data as .RDS file
# a) Current dataSet object
# Note, the memory will only contain one dataSet object. By default, the last one will be the current dataSet object;
# Users can switch this (from the interface) to specify which data is load into memory (dataSet object)
# b) Include for certain analysis
# For chord and heatmap analysis, users can do multiple selection (include)
# All datasets are selected by default (1 for selected, 0 for unselected)
# note, dataSet need to have "name" property
# for now only most current will be selected
RegisterData <- function(dataSet){
  dataName <- dataSet$name;
  saveRDS(dataSet, file=dataName);
  mdata.all <<- lapply(mdata.all, function(x){ x <- 0;});
  mdata.all[[dataName]] <<- 1;
  print(paste("Sucessfully registered data:", dataName));
  return(1);
}

# only for switching single expression data results
SetCurrentData <- function(nm){
  #if(dataSet$name != nm){
    dataSet <- readRDS(nm);
  #}
  return(1);
}

# remove data object, the current dataSet will be the last one by default 
RemoveData <- function(dataName){
  if(!is.null(mdata.all[[dataName]])){
    mdata.all[[dataName]] <<- NULL;
  }
}

# users can select one or more data for analysis
# note, we use 1 to indicate this is selected
# and by default is all selected. 
SelectData <- function(){
  if(!exists('nm.vec')){
    current.msg <<-"No dataset is selected for analysis!";
    return(0);
  }
  
  all.nms <- names(mdata.all);
  for(nm in all.nms){
    if(nm %in% nm.vec){
      mdata.all[[nm]] <<- 1;
    }else{
      mdata.all[[nm]] <<- 0;
    }
  }
  
  if("merged" %in% nm.vec){
    merged.reduction <<- TRUE;
  }else{
    merged.reduction <<- FALSE;
  }
  
  rm('nm.vec', envir = .GlobalEnv);
  return(1);
}

GetAllDataNames <- function(){
  names(mdata.all);
}

# given a gene id, plot its expression profile as box plot
PlotSelectedGene<-function(filenm, gene.id, meta, selected_meta){
  library(Cairo)
  selected_meta <- strsplit(selected_meta, "; ")
  selected_meta = selected_meta[[1]]  
  require(lattice);
  if(selected_meta == "null"){
    metadf =dataSet$meta[,meta]
    trimmed.data.proc = dataSet$data.proc
  }else{
    metadf = dataSet$meta[,meta][dataSet$meta[,meta] %in% selected_meta]
    trimmed.data.proc = dataSet$data.proc[dataSet$meta[,meta] %in% selected_meta,]
  }
  
  Cairo(file = filenm, width=280, height=320, type="png", bg="white");
  myplot <- bwplot(trimmed.data.proc[gene.id,] ~ as.character(metadf), fill="#0000ff22", scales=list(x=list(rot=30)),
                   xlab="Class", ylab="Expression Pattern", main=gene.id);
  print(myplot); 
  dev.off();
}

UpdateSampleBasedOnLoading<-function(filenm, gene.id, omicstype){
  if(omicstype != "NA"){
  sel.nms <- names(mdata.all)[mdata.all==1];
    for(i in 1:length(sel.nms)){
      dat = readRDS(sel.nms[i])
      if(dat$type == omicstype){
        dataSet <- dat;
      }
    }
  }else{
    dataSet <- .get.rdt.set();
  }
  
  inx <- which(dataSet$enrich_ids == gene.id)
  print(head(dataSet$enrich_ids));
  id <- unname(dataSet$enrich_ids[inx])
  vec = as.vector(dataSet$data.proc[rownames(dataSet$data.proc) == gene.id,])
  colors<- ComputeColorGradient(as.numeric(vec), "black", F, F);
  sink(filenm);
  cat(toJSON(colors));
  sink();
}

DoDimensionReductionIntegrative <- function(reductionOpt){

    if(!exists("reduce.dimension")){ # public web on same user dir
        compiler::loadcmp("../../rscripts/OmicsAnalystR/R/_util_dimreduction.Rc");    
    }
    dr.res <- reduce.dimension(reductionOpt);

    return(dr.res)
}


PerformClustering <- function(init, nclust, type){
  dataSet <- .get.rdt.set();
  if(nrow(dataSet$pos.xyz)<21){
    return(1);
  }
  nclust = as.numeric(nclust)
  library(cluster)
    pos.xyz = dataSet$pos.xyz
  if(type == "density"){
    library(ADPclust)
    if(init == "true"){
      ans = ADPclust::adpclust(pos.xyz)
    }else{
      nclust = as.numeric(nclust)
      ans = ADPclust::adpclust(pos.xyz, nclust=nclust)
    }
    cluster = ans$clusters
    dataSet$clusterObject = ans
  }else if (type == "kmeans"){
    ans = kmeans(pos.xyz, nclust)
    cluster = ans$cluster
  }else{
    distMat = dist(pos.xyz)
    hc = hclust(distMat, "complete")
    cluster <- cutree(hc, k = nclust)
  }
  dataSet$newmeta$cluster=cluster
  
    loading.pos.xyz = dataSet$loading.pos.xyz
    nclust = as.numeric(nclust)
    ans = ADPclust::adpclust(pos.xyz, nclust=nclust)
    cluster = ans$clusters
    dataSet$loadingCluster=cluster
  
  .set.rdt.set(dataSet);
  return(1)
}

doScatterJson <- function(filenm){
    if(!exists("my.json.scatter")){ # public web on same user dir
        compiler::loadcmp("../../rscripts/OmicsAnalystR/R/_util_scatter_json.Rc");    
    }
    return(my.json.scatter(filenm));
}

PlotDataProfile<-function(dataName,type, boxplotName, pcaName){
  dataSet <- readRDS(dataName);
  if(type=="normalize"){
    qc.boxplot2(as.matrix(dataSet$data.proc), boxplotName);
    qc.pcaplot2(as.matrix(dataSet$data.proc), pcaName);
  }else{
    qc.boxplot2(as.matrix(dataSet$data.raw), boxplotName);
    qc.pcaplot2(as.matrix(dataSet$data.raw), pcaName);
  }
}


qc.boxplot2 <- function(dat, imgNm){
  require('lattice');
  require('Cairo');
  imgNm = paste(imgNm, "dpi", "72", ".png", sep="");
  subgene=10000;
  if (nrow(dat)>subgene) {
    set.seed(28051968);
    sg  = sample(nrow(dat), subgene)
    Mss = dat[sg,,drop=FALSE]
  } else {
    Mss = dat
  }
  
  subsmpl=100;
  if (ncol(Mss)>subsmpl) {
    set.seed(28051968);
    ss  = sample(ncol(Mss), subsmpl)
    Mss = Mss[,ss,drop=FALSE]
  } else {
    Mss = Mss
  }
  
  sample_id = rep(seq_len(ncol(Mss)), each = nrow(Mss));
  values  = as.numeric(Mss)
  formula = sample_id ~ values
  
  box = bwplot(formula, groups = sample_id, layout = c(1,1), as.table = TRUE,
               strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
               horizontal = TRUE,
               pch = "|",  col = "black", do.out = FALSE, box.ratio = 2,
               xlab = "", ylab = "Samples",
               fill = "#1c61b6AA",
               panel = panel.superpose,
               scales = list(x=list(relation="free"), y=list(axs="i")),
               ylim = c(ncol(Mss)+0.7,0.3),
               prepanel = function(x, y) {
                 list(xlim = quantile(x, probs = c(0.01, 0.95), na.rm=TRUE))
               },
               panel.groups = function(x, y, ...) {
                 panel.bwplot(x, y, ...)
               })
  
  Cairo(file=imgNm, width=460, height=420, type="png", bg="white");
  print(box);
  dev.off();
}

qc.pcaplot2 <- function(x, imgNm){
  imgNm = paste(imgNm, "dpi", "72", ".png", sep="");
  require('lattice');
  pca <- prcomp(t(na.omit(x)));
  names <- colnames(x);
  pca.res <- as.data.frame(pca$x);
  
  xlim <- GetExtendRange(pca.res$PC1);
  ylim <- GetExtendRange(pca.res$PC2);
  pcafig = xyplot(PC2~PC1, data=pca.res, pch=19, cex=1,aspect = "iso", xlim = xlim, ylim=ylim,
                  panel=function(x, y, ...) {
                    panel.xyplot(x, y, ...);
                    ltext(x=x, y=y, labels=names, pos=1, offset=1, cex=0.8, col="magenta")
                  })
  
  Cairo(file=imgNm, width=480, height=480, type="png", bg="white");
  print(pcafig);
  dev.off();
}

ScalingData <-function (nm,opt){
  #if(dataSet$name != nm){
    dataSet <- readRDS(nm);
  #}
  return(ScalingDataOmics(dataSet, opt))
}

ScalingDataOmics <-function (dataSet, norm.opt){
  return(1) 
}

PlotCluster <-function(opt, name, filenm, dpi, type){
  library(purrr)
  library(cluster)
  pos.xyz = dataSet$pos.xyz
  dpi = as.numeric(dpi)
  jsonnm = paste0(filenm, ".json");
  imgNm = paste(filenm, "dpi", dpi, ".", type, sep="");
  
  #if(opt == "kmeans"){
  mss <- (nrow(pos.xyz)-1)*sum(apply(pos.xyz,2,var))
  for (i in 2:16) mss[i] <- sum(kmeans(pos.xyz,centers=i)$withinss)
  
  
  # function to compute average silhouette for k clusters
  avg_sil <- function(k) {
    km.res <- kmeans(pos.xyz, centers = k, nstart = 25)
    ss <- silhouette(km.res$cluster, dist(pos.xyz))
    mean(ss[, 3])
  }
  
  # Compute and plot wss for k = 2 to k = 15
  k.values <- 2:10;
  # extract avg silhouette for 2-15 clusters
  avg_sil_values <- map_dbl(k.values, avg_sil);

  library(Cairo)
  Cairo(file=imgNm, width=800, height=600, type=type, bg="white", dpi=dpi, unit="px");
  
  g = plot(1:16, mss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares", main="Elbow plot")
  print(g)
  dev.off();
}

###
### Data utilities
###

ReadDataForMetaInfo<-function(dataName){
  dataSet <- readRDS(dataName)
  return(colnames(dataSet$meta));
}

GetFeatureNum <-function(dataName){
  #if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  #}
  return(nrow(dataSet$data.proc));
}

ComputeEncasing <- function(filenm, type, names.vec, level=0.95, omics="NA"){

  level <- as.numeric(level)
  names = strsplit(names.vec, "; ")[[1]]
  reductionSet <- .get.rdt.set();
  if(reductionOptGlobal %in% c("diablo") || omics != "NA"){
    if(grepl("pca_", omics, fixed=TRUE)){
        pca.scatter <- qs::qread("pca.scatter.qs");
        pos.xyz<-pca.scatter[[ omics ]]$score/1000
        print(head(pos.xyz))
    }else{
        omics.inx = 1;
        sel.nms <- names(mdata.all)[mdata.all==1];
        for(i in 1:length(sel.nms)){
        dataSet <- readRDS(sel.nms[i]);
            if(omics == dataSet$type){
                omics.inx = i;
            }
        }
        if(omics.inx == 1){
            pos.xyz = reductionSet$pos.xyz
        }else{
            pos.xyz = reductionSet$pos.xyz2
        }
    }

  }else{
  pos.xyz = reductionSet$pos.xyz
  }

  inx = rownames(pos.xyz) %in% names;
  coords = as.matrix(pos.xyz[inx,c(1:3)])
  mesh = list()
  if(type == "alpha"){
    library(alphashape3d)
    library(rgl)
    sh=ashape3d(coords, 1.0, pert = FALSE, eps = 1e-09);
    mesh[[1]] = as.mesh3d(sh, triangles=T);
  }else if(type == "ellipse"){
    library(rgl);
    pos=cov(coords, y = NULL, use = "everything");
    mesh[[1]] = ellipse3d(x=as.matrix(pos), level=level);
  }else{
    library(ks);
    res=kde(coords);
    r = plot(res, cont=level*100, display="rgl");
    sc = scene3d();
    mesh = sc$objects;
  }
  library(RJSONIO);
  sink(filenm);
  cat(toJSON(mesh));
  sink();
  return(filenm);
}

ClearFactorStrings<-function(cls.nm, query){
  # remove leading and trailing space
  query<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", query, perl=TRUE);
  
  # kill multiple white space
  query <- gsub(" +","_",query);
  # remove non alphabets and non numbers 
  query <- gsub("[^[:alnum:] ]", "_", query);
  
  # test all numbers (i.e. Time points)
  chars <- substr(query, 0, 1);
  num.inx<- chars >= '0' & chars <= '9';
  if(all(num.inx)){
    query = as.numeric(query);
    nquery <- paste(cls.nm, query, sep="_");
    query <- factor(nquery, levels=paste(cls.nm, sort(unique(query)), sep="_"));
  }else{
    query[num.inx] <- paste(cls.nm, query[num.inx], sep="_");
    query <- factor(query);
  }
  return (query);
}

GetCurrentJson <-function(type){
  return(jsonNms[[type]]);
}

.set.dataSet <- function(dataSetObj=NA){
  RegisterData(dataSetObj);
  return (1);
}

GetCurrentDatasets <-function(type){
  sel.nms <- names(mdata.all)[mdata.all==1];
  return(sel.nms);
}

SetGroupContrast <- function(dataName, grps, meta="NA"){
  #if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  #}

    if(meta == "NA"){
        meta <- 1;
    }

  if(length(levels(dataSet$meta[,meta]))>2){ 
    cls <- dataSet$meta[,meta]
    print("Updating group contrasts .....");
    grp.nms <- strsplit(grps, " vs. ")[[1]];
    sel.inx <- as.character(cls) %in% grp.nms;
    
    # regenerate factor to drop levels, force the levels order
    group <- factor(cls[sel.inx], levels=grp.nms);  
    data <- dataSet$data.proc[, sel.inx];
    dataSet$cls <- group;
    dataSet$data <- data;
    RegisterData(dataSet);  
  }

}


# here should first try to load the original data
# the data in the memory could be changed
GetGroupNames <- function(dataName, meta="NA"){
    #if(dataSet$name != dataName){
        dataSet <- readRDS(dataName);
    #}
    if(meta == "NA"){
        return(levels(factor(dataSet$meta[,1])));
    }else{
        return(levels(factor(dataSet$meta[,meta])));
    }

}

ReadMetaData <- function(metafilename){
  metadata <- .readDataTable(metafilename);
    metadata[is.na(metadata)] = "NA"
    if(class(metadata) == "try-error"){
      AddErrMsg("Failed to read in the metadata file! Please make sure that the metadata file is in the right format and does not have empty cells or contains NA.");
      return(0);
    }

    # need to add metadata sanity check
    # are sample names identical to data$orig
    # order samples in same way as in abundance table
    smpl.nms <- metadata[,1];
    smpl.var <- colnames(metadata)[-1];  
    sel.nms <- names(mdata.all)

   for(i in 1:length(sel.nms)){
    dataSet <- readRDS(sel.nms[i]);
    data.smpl.nms <- colnames(dataSet$data.proc)
    nm.hits <- data.smpl.nms %in% smpl.nms;
    if(!all(nm.hits)){
      AddErrMsg("Some sample names in your data are not in the metadata file!");
      mis.nms <- data.smpl.nms[!nm.hits];
      AddErrMsg(paste(mis.nms, collapse="; "));
      return(0);
    }

    # now remove extra meta if present, and order them
    nm.hits2 <- which(smpl.nms %in% data.smpl.nms);
    metadata1 <- metadata[nm.hits2,];
    metadata1 <- metadata1[,-1];

    if(!is.data.frame(metadata1)){
      metadata1 <- data.frame(metadata1, stringsAsFactors=T);
      colnames(metadata1) <- colnames(metadata)[2]
    }else{
      metadata1[] <- lapply( metadata1, factor)
    }
    rownames(metadata1) <- data.smpl.nms;
    dataSet$meta <- metadata1;
    RegisterData(dataSet);
  }
  return(1);
}

CheckDataType <- function(dataName, type){
  dataSet <- readRDS(dataName);
  isOk <- T;
  data <- dataSet$data.raw
  containsNeg <- "TRUE" %in% names(table(data < 0)) 
  msg.vec <- ""
  negativeBool <- F
  logBool <- F
  countBool <- T
  if(containsNeg){
    data <- data + abs(min(data))
    negativeBool <- T
  }
  
  # check percentage of values < 30 to guess if log transformed
  tbl <- table(data < 30)
  if("TRUE" %in% names(tbl)){
    num2 <- tbl[["TRUE"]]
    total <- dim(data)[1]*dim(data)[2]
    pct <- num2/total
    if(pct > 0.85){
      logBool <- T  
    }
  }else{
    logBool <- F  
  }
  
  # check if any decimal values
  if(any(data %% 1 != 0)){
    countBool <- F;
  }
  
  
  if(type == "false"){
    if(!countBool){
      msg.vec <- paste(msg.vec, "Decimal values detected;")
      isOk <-  F;
    }
    if(negativeBool){
      msg.vec <- paste(msg.vec, "Negative values detected;")
      isOk <-  F;
    }
  }

  dataSet$isValueNormalized <- type
  RegisterData(dataSet);

  if(!isOk){
    msg.vec <<- msg.vec
    return(0)
  }else{

    return(1)
  }
}

CheckNormalizedData <- function(dataName, omicsType){
  dataSet <- readRDS(dataName);
  dataSet$isValueNormalized <- "true";

  dataSet$type <- omicsType
  if(omicsType == "rna_b"){
    readableType <- "Transcriptomics";
  }else if (omicsType == "met_t"){
    readableType <- "Metabolomics";
  }else if (omicsType == "mic_m"){
    readableType <- "Microbiome";
  }else if (omicsType == "prot"){
    readableType <- "Proteomics";
  }else if (omicsType == "mirna"){
    readableType <- "miRNA";
  }else{
    readableType <-  omicsType;
  }
  dataSet$readableType <- readableType;

  RegisterData(dataSet);
}

SetParamsNormalizedData <- function(dataName){
    dataSet <- readRDS(dataName);

    int.mat <- dataSet$data.annotated;
    msg.vec <- "";
  
    if(sum(is.na(int.mat)) == 0){ # check if any missing values
      dataSet$data.proc <- int.mat;
      msg.vec <<- msg.vec;
      res <- 1;
    } else {
      msg.vec <<- c(msg.vec, "Missing values detected! Please upload normalized data with no missing values.");
      res <- 0;
    }

    RegisterData(dataSet);
    return(res)
}

RemoveMissingPercent <- function(dataName="", percent=0.5){

  dataSet <- readRDS(dataName);
  
  int.mat <- dataSet$data.annotated;
  good.inx1 <- apply(is.na(int.mat), 1, sum)/ncol(int.mat) < percent; # check less than 50% NA for each feature
  good.inx2 <- apply(R.utils:::isZero(int.mat), 1, sum)/ncol(int.mat) < percent; # check less than 50% 0 for each feature
  
  dataSet$data.annotated <- as.data.frame(int.mat[good.inx1, , drop=FALSE]);
  dataSet$data.annotated <- as.data.frame(int.mat[good.inx2, , drop=FALSE]);
  
  good.inx1 <- good.inx1[!is.na(good.inx1)]
  good.inx2 <- good.inx2[!is.na(good.inx2)]
  
  if(sum(!good.inx1) > 0 || sum(!good.inx2) > 0){
    msg.vec <<- paste(sum(!good.inx1) + sum(!good.inx2), " variables were removed for containing missing values over threshold", round(100*percent, 2), "percent.");
  }
  
  RegisterData(dataSet);
  return(1);
}

ImputeMissingVar <- function(dataName="", method="min"){

  # get parameters
  dataSet <- readRDS(dataName);
  int.mat <- dataSet$data.annotated;
  new.mat <- NULL;
  msg.vec <- "";
  
  if(sum(is.na(int.mat)) == 0){ # check if any missing values
    new.mat <- int.mat;
    msg.vec <<- msg.vec;
  } else {
    if(method=="exclude"){
      good.inx<-apply(is.na(int.mat), 1, sum)==0
      new.mat<-int.mat[,good.inx, drop=FALSE];
      msg.vec <<- c(msg.vec ,"Variables with missing values were excluded.");
      
    }else if(method=="min"){
      new.mat<- ReplaceMissingByLoD(int.mat);
      msg.vec <<- c(msg.vec, "Missing variables were replaced by LoDs (1/5 of the min positive value for each variable)");
    }else if(method=="colmin"){
      new.mat<-apply(int.mat, 1, function(x){
        if(sum(is.na(x))>0){
          x[is.na(x)]<-min(x,na.rm=T)/2;
        }
        x;
      });
      msg.vec <<- c(msg.vec,"Missing variables were replaced by 1/2 of min values for each feature column.");
    }else if (method=="mean"){
      new.mat<-apply(int.mat, 1, function(x){
        if(sum(is.na(x))>0){
          x[is.na(x)]<-mean(x,na.rm=T);
        }
        x;
      });
      msg.vec <<- c(msg.vec,"Missing variables were replaced with the mean value for each feature column.");
    }else if (method == "median"){
      new.mat<-apply(int.mat, 1, function(x){
        if(sum(is.na(x))>0){
          x[is.na(x)]<-median(x,na.rm=T);
        }
        x;
      });
      msg.vec <<- c(msg.vec,"Missing variables were replaced with the median for each feature column.");
    }else{
      if(method == "knn_var"){
        new.mat<-t(impute::impute.knn(int.mat)$data);
      }else if(method == "knn_smp"){
        new.mat<-impute::impute.knn(data.matrix(t(int.mat)))$data;
      }else{
        if(method == "bpca"){
          new.mat<-pcaMethods::pca(int.mat, nPcs =5, method="bpca", center=T)@completeObs;
        }else if(method == "ppca"){
          new.mat<-pcaMethods::pca(int.mat, nPcs =5, method="ppca", center=T)@completeObs;
        }else if(method == "svdImpute"){
          new.mat<-pcaMethods::pca(int.mat, nPcs =5, method="svdImpute", center=T)@completeObs;
        }
      }
      msg.vec <<- c(msg.vec, paste("Missing variables were imputated using", toupper(method)));
    }
  }
  
  if(!is.null(new.mat)){
    dataSet$data.missed <- as.data.frame(new.mat);
  }
  
  print("Impute Missing")
  RegisterData(dataSet);
  return(1)
}


FilteringData <- function(nm, countOpt="pct",count, var){
  dataSet <- readRDS(nm);
  return(FilteringDataOmics(dataSet,countOpt, count,  var))
}

FilteringDataOmics <- function(dataSet, countOpt="pct",count="2", var="15"){
  data = dataSet$data.missed
  msg <- ""
  
  count.thresh = as.numeric(count);
  var.thresh = as.numeric(var);
  if((count.thresh + var.thresh) >0){
  sum.counts <- apply(data, 1, sum, na.rm=TRUE);
  if(countOpt == "pct"){
  inx <- order(sum.counts)
  data <- data[inx,]
  rm.inx <- round(count.thresh/100 * nrow(data))
  data <- data[-c(1:rm.inx),];
  rmSum <- rm.inx
}else{
  rm.inx <- sum.counts < count.thresh;
  data <- data[!rm.inx,];
  rmSum <- sum(rm.inx)
}
  msg <- paste(msg, "Filtered ",rmSum, " genes with low counts.", collapse=" ");
  filter.val <- apply(data, 1, IQR, na.rm=T);
  nm <- "Interquantile Range";
  rk <- rank(-filter.val, ties.method='random');
  kp.pct <- (100 - var.thresh)/100;
  remain <- rk < nrow(data)*kp.pct;
  data <- data[remain,];
  msg <- paste("Filtered ", sum(!remain), " low variance genes based on IQR");
  
  if((sum(!remain) + rmSum)/nrow(data) > 0.8){
    msg <- paste("Over 80% of features are removed. Please readjust filtering thresholds." );
    msg.vec <<- msg;
    return(0);
  }
}

  dataSet$data.filtered = data
  dataSet$data.proc = data
if(exists("m2m",dataSet)){
  dataSet$data.filt.taxa <- list()
  mic.vec <- dataSet$taxa_table
  for(i in 1:length(colnames(mic.vec))){
    dataSet$data.filt.taxa[[i]]<- as.matrix(dataSet$data.filtered);
    rownames(dataSet$data.filt.taxa[[i]]) <- mic.vec[match(rownames(dataSet$data.filtered),rownames(mic.vec)),i];
    dataSet$data.filt.taxa[[i]] <- RemoveDuplicates(dataSet$data.filt.taxa[[i]], "sum", quiet=T); # remove duplicates
    dataSet$data.filt.taxa[[i]] <- as.data.frame(dataSet$data.filt.taxa[[i]]);
  }
names(dataSet$data.filt.taxa) <- colnames(dataSet$taxa_table);
  
}

dataSet$data.proc.taxa <- dataSet$data.filt.taxa;


  RegisterData(dataSet);
  msg.vec <<- msg

  return(1)
}


# Limit of detection (1/5 of min for each var)
.replace.by.lod <- function(x){
    #print(min(x[x>0], na.rm=T))
    lod <- min(x[x>0], na.rm=T)/5;
    orig.x <- x;
    x[x==0|is.na(x)] <- lod;
    x[is.infinite(x)] <- orig.x[is.infinite(x)] #replace Infinite values with original values
    return(x);
}

ReplaceMissingByLoD <- function(int.mat){
    int.mat <- as.matrix(int.mat);

    rowNms <- rownames(int.mat);
    colNms <- colnames(int.mat);
    int.mat <- t(apply(int.mat, 1, .replace.by.lod));
    rownames(int.mat) <- rowNms;
    colnames(int.mat) <- colNms;
    return (int.mat);
}


SetCustomSig <- function(dataName, ids){
    dataSet <- readRDS(dataName);
    lines <- strsplit(ids, "\r|\n|\r\n")[[1]];
    lines<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", lines, perl=TRUE);

    orig.nms.vec = names(dataSet$rawToEntrez)
    # need to first convert to correct id used in the graph
    hit.inx <- orig.nms.vec %in% lines;
    id.vec <- dataSet$rawToEntrez
    id.vec <- unname(id.vec[hit.inx])
    id.vec <- id.vec[!is.na(id.vec)]
    hit.inx2 <- rownames(dataSet$comp.res) %in% id.vec;
    dataSet$custom.sig.mat <- dataSet$comp.res[hit.inx2,]
    #dataSet$custom.sig.vec <- id.vec;
    print(nrow(dataSet$comp.res));
    sigNum <- nrow(dataSet$custom.sig.mat)

    if(sigNum > 0){
      RegisterData(dataSet);
    }
    matched.vec <- rownames(dataSet$custom.sig.mat)
    matched.vec  <- dataSet$rawToEntrez[unname(dataSet$rawToEntrez) %in% matched.vec]
    matched.vec <- names(matched.vec)
    return(matched.vec);
}


#'Record R Commands
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param cmd Commands 
#'@export
RecordRCommand <- function(cmd){
  cmdSet <- readSet(cmdSet, "cmdSet"); 
  cmdSet$cmdVec <- c(cmdSet$cmdVec, cmd);
  saveSet(cmdSet, "cmdSet");
  return(1);
}

SaveRCommands <- function(){
  cmdSet <- readSet(cmdSet, "cmdSet"); 
  cmds <- paste(cmdSet$cmdVec, collapse="\n");
  pid.info <- paste0("# PID of current job: ", Sys.getpid());
  cmds <- c(pid.info, cmds);
  write(cmds, file = "Rhistory.R", append = FALSE);
}

#'Export R Command History
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@export
GetRCommandHistory <- function(){
  cmdSet <- readSet(cmdSet, "cmdSet"); 
  if(length(cmdSet$cmdVec) == 0){
    return("No commands found");
  }
  return(cmdSet$cmdVec);
}

