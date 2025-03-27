
##############################################
##############################################
########## Utilities for web-server ##########
##############################################
##############################################

GetCovSigFileName <-function(dataName){
  dataSet <- readDataset(dataName);
  dataSet$analSet$cov$sig.nm;
}

GetCovSigMat<-function(dataName){
  dataSet <- readDataset(dataName);
  drops <- c("ids","label")
  return(CleanNumber(as.matrix(dataSet$analSet$cov$sig.mat[, !(names(dataSet$analSet$cov$sig.mat) %in% drops)])));
}

GetCovSigRowNames<-function(dataName){
  dataSet <- readDataset(dataName);
  rownames(dataSet$analSet$cov$sig.mat);
}

GetCovSigSymbols<-function(dataName){
  dataSet <- readDataset(dataName);
  dataSet$analSet$cov$sig.mat$label
}

GetCovSigColNames<-function(dataName){
  dataSet <- readDataset(dataName);
  drops <- c("ids","label");
  colnames(dataSet$analSet$cov$sig.mat[,!(names(dataSet$analSet$cov$sig.mat) %in% drops)]);
}

GetCovDENums <- function(dataName){
    deNum <- nrow(dataSet$analSet$cov$sig.mat);
    nonDeNum <- nrow(dataSet$comp.res) - deNum;
    return(c(deNum, nonDeNum));
}



# users can manually update sample names
UpdateFeatureName <-function(dataName, old.nm, new.nm){
  dataSet <- readDataset(dataName);
  inx <- names(dataSet$enrich_ids) == old.nm;
  names(dataSet$enrich_ids)[inx] <- new.nm;
  return(RegisterData(dataSet));
}

SetDimMethod <- function(method){
  reductionSet<-.get.rdt.set();
  reductionSet$reductionOpt <- method;
  .set.rdt.set(reductionSet);
}


GetLoadingFileName <-function(dataName){
  reductionSet<-.get.rdt.set();
  dataSet <- readDataset(dataName);
  reductionSet[[reductionSet$reductionOpt]]$loading.file.nm;
}

GetLoadingMat<-function(dataName){
  reductionSet<-.get.rdt.set();
  dataSet <- readDataset(dataName);
  omicstype <- dataSet$type
  loading.pos.xyz <- reductionSet[[reductionSet$reductionOpt]]$loading.pos.xyz
  inx <- loading.pos.xyz$type %in% omicstype;
  drops <- c("ids","label", "type")
  return(CleanNumber(as.matrix(loading.pos.xyz[inx,!(names(loading.pos.xyz) %in% drops)])));
}

GetLoadingIds<-function(dataName){
  reductionSet<-.get.rdt.set();
  dataSet <- readDataset(dataName);
  omicstype <- dataSet$type;
  loading.pos.xyz <- reductionSet[[reductionSet$reductionOpt]]$loading.pos.xyz
  inx <- loading.pos.xyz$type %in% omicstype;
  loading.pos.xyz$ids[inx];
}

GetLoadingSymbols<-function(dataName){
  reductionSet<-.get.rdt.set();
  dataSet <- readDataset(dataName);
  omicstype <- dataSet$type
  loading.pos.xyz <- reductionSet[[reductionSet$reductionOpt]]$loading.pos.xyz
  inx <- loading.pos.xyz$type %in% omicstype;
  loading.pos.xyz$label[inx];
}

GetLoadingColNames<-function(dataName){
  reductionSet<-.get.rdt.set();
  dataSet <- readDataset(dataName);
  drops <- c("ids","label", "type")
  loading.pos.xyz <- reductionSet[[reductionSet$reductionOpt]]$loading.pos.xyz
  colnames(loading.pos.xyz[!(names(loading.pos.xyz) %in% drops)]);
}

GetVarianceArr<-function(dataName){
  reductionSet <- .get.rdt.set();
  dataSet <- readDataset(dataName); 
  df <- reductionSet[[reductionSet$reductionOpt]]$var.exp;
  varArr <- df[,dataSet$type];
  varArr <- signif(varArr,4)*100;
  return(varArr);
}

GetOmicsDataDims <- function(dataName){
  dataSet <- readDataset(dataName);
  dm <- dim(dataSet$data.proc);
  naNum <- sum(is.na(dataSet$data.proc));
  return(c(dm, naNum));
} 

GetMultiSummary <- function(){
  sel.nms <- names(mdata.all);
  featureNumAnn <- "";
  featureNumFilter <- "";
  dat.nms <- "";
  for(i in 1:length(sel.nms)){
    dataSet = readDataset(sel.nms[i]);
    datAnn <- qs::qread(dataSet$data.annotated.path);
    datProc <- dataSet$data.proc;
    if(i == 1){
      cls.lbls <- dataSet$meta[,1]
      featureNumAnn <- nrow(datAnn);
      featureNumFilter <- nrow(datProc);
      sampleNum <- ncol(datProc);
      dat.nms <- sel.nms[i];
    }else{
      featureNumAnn <- c(featureNumAnn, nrow(datAnn));
      featureNumFilter <- c(featureNumFilter, nrow(datProc));
      dat.nms <- c(dat.nms, sel.nms[i]);
    }
  }
  featureNumAnn <- paste(featureNumAnn, collapse="; ");
  featureNumFilter <- paste(featureNumFilter, collapse="; ");
  dat.nms <- paste(dat.nms, collapse="; ");
  cls.lbls <- unique(cls.lbls);
  cls.lbls <- paste(cls.lbls, collapse="; ");
  res <- list(
    sampleNum = sampleNum,
    featureNumAnn = featureNumAnn,
    dat.nms = dat.nms,
    cls.lbls = cls.lbls,   
    featureNumFilter = featureNumFilter
  )
  infoSet <- readSet(infoSet, "infoSet");
  infoSet$paramSet$summaryUploadedData <- res;
  saveSet(infoSet);
  return(unlist(res))
}

SetReductionOpt <- function(opt){
  reductionSet<-.get.rdt.set();
  reductionSet$reductionOpt <- opt;
}


CheckDetailsTablePerformed <-function(type, dataName){
  reductionSet<-.get.rdt.set();
  dataSet <- readDataset(dataName);
  omicstype <- dataSet$type;

  performed <- T;
  if(grepl("loading", type)){
    reductOpt <- gsub("loading_","", type);
    performed <- !is.null(reductionSet[[reductOpt]]$loading.pos.xyz);
  }else if(startsWith(type, "OmicsData #")){
    performed <- !is.null(dataSet$analSet$cov$sig.mat);
  }

  print(paste("checkPerformed=", type, "====",performed));

  return(performed)
}