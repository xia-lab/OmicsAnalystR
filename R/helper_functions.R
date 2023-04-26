
##############################################
##############################################
########## Utilities for web-server ##########
##############################################
##############################################

GetCovSigFileName <-function(dataName){
  dataSet <- qs::qread(dataName);
  dataSet$analSet$cov$sig.nm;
}

GetCovSigMat<-function(dataName){
  dataSet <- qs::qread(dataName);
  drops <- c("ids","label")
  return(CleanNumber(as.matrix(dataSet$analSet$cov$sig.mat[, !(names(dataSet$analSet$cov$sig.mat) %in% drops)])));
}

GetCovSigRowNames<-function(dataName){
  dataSet <- qs::qread(dataName);
  rownames(dataSet$analSet$cov$sig.mat);
}

GetCovSigSymbols<-function(dataName){
  dataSet <- qs::qread(dataName);
  dataSet$analSet$cov$sig.mat$label
}

GetCovSigColNames<-function(dataName){
  dataSet <- qs::qread(dataName);
  drops <- c("ids","label");
  colnames(dataSet$analSet$cov$sig.mat[,!(names(dataSet$analSet$cov$sig.mat) %in% drops)]);
}

GetPrimaryType <- function(analysis.var){
    rdtSet <- .get.rdt.set();
    primary.type <- unname(rdtSet$dataSet$meta.types[analysis.var]);
    return(primary.type);
}

GetCovDENums <- function(dataName){
    deNum <- nrow(dataSet$analSet$cov$sig.mat);
    nonDeNum <- nrow(dataSet$comp.res) - deNum;
    return(c(deNum, nonDeNum));
}



# users can manually update sample names
UpdateFeatureName <-function(dataName, old.nm, new.nm){
  dataSet <- qs::qread(dataName);
  inx <- names(dataSet$enrich_ids) == old.nm;
  names(dataSet$enrich_ids)[inx] <- new.nm;
  return(RegisterData(dataSet));
}


GetLoadingFileName <-function(dataName){
  reductionSet<-.get.rdt.set();
  dataSet <- qs::qread(dataName);
  reductionSet$loading.file.nm;
}

GetLoadingMat<-function(dataName){
  reductionSet<-.get.rdt.set();
  dataSet <- qs::qread(dataName);
  omicstype <- dataSet$type
  inx <- reductionSet$loading.pos.xyz$type %in% omicstype;
  drops <- c("ids","label", "type")
  return(CleanNumber(as.matrix(reductionSet$loading.pos.xyz[inx,!(names(reductionSet$loading.pos.xyz) %in% drops)])));
}

GetLoadingIds<-function(dataName){
  reductionSet<-.get.rdt.set();
  dataSet <- qs::qread(dataName);
  omicstype <- dataSet$type;
  print(head(reductionSet$loading.pos.xyz));
  inx <- reductionSet$loading.pos.xyz$type %in% omicstype;
  reductionSet$loading.pos.xyz$ids[inx];
}

GetLoadingSymbols<-function(dataName){
  reductionSet<-.get.rdt.set();
  dataSet <- qs::qread(dataName);
  omicstype <- dataSet$type
  inx <- reductionSet$loading.pos.xyz$type %in% omicstype;
  reductionSet$loading.pos.xyz$label[inx];
}

GetLoadingColNames<-function(dataName){
  reductionSet<-.get.rdt.set();
  dataSet <- qs::qread(dataName);
  drops <- c("ids","label", "type")
  colnames(reductionSet$loading.pos.xyz[!(names(reductionSet$loading.pos.xyz) %in% drops)]);
}

GetVarianceArr<-function(omicsType){
  reductionSet <- .get.rdt.set();
  df <- reductionSet$var.exp;
  varArr <- df[,omicsType];
  varArr <- signif(varArr,4)*100;
  print(varArr);
  return(varArr);
}

