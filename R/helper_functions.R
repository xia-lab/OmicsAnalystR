
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


