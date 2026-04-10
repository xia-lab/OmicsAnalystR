
##############################################
##############################################
########## Utilities for web-server ##########
##############################################
##############################################


GetNodeMat <- function(){
  #reductionSet <- .get.rdt.set();

  #if(is.null(reductionSet$imgSet$node_table)){
    df <- .readDataTable('node_table.csv')
    df[,-c(1:2)] <- lapply(df[,-c(1:2)], function(col) as.numeric(as.character(col)))
    #reductionSet$imgSet$node_table <<- df;
  #}
  return(as.matrix(df[,-c(1:2)]))  # ensure matrix of numerics
}

GetNodeRowNames <- function(){
  #reductionSet <- .get.rdt.set();

  #if(is.null(reductionSet$imgSet$node_table)){
  df <- .readDataTable('node_table.csv')
    #reductionSet$imgSet$node_table <<- df;

  #}
  df$Id;
}

GetNodeGeneSymbols <- function(){
  #reductionSet <- .get.rdt.set();

  #if(is.null(reductionSet$imgSet$node_table)){
  df <- .readDataTable('node_table.csv')
    #reductionSet$imgSet$node_table <<- df;

  #}
  df$Label;
}

GetNodeColNames <- function(){
  #reductionSet <- .get.rdt.set();

  #if(is.null(reductionSet$imgSet$node_table)){
  df <- .readDataTable('node_table.csv')
    #reductionSet$imgSet$node_table <<- df;

  #}
  return(colnames(df[,-c(1:2)]));

}

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


GetCovSigIDs<-function(dataName){
  dataSet <- readDataset(dataName);
  dataSet$analSet$cov$sig.mat$ids
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
  if (is.null(reductionSet$reductionOpt) || is.null(reductionSet[[reductionSet$reductionOpt]])) return("");
  dataSet <- readDataset(dataName);
  reductionSet[[reductionSet$reductionOpt]]$loading.file.nm;
}

GetLoadingMat<-function(dataName){
  reductionSet<-.get.rdt.set();
  if (is.null(reductionSet$reductionOpt) || is.null(reductionSet[[reductionSet$reductionOpt]])) return(matrix(0));
  dataSet <- readDataset(dataName);
  omicstype <- dataSet$type
  loading.pos.xyz <- reductionSet[[reductionSet$reductionOpt]]$loading.pos.xyz
  inx <- loading.pos.xyz$type %in% omicstype;
  drops <- c("ids","label", "type")
  return(CleanNumber(as.matrix(loading.pos.xyz[inx,!(names(loading.pos.xyz) %in% drops)])));
}

GetLoadingIds<-function(dataName){
  reductionSet<-.get.rdt.set();
  if (is.null(reductionSet$reductionOpt) || is.null(reductionSet[[reductionSet$reductionOpt]])) return(character(0));
  dataSet <- readDataset(dataName);
  omicstype <- dataSet$type;
  loading.pos.xyz <- reductionSet[[reductionSet$reductionOpt]]$loading.pos.xyz
  inx <- loading.pos.xyz$type %in% omicstype;
  loading.pos.xyz$ids[inx];
}

GetLoadingSymbols<-function(dataName){
  reductionSet<-.get.rdt.set();
  if (is.null(reductionSet$reductionOpt) || is.null(reductionSet[[reductionSet$reductionOpt]])) return(character(0));
  dataSet <- readDataset(dataName);
  omicstype <- dataSet$type
  loading.pos.xyz <- reductionSet[[reductionSet$reductionOpt]]$loading.pos.xyz
  inx <- loading.pos.xyz$type %in% omicstype;
  loading.pos.xyz$label[inx];
}

GetLoadingColNames<-function(dataName){
  reductionSet<-.get.rdt.set();
  if (is.null(reductionSet$reductionOpt) || is.null(reductionSet[[reductionSet$reductionOpt]])) return(character(0));
  dataSet <- readDataset(dataName);
  drops <- c("ids","label", "type")
  loading.pos.xyz <- reductionSet[[reductionSet$reductionOpt]]$loading.pos.xyz
  colnames(loading.pos.xyz[!(names(loading.pos.xyz) %in% drops)]);
}

GetVarianceArr<-function(dataName){
  reductionSet <- .get.rdt.set();
  if (is.null(reductionSet$reductionOpt) || is.null(reductionSet[[reductionSet$reductionOpt]])) return(c(0, 0, 0));
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
  }else if(type == "node"){
    performed <- file.exists("node_table.csv");
  }

  print(paste("checkPerformed=", type, "====",performed));

  return(performed)
}

# =============================================================================
# RSclient Subprocess Execution (Rserve fork isolation)
# =============================================================================
# Core infrastructure for isolating heavy R packages in forked Rserve children.
# RSclient comes with Rserve and is available on all deployments.
# pro_report_utils.R provides report rendering isolation.

#' Execute R function in forked Rserve child via RSclient
#' @param func Function to run
#' @param args List of arguments passed via do.call
#' @param timeout_sec Hard timeout in seconds
#' @return Result of do.call(func, args)
run_func_via_rsclient <- function(func, args = list(), timeout_sec = 60) {
  conn <- RSclient::RS.connect(host = "localhost", port = 6311)
  on.exit(try(RSclient::RS.close(conn), silent = TRUE))
  RSclient::RS.assign(conn, ".exec_wd", getwd())
  RSclient::RS.assign(conn, ".exec_func", func)
  RSclient::RS.assign(conn, ".exec_args", args)
  RSclient::RS.assign(conn, ".exec_timeout", timeout_sec)
  RSclient::RS.eval(conn, quote({
    setwd(.exec_wd)
    setTimeLimit(elapsed = .exec_timeout, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf))
    do.call(.exec_func, .exec_args)
  }))
}

#' Execute heavy package function in isolated RSclient fork
#' @param func_body Function(input_data) to run in child
#' @param input_data List serialized via qs to child
#' @param packages Packages to load in child before execution
#' @param timeout Timeout in seconds
#' @param output_type "qs" for complex objects
#' @return Result from child process
rsclient_isolated_exec <- function(func_body, input_data, packages = character(0),
                                   timeout = 180, output_type = "qs") {
  bridge_tmp <- file.path(tempdir(), "rsclient_bridge")
  if (!dir.exists(bridge_tmp)) dir.create(bridge_tmp, recursive = TRUE)
  uid <- paste0(sample(letters, 6), collapse = "")
  input_path <- file.path(bridge_tmp, paste0(uid, "_in.qs"))
  output_path <- file.path(bridge_tmp, paste0(uid, "_out.qs"))
  qs::qsave(input_data, input_path, preset = "fast")
  Sys.sleep(0.02)
  on.exit({ for (p in c(input_path, output_path)) if (file.exists(p)) unlink(p) }, add = TRUE)
  result <- run_func_via_rsclient(
    func = function(input_path, output_path, func_body, pkgs) {
      tryCatch({
        # Suppress quartz/X11 popups from rgl
        Sys.setenv(RGL_USE_NULL = TRUE)
        # HDF5Array .onLoad calls init_global_counter() which fails if counter
        # files already exist in tempdir(). Remove stale counters before loading.
        if ("MOFA2" %in% pkgs || "HDF5Array" %in% pkgs) {
          stale <- list.files(tempdir(), pattern = "^HDF5Array_", full.names = TRUE)
          if (length(stale) > 0) unlink(stale)
        }
        for (pkg in pkgs) suppressPackageStartupMessages(library(pkg, character.only = TRUE))
        input_data <- qs::qread(input_path)
        res <- func_body(input_data)
        qs::qsave(res, output_path, preset = "fast")
        Sys.sleep(0.02)
        list(success = TRUE)
      }, error = function(e) {
        list(success = FALSE, message = e$message)
      })
    },
    args = list(input_path = input_path, output_path = output_path,
                func_body = func_body, pkgs = packages),
    timeout_sec = timeout
  )
  if (isTRUE(result$success) && file.exists(output_path)) {
    return(qs::qread(output_path))
  }
  msg <- if (!is.null(result$message)) result$message else "RSclient subprocess failed"
  message("[rsclient_isolated_exec] ", msg)
  return(list(success = FALSE, message = msg))
}

# Null-coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

