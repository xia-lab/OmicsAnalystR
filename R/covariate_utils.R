##################################################
## R scripts for OmicsAnalyst
## Description: Related to linear modeling
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

GetCovUpMat<-function(){
  rdtSet <- .get.rdt.set();
  lod <- rdtSet$analSet$cov$p.log;
  pval.no <- rdtSet$analSet$cov$p.value.no;
  red.inx<- which(rdtSet$analSet$cov$inx.imp);
  if(sum(red.inx) > 0){
    return(as.matrix(cbind(pval.no[red.inx], lod[red.inx])));
  }else{
    return(as.matrix(cbind(-1, -1)));
  }
}

GetCovUpIDs <- function(){
  rdtSet <- .get.rdt.set();
  red.inx<- which(rdtSet$analSet$cov$inx.imp);
  if(sum(red.inx) > 0){
    return(names(rdtSet$analSet$cov$p.log)[red.inx]);
  }else{
    return("NA");
  }
}

GetCovDnMat<-function(){
  rdtSet <- .get.rdt.set();
  lod <- rdtSet$analSet$cov$p.log;
  pval.no <- rdtSet$analSet$cov$p.value.no;
  blue.inx <- which(!rdtSet$analSet$cov$inx.imp);
  if(sum(blue.inx) > 0){
    return(as.matrix(cbind(pval.no[blue.inx], lod[blue.inx])));
  }else{
    return(as.matrix(cbind(-1, -1)));
  }
}

GetCovDnIDs <- function(){
  rdtSet <- .get.rdt.set();
  blue.inx<- which(!rdtSet$analSet$cov$inx.imp);
  if(sum(blue.inx) > 0){
    return(names(rdtSet$analSet$cov$p.log)[blue.inx]);
  }else{
    return("NA");
  }
}


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
  return(CleanNumber(as.matrix(dataSet$analSet$cov$sig.mat)));
}

GetCovSigRowNames<-function(dataName){
  dataSet <- qs::qread(dataName);
  rownames(dataSet$analSet$cov$sig.mat);
}

GetCovSigColNames<-function(dataName){
  dataSet <- qs::qread(dataName);
  colnames(dataSet$analSet$cov$sig.mat);
}

GetPrimaryType <- function(analysis.var){
    rdtSet <- .get.rdt.set();
    primary.type <- unname(rdtSet$dataSet$meta.types[analysis.var]);
    return(primary.type);
}


#' CovariateScatter.Anal
#' @param imgName image name
#' @param format image format
#' @param analysis.var variable of analysis
#' @param ref reference group
#' @param block block name
#' @param thresh threshold
#' @param contrast.cls contrast group
#' @export
CovariateScatter.Anal <- function(dataName, 
                                  imgName="NA", 
                                  format="png", 
                                  analysis.var, 
                                  ref = NULL, 
                                  block = "NA", 
                                  thresh=0.05,
                                  contrast.cls = "anova"){
  save.image("scatter.RData");
  dataSet <- qs::qread(dataName);
  rdtSet <- .get.rdt.set();
  
  # load libraries
  library(limma)
  library(dplyr)
  
  # get inputs
  if(!exists('adj.vec')){
    adj.bool = F;
    vars <- analysis.var;
  }else{
    if(length(adj.vec) > 0){
      adj.bool = T;
      vars <- c(analysis.var, adj.vec)
    }else{
      adj.bool = F;
      vars <- analysis.var;
    }
  }
  
  covariates <- rdtSet$dataSet$meta.info
  var.types <- rdtSet$dataSet[["meta.types"]]
  feature_table <- dataSet$data.proc
  
  # process inputs
  thresh <- as.numeric(thresh)
  ref <- make.names(ref)
  analysis.type <- unname(rdtSet$dataSet$meta.types[analysis.var]);
  
  # process metadata table (covariates)
  for(i in c(1:length(var.types))){ # ensure all columns are the right type
    if(var.types[i] == "disc"){
      covariates[,i] <- covariates[,i] %>% make.names() %>% factor()
    } else {
      covariates[,i] <- covariates[,i] %>% as.character() %>% as.numeric()
    }
  }
  #subset to samples contained in dataset
  covariates <- covariates[match(colnames(feature_table), rownames(covariates)),]
  if (block != "NA"){    
    if(dataSet$meta.types[block] == "cont"){
      AddErrMsg("Blocking factor can not be continuous data type.")
      return(c(-1,-1));
    }
    # recent update: remove check for unbalanced design. Limma can handle.
  }
  
  sig.num <- 0;
  
  if(analysis.type == "disc"){
    # build design and contrast matrix
    covariates[, analysis.var] <- covariates[, analysis.var] %>% make.names() %>% factor();
    grp.nms <- levels(covariates[, analysis.var]);
    design <- model.matrix(formula(paste0("~ 0", paste0(" + ", vars, collapse = ""))), data = covariates);
    colnames(design)[1:length(grp.nms)] <- grp.nms;
    myargs <- list();
    
    # perform specified contrast
    if(contrast.cls == "anova"){
      cntr.cls <- grp.nms[grp.nms != ref];
      myargs <- as.list(paste(cntr.cls, "-", ref, sep = ""));
    } else {
      myargs <- as.list(paste(contrast.cls, "-", ref, sep = ""));
    }
    myargs[["levels"]] <- design;
    contrast.matrix <- do.call(makeContrasts, myargs);
    feature_table <- feature_table[,which(colnames(feature_table) %in% rownames(design)) ]
    # handle blocking factor
    if (block == "NA") {
      fit <- lmFit(feature_table, design)
    } else {
      block.vec <- covariates[,block];
      corfit <- duplicateCorrelation(feature_table, design, block = block.vec)
      fit <- lmFit(feature_table, design, block = block.vec, correlation = corfit$consensus)
    }
    
    fit <- contrasts.fit(fit, contrast.matrix);
    fit <- eBayes(fit);
    rest <- topTable(fit, number = Inf);
    
    ### get results with no adjustment
    design <- model.matrix(formula(paste0("~ 0", paste0(" + ", analysis.var, collapse = ""))), data = covariates);
    colnames(design)[1:length(grp.nms)] <- grp.nms;
    myargs[["levels"]] <- design;
    contrast.matrix <- do.call(makeContrasts, myargs);
    fit <- lmFit(feature_table, design)
    fit <- contrasts.fit(fit, contrast.matrix);
    fit <- eBayes(fit);
    res.noadj <- topTable(fit, number = Inf);
    
  } else { 
    covariates[, analysis.var] <- covariates[, analysis.var] %>% as.numeric();
    types <- unname(rdtSet$dataSet$meta.types[vars])
    if(sum(types == "cont") == length(vars)){ #in case of single, continuous variable, must use different intercept or limma will give unreasonable results
      design <- model.matrix(formula(paste0("~", paste0(" + ", vars, collapse = ""))), data = covariates);
    } else {
      design <- model.matrix(formula(paste0("~ 0", paste0(" + ", vars, collapse = ""))), data = covariates);
    }
    
    feature_table <- feature_table[,which(colnames(feature_table) %in% rownames(design)) ]
    # recent update: enable blocking factor for continuous primary metadata
    if (block == "NA") {
      fit <- lmFit(feature_table, design)
    } else {
      block.vec <- covariates[,block];
      corfit <- duplicateCorrelation(feature_table, design, block = block.vec)
      fit <- lmFit(feature_table, design, block = block.vec, correlation = corfit$consensus)
    }
    
    fit <- eBayes(fit);
    rest <- topTable(fit, number = Inf, coef = analysis.var);
    colnames(rest)[1] <- analysis.var;
    
    ### get results with no adjustment
    design <- model.matrix(formula(paste0("~", analysis.var)), data = covariates);
    
    fit <- eBayes(lmFit(feature_table, design));
    res.noadj <- topTable(fit, number = Inf);
  }
  
  
  # make visualization
  adj.mat <- rest[, c("P.Value", "adj.P.Val")]
  noadj.mat <- res.noadj[, c("P.Value", "adj.P.Val")]
  
  colnames(adj.mat) <- c("pval.adj", "fdr.adj")
  colnames(noadj.mat) <- c("pval.no", "fdr.no")
  
  both.mat <- merge(adj.mat, noadj.mat, by = "row.names")
  both.mat$pval.adj <- -log10(both.mat$pval.adj)
  both.mat$fdr.adj <- -log10(both.mat$fdr.adj)
  both.mat$pval.no <- -log10(both.mat$pval.no)
  both.mat$fdr.no <- -log10(both.mat$fdr.no)
  
  # make plot
  if( "F" %in% colnames(rest)){
    fstat <- rest[, "F"];
  }else{
    fstat <- rest[, "t"];
  }    
  
  p.value <- rest[,"P.Value"];
  
  names(fstat) <- names(p.value) <- colnames(dataSet$data.proc);
  fdr.p <- rest[,"adj.P.Val"];
  inx.imp <- p.value <= thresh;
  sig.num <- sum(inx.imp);
  
  if(sig.num > 0){ 
    sig.p <- p.value[inx.imp];
    sig.mat <- rest[inx.imp,];
    sig.mat <- sapply(sig.mat, function(x) signif(x, 5));
    rownames(sig.mat) <- rownames(rest)[inx.imp]
    # order the result simultaneously
    ord.inx <- order(sig.p, decreasing = FALSE);
    sig.mat <- sig.mat[ord.inx,,drop=F];
  }
  
  AddMsg(paste(c("A total of", sum(inx.imp), "significant features were found."), collapse=" "));
  rownames(both.mat) = both.mat[,1]
  both.mat <- both.mat[rownames(rest),]
  if(sig.num> 0){
    res <- 1;
    fileName <- "covariate_result.csv"
    fast.write.csv(sig.mat,file=fileName);
    cov<-list (
      sig.num = sig.num,
      sig.nm = fileName,
      raw.thresh = thresh,
      thresh = -log10(thresh), # only used for plot threshold line
      p.value = p.value,
      p.value.no = both.mat$pval.no,
      p.log = -log10(p.value),
      inx.imp = inx.imp,
      sig.mat = sig.mat
    );
  }else{
    res <- 0;
    cov<-list (
      sig.num = sig.num,
      raw.thresh = thresh,
      thresh = -log10(thresh), # only used for plot threshold line
      p.value = p.value,
      p.value.no = both.mat$pval.no,
      p.log = -log10(p.value),
      inx.imp = inx.imp
    );
  }
  
  # for detail table
  dataSet$analSet$cov <- cov; 
  # for plotting adjp vs p
  dataSet$analSet$cov.mat <- both.mat; 
  

  jsonNm <- gsub(paste0(".", format), ".json", imgName);
  jsonObj <- RJSONIO::toJSON(both.mat);
  sink(jsonNm);
  cat(jsonObj);
  sink();

  #reformat for comp.res
  colnames(both.mat)[1] <- c("ids");
  dataSet$comp.res <- both.mat;
  dataSet$sig.mat <- sig.mat

  return(RegisterData(dataSet));
}


PlotCovariateMap <- function(dataName, theme="default", imgName="NA", format="png", dpi=72){
  dataSet <- qs::qread(dataName);
  both.mat <- dataSet$cov.mat
  both.mat <- both.mat[order(-both.mat[,"pval.adj"]),]
  logp_val <- dataSet$cov$thresh
  load_ggplot();
  library(ggrepel);
  topFeature <- 5;
  if(nrow(both.mat) < topFeature){
    topFeature <- nrow(both.mat);
  }
  if(theme == "default"){
    p <- ggplot(both.mat, mapping = aes(x = pval.no, y = pval.adj, label = Row.names)) +
      geom_rect(mapping = aes(xmin = logp_val, xmax = Inf, 
                              ymin = logp_val, ymax = Inf),
                fill = "#6699CC") +
      geom_rect(mapping = aes(xmin = -Inf, xmax = logp_val, 
                              ymin = -Inf, ymax = logp_val),
                fill = "grey") +
      geom_rect(mapping = aes(xmin = logp_val, xmax = Inf, 
                              ymin = -Inf, ymax = logp_val),
                fill = "#E2808A") +
      geom_rect(mapping = aes(xmin = -Inf, xmax = logp_val, 
                              ymin = logp_val, ymax = Inf),
                fill = "#94C973") +
      guides(size="none") +
      #annotate("text", x = 0.8, y = 0, label = "Never significant", size = 3) +
      #annotate("text", x = 2, y = 0, label = "Significant without adjustment", size = 3) +
      #annotate("text", x = 0.4, y = 1.5, label = "Significant with adjustment", size = 3) +
      #annotate("text", x = 2.25, y = 1.5, label = "Always significant", size = 3) +
      geom_point(aes(size=pval.adj), alpha=0.5) +
      geom_abline(slope=1, intercept = 0, linetype="dashed", color = "red", size = 1) +
      xlab("-log10(P-value): no covariate adjustment") +
      ylab("-log10(P-value): adjusted") +
      geom_text_repel(data = both.mat[c(1:topFeature),], 
                  aes(x=pval.no,y=pval.adj,label=Row.names)) +
      theme_bw();
  }else{
    p <- ggplot(both.mat, mapping = aes(x = pval.no, y = pval.adj, label = Row.names)) +
      guides(size="none") +
      geom_point(aes(size=pval.adj), alpha=0.5) +
      geom_abline(slope=1, intercept = 0, linetype="dashed", color = "red", size = 1) +
      geom_vline(xintercept = logp_val) +
      geom_hline(yintercept = logp_val) +
      xlab("-log10(P-value): no covariate adjustment") +
      ylab("-log10(P-value): adjusted") +
      geom_text_repel(data = both.mat[c(1:topFeature),], 
                  aes(x=pval.no,y=pval.adj,label=Row.names))
  }
  
  dataSet$covAdj <- imgName;

  width <- 8;
  height <- 8.18;
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=width, height=height, type=format, bg="white");
  print(p)
  dev.off()
  
  return(RegisterData(dataSet));
}


# general message only print when running local
AddMsg <- function(msg){
  if(!exists("msg.vec")){
    msg.vec <<- "";
  }
  msg.vec <<- c(msg.vec, msg);
  if(!.on.public.web){
    print(msg);
  }
}


#'Plot compound summary for multi-linear regression tool
#'@param cmpdNm Input the name of the compound to plot
#'@param format Input the format of the image to create
#'@param dpi Input the dpi of the image to create
#'@param width Input the width of the image to create
#'@param meta Input the metadata to visualize
#'@param version version
#'@author Jessica Ewald\email{jessica.ewald@mcgill.ca}
#'McGill University, Canada
#'License: GPL-3 License
#'@export
#'
PlotMultiFacCmpdSummary <- function(dataName, cmpdNm, meta, version, format="png", dpi=72, width=NA){
  dataSet <- qs::qread(dataName);
  rdtSet <- .get.rdt.set();
  
  if(.on.public.web){
    load_ggplot()
  }
  
  if(is.na(width)){
    w <- 7.5;
  }else{
    w <- width;
  }
  
  meta.info <- rdtSet$dataSet$meta.info
  sel.cls <- meta.info[which(rownames(meta.info) %in% colnames(dataSet$data.proc)),meta]
  cls.type <- unname(rdtSet$dataSet$meta.types[meta])
  xlab = meta;
  h <- 6;
  imgName <- rdtSet$dataSet$url.var.nms[cmpdNm];
  imgName <- paste(imgName, "_", meta, "_", version, "_summary_dpi", dpi, ".", format, sep="");
  
  inx <- which(rownames(dataSet$data.proc) == cmpdNm)
  df.norm <- data.frame(value=as.vector(t(dataSet$data.proc)[, inx]), name = sel.cls)
  col <- unique(GetColorSchema(sel.cls));
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  if(cls.type == "disc"){
    p <- ggplot2::ggplot(df.norm, aes(x=name, y=value, fill=name)) + geom_boxplot(outlier.shape = NA, outlier.colour=NA) + theme_bw() + geom_jitter(size=1) 
    p <- p + scale_fill_manual(values=col) + theme(axis.text.x = element_text(angle=90, hjust=1))
    p <- p + ggtitle(cmpdNm) + theme(plot.title = element_text(size = 11, hjust=0.5, face = "bold")) + ylab("Abundance") + xlab(meta)
    p <- p + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) # remove gridlines
    p <- p + theme(plot.margin = margin(t=0.15, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10)) 
  }else{
    p <- ggplot2::ggplot(df.norm, aes(x=name, y=value)) 
    p <- p + geom_point(size=2) + theme_bw()  + geom_smooth(method=lm,se=T)     
    p <- p + theme(axis.text.x = element_text(angle=90, hjust=1)) + guides(size="none")
    p <- p + ggtitle(cmpdNm) + theme(plot.title = element_text(size = 11, hjust=0.5, face = "bold")) + ylab("Abundance") + xlab(meta)
    p <- p + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) # remove gridlines
    p <- p + theme(plot.margin = margin(t=0.15, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10)) 
  }
  print(p)
  dev.off()

  if(.on.public.web){
    return(imgName);
  }else{
    return(.set.rdt.set(rdtSet));
  }
}
