##################################################
## R script for ExpressAnalyst
## Description: GSEA functions
## Author: G. Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#'Perform Gene Set Enrichment Analysis test on single or multiple gene expression matrices
#'@param dataSetObj file name of the data, .txt format
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
PerformGSEA<- function(dataName, file.nm, fun.type,omics.type="", input.type="loading",loading.comp=1, mode = "multi"){
  rdtSet <- .get.rdt.set();
  setres <- .loadEnrichLib(fun.type, data.org);
  current.geneset <- setres$current.geneset;
  require("fgsea");
  
  if(input.type == "loading"){
    loading.pos.xyz <- rdtSet$loading.pos.xyz.orig;
    loading.pos.xyz <- loading.pos.xyz[loading.pos.xyz$omicstype == omics.type,]
    rankedVec <- loading.pos.xyz[,loading.comp];
    names(rankedVec) <- rownames(loading.pos.xyz);
  }else{
    dataSet <- qs::qread(dataName);
    rankedVec <- dataSet$comp.res[,"coefficient"];
    names(rankedVec) <- dataSet$comp.res$ids;
  }
  
  if(mode == "simple"){
    fgseaRes <- fgsea(pathways = current.geneset, 
                      stats = rankedVec,
                      minSize=15,
                      maxSize=500,
                      nperm=10000);
  } else {
    fgseaRes <- fgsea(pathways = current.geneset, 
                      stats = rankedVec,
                      minSize=15,
                      maxSize=500);
  }
  
  fgseaRes <- fgseaRes[!duplicated(fgseaRes$pathway),];
  rownames(fgseaRes) <- make.names(fgseaRes$pathway, unique=TRUE)
  fgseaRes <- fgseaRes[,c("size","ES", "pval", "pathway", "padj")]
  fgseaRes <- fgseaRes[order(fgseaRes$pval),]
  if(nrow(fgseaRes[which(fgseaRes$pval < 0.05),])<20 ){
    if(nrow(fgseaRes)>20){
      fgseaRes <- fgseaRes[c(1:20),]
    }
  }else{
    fgseaRes <- fgseaRes[which(fgseaRes$pval < 0.05),]
  } 
  
  current.mset <- current.geneset[fgseaRes$pathway]
  current.mset <- current.mset[!duplicated(names(current.mset))]
  
  ora.vec <- names(rankedVec)
  ora.nms <- ora.vec #doEntrez2SymbolMapping(ora.vec)
  hits.query <- lapply(current.mset, 
                       function(x) {
                         unique(ora.nms[ora.vec%in%unlist(x)]);
                       });
  
  set.num <- unlist(lapply(current.mset, function(x){length(unique(x))}), use.names=TRUE);
  names(hits.query) <- names(current.mset);
  hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=TRUE);
  
  qs::qsave(hits.query, "hits_query.qs");
  fgseaRes$hits <- hit.num[which(fgseaRes$pathway %in% names(hit.num))] 
  fgseaRes$total <- set.num[which(fgseaRes$pathway %in% names(set.num))]
  
  fgseaRes <- fgseaRes[which(fgseaRes$hits>1),];
  fgseaRes <- fgseaRes[which(fgseaRes$hits<500),];
  fgseaRes <- fgseaRes[which(fgseaRes$total<2000),];
  
  fgseaRes=fgseaRes[order(fgseaRes$pval),];
  if(nrow(fgseaRes[which(fgseaRes$pval < 0.05),])<20 ){
    if(nrow(fgseaRes)>20){
      fgseaRes <- fgseaRes[c(1:20),];
    }
  }else{
    fgseaRes <- fgseaRes[which(fgseaRes$pval < 0.05),];
  } 
  
  fgseaRes <- data.frame(fgseaRes, stringsAsFactors=FALSE);
  
  #get gene symbols
  current.msg <<- "Functional analysis was completed";
  
  # write json
  fun.anot <- hits.query; 
  fun.pval <- fgseaRes[,3]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  fun.padj <- fgseaRes[,5]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  es.num <- fgseaRes[,2]; if(length(es.num) ==1) { es.num <- matrix(es.num) };
  fun.ids <- as.vector(setres$current.setids[names(fun.anot)]); 
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  
  json.res <- list(
    fun.link = setres$current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    pathname = fgseaRes[,"pathway"],
    es.num = es.num,
    hit.num = fgseaRes[,"hits"],
    total = fgseaRes[,"total"]
  );
  

  json.res$org <- data.org
  json.res$analType <- anal.type
  json.res$naviString <- "GSEA";
  
  json.mat <- RJSONIO::toJSON(json.res);
  json.nm <- paste(file.nm, ".json", sep="");
  partialToBeSaved <- c(partialToBeSaved, c(json.nm, "current_geneset.qs"))
  sink(json.nm)
  cat(json.mat);
  sink();
  
  # write csv
  
  ftype <- fun.type
  if(fun.type %in% c("bp", "mf", "cc")){
    ftype <- paste0("go_", fun.type);
  }
  
  csvDf <- data.frame(Name=fgseaRes$pathway, Total=fgseaRes$total, Hits=fgseaRes$hits, EnrichmentScore=fgseaRes$ES, Pval=fgseaRes$pval, Padj=fgseaRes$padj);
  fast.write(csvDf, file=paste0(file.nm, ".csv"));
  
  return(1);
}