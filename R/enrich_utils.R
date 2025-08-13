
PerformHeatmapEnrichment <- function(file.nm, fun.type, IDs, type="heatmap0"){
  if(IDs=="NA"){

  }else{
    gene.vec <- unlist(strsplit(IDs, "; "));
  }
  sym.vec <- doEntrez2SymbolMapping(gene.vec);
  names(gene.vec) <- sym.vec;
  res <- PerformEnrichAnalysis(file.nm, fun.type, gene.vec, type);
  return(res);
}


PerformScatterOrNetEnrichment <- function(file.nm, fun.type, IDs, type){
  # prepare query
  print("PerformScatterOrNetEnrichment")
  id_type <<- "entrez";
  ora.vec <- unlist(strsplit(IDs, "; "));
  names(ora.vec) <- ora.vec;
  res <- PerformEnrichAnalysis(file.nm, fun.type, ora.vec, type);
  return(res);   
}


PerformKeggEnrichment <- function(file.nm, fun.type, ids){
  if(ids=="Not_applicable"){
    ora.vec <- keggp.allfeatures
    names(ora.vec) <- ora.vec;    
  }else{
    ora.vec <- unlist(strsplit(ids, "; "));
    names(ora.vec) <- ora.vec;
    ora.vecu <<- ora.vec
  }
  noGene <- F
  noMet <- F
  if(length(rownames(dataSet$ko.seed))== 0 && length(rownames(dataSet$gene.seed))== 0){
    noGene <- T
  }else if(length(rownames(dataSet$met.seed))== 0){
    noMet <- T
  }
  
  if(fun.type == "keggm" && noMet){
    return(2)
  }else if(fun.type != "keggm" && fun.type != "integ" && noGene){
    return(3)
  }
  
  if(fun.type == "integ"){
    res1 = PerformEnrichAnalysisKegg(file.nm, "kegg", ora.vec)
    res2 = PerformEnrichAnalysisKegg(file.nm, "keggm", ora.vec)
    res3 = PerformEnrichAnalysisKegg(file.nm, "integ", ora.vec)
  }else{
    res1 = PerformEnrichAnalysisKegg(file.nm, fun.type, ora.vec)
    res.mat <- as.data.frame(res1)
  }
  if(fun.type == "integ" || fun.type == "ginteg"){
    inx = which(rownames(res1) %in% rownames(res2))
    subres1 = as.data.frame(res1[inx,])
    inx = which(rownames(res2) %in% rownames(res1))
    subres2 = as.data.frame(res2[inx,])
    inx = which(rownames(res3) %in% rownames(subres2))
    subres3 = as.data.frame(res3[inx,])
    ord = order(rownames(subres1));
    subres1 = subres1[ord,]
    ord = order(rownames(subres2));
    subres2 = subres2[ord,]
    ord = order(rownames(subres3));
    subres3 = subres3[ord,]
    integ=data.frame(hitsG = subres1$Hits,hitsM = subres2$Hits,hitsTotal = subres3$Hits, P.ValueG=subres1$P.Value, P.ValueM=subres2$P.Value, P.ValueMerge=subres3$P.Value, P.ValueJoint=subres2$P.Value)
    
    rownames(integ) = rownames(subres1)
    for(i in 1:nrow(integ)){
      if(integ$P.ValueG[i] != 1 && integ$P.ValueM[i] != 1){
        #integ$P.ValueJoint[i] = metap::sumlog(c(integ$P.ValueG[i], integ$P.ValueM[i]))$p
        integ$P.ValueJoint[i] = metap::sumz(p=c(integ$P.ValueG[i], integ$P.ValueM[i]), weight=c(stouffer_gene_percent, stouffer_compound_percent))$p
      }else{
        integ$P.ValueJoint[i]=1
      }
    }
    
    inxM = integ[,"hitsG"] == 0
    inxG = integ[,"hitsM"] == 0
    inxJ = integ[,"hitsM"] != 0 & integ[,"hitsM"] != 0 
    integP = integ[,"P.ValueJoint"]
    integP[inxM] = integ[,"P.ValueM"][inxM]
    integP[inxG] = integ[,"P.ValueG"][inxG]
    integ$integP = integP
    res.mat <- integ
    hits.query<<-hits.query
    SaveIntegEnr(file.nm,res.mat);
  }else{
    SaveSingleOmicsEnr(file.nm,res.mat);
    if(fun.type == "kegg"){
      kg.query <<-hits.query
    }else if(fun.type == "keggm"){
      km.query <<-hits.query
    }
  }
  
  containsKo = length(rownames(dataSet$ko.seed))>0
  containsGene = length(rownames(dataSet$gene.seed))>0
  containsMet = length(rownames(dataSet$met.seed))>0
  if(containsMet && (containsGene || containsKo)){
    met.vec <- dataSet$met.seed
    #if(isKo){
    #  gene.vec <- dataSet$ko.seed
    #}else{
      gene.vec <- dataSet$gene.seed
    #}
    res= pathwayMetGenePair(met.vec, gene.vec)
    nohits.vec = vector()
    
    for(i in 1:length(hits.query)){
      query.vec = hits.query[[i]]
      table.name = names(hits.query)[i]
      res1 = res[which(res[,2] %in% query.vec),]
      res2 = res[which(res[,5] %in% query.vec),]
      combinedres = rbind(res1, res2)
      combinedres = unique(combinedres)
      if(length(rownames(combinedres)) >0){
        combinedres$pathway = table.name
        if(i==1){
          alltables = combinedres
        }else{
          alltables = rbind(alltables, combinedres)
        }
      }
    }
    fast.write.csv(alltables, paste0("pairs_", file.nm, ".csv"), row.names = F)
  }
  return(1)
}


# note: hit.query, resTable must synchronize
# ora.vec should contains entrez ids, named by their gene symbols
PerformEnrichAnalysis <- function(file.nm, fun.type, ora.vec,type,ifNet=F){
  # prepare lib
  print(c("type",type))
  setres <- .loadEnrichLib(fun.type, data.org)
   current.geneset <- setres$current.geneset;
  current.setids <- setres$current.setids;
  current.setlink <- setres$current.setlink;
  current.universe <- unique(unlist(current.geneset));
  ###remove the NA pathway
  rm <- which(is.na(names(current.geneset)))
  if(length(rm)>0){
  current.geneset <- current.geneset[-rm]
  current.setids <- current.setids[-rm]

  } 
  # prepare query
  ora.nms <- names(ora.vec);
  
  # prepare for the result table
  set.size<-length(current.geneset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.geneset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "P.Value", "FDR");
  
  # need to cut to the universe covered by the pathways, not all genes 
  if(tolower(fun.type) %in% c("chea", "jaspar", "encode", "mir")){ 
    current.universe = toupper(current.universe)
  }
  hits.inx <- ora.vec %in% current.universe;
  ora.vec <- ora.vec[hits.inx];
  ora.nms <- ora.nms[hits.inx];
 
  if(length(ora.vec)==0){
  return(0)
  }

  q.size<-length(ora.vec);
  
  # get the matched query for each pathway
  hits.query <- lapply(current.geneset, 
                       function(x) {
                         ora.nms[ora.vec%in%unlist(x)];
                       }
  );
  hits.query <- lapply(hits.query, function(x){unique(x)});
  names(hits.query) <- names(current.geneset);
  hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=FALSE);
   lib.num<-unlist(lapply(current.geneset, function(x){length(x)}), use.names=FALSE);
  # total unique gene number
  uniq.count <- length(current.universe);
  
  # unique gene count in each pathway
  set.size <- unlist(lapply(current.geneset, length));

  res.mat[,1]<-set.size;
  res.mat[,2]<-q.size*(set.size/uniq.count);
  res.mat[,3]<-paste0(hit.num,"/", lib.num);
  
  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
  res.mat[,4]<- raw.pvals;
  res.mat[,5] <- p.adjust(raw.pvals, "fdr");
 
  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];
  
  if(nrow(res.mat)> 1){
    # order by p value
    ord.inx<-order(res.mat[,4]);
    res.mat <- res.mat[ord.inx,]
  #  res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];
    
    imp.inx <- res.mat[,4] <= 0.05;
    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      if(sum(imp.inx) > 120){
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }
  
  #get gene symbols
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
  current.msg <<- "Functional enrichment analysis was completed";
  if(!grepl("pathway",file.nm )){
  # write json
  require(RJSONIO);
  fun.anot = hits.query; 
  fun.padj = resTable$FDR; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  fun.pval = resTable$P.Value; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  hit.num = resTable$Hits; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
  fun.ids <- as.vector(current.setids[current.setids %in% names(fun.anot)]);
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
   print( fun.pval)
  json.res <- list(
    fun.link = current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    hit.num = hit.num
  );
  json.mat <- toJSON(json.res);
  json.nm <- paste(file.nm, ".json", sep="");
  sink(json.nm)
  cat(json.mat);
  sink();
 }
  # write csv
  csv.nm <- paste(file.nm, ".csv", sep="");
  fast.write.csv(resTable, file=csv.nm, row.names=F);
  #record table for report
  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$enrTables[[type]] <- list()
  infoSet$imgSet$enrTables[[type]]$table <- resTable;
  infoSet$imgSet$enrTables[[type]]$library <- fun.type
  infoSet$imgSet$enrTables[[type]]$algo <- "Overrepresentation Analysis"
  infoSet$imgSet$enrTables[[type]]$hits.query <- hits.query
  infoSet$imgSet$enrTables[[type]]$current.geneset <- current.geneset
  saveSet(infoSet);
  
  return(1);
}


.loadEnrichLib <- function(fun.type, data.org){
  folderNm <- data.org;
  my.path <- paste(lib.path, folderNm, "/", fun.type, ".rds", sep="");
  print(my.path)
  my.lib <- readRDS(my.path);
 
  if(substr(fun.type, 0, 2)=="go"){
    if(is.null(names(my.lib))){ # some go lib does not give names
      names(my.lib) <- c("link", "term", "sets");
    }
  }
  res <- list();
  if(fun.type=="integ"){
   current.geneset <- my.lib
  qs::qsave(current.geneset, "current_geneset.qs");
 res$current.setlink <- "";
    res$current.setids <- names(current.geneset);
    res$current.geneset <- current.geneset;

  return(res);
  }
  current.geneset <- my.lib$sets;
  
  #remove empty pathways
  keep.inx <- lapply(current.geneset,length)>0
  current.geneset <- current.geneset[keep.inx]
  my.lib$term <- my.lib$term[keep.inx]
  set.ids<- names(current.geneset); 
  names(set.ids) <- names(current.geneset) <- my.lib$term;
  
  qs::qsave(current.geneset, "current_geneset.qs");

  
  if(fun.type == "keggm"){
    res$current.setlink <- "";
    res$current.setids <- names(my.lib);
    res$current.geneset <- my.lib
  }else{
    res$current.setlink <- my.lib$link;
    res$current.setids <- set.ids;
    res$current.geneset <- current.geneset;
  }
  return(res);
}


PerformEnrichAnalysisKegg <- function(file.nm, fun.type, ora.vec){
  toremove <- c("Metabolic pathways",
                "Biosynthesis of secondary metabolites",
                "Microbial metabolism in diverse environments",
                "Biosynthesis of antibiotics",
                "Carbon metabolism",
                "2-Oxocarboxylic acid metabolism",
                "Fatty acid metabolism",
                "Biosynthesis of amino acids",
                "Degradation of aromatic compounds");
  # prepare lib
  setres <- .loadEnrichLib(fun.type, data.org)
  current.geneset <- setres$current.geneset;
  current.setids <- setres$current.setids;
  current.setlink <- setres$current.setlink;
  current.universe <- unique(unlist(current.geneset));
  
  # prepare query
  ora.nms <- names(ora.vec);
  lens <- lapply(current.geneset, 
                 function(x) {
                   length(unique(x))
                 }
  );
  inx = lens > 2
  
  current.geneset = current.geneset[inx]
  current.geneset = current.geneset[which(!names(current.geneset) %in% toremove)]
  
  # prepare for the result table
  set.size<-length(current.geneset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.geneset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "Pval", "FDR");
  
  # need to cut to the universe covered by the pathways, not all genes 
  
  hits.inx <- ora.vec %in% current.universe;
  
  #calculate weight for stouffer
  if(fun.type == "kegg"){
    stouffer_gene_percent <<- length(hits.inx)/length(current.universe)
  }else if(fun.type == "keggm"){
    stouffer_compound_percent <<- length(hits.inx)/length(current.universe)
  }
  
  ora.vec <- ora.vec[hits.inx];
  ora.nms <- ora.nms[hits.inx];
  
  q.size<-length(ora.vec);
  
  # get the matched query for each pathway
  hits.query <- lapply(current.geneset, 
                       function(x) {
                         ora.nms[ora.vec%in%unlist(x)];
                       }
  );
  hits.query <- lapply(hits.query, function(x){unique(x)})
  
  names(hits.query) <- names(current.geneset);
  hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=FALSE);
  
  # total unique gene number
  uniq.count <- length(current.universe);
  
  # unique gene count in each pathway
  set.size <- unlist(lapply(current.geneset, length));
  
  res.mat[,1]<-set.size;
  res.mat[,2]<-q.size*(set.size/uniq.count);
  res.mat[,3]<-hit.num;
  
  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
  res.mat[,4]<- raw.pvals;
  res.mat[,5] <- p.adjust(raw.pvals, "fdr");
  all.res.mat <<- res.mat
  hits.query <<- hits.query
  if(fun.type == "integ"){
    integ.query<- list()
    integ.query<- hits.query
  }
  # now, clean up result, synchronize with hit.query
  
  return(all.res.mat);
}
 

SaveSingleOmicsEnr <- function(file.nm,res.mat){
  inx = res.mat[,3]>0
  res.mat <- res.mat[inx,];
  hits.query <- hits.query[inx];
  if(nrow(res.mat)> 1){
    # order by p value
    ord.inx<-order(res.mat[,4]);
    res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];
    
    imp.inx <- res.mat[,4] < 0.05
    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      if(sum(imp.inx) > 120){
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }
  
  #get gene symbols
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
  current.msg <<- "Functional enrichment analysis was completed";
  
  # write json
  require(RJSONIO);
  hi= hits.query
  fun.anot = hi; 
  fun.padj = resTable$FDR; if(length(fun.padj) ==1) { fun.pval <- matrix(fun.padj)};
  fun.pval = resTable$Pval; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval)};
  hit.num = paste0(resTable$Hits,"/",resTable$Total); if(length(hit.num) ==1) { hit.num <- matrix(hit.num)};
  fun.ids <- as.vector(current.setids[names(fun.anot)]);
  if(length(fun.ids) ==1) {fun.ids <- matrix(fun.ids)};
  
  json.res <- list(
    fun.link = current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    hit.num = hit.num
  );
  json.mat <- toJSON(json.res);
  json.nm <- paste(file.nm, ".json", sep="");
  
  sink(json.nm)
  cat(json.mat);
  sink();
  hitss = lapply(hits.query, function(x){paste(x, collapse=" ")})
  hitss = unlist(hitss)
  resTable$Features = hitss
  csv.nm <- paste(file.nm, ".csv", sep="");
  fast.write.csv(resTable, file=csv.nm, row.names=F);
}

InitEnrichmentNetwork <- function(file.nm, fun.type,type){
   print(file.nm)
    sel.inx <- mdata.all==1; 
  sel.nms <- names(mdata.all)[sel.inx];
  dataSetList <- lapply(sel.nms, readDataset);
  
  if(fun.type=="kegg"){
   kp = unlist(lapply(dataSetList,function(x) x[["type"]] %in% c("rna_b","prot")))

  }else if(fun.type=="keggm"){
   kp = lapply(dataSetList,function(x) x[["type"]] %in% c("met_t","met_u"))

  }else if(fun.type=="tam_hmdd"|fun.type=="tam_func"|fun.type=="mirfamily"){
   kp = lapply(dataSetList,function(x) x[["type"]] %in% c("mirna"))
 
  }else if(fun.type=="integ"){
  kp = lapply(dataSetList,function(x) x[["type"]] %in% c("met_t","met_u","rna_b","prot"))
 }

  if(type=="limma"){
    idList <- lapply(dataSetList[kp],function(x){
       return(x[["sig.mat"]]$ids)
     })
   }
  sig.mat <- do.call(rbind,lapply(dataSetList[kp],function(x) x[["sig.mat"]]))
   
  
  id_type <<- "entrez";
  ora.vec <- unique(unlist(idList));
  names(ora.vec) <- ora.vec;
  res <- PerformEnrichAnalysis(file.nm, fun.type, ora.vec, type,T);
  if(res==0){
  return(0)
  }
 infoSet <- readSet(infoSet, "infoSet");
   infoSet$imgSet$enrTables[[type]]$sig.mat <- sig.mat
   infoSet$imgSet$enrTables[[type]]$selDataNm <-  sel.nms[kp];
     saveSet(infoSet);
  .prepareEnrichNet(netNm=file.nm,type, "mixed")
  return(1)
}

.prepareEnrichNet<-function( netNm, type, overlapType){
    if(!exists("my.enrich.net")){ 
        compiler::loadcmp("../../rscripts/OmicsAnalystR/R/utils_enrichnet.Rc");      
    }
    return(my.enrich.net(netNm, type, overlapType));
}


GetEnrichList <- function(dataNm, type,fileNm){
  print(dataNm)
  dataSet<- readDataset(dataNm);
  if(type=="limma"){
   if(nrow(dataSet[["sig.mat"]])==0){
    return(0)

    }
  print(head(dataSet[["sig.mat"]]))


   if(dataSet$idType=="name"|dataSet$idType=="symbol"){
     #all_str <- dataSet[["sig.mat"]]$label
      if("logFC" %in%  colnames(dataSet[["sig.mat"]])){
 df <- data.frame("#id" = dataSet[["sig.mat"]]$label, logFC = dataSet[["sig.mat"]]$logFC,check.names = FALSE)
   }else{
df <- data.frame("#id" = dataSet[["sig.mat"]]$label, F = dataSet[["sig.mat"]]$F,check.names = FALSE)
 
}
     }else{
     
      if("logFC" %in% colnames(dataSet[["sig.mat"]])){
     df <- data.frame("#id" = dataSet[["sig.mat"]]$ids, logFC = dataSet[["sig.mat"]]$logFC,check.names = FALSE)
   }else{
df <- data.frame("#id" = dataSet[["sig.mat"]]$ids, F = dataSet[["sig.mat"]]$F,check.names = FALSE)
 
}
  } 

   }else if(type=="var"){
    if(is.null(dataSet[["varPar"]])){
    return(0)
    }
if(dataSet$idType=="name"|dataSet$idType=="symbol"){
 
     df <- data.frame("#id" = dataSet[["varPar"]]$label, percentage = round(dataSet[["varPar"]]$varPart.df[[2]],4),check.names = FALSE)

     }else{
     df <- data.frame("#id" = dataSet[["varPar"]]$ids, percentage = round(dataSet[["varPar"]]$varPart.df[[2]],4),check.names = FALSE)

  } 
df <- df[1:min(nrow(df),dataSet$varPar$topNum),]
}

write.table(df, file = fileNm, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  #writeLines(all_str, fileNm)
  return(1)
}

GetNetEnrichlist <- function(netNm, mode, dataName="",enrNm="",queryFrom="dr",query=""){
   print(queryFrom)

   if(mode=="fileNms"){
     netFile <- RJSONIO::fromJSON(netNm)
     netEnrTypes <<- netFile$omicstype
     netEnridtypes <<- netFile$idTypes

    return(netFile$fileNms)
   }else if(mode=="types"){
      type_map <- c(
  "met_t" = "met",
  "rna_b" = "gene",
  "prot"  = "protein",
  "mic_m" = "mic",
  "mirna" = "mir"
)

      return(unname(type_map[netEnrTypes]))
   }else if(mode=="idtypes"){
     return(netEnridtypes)
  }else if(mode=="data"){
    if(queryFrom=="corr"){
netFile <- RJSONIO::fromJSON(netNm)
   idx = which(netFile$fileNms==dataName)
   nodes <- netFile$nodes[unlist(lapply(netFile$nodes, function(x) x$molType==netFile$nodeTypes[idx]))] 
  df = data.frame("#id" = unlist(lapply(nodes,function(x) x$featureId)), btw =  unlist(lapply(nodes,function(x) x$size)),check.names = FALSE)

  }else if(queryFrom=="dr"){
loading = qs::qread("loading.qs")[[dataName]]
    print(head(loading))
df = data.frame("#id" = loading$ids, loading=round(loading$Loading,4),check.names = FALSE)

}
     
  write.table(df, file = enrNm, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
return(1)
 }
   

}


GetDREnrichList <- function(query,tupNum){
  rdtSet <- .get.rdt.set();
 tupNum <- as.numeric(tupNum)
  loading <- rdtSet[[rdtSet$reductionOpt]]$loading.pos.xyz
  loading <- split(loading[,c(1:3,6)],loading$filenm)
 
  if(query=="1"){ 
    loading <- lapply(loading,function(x){
       z= x[order(abs(x$Loading1),decreasing = T),c("Loading1","ids")]
       names(z)[1] = "Loading"
       z=z[1:(tupNum*nrow(z)/100),]
       return(z)
    })
  }else if(query=="2"){
    loading <- lapply(loading,function(x){
      z= x[order(abs(x$Loading1),decreasing = T),c("Loading2","ids")]
      names(z)[1] = "Loading"
      z=z[1:(tupNum*nrow(z)/100),]
      return(z)
    })
    
  }else if(query=="3"){
    loading <- lapply(loading,function(x){
      z= x[order(abs(x$Loading1),decreasing = T),c("Loading3","ids")]
      names(z)[1] = "Loading"
      z=z[1:(tupNum*nrow(z)/100),]
      return(z)
    })
  }

 qs::qsave(loading,"loading.qs") 

}

PrepDashboardList <- function(mode,jsonnm){
 
  rdtSet <- .get.rdt.set();
 # sel.inx <- mdata.all==1; 
  sel.nms <- names(mdata.all)#[sel.inx];
  dataSetList <- lapply(sel.nms, readDataset);
 
listnms = lapply(dataSetList,function(x) return(x[["readableType"]]))

featls <- lapply(dataSetList, function(x) x[["sig.mat"]])
names(featls) <- listnms
jsonls <- list(featls=featls)
jsonls$filenms <- lapply(dataSetList,function(x) return(x[["name"]]))
jsonls$types <- lapply(dataSetList,function(x) return(x[["type"]]))
jsonls$idTypes <- lapply(dataSetList,function(x) return(x[["idType"]]))
if("logFC" %in% colnames(dataSetList[1][["sig.mat"]])){
jsonls$mode<-"comp"
}else{
jsonls$mode<-"anova"
}

 jsonfile <- RJSONIO::toJSON(jsonls);
  sink(jsonnm)
  cat(jsonfile);
  sink();
return(1)
}

PrepDashboardListDR <- function(mode,jsonnm,featNumPerComp,components){
  library(dplyr)
  rdtSet <- .get.rdt.set();
  sel.nms <- names(mdata.all)#[sel.inx];
  print("PrepDashboardListDR")
  dataSetList <- lapply(sel.nms, readDataset); 
  listnms = lapply(dataSetList,function(x) return(x[["readableType"]]))
  loadings = rdtSet[[mode]][["loading.pos.xyz"]] 
  components <- as.numeric(unlist(strsplit(components,split = ",")))
  loadings <- loadings[,c(components,6:8)] 
  types <- unlist(lapply(dataSetList,function(x) x[["readableType"]]))
  names(types) <- unlist(lapply(dataSetList,function(x) x[["type"]]))
 
  tops <- reshape2::melt(loadings,id.vars=c("ids","type","label"  )) %>%
    group_by(type, variable) %>%
    arrange(desc(abs(value)), .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    dplyr::select(label, ids, value, type, variable)  %>% ungroup()
  loadings <- loadings[loadings$ids %in% tops$ids,]
  loadings <- loadings %>%
    rowwise() %>%
    mutate(
       loading = c(Loading1, Loading2, Loading3)[which.max(abs(c(Loading1, Loading2, Loading3)))]
    ) %>%
    ungroup()

 
  featls <- split(loadings,list(loadings$type))
  indx <- match(names(featls) ,names(types) )
  featls <- featls[indx]
    names(featls) <- types[ names(featls) ]
    rdtSet[[mode]]$featls <- featls
    .set.rdt.set(rdtSet);
  jsonls <- list(featls=featls)
  jsonls$filenms <- lapply(dataSetList,function(x) return(x[["name"]]))
  jsonls$types <- lapply(dataSetList,function(x) return(x[["type"]]))
  jsonls$idTypes <- lapply(dataSetList,function(x) return(x[["idType"]]))
 
  jsonls$mode<-"dr"
  
  jsonfile <- RJSONIO::toJSON(jsonls);
  sink(jsonnm)
  cat(jsonfile);
  sink();
  return(1)
}

#############################
##############GETTERS########
#############################

GetFeatureRownames <- function(selectedBar,mode){
 
  sel.inx <- mdata.all==1; 
  sel.nms <- names(mdata.all)[sel.inx];
selectedBar <- as.numeric(selectedBar)+1
dataSet <- readDataset(sel.nms[selectedBar]) 
if(mode=="default"){ 

 return(rownames(dataSet[["sig.mat"]]))
}else if(mode %in% c("mcia","mofa","diablo")){
  rdtSet <- .get.rdt.set();
  mat <- rdtSet[[mode]]$featls[[dataSet$readableType]] 
 return(rownames( mat))
}
}

GetFeatureIds <- function(selectedBar,mode){

 sel.inx <- mdata.all==1; 
  sel.nms <- names(mdata.all)[sel.inx];
selectedBar <- as.numeric(selectedBar)+1

dataSet <- readDataset(sel.nms[selectedBar])
  
if(mode=="default"){
 mat <- dataSet[["sig.mat"]]
 
}else if(mode %in% c("mcia","mofa","diablo")){
  rdtSet <- .get.rdt.set();
 mat <-  rdtSet[[mode]]$featls[[dataSet$readableType]] 
}

 print(head(mat$ids))
 return(mat$ids)
}

GetFeatureLabels <- function(selectedBar,mode){
 sel.inx <- mdata.all==1; 
  sel.nms <- names(mdata.all)[sel.inx];
selectedBar <- as.numeric(selectedBar)+1

dataSet <- readDataset(sel.nms[selectedBar])
  
if(mode=="default"){

 
 mat <- dataSet[["sig.mat"]] 

}else if(mode %in% c("mcia","mofa","diablo")){
  rdtSet <- .get.rdt.set(); 
 mat <-  rdtSet[[mode]]$featls[[dataSet$readableType]] 

}
 return(mat$label)
}


GetFeatureColnms <- function(selectedBar,mode){
 
  sel.inx <- mdata.all==1; 
  sel.nms <- names(mdata.all)[sel.inx];
selectedBar <- as.numeric(selectedBar)+1

dataSet <- readDataset(sel.nms[selectedBar])
 
if(mode=="default"){
 mat <- dataSet[["sig.mat"]] 
 return( setdiff(colnames(mat),c("ids","label")))
}else if(mode %in% c("mcia","mofa","diablo")){
  rdtSet <- .get.rdt.set();
 mat <-  rdtSet[[mode]]$featls[[dataSet$readableType]] 
 return( setdiff(colnames(mat),c("ids","label","type","loading")))
}

}

GetFeatureSigMat <- function(selectedBar,mode) {
 sel.inx <- mdata.all==1; 
  sel.nms <- names(mdata.all)[sel.inx];
selectedBar <- as.numeric(selectedBar)+1

dataSet <- readDataset(sel.nms[selectedBar])
 
if(mode=="default"){
 mat <- dataSet[["sig.mat"]] 
 mat <- as.matrix(mat[,setdiff(colnames(mat),c("ids","label"))])

}else if(mode %in% c("mcia","mofa","diablo")){
  rdtSet <- .get.rdt.set();
 mat <-  rdtSet[[mode]]$featls[[dataSet$readableType]] 
 mat <- as.matrix(mat[,setdiff(colnames(mat),c("ids","label","type","loading"))])
}
  return( mat)
}

PerformDashEnrichment <- function(file.nm, fun.type, type, dtInx){
  # prepare query
  print(mdata.all)
print(c(fun.type, type, dtInx))
 sel.inx <- mdata.all==1; 
  sel.nms <- names(mdata.all)[sel.inx];
dtInx <- as.numeric(unlist(strsplit(dtInx,split=",")))+1
 print(dtInx)
 id_type <<- "entrez";
if(type=="default"){
 dataSetList <- lapply(sel.nms[dtInx], readDataset);
 ora.vec<- unique(unlist(unique(lapply(dataSetList,function(x) x[["sig.mat"]][["ids"]]))))
}else if(type %in% c("mcia","mofa","diablo")){
  rdtSet <- .get.rdt.set();

  featls <- rdtSet[[type]]$featls
  ora.vec<- unique(unlist(unique(lapply(featls,function(x) x[["ids"]]))))
}

  names(ora.vec) <- ora.vec;
  print(head(ora.vec))
  res <- PerformEnrichAnalysis(file.nm, fun.type, ora.vec, type);
  return(res);   
}


GetPathRownames <- function(type){

infoSet <- readSet(infoSet, "infoSet");
 return(rownames(infoSet[["imgSet"]][["enrTables"]][[type]][["table"]]))

}



GetPathIds <- function(type){

infoSet <- readSet(infoSet, "infoSet"); 

pathIds <- infoSet$imgSet$enrTables[[type]]$table$Pathway
return(pathIds)

}



GetPathHits <- function(type){

infoSet <- readSet(infoSet, "infoSet"); 

hits <- infoSet$imgSet$enrTables[[type]]$table$Hits
return(hits)

}


GetPathPval <- function(type){

infoSet <- readSet(infoSet, "infoSet"); 

pvals <- infoSet$imgSet$enrTables[[type]]$table$P.Value
return(round(as.numeric(pvals),4))

}


GetPathAdjPval <- function(type){

infoSet <- readSet(infoSet, "infoSet"); 

fdrs <- infoSet$imgSet$enrTables[[type]]$table$FDR
return(round(as.numeric(fdrs),4))

}



GetHitInfo <- function(pathnm,mode,dtInx){
  if(mode=="default"){
    sel.inx <- mdata.all==1; 
    sel.nms <- names(mdata.all)[sel.inx];
    dtInx <- as.numeric(unlist(strsplit(dtInx,split=",")))+1
    infoSet <- readSet(infoSet, "infoSet");
    hit.query<-infoSet$imgSet$enrTables[[mode]]$hits.query
    hit.query <- hit.query[[pathnm]]
    dataSetList <- lapply(sel.nms[dtInx], readDataset);
    hitList <- lapply(dataSetList, function(x) x$sig.mat)
    hitList <- lapply(hitList, function(x) x[x$ids %in% hit.query,])
    
    hitList <- lapply(hitList, function(x){
      colors <- ifelse(x[,1]>0,"red","blue")
      nms <-paste("<font color=",colors,">", "<b>", x$label, "</b>", "</font>",sep="")

      return( paste(unique(nms), collapse="; "))
    } )
    
    hitList <- unlist(hitList)
    
   return(hitList)
  }

  
  
  
}
