##################################################
## R script for OmicsAnalyst
## Description: functions for id annotation
##
## Authors: 
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#'Annotate microbiome data to internal database
#'@param fileName File name of data file.
#'@param org Organism three-letters id
#'@param idtype ID type
#'@author Guangyan Zhou \email{guangyan.zhou@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
AnnotateMicrobiomeData <- function(dataName,org){
  library("tidyverse")
  dataSet <- qs::qread(dataName);
  
  data <- dataSet$data.raw;
  #  mic.vec <- rownames(data);
  
  data.org <<- org
  
  mic.vec <- data.frame(dataSet$orig.var.nms,check.names=FALSE);
  rownames(mic.vec) <- NULL
  mic.vec <- data.frame(apply(mic.vec, 1, function(x) gsub(";\\s;", ";", x)),check.names=FALSE) # remove empty taxa
  
  names(mic.vec) <- "Rank"
  if(all(grepl("Bacteria",mic.vec$Rank))){
    
    mic.vec <- separate(mic.vec,Rank,into=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species") , sep = ";")
    mic.vec <- mic.vec[,-1]
  }else if(all(!(grepl("Bacteria",mic.vec$Rank)))){
    mic.vec <- separate(mic.vec,Rank,into=c("Phylum", "Class", "Order", "Family", "Genus", "Species") , sep = ";")
    
  }else{
    msg.vec <<- "Please check the taxonomy labels!"
    return(0)
    
  }
  
  idx <- names(which(apply(mic.vec, 2, function(x) length(which(is.na(x))))==nrow(mic.vec)))
  mic.vec[,idx] <- NULL
  
  mic.vec <- apply(mic.vec,2, function(x) gsub("[a-z]__","",x))
  mic.vec <- apply(mic.vec,2, function(x) gsub("[a-z]_","",x))
  
  mic.vec[which(is.na(mic.vec)|mic.vec=="NA"|mic.vec=="")] <- "Not_Assigned"
  
  rownames(data) <- apply(mic.vec,1,function(x) paste(x,collapse=";"))
  rownames(mic.vec)<- rownames(data)
  
  dataSet$taxa_table <- mic.vec;
  dataSet$taxaLev <- names(mic.vec)
  dataSet$data.taxa <- list()

  
  for(i in 1:length(colnames(mic.vec))){
    dataSet$data.taxa[[colnames(mic.vec)[i]]] <- as.matrix(data)
    rownames(dataSet$data.taxa[[i]]) <- mic.vec[,i]
    dataSet$data.taxa[[i]] <- RemoveDuplicates(dataSet$data.taxa[[i]], "sum", quiet=T); # remove duplicates
    dataSet$data.taxa[[i]] <- as.data.frame(dataSet$data.taxa[[i]])
  }

  dataSet$data.annotated <- data
  dataSet$enrich_ids <- rownames(data);
  names(dataSet$enrich_ids) = rownames(data);
  dataSet$m2m <- 1
  RegisterData(dataSet);
  return(1)
}


#'Annotate gene data to internal database
#'@param fileName File name of data file.
#'@param org Organism three-letters id
#'@param idtype ID type
#'@author Guangyan Zhou \email{guangyan.zhou@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
AnnotateGeneData <- function(dataName, org, idtype){
  
  #if(dataSet$name != dataName){
  dataSet <- qs::qread(dataName);
  #}
  
  if(org == "NA"){
    msg.vec <<- "Invalid organism!"
    return(1)
  }
  
  data <- dataSet$data.raw;
  gene.vec <- rownames(data);
  
  #record the info
  data.org <<- org
  dataSet$q.type.gene <- idtype;
  dataSet$gene.org <- org;
  dataSet$gene <- gene.vec;
  
  if(idtype %in% c("mir_id", "mir_acc", "mirnet")){
    enMat <- doIdMapping(gene.vec, idtype);
    if(enMat == 0){
      msg.vec <<- "Please make sure the right ID type and organism are selected.";
      return(0)
    }
    enIDs <- gene.vec
    enIDs[enIDs %in% enMat[,2]] <- enMat[,1]
  }else{
    enIDs <- doGeneIDMapping(gene.vec, org, idtype);
  }
  
  dataSet$rawToEntrez <- enIDs
  names(dataSet$rawToEntrez) <- gene.vec;
  
  if(idtype == "kos"){
    kos <- enIDs$kos;
    enIDs <- enIDs$entrezs;
    dataSet$kos.name.map <- kos
  }
  
  # Handle case when only KOs are mapped with no corresponding entrez id
  na.inx <- is.na(enIDs);
  
  if(sum(!na.inx) == 0 && idtype == "kos"){
    na.inx <- is.na(kos);
  }
  
  dataSet$gene.name.map <- list(
    hit.values=enIDs,
    match.state = ifelse(is.na(enIDs), 0, 1)
  );
  
  hit.inx <- which(!is.na(enIDs));
  matched.len <- length(hit.inx);
  if(matched.len > 1){
    data.proc <- dataSet$data.raw[hit.inx,];
    matched.entrez <- enIDs[hit.inx];
    
    
    # now, deal with duplicated entrez id
    # first, average duplicate rows
    
    myave <- function (x, ...) {
      n <- length(list(...))
      if (n) {
        g <- interaction(...)
        split(x, g) <- lapply(split(x, g), mean, na.rm=T)
      }
      else x[] <- FUN(x, na.rm=T)
      return(x);
    }
    ave.data <- apply(data.proc, 2, myave, matched.entrez); 
    # then removed duplicated entries
    dup.inx <- duplicated(matched.entrez);
    matched.entrez <- matched.entrez[!dup.inx]
    int.mat <- ave.data[!dup.inx,];
    # update
    dataSet$data.annotated <-int.mat;
    rownames(int.mat) <- matched.entrez
    if(idtype %in% c("mir_id", "mir_acc", "mirnet")){
      rownames(dataSet$data.annotated) <- rownames(int.mat);
    }else{
      rownames(dataSet$data.annotated) <- matched.entrez
    }
    dataSet$enrich_ids = rownames(int.mat);
    names(dataSet$enrich_ids) = doEntrez2SymbolMapping(rownames(int.mat))
    
    dataSet$id.type <- "entrez";
  }else{
    dataSet$data.annotated <- dataSet$data.raw;
    dataSet$enrich_ids = rownames(dataSet$data.annotated)
    dataSet$id.type <- "none";
  }
  
  if(idtype != "NA"){

    if(length(unique(enIDs))/length(gene.vec) < 0.3){
      msg <- paste("Less than ", round( length(unique(enIDs))/length(gene.vec) * 100, 2), "% features were mapped");
      msg.vec <<- msg
      return(0)
    }else{
      msg <- paste("A total of ", length(unique(enIDs)), "unique features were mapped");
    }
  }else{
    msg <- paste("There is a total of ", length(unique(gene.vec)), "unique features.");
  }
  msg.vec <<- msg
  
  RegisterData(dataSet);
  return(1)
}


#'Annotate metabolite data to internal database
#'@param fileName File name of data file.
#'@param idtype ID type
#'@author Guangyan Zhou \email{guangyan.zhou@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
AnnotateMetaboliteData <- function(dataName, idtype){
  dataSet <- qs::qread(dataName);
  dataSet$name <- dataName
  
  data <- dataSet$data.raw;
  qvec <- rownames(data);
  
  # record all the data
  if(!exists("name.map", where = dataSet)){
    dataSet$name.map <- list();
  }
  
  # map to cpd db from metaboanalyst
  dataSet <- MetaboliteMappingExact(dataSet, qvec, idtype)
  
  # do some sanity check
  todo.inx <- which(is.na(dataSet$name.map$hit.inx));
  resint <- 1;
  print(length(todo.inx)/length(dataSet$name.map$hit.inx));
  if(length(todo.inx)/length(dataSet$name.map$hit.inx) > 0.5){
    msg <- c("Over half of the compound IDs could not be matched to our database. Please make 
             sure that correct compound IDs or common compound names are used.");
    resint <- 2;
  }else if (length(todo.inx) > 15){
    msg <- c("There are >15 compounds without matches. You can either proceed or if necessary, update these compound IDs and upload again.");   
    resint <- 1;     
  }else{
    msg <- paste0("Name matching OK!", " A total of ", nrow(data), " are mapped.");  
    resint <- 1; 
  }
  
  msg <- c(msg, paste0("A total of ", length(na.omit(dataSet$name.map$hit.inx)), " are mapped."))

  data <- as.matrix(dataSet$data.raw);
  rownames(data) <- unname(dataSet$enrich_ids);
  data <- RemoveDuplicates(data, "mean", quiet=T); # remove duplicates
  data <- as.data.frame(data)
  dataSet$data.annotated <- data
  RegisterData(dataSet);
  msg.vec <<- msg
  
  return(resint)
}

#'Skipping data annotation
#'@param fileName File name of data file.
#'@param idtype ID type (not used);
#'@author Guangyan Zhou \email{guangyan.zhou@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
SkippingAnnotation <- function(dataName, idtype){
  
  dataSet <- qs::qread(dataName);
  
  data <- dataSet$data.raw;
  qvec <- rownames(data);
  
  dataSet$enrich_ids <- rownames(dataSet$data.raw)
  names(dataSet$enrich_ids) <- rownames(dataSet$data.raw)
  dataSet$data.annotated <- dataSet$data.raw
  
  RegisterData(dataSet);
  
  
  return(1)
}
