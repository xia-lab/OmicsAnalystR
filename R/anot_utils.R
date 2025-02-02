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
AnnotateMicrobiomeData <- function(dataName,org,feattype){
  library("tidyverse")
  dataSet <- readDataset(dataName);
  data <- qs::qread(dataSet$data.raw.path);
  #  mic.vec <- rownames(data);
  
  data.org <<- org
  
  if(feattype=="taxa"){
  mic.vec <- data.frame(dataSet$orig.var.nms,check.names=FALSE);
  rownames(mic.vec) <- NULL
  mic.vec <- data.frame(apply(mic.vec, 1, function(x) gsub(";\\s;", ";", x)),check.names=FALSE) # remove empty taxa
  names(mic.vec) <- "Rank"
  mic.vec <- separate(mic.vec,Rank,into=c("R1","R2", "R3", "R4", "R5", "R6", "R7") , sep = ";")

  if(any(grepl("Bacteria",mic.vec$R1))|any(grepl("Archaea",mic.vec$R1))){
   if("Firmicutes" %in% mic.vec$R1 | "Bacteroidota" %in% mic.vec$R1){
    msg.vec <<- "Please check the taxonomy labels!"
    return(0)
   }
   if(length(unique(mic.vec$R1))==1){
     mic.vec = mic.vec[,-1]
     colnames(mic.vec) <- c("Phylum", "Class", "Order", "Family", "Genus", "Species") 
      }else{
     colnames(mic.vec) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species") 
     }
  }else if(all(mic.vec$R1!="Bacteria"  & mic.vec$R1!="Archaea" & mic.vec$R1!="Firmicutes" & mic.vec$R1!="Bacteroidota" )){
    msg.vec <<- "Please check the taxonomy labels!"
   return(0)
  }else{
   colnames(mic.vec) <- c("Phylum", "Class", "Order", "Family", "Genus", "Species") 
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
  }

  qs::qsave(data, dataSet$data.annotated.path);
  fast.write.csv(data,file=paste0(dataSet$folderName, "/data.annotated.csv"));
  dataSet$enrich_ids <- rownames(data);
  names(dataSet$enrich_ids) = rownames(data);
  if(feattype != "otu"){
    dataSet$m2m <- 1
  }
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
  
  if(org == "NA"){
    msg.vec <<- "Invalid organism!"
    return(0)
  }
  
  dataSet <- readDataset(dataName);
  data.raw <- qs::qread(dataSet$data.raw.path);
  gene.vec <- rownames(data.raw);
  #saveRDS(data.raw,"/Users/lzy/Documents/OmicsAnalystR/data.raw.rds")
  #print(idtype)
  #record the info
  data.org <<- org
  dataSet$idType <- idtype;
  dataSet$gene.org <- org;
  dataSet$gene <- gene.vec;
  
  if(idtype %in% c("mir_id", "mir_acc", "mirnet")){
    enMat <- doIdMapping(gene.vec, idtype);
    if(nrow(enMat) == 0){
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
    data.raw <- data.raw[hit.inx,];
    matched.entrez <- enIDs[hit.inx];

    # now, deal with duplicated entrez id
    # first, average duplicate rows

    ave.data <- RemoveDuplicates(data.raw, "max");

    # then removed duplicated entries
    dup.inx <- duplicated(matched.entrez);
    matched.entrez <- matched.entrez[!dup.inx]
    int.mat <- ave.data[!dup.inx,];
    # update
    data.annotated <-int.mat;
    rownames(int.mat) <- matched.entrez
    if(idtype %in% c("mir_id", "mir_acc", "mirnet")){
      rownames(data.annotated) <- rownames(int.mat);
    }else{
      rownames(data.annotated) <- matched.entrez
    }
    dataSet$enrich_ids <- rownames(int.mat);
    names(dataSet$enrich_ids) <- doEntrez2SymbolMapping(rownames(int.mat))
    
    dataSet$id.type <- "entrez";
  }else{
    data.annotated <- data.raw;
    dataSet$enrich_ids = rownames(data.annotated)
    dataSet$id.type <- "none";
  }
  
  if(idtype != "NA"){
    if(length(unique(enIDs))/length(gene.vec) < 0.3){
      msg <- paste("Less than ", round( length(unique(enIDs))/length(gene.vec) * 100, 2), "% features were mapped in ", dataSet$name);
      msg.vec <<- msg
      return(0)
    }else{
      msg <- paste("A total of ", length(unique(enIDs)), "unique features were mapped");
    }
  }else{
    msg <- paste("There is a total of ", length(unique(gene.vec)), "unique features.");
  }
  msg.vec <<- msg;
  qs::qsave(data.annotated, dataSet$data.annotated.path);
  fast.write.csv(data.annotated,file=paste0(dataSet$folderName, "/data.annotated.csv"));
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
  dataSet <- readDataset(dataName);
  dataSet$name <- dataName
  dataSet$idType <- idtype;
  data <- qs::qread(dataSet$data.raw.path);
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
  #print(length(todo.inx)/length(dataSet$name.map$hit.inx));
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

  data <- as.matrix(data);
  rownames(data) <- unname(dataSet$enrich_ids);
  data <- RemoveDuplicates(data, "mean", quiet=T); # remove duplicates
  data <- as.data.frame(data)
  qs::qsave(data, dataSet$data.annotated.path);
  fast.write.csv(data,file=paste0(dataSet$folderName, "/data.annotated.csv"));
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
  
  dataSet <- readDataset(dataName);
  
  data <- qs::qread(dataSet$data.raw.path);
  qvec <- rownames(data);
  
  dataSet$enrich_ids <- rownames(data)
  names(dataSet$enrich_ids) <- rownames(data)
  qs::qsave(data, dataSet$data.annotated.path);
  fast.write.csv(data,file=paste0(dataSet$folderName, "/data.annotated.csv"));
  RegisterData(dataSet);
  
  
  return(1)
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

queryGeneDB <- function(table.nm, data.org){
  require('RSQLite')
  
  conv.db <- dbConnect(SQLite(), paste(sqlite.path, data.org, "_genes.sqlite", sep="")); 
  db.map <- dbReadTable(conv.db, table.nm)
  dbDisconnect(conv.db); cleanMem();
  
  return(db.map)
}

.query.sqlite <- function(db.con, statement, offline=TRUE){
  rs <- dbSendQuery(db.con, statement);
  res <- fetch(rs, n=-1); # get all records
  dbClearResult(rs);
  if(offline){
    dbDisconnect(db.con);
  }
  cleanMem();
  return(res);
}