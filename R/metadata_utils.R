##################################################
## R scripts for OmicsAnalyst
## Description: Related to meta-data
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

ReadMetaDataFile <- function(metafilename){
  save.image("meta.RData");
  reductionSet <- .get.rdt.set();
  res <- .readMetaData(metafilename,"", F)
  meta.types <- rep("disc", ncol(res$meta.info));
  meta.types[res$cont.inx] <- "cont";
  names(meta.types) <- colnames(res$meta.info);
  print(meta.types);
  reductionSet$dataSet$meta.types <- meta.types;
  #reductionSet$dataSet$meta.status <- rep("OK", ncol(res$meta.info) - 1);
  reductionSet$dataSet$meta.status <- rep("OK", ncol(res$meta.info) );
  reductionSet$dataSet$cont.inx <- res$cont.inx;
  reductionSet$dataSet$disc.inx <- res$disc.inx;
  reductionSet$dataSet$meta.info <- res$meta.info;

  return(.set.rdt.set(reductionSet));
}


GetUniqueMetaNames <-function(metadata){
  rdtSet <- .get.rdt.set();
  data.type <- rdtSet$dataSet[["meta.types"]][metadata];
  if(data.type == "cont"){
    return("--- NA ---");
  } else {
    return(levels(as.factor(rdtSet$dataSet$meta.info[,metadata])));
  }
}

.readMetaData <- function(metafileName,datOrig, metaContain) {
  msgSet <- readSet(msgSet, "msgSet");
  na.msg = ""
  if(is.null(msgSet$current.msg)){
    msg <-""
  }else{
    msg <- msgSet$current.msg
  }
  if(metaContain){
    meta.info <- list();
    # look for #CLASS, could have more than 1 class labels, store in a list
    cls.inx <- grep("^#CLASS", datOrig[,1]);
    if(length(cls.inx) > 0){
      for(i in 1:length(cls.inx)){
        inx <- cls.inx[i];
        cls.nm <- substring(datOrig[inx, 1],2); # discard the first char #
        if(nchar(cls.nm) > 6){
          cls.nm <- substring(cls.nm, 7); # remove class
        }
        if(grepl("[[:blank:]]", cls.nm)){
          cls.nm<- gsub("\\s+","_", cls.nm);
          msg <- c(msg, " Blank spaces in group names are replaced with underscore '_'! ");
        }
        cls.lbls <- datOrig[inx, -1];
        # test NA
        na.inx <- is.na(cls.lbls);
        cls.lbls[na.inx] <- "NA";
        cls.lbls <- ClearFactorStrings(cls.lbls,cls.nm);
        meta.info[[cls.nm]] <- cls.lbls;
      }
    }else{
      msgSet$current.msg <- "No metadata labels #CLASS found in your data!";
      saveSet(msgSet, "msgSet");
      return(NULL);
    }
    meta.info <- data.frame(meta.info);
  }else{ # metadata input as an individual table
    mydata <- try(data.table::fread(metafileName, header=TRUE, check.names=FALSE, data.table=FALSE));
    mydata[is.na(mydata)] <- "NA";
    # look for #NAME, store in a list
    sam.inx <- grep("^#NAME", colnames(mydata)[1]);
    if(length(sam.inx) > 0){
      smpl_nm<-mydata[,1];
      smpl_var<-colnames(mydata[-1]);
    }else{
      msgSet$current.msg <- "Please make sure you have the label #NAME in your sample data file!";
      saveSet(msgSet, "msgSet");
      return(NULL);
    }

    # covert to factor
    #mydata <-data.frame(lapply(1:ncol(mydata),function(x){
    #  mydata[,x]=unlist(ClearFactorStrings(mydata[,x], colnames(mydata)[x]))
    #}))

    mydata <- mydata[,-1]; # converting to character matrix as duplicate row names not allowed in data frame.
    if(nrow(mydata)==1){
      current.msg <- "Only one sample in the dataset or the metadata file must be transposed!";
      saveSet(msgSet, "msgSet");
      return(NULL);
    }
    rownames(mydata) <- smpl_nm;
    colnames(mydata) <- smpl_var;

    # empty cell or NA cannot be tolerated in metadata
    na.inx  <- is.na(mydata);
    na.msg <- na.msg1 <- NULL;
    if(sum(na.inx) > 0){
      na.msg1 <- paste("A total of", sum(na.inx), "empty or NA values were detected. Please update in using metadata editor");
    }
    #Check group label names for spaces and replace with underscore
    meta.info <- data.frame(mydata,check.names=FALSE);
    if(any(grepl("[[:blank:]]", names(meta.info)))){
      names(meta.info) <- gsub("\\s+","_", names(meta.info));
      na.msg1 <- c(na.msg1, "Blank spaces in group names are replaced with underscore '_'");
    }
  }
  disc.inx <- GetDiscreteInx(meta.info);
  if(sum(disc.inx) == length(disc.inx)){
    na.msg <- c(na.msg,"All metadata columns are OK!")
  }else{
    bad.meta<- paste(names(disc.inx)[!disc.inx], collapse="; ");
    na.msg <- c(na.msg, paste0("<font style=\"color:red\">Detected presence of unique values in the following columns: <b>", bad.meta, "</b></font>","Please use meta editor to update or these columns will be removed in the following analysis!"));
  }
  cont.inx <- GetNumbericalInx(meta.info);
  cont.inx <- !disc.inx & cont.inx; # discrete is first
  #if(sum(cont.inx)>0){
    # make sure the discrete data is on the left side
  #  meta.info <- cbind(meta.info[,disc.inx, drop=FALSE], meta.info[,cont.inx, drop=FALSE]);
  #}
  msgSet$na.msg <- na.msg
  saveSet(msgSet, "msgSet");
  print(head(meta.info));
  return(list(meta.info=meta.info,disc.inx=disc.inx,cont.inx=cont.inx))
}

GetPrimaryType <- function(analysis.var){
  rdtSet <- .get.rdt.set();
  primary.type <- unname(rdtSet$dataSet$meta.types[analysis.var]);
  return(primary.type);
}

GetMetaDataGroups <- function(){
  rdtSet <- .get.rdt.set();
  groups <- colnames(rdtSet$dataSet$meta.info);
  return(groups);
}

GetMetaDataStatus <- function(){
  rdtSet <- .get.rdt.set();
  res <- unname(rdtSet$dataSet$meta.status);
  return(res);
}

UpdateMetaType <-function(metadata="NA", type="disc"){
  rdtSet <- .get.rdt.set();
  rdtSet$dataSet$meta.types[metadata] = type;
  return(.set.rdt.set(rdtSet));
}

GetMetaTypes <- function(){
  rdtSet <- .get.rdt.set();
  return(unname(rdtSet$dataSet$meta.types));
}

SetMetaTypes <- function(metaTypes.vec){
  rdtSet <- .get.rdt.set();
  names(metaTypes.vec) <- colnames(rdtSet$dataSet$meta.info)
  rdtSet$dataSet$meta.types <- metaTypes.vec;
  return(.set.rdt.set(rdtSet));
}

UpdateMetaOrder <- function(metaName){
  rdtSet <- .get.rdt.set();
  if(exists('meta.ord.vec')){
    metadata <- rdtSet$dataSet$meta.info[,metaName];
    rdtSet$dataSet$meta.info[,metaName] <- factor(as.character(metadata), levels=meta.ord.vec)
  }
  return(.set.rdt.set(rdtSet));
}
