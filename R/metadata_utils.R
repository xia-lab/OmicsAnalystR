##################################################
## R scripts for OmicsAnalyst
## Description: Related to meta-data
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

ReadMetaDataFile <- function(metafilename){
  reductionSet <- .get.rdt.set();
  res <- .readMetaData(metafilename,"", "false");
  res$meta.info <- removeXPrefix(res$meta.info);
  meta.types <- rep("disc", ncol(res$meta.info));
  meta.types[res$cont.inx] <- "cont";
  names(meta.types) <- colnames(res$meta.info);
  reductionSet$dataSet$meta.types <- meta.types;
  reductionSet$dataSet$meta.status <- rep("OK", ncol(res$meta.info) );
  reductionSet$dataSet$cont.inx <- res$cont.inx;
  reductionSet$dataSet$disc.inx <- res$disc.inx;
  reductionSet$dataSet$meta.info <- res$meta.info;
  .set.rdt.set(reductionSet)
  return(res$check.bool);
}

GetPrimaryMeta <- function(){
    rdtSet <- .get.rdt.set();
    return(colnames(rdtSet$dataSet$meta.info)[1]);
}

GetMetaDims <- function(){
  rdtSet <- .get.rdt.set();
  dm <- dim(rdtSet$dataSet$meta.info);
  return(dm);
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

# note, try to use the fread, however, it has issues with 
# some windows 10 files "Line ending is \r\r\n. .... appears to add the extra \r in text mode on Windows"
# in such as, use the slower read.table method
.readDataTable <- function(fileName){
  infoSet <- readSet(infoSet, "infoSet");
  msgSet <- infoSet$msgSet; 
  if(length(grep('\\.zip$',fileName,perl=TRUE))>0){
    fileName <- unzip(fileName);
    if(length(fileName) > 1){
      # test if "__MACOSX" or ".DS_Store"
      osInx <- grep('MACOSX',fileName,perl=TRUE);
      if(length(osInx) > 0){
        fileName <- fileName[-osInx];
      }
      dsInx <- grep('DS_Store',fileName,perl=TRUE);
      if(length(dsInx) > 0){
        fileName <- fileName[-dsInx];
      }
      dat.inx <- grep(".[Tt][Xx][Tt]$", fileName);
      if(length(dat.inx) != 1){
        msg.vec <<- "More than one text files (.txt) found in the zip file.";
        return(NULL);
      }
    }
  }
  dat <- try(data.table::fread(fileName, header=TRUE, check.names=FALSE, data.table=FALSE));
  rm.inx <- apply(dat,2,function(x){all(is.na(x))});
  dat <- dat[,!rm.inx];
  if(class(dat) == "try-error"){
    #try to use "tr" to remove double return characters
    trFileName <- paste("tr -d \'\\r\' <", fileName);
    dat <- try(data.table::fread(trFileName, header=TRUE, check.names=FALSE, data.table=FALSE));
    if(class(dat) == "try-error"){
      print("Using slower file reader ...");
      formatStr <- substr(fileName, nchar(fileName)-2, nchar(fileName))
      if(formatStr == "txt"){
        dat <-try(read.table(fileName,header=TRUE,comment.char = "", check.names=F, as.is=T));
      }else{ # note, read.csv is more than read.table with sep=","
        dat <-try(read.csv(fileName,header=TRUE,comment.char = "", check.names=F, as.is=T));
      }  
    }
  }
  if(class(dat) == "try-error"){
    msg.vec <<- "Failed to read the data table! Please check your data format.";
    return(NULL);
  }
  
  # need to remove potential empty columns
  dat <- dat[!sapply(dat, function(x) all(x == "" | is.na(x)))];
  return(dat);
}

.readMetaData <- function(metafileName,datOrig,metaContain) {
  infoSet <- readSet(infoSet, "infoSet");
  msgSet <- infoSet$msgSet; 
  na.msg = ""
  if(is.null(msg.vec)){
    msg <-""
  }else{
    msg <- msg.vec;
  }

  #any warning or error, 0 error, 1 success, 2 warning
  check.bool <- 1;
  if(metaContain=="true"){
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
        cls.lbls <- ClearFactorStrings(cls.nm, cls.lbls);
        meta.info[[cls.nm]] <- cls.lbls;
      }
    }else{
      msg.vec <<-"No metadata labels #CLASS found in your data!";
      return(NULL);
    }
    
    meta.info <- data.frame(meta.info);
    rownames(meta.info) = colnames(datOrig)[-1]
  }else{ # metadata input as an individual table
    mydata <- try(data.table::fread(metafileName, header=TRUE, check.names=FALSE, data.table=FALSE));
    if(class(mydata) == "try-error"){
      msg.vec <<- "Failed to read the metadata table! Please check your data format.";
      return(NULL);
    }
    mydata[is.na(mydata)] <- "NA";
    # look for #NAME, store in a list
    sam.inx <- grep("^#NAME", colnames(mydata)[1]);
    if(length(sam.inx) > 0){
      smpl_nm<-mydata[,1];
      smpl_var<-colnames(mydata[-1]);
    }else{
      msg.vec <<- "Please make sure you have the label #NAME in your sample data file!";
      return(NULL);
    }
    
    # covert to factor
    mydata <-data.frame(lapply(1:ncol(mydata),function(x){
      mydata[,x]=unlist(ClearFactorStrings(mydata[,x]))
    }));
    
    mydata <- mydata[,-1,drop=F]; # converting to character matrix as duplicate row names not allowed in data frame.
    if(nrow(mydata)==1){
      msg.vec <<- "Only one sample in the dataset or the metadata file must be transposed!";

      return(NULL);
    }
    rownames(mydata) <- smpl_nm;
    colnames(mydata) <- smpl_var;
    
    #Check group label names for spaces and replace with underscore
    meta.info <- data.frame(mydata,check.names=FALSE);
    if(any(grepl("[[:blank:]]", names(meta.info)))){
      names(meta.info) <- gsub("\\s+","_", names(meta.info));
      na.msg1 <- c(na.msg1, "Blank spaces in group names are replaced with underscore '_'");
    }
    
  }
    
  disc.inx <- GetDiscreteInx(meta.info);

  # make sure categorical metadata are valid names
  if(class(meta.info[,disc.inx]) == "data.frame"){
    meta.info[,disc.inx] <- apply(meta.info[,disc.inx], 2, function(x){x[x != "NA"] = make.names(x[x != "NA"]); return(x)});
    meta.info[,disc.inx] <- lapply(meta.info[,disc.inx], function(x) factor(x));
  }else{
    x <- meta.info[,disc.inx];
    x[x != "NA"] = make.names(x[x != "NA"])
    x <- factor(x);
    meta.info[,disc.inx] <- x;
  }

  if(sum(disc.inx) == length(disc.inx)){
    na.msg <- c(na.msg,"All metadata columns are OK!")
  }else{
    bad.meta<- paste(names(disc.inx)[!disc.inx], collapse="; ");
    na.msg <- c(na.msg, paste0("<font style=\"color:red\">Detected presence of unique values in the following columns: <b>", bad.meta, "</b></font>","Please make sure the metadata is in right format! You can use meta editor to update the information !"));
  }
  
  cont.inx <- GetNumbericalInx(meta.info);
  cont.inx <- !disc.inx & cont.inx; # discrete is first
  
  rmcol <- intersect(which(!disc.inx),which(!cont.inx ))
  
  if(sum(cont.inx)>0){
    # make sure the discrete data is on the left side
    meta.info <- cbind(meta.info[,disc.inx, drop=FALSE], meta.info[,cont.inx, drop=FALSE]);
  }


  disc.inx <- disc.inx[colnames(meta.info)]
  cont.inx <- cont.inx[colnames(meta.info)]

  meta.info <- as.data.frame(meta.info);

  check.inx <-apply(meta.info , 2, function(x){ ( sum(is.na(x))/length(x) + sum(x=="NA")/length(x) + sum(x=="")/length(x) ) >0})
  
  init <- 1;

  cls.vec <- vector()
  lowrep.vec <- vector()
  toolow.vec <- vector();

  for(i in 1:ncol(meta.info)){
      cls.lbl <- meta.info[,i];
      qb.inx <- tolower(cls.lbl) %in% c("qc", "blank");
      if(sum(qb.inx) > 0){
        cls.Clean <- as.factor(as.character(cls.lbl[!qb.inx])); # make sure drop level
      } else {
        cls.Clean <- cls.lbl;
      }
      meta.name <- colnames(meta.info)[i]
      min.grp.size <- min(table(cls.Clean));
      cls.num <- length(levels(cls.Clean));


    # checking if too many groups but a few samples in each group
      if(cls.num/min.grp.size > 3 && !tolower(meta.name) %in% c("subject", "time")){
        if(init == 1){
           isNum <- grepl("^-?[0-9.]+$", cls.Clean);
           if(all(isNum)){
             cls.vec <- c(cls.vec, meta.name)
           }else{
             if(!check.inx[i]){
             toolow.vec <- c(toolow.vec, meta.name)
             }
           }
        }
    # checking if some groups have low replicates
      } else if(min.grp.size < 3 | cls.num < 2){
        if(init == 1){
           isNum <- grepl("^-?[0-9.]+$", cls.Clean);
           if(all(isNum)){
             cls.vec <- c(cls.vec, meta.name)
           }else{
             if(!check.inx[i] && !meta.name %in% toolow.vec){
             lowrep.vec <- c(lowrep.vec, meta.name)
             }
           }
        }
      }
    
  }

  if(length(toolow.vec)>0 && init == 1){
    msg <- paste0( "<b>",paste0(toolow.vec, collapse=", "),"</b>", " meta-data factors have too many groups with low replicates (less than 3) per group.");
    check.bool = 2;
   }

  if(length(lowrep.vec)>0 && init == 1){
    msg <-paste0( "<b>",paste0(lowrep.vec, collapse=", "),"</b>", " meta-data factors have some groups with low replicates (less than 3) per group.");
    check.bool = 2;
}

  if(nrow(meta.info) < 4){
    msg <-  paste0("Less than 4 samples are detected. More samples are required!");
    check.bool = 0;
}

  msg.vec <<- paste(na.msg, msg);
  return(list(meta.info=meta.info,disc.inx=disc.inx,cont.inx=cont.inx,check.bool=check.bool))
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

GetMetaTypes <- function(colNm="NA"){
  rdtSet <- .get.rdt.set();
  if(colNm=="NA"){
    meta.types <- rdtSet$dataSet$meta.types
  }else{
    meta.types <- rdtSet$dataSet$meta.types[colNm]
  }
  return(unname(meta.types));
}

SetMetaTypes <- function(metaTypes.vec){
  rdtSet <- .get.rdt.set();
  names(metaTypes.vec) <- colnames(rdtSet$dataSet$meta.info)
  rdtSet$dataSet$meta.types <- metaTypes.vec;
  return(.set.rdt.set(rdtSet));
}

UpdateMetaOrder <- function(metacol){
  rdtSet <- .get.rdt.set();

  meta <- rdtSet$dataSet$meta.info
  if(length(metaVec)>0 & metacol %in% colnames(meta)){
   meta[,metacol] <- factor(as.character(meta[,metacol]),levels = metaVec)
   rdtSet$dataSet$meta.info <- meta
  }else{
  msg.vec <- "The metadata column is empty! Please check your selection!"
    return(0)
  }
  .set.rdt.set(rdtSet);
  return(1)
}

GetMetaDataSmpl <- function(){
  rdtSet <- .get.rdt.set();
  return(rownames(rdtSet$dataSet$meta.info));
}

GetMetaCell <- function(ridx=1,cidx=1){
  rdtSet <- .get.rdt.set();
  return(rdtSet$dataSet$meta.info[ridx,cidx]);
}

removeXPrefix <- function(df) {
  for (col in 1:ncol(df)) {
    values <- df[[col]]
    # Check if all values start with "X"
    if (all(grepl("^X", values))) {
      df[[col]] <- sub("^X", "", values)  # Remove "X" prefix
    }
  }
  return(df)
}


UpdateMetaStatus <- function(colNm){
  rdtSet <- .get.rdt.set();
  cidx <- which(colnames(rdtSet$dataSet$meta.info)==colNm)
  old = rdtSet$dataSet$meta.types[cidx];
  if(old=="disc"){
    if(all(is.na( as.numeric(as.character(rdtSet$dataSet$meta.info[,cidx]))))){
      msg.vec <<- "Category metadata cannot be continuous data!"
      return(0)
    }
    rdtSet$dataSet$meta.types[cidx] = "cont"
    rdtSet$dataSet$disc.inx[cidx]=FALSE;
    rdtSet$dataSet$cont.inx[cidx]=TRUE;
  }else{
    if(all(!duplicated(as.character(rdtSet$dataSet$meta.info[,cidx])))){
      msg.vec <<- "No duplicates were detected! The metadata cannot be discrete!"
      return(0)
    }
    rdtSet$dataSet$meta.types[cidx] = "disc"
    rdtSet$dataSet$disc.inx[cidx]=TRUE;
    rdtSet$dataSet$cont.inx[cidx]=FALSE;
  }
  new = rdtSet$dataSet$meta.types[cidx]
  msg.vec <<- paste0("Metadata type of ",colnames(rdtSet$dataSet$meta.info)[cidx]," has been changed to ", new, " !")
  .set.rdt.set(rdtSet)
  return(1);
}

GetSampleNm <- function(ridx=1){
  rdtSet <- .get.rdt.set();
  return( rownames(rdtSet$dataSet$meta.info)[ridx]);
}

DeleteSample <- function(samplNm){
  rdtSet <- .get.rdt.set();
  rdtSet$dataSet$meta.info <- rdtSet$dataSet$meta.info[rownames(rdtSet$dataSet$meta.info)!=samplNm,,drop=F]
  sel.nms <- names(mdata.all)
  for(nm in sel.nms){
    dataSet <- qs::qread(nm);
    dataSet$meta <-  rdtSet$dataSet$meta.info;
    dataSet$data.proc <- dataSet$data.proc[,colnames(dataSet$data.proc)!=samplNm]
  }
 
  .set.rdt.set(rdtSet);
  return(1);
}

ResetMetaTab <- function(){
  rdtSet <- .get.rdt.set();
  rdtSet$dataSet <- rdtSet$dataSet.origin

  sel.nms <- names(mdata.all)
  for(nm in sel.nms){
    dataSet <- qs::qread(nm);
    dataSet$data.proc <- dataSet$data.proc.origin;
    dataSet$meta <- rdtSet$dataSet.origin$meta.info;
    RegisterData(dataSet)
  }
  .set.rdt.set(rdtSet);
  return(1);
}

GetDiscMetas <- function(){
  keepVec<-keepVec
  rdtSet <- .get.rdt.set();
  meta <- rdtSet$dataSet$meta.info
  if(length(keepVec)>0){
    keepidx <- which(keepVec %in% colnames(meta) ) 
    keepidx <- intersect(keepidx,which(rdtSet$dataSet$disc.inx))
  }else{
    keepidx <-  which(rdtSet$dataSet$disc.inx)
  }
  colnms<- colnames(meta)[keepidx]
  return(colnms);
}

GetMetaDataCol <- function(colnm){
  rdtSet <- .get.rdt.set();
  cls<-levels(rdtSet$dataSet$meta.info[,colnm])
  return(cls[cls!="NA"]);
}

DeleteMetaCol <- function(metaCol){
  rdtSet <- .get.rdt.set();
  meta <- rdtSet$dataSet$meta.info
  if(metaCol %in% colnames(meta)){
    idx = which(colnames(meta)==metaCol)
    rdtSet$dataSet$meta.info <- meta[,-idx,drop=F]
    rdtSet$dataSet$meta.types <- rdtSet$dataSet$meta.types[-idx]
    rdtSet$dataSet$meta.status <- rdtSet$dataSet$meta.status[-idx]
    rdtSet$dataSet$disc.inx <- rdtSet$dataSet$disc.inx[-idx]
    rdtSet$dataSet$cont.inx <- rdtSet$dataSet$cont.inx[-idx]
    if(!exists("rmMetaCol",dataSet)){
      dataSet$rmMetaCol <- vector()
    }
    dataSet$rmMetaCol <- unique(c(dataSet$rmMetaCol,metaCol))
  }
  .set.rdt.set(rdtSet);
  return(1);
}

UpdateMetaName <-  function(oldvec,newvec){
  rdtSet <- .get.rdt.set();
  idx <- which(colnames(rdtSet$dataSet$meta.info)==oldvec)
  if(length(idx)==1){
    colnames(rdtSet$dataSet$meta.info)[idx] <- 
      names(rdtSet$dataSet$disc.inx)[idx] <- names(rdtSet$dataSet$cont.inx)[idx] <-
      names(rdtSet$dataSet$meta.types)[idx] <- names(rdtSet$dataSet$meta.status)[idx] <- newvec
  }else{
    return(0)
  }
  .set.rdt.set(rdtSet);
  return(1);
}

CheckMetaNAs <- function(){
  rdtSet <- .get.rdt.set();
  meta <- rdtSet$dataSet$meta.info
  if(any(is.na(meta))|any(meta=="")|any(meta=="NA")){
    return(2)
  }else{
   return(1)
  }

}

CheckEditRes <- function(){
  
  rdtSet <- .get.rdt.set();
  meta <- rdtSet$dataSet$meta.info
  # use first column by default
  cls <- droplevels(meta[meta[,1]!="NA",1])
  # check class info
  min.grp.size <- min(table(cls));
  cls.num <- length(levels(cls));
  if(min.grp.size<2){
    msg.vec <<- paste0( "No replicates were detected for group  ",as.character(cls[which( table(cls)<2)])," in  ",colnames(meta)[1])
    return(0)
  }
  sel.nms <- names(mdata.all)
  data.list = list();
  for(i in 1:length(sel.nms)){
    dataSet <- qs::qread(sel.nms[i])
    dataSet$meta <- rdtSet$dataSet$meta.info;
    dataSet$data.proc <- dataSet$data.proc[,colnames(dataSet$data.proc) %in% rownames(meta)]
    dataSet$data.proc <- dataSet$data.proc[,match(rownames(meta),colnames(dataSet$data.proc))]
    dataSet$cls <- cls
    dataSet$rmidx <- which(meta[,1]=="NA")
    dataSet$min.grp.size <- min.grp.size;
    dataSet$cls.num <- cls.num;
    RegisterData(dataSet);
  }
  
  .set.rdt.set(rdtSet)
  return(1)
  }


UpdateSampInfo <-  function(rowNm,colNm,cell){
  rdtSet <- .get.rdt.set();
  meta <- rdtSet$dataSet$meta.info
  ridx <- which(rownames(meta)==rowNm)
  if(colNm==""){
    if(rowNm !=cell){
    rownames(meta)[ridx]=cell
    sel.nms <- names(mdata.all)
    data.list = list();
    for(i in 1:length(sel.nms)){
    dataSet <- qs::qread(sel.nms[i])
    nmidx<-which(colnames(dataSet$data.proc)==rowNm)
    colnames(dataSet$data.proc)[nmidx] <- cell
    RegisterData(dataSet);
  }
    }
  
  }else{  
    cidx<- which(colnames(meta)==colNm)
    if(cell!= as.character(meta[ridx,cidx])){
      if(cell %in% levels(meta[,cidx])){
        meta[ridx,cidx] = cell
      }else{
        levels(meta[,cidx]) <- c(levels(meta[,cidx]), cell)
        meta[ridx,cidx] = cell
      }
      meta[,cidx] <- droplevels(meta[,cidx])
    }
  }
  rdtSet$dataSet$meta.info = meta
  .set.rdt.set(rdtSet)
  return(1);
}

UpdatePrimaryMeta <- function(primaryMeta){
  rdtSet <- .get.rdt.set(); 
  meta <- rdtSet$dataSet$meta.info
  if(primaryMeta %in% colnames(meta)){
    cidx <- which(colnames(meta)==primaryMeta)
    rdtSet$dataSet$meta.info<-cbind(meta[,cidx,drop=F],meta[,-cidx,drop=F])
    rdtSet$dataSet$meta.types = c(rdtSet$dataSet$meta.types[cidx],rdtSet$dataSet$meta.types[-cidx])
    rdtSet$dataSet$disc.inx=c(rdtSet$dataSet$disc.inx[cidx],rdtSet$dataSet$disc.inx[-cidx])
    rdtSet$dataSet$cont.inx=c(rdtSet$dataSet$cont.inx[cidx],rdtSet$dataSet$cont.inx[-cidx])
  }else{
    msg.vec <- "The metadata column is empty! Please check your selection!"
    return(0)
  }
  .set.rdt.set(rdtSet);
  return(1)
}

#'Generate correlation heatmap for metadata
#'@description Plot correlation coefficients between metadata
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param cor.opt Meethod for computing correlation coefficient
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotMetaCorrHeatmap <- function(cor.opt="pearson", imgName, dpi=96, imgFormat="png"){
  
  imgName <- paste(imgName, "dpi", dpi, ".", imgFormat, sep="");
  dpi <- as.numeric(dpi);
  rdtSet <- .get.rdt.set();
  metaData <- rdtSet$dataSet$meta.info
  meta.types <- rdtSet$dataSet$meta.types
  disc.inx <- which(meta.types == "disc")
  cont.inx <- which(meta.types == "cont")
  meta.num <- ncol(metaData)
  
  textSize = 4;
  if(meta.num > 25){
    w <- 24
    h <- 18
    textSize = 3.5;
  }else if(meta.num > 10){
    w <- 16
    h <- 12
  } else {
    w <- 10
    h <- 7.5
  }
  
  library(reshape2)
  load_ggplot();
  library(scales);
  
  metaData[metaData == "NA"] <- NA;
  for(i in c(1:length(disc.inx))){
    metaData[,disc.inx[i]] <- as.integer(metaData[,disc.inx[i]], na.rm = TRUE);
  }
  for(i in c(1:length(cont.inx))){
    metaData[,cont.inx[i]] <- as.numeric(as.character(metaData[,cont.inx[i]], na.rm = TRUE));
  }
  

  cormat <- round(cor(metaData, method=cor.opt, use="pairwise.complete.obs"),3);
  upper_tri <- get_upper_tri(cormat);
  melted_cormat <- melt(upper_tri, na.rm = TRUE);
  
  ggheatmap <- ggplot2::ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white")+ scale_y_discrete("Var1", position="right") +
    scale_fill_gradient2(low = muted("blue"), mid="white", high = muted("red"), midpoint = 0,
                         limit = c(-1,1), space = "Lab", name="Correlation") + theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.y = element_blank(), axis.text.y.right = element_text(),
          legend.direction = "vertical", legend.position="left")+ coord_fixed();
  
  ggheatmap <- ggheatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = textSize);
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=imgFormat, bg="white");
  print(ggheatmap);
  dev.off();
  return(1);
}

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
