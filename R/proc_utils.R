
#'Mapping from different gene IDs
#'@param q.vec Input the genes to be mapped.
#'@param org Input the name of the organism.
#'@param q.type Input the query-type.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
doGeneIDMapping <- function(q.vec, org, type){
  org <- data.org
  library(RSQLite)  
  db.path <- paste0(sqlite.path, org, "_genes.sqlite");
  con <- dbConnect(SQLite(), db.path ); 
  
  if(type == "symbol"){
    db.map = dbReadTable(con, "entrez")
    hit.inx <- match(q.vec, db.map[, "symbol"]);
  }else if(type == "entrez"){
    db.map = dbReadTable(con, "entrez")
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
  }else{
    if(type == "gb"){
      db.map = dbReadTable(con, "entrez_gb");
    }else if(type == "embl_gene" || type == "embl"){
      db.map = dbReadTable(con, "entrez_embl_gene");
    }else if(type == "uniprot"){
      db.map = dbReadTable(con, "entrez_uniprot");
    }else if(type == "refseq"){
      db.map = dbReadTable(con, "entrez_refseq");
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
    }else if(type == "kos"){
      db.map = dbReadTable(con, "entrez_ortholog");
    }
    hit.inx <- match(q.vec, db.map[, "accession"]);
  }
  
  entrezs=db.map[hit.inx, "gene_id"];
  rm(db.map, q.vec); 
  gc();
  dbDisconnect(con);
  return(entrezs);
}


cleanMem <- function(n=10) { for (i in 1:n) gc() }

doEntrez2SymbolMapping <- function(entrez.vec){
    gene.map <-  queryGeneDB("entrez", data.org);
    gene.map[] <- lapply(gene.map, as.character)

    hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
    symbols <- gene.map[hit.inx, "symbol"];

    # if not gene symbol, use id by itself
    na.inx <- is.na(symbols);
    symbols[na.inx] <- entrez.vec[na.inx];
    return(symbols);
}

#'Mapping from different metabolite IDs
#'@description For compound names to other ids, can do exact or approximate matches
#'For other IDs, except HMDB ID, all others may return multiple/non-unique hits
#'Multiple hits or non-unique hits will allow users to manually select
#'@param qvec Input the metabolites to be mapped.
#'@param q.type Input the query-type, "name" for compound names, "hmdb" for HMDB IDs, "kegg" for KEGG IDs, "pubchem"
#'for PubChem CIDs, "chebi" for ChEBI IDs, "metlin" for METLIN IDs, and "hmdb_kegg" for a both KEGG and HMDB IDs.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
MetaboliteMappingExact <- function(dataSet, qvec, q.type){
  # variables to record results
  hit.inx <- vector(mode='numeric', length=length(qvec)); # record hit index, initial 0
  names(hit.inx) <- qvec;
  match.values <- vector(mode='character', length=length(qvec)); # the best matched values (hit names), initial ""
  match.state <- vector(mode='numeric', length=length(qvec));  # match status - 0, no match; 1, exact match; initial 0 
  
  cmpd.db <- readRDS(paste0(lib.path, "compound_db.rds")); 
  
  if(q.type == "hmdb"){
    n <- 5 # Number of digits for V3 of HMDB
    hmdb.digits <- as.vector(sapply(cmpd.db$hmdb, function(x) strsplit(x, "HMDB")[[1]][2]))
    hmdb.v3.ids <- paste0("HMDB", substr(hmdb.digits, nchar(hmdb.digits)-n+1, nchar(hmdb.digits)))
    hit.inx.v3 <- match(tolower(qvec), tolower(hmdb.v3.ids));
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$hmdb));
    hit.inx[is.na(hit.inx)] <- hit.inx.v3[is.na(hit.inx)]
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "pubchem"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$pubchem));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "chebi"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$chebi));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "metlin"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$metlin));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "kegg"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$kegg));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "name"){
    # first find exact match to the common compound names
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$name));
    match.values <- cmpd.db$name[hit.inx];
    
    match.state[!is.na(hit.inx)] <- 1;
    
    # then try to find exact match to synanyms for the remaining unmatched query names one by one
    syn.db <- readRDS(paste0(lib.path, "syn_nms.rds")); 
    syns.list <-  syn.db$syns.list;
    todo.inx <-which(is.na(hit.inx));
    if(length(todo.inx) > 0){
      for(i in 1:length(syns.list)){
        syns <-  syns.list[[i]];
        hitInx <- match(tolower(qvec[todo.inx]), tolower(syns));
        
        hitPos <- which(!is.na(hitInx));
        if(length(hitPos)>0){
          # record matched ones
          orig.inx<-todo.inx[hitPos];
          hit.inx[orig.inx] <- i;                  
          # match.values[orig.inx] <- syns[hitInx[hitPos]];  # show matched synnames
          match.values[orig.inx] <- cmpd.db$name[i];    # show common name
          match.state[orig.inx] <- 1;
          
          # update unmatched list
          todo.inx<-todo.inx[is.na(hitInx)];
        }
        if(length(todo.inx) == 0) break;
      }
    }
    
  }else{
    print(paste("Unknown compound ID type:", q.type));
    # guess a mix of kegg and hmdb ids
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$hmdb_id));
    hit.inx2 <- match(tolower(qvec), tolower(cmpd.db$kegg_id));
    nohmdbInx <- is.na(hit.inx);
    hit.inx[nohmdbInx]<-hit.inx2[nohmdbInx]
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
    
  }
  match.values[!is.na(match.values)] <- qvec[!is.na(hit.inx)]
  # empty memory
  gc();
  
  dataSet$query.vec <- qvec; 
  dataSet$name.map$hit.inx <- hit.inx;
  dataSet$name.map$hit.values <- match.values;
  dataSet$name.map$match.state <- match.state;
  ids <- cmpd.db$kegg_id[hit.inx];
  ids[is.na(ids) | ids == ""] <- qvec[is.na(ids) | ids == ""]
  nms <- match.values
  nms[is.na(nms)]<- qvec[is.na(nms)]
  dataSet$enrich_ids <- ids;
  names(dataSet$enrich_ids) = nms;
  print(head(dataSet$enrich_ids));
  dataSet$rawToEntrez <- ids
  names(dataSet$rawToEntrez) <- qvec;
  return(dataSet);
}

  bindList <- function(A,B){
    if(length(A)==length(B)){
      res <- list()
      for(i in 1:length(A)){
       res[[i]] <-  rbind(A[[i]],B[[i]])
      }
      names(res) <- names(A)
      return(res)
    }else{
      return("Please check the length of input data")
    }
  
    }

  filtList <- function(A,B){
    if(length(A)==length(B)){
      res <- list()
      for(i in 1:length(A)){
        res[[i]] <-  A[[i]][B[[i]],]
      }
      names(res) <- names(A)
      return(res)
    }else{
      return("Please check the length of input data")
    }
    
  }


DoMetMapping <- function(mvec){
  
  hit.inx <- vector(mode='numeric', length=length(mvec)); # record hit index, initial 0
  names(hit.inx) <- mvec;
  match.values <- vector(mode='character', length=length(mvec)); # the best matched values (hit names), initial ""
  match.state <- vector(mode='numeric', length=length(mvec));  # match status - 0, no match; 1, exact match; initial 0 
  
  cmpd.db <- readRDS(paste0(lib.path, "compound_db.rds")); 
  if(all(grepl("HMDB",mvec))){
    n <- 5 # Number of digits for V3 of HMDB
    hmdb.digits <- as.vector(sapply(cmpd.db$hmdb, function(x) strsplit(x, "HMDB")[[1]][2]))
    hmdb.v3.ids <- paste0("HMDB", substr(hmdb.digits, nchar(hmdb.digits)-n+1, nchar(hmdb.digits)))
    hit.inx.v3 <- match(tolower(qvec), tolower(hmdb.v3.ids));
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$hmdb));
    hit.inx[is.na(hit.inx)] <- hit.inx.v3[is.na(hit.inx)]
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(all(grepl("^C",mvec))){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$kegg));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;

  }else if(all(grepl("^[0-9]",mvec))){
    return("Unknown compound ID type")
  }else{
    hit.inx <- match(tolower(mvec), tolower(cmpd.db$name));
    match.values <- cmpd.db$name[hit.inx];
    
    match.state[!is.na(hit.inx)] <- 1;
    
    # then try to find exact match to synanyms for the remaining unmatched query names one by one
    syn.db <- readRDS(paste0(lib.path, "syn_nms.rds")); 
    syns.list <-  syn.db$syns.list;
    todo.inx <-which(is.na(hit.inx));
    if(length(todo.inx) > 0){
      for(i in 1:length(syns.list)){
        syns <-  syns.list[[i]];
        hitInx <- match(tolower(mvec[todo.inx]), tolower(syns));
        
        hitPos <- which(!is.na(hitInx));
        if(length(hitPos)>0){
          # record matched ones
          orig.inx<-todo.inx[hitPos];
          hit.inx[orig.inx] <- i;                  
          # match.values[orig.inx] <- syns[hitInx[hitPos]];  # show matched synnames
          match.values[orig.inx] <- cmpd.db$name[i];    # show common name
          match.state[orig.inx] <- 1;
          
          # update unmatched list
          todo.inx<-todo.inx[is.na(hitInx)];
        }
        if(length(todo.inx) == 0) break;
      }
    }

  }
  match.values[!is.na(match.values)] <- mvec[!is.na(hit.inx)]
  # empty memory
  gc();

  query.vec <- mvec; 
  kegg.ids <- cmpd.db$kegg_id[hit.inx];
  hmdb.ids <- cmpd.db$hmdb_id[hit.inx]
 # ids[is.na(ids) | ids == ""] <- mvec[is.na(ids) | ids == ""]
  nms <- match.values
#  nms[is.na(nms)]<- mvec[is.na(nms)]
 
  res=data.frame(mvec=mvec,name=match.values,KEGG=kegg.ids,HMDB=hmdb.ids,stringsAsFactors = F)
  
  return(res)
  
}

M2Mscore <- function(qvec,mvec,taxlvl="Genus",dataGem="agora"){
  require('RSQLite'); 

  path <- paste0(sqlite.path, "omicsnet_",dataGem,".sqlite")
  
  m2m.db <- dbConnect(SQLite(), path);
  q.vec <- qvec
  
  query <- paste (shQuote(q.vec),collapse=",");
  
  taxlvl <- tolower(taxlvl)
  statement <- paste("SELECT * FROM ", taxlvl, " WHERE (( ",taxlvl," IN (",query,")))", sep="");
  
  m2m.dic <- .query.sqlite(m2m.db, statement);   
  metInfo <- readRDS(paste0(lib.path, "metInfo.rds")); 
  
  m2m.dic$KEGG <- metInfo$KEGG[match(m2m.dic$metID,metInfo$metID)]
  m2m.dic$HMDB <- metInfo$HMDB[match(m2m.dic$metID,metInfo$metID)]
  rm(metInfo)
  m.vec <- DoMetMapping(mvec)
  
  
  keep.met.idx <- unique(which(m2m.dic$metabolite %in% m.vec$mvec | m2m.dic$metID %in% m.vec$metabolite[!(is.na(m.vec$name))] |
                                 m2m.dic$KEGG %in% m.vec$KEGG[!(is.na(m.vec$KEGG))] |  m2m.dic$HMDB %in% m.vec$HMDB[!(is.na(m.vec$HMDB))] ))
    
 
  m2m.dic <- m2m.dic[keep.met.idx,] 
  m2m.dic$KEGG[which(!(m2m.dic$KEGG %in% m.vec$KEGG))] <- NA
  m2m.dic$HMDB[which(!(m2m.dic$HMDB %in% m.vec$HMDB))] <- NA
  m2m.dic$metabolite[which(!(m2m.dic$metabolite %in% m.vec$mvec) & !(is.na(m2m.dic$KEGG)))] <- m.vec$mvec[match(  m2m.dic$KEGG[which(!(m2m.dic$metabolite %in% m.vec$mvec) & !(is.na(m2m.dic$KEGG)))],m.vec$KEGG)]
  m2m.dic$metabolite[which(!(m2m.dic$metabolite %in% m.vec$mvec) & !(is.na(m2m.dic$HMDB)))] <- m.vec$mvec[match(  m2m.dic$HMDB[which(!(m2m.dic$metabolite %in% m.vec$mvec) & !(is.na(m2m.dic$HMDB)))],m.vec$HMDB)]
  
  m2m.dic <- m2m.dic[,c(1,3,4)]
  names(m2m.dic)[1:2] <- c("from","to")
  
  if(nrow(m2m.dic) == 0){
    return(0)
  }
  
  return(m2m.dic)
  
}



ReadOmicsDataFile <- function(fileName, omics.type=NA) {
  # need to handle reading .csv files too!
  rdtSet <- .get.rdt.set();

  data <- .readDataTable(fileName)
  dataSet <- list();
  
  meta.info <- rdtSet$dataSet$meta.info;

  if(class(data) == "try-error" || ncol(data) == 1){
    AddErrMsg("Data format error. Failed to read in the data!");
    AddErrMsg("Make sure the data table is saved as comma separated values (.csv) format!");
    AddErrMsg("Please also check the followings: ");
    AddErrMsg("Either sample or feature names must in UTF-8 encoding; Latin, Greek letters are not allowed.");
    AddErrMsg("We recommend to use a combination of English letters, underscore, and numbers for naming purpose.");
    AddErrMsg("Make sure sample names and feature (peak, compound) names are unique.");
    AddErrMsg("Missing values should be blank or NA without quote.");
    AddErrMsg("Make sure the file delimeters are commas.");
    return(0);
  }
  
  var.nms <- data[,1];
  data[,1] <- NULL;
  smpl.nms <- colnames(data);
  data <- as.matrix(data);
  rownames(data) <- var.nms;
  
  data <- RemoveDuplicates(data, "mean", quiet=T); # remove duplicates
  data <- as.data.frame(data)
  var.nms <- rownames(data)
  
  msg <- paste("A total of ", ncol(data), " samples and ", nrow(data), " features were found")
  
  # Basic checks - no duplicate samples names
  # Check for uniqueness of sample name
  if(length(unique(smpl.nms))!=length(smpl.nms)){
    dup.nm <- paste(smpl.nms[duplicated(smpl.nms)], collapse=" ");
    AddErrMsg("Duplicate sample names are not allowed!");
    AddErrMsg(dup.nm);
    return(0);
  }
  
  # checking variable names - no duplicate variables for metabolites and microbiome?
  #if(length(unique(var.nms))!=length(var.nms)){
  #  dup.nm <- paste(var.nms[duplicated(var.nms)], collapse=" ");
  #  AddErrMsg("Duplicate feature names are not allowed!");
  #  AddErrMsg(dup.nm);
  #  return(0);
  #}
  
  # now check for special characters in the data labels
  if(sum(is.na(iconv(smpl.nms)))>0){
    na.inx <- is.na(iconv(smpl.nms));
    nms <- paste(smpl.nms[na.inx], collapse="; ");
    AddErrMsg(paste("No special letters (i.e. Latin, Greek) are allowed in sample names!", nms, collapse=" "));
    return(0);
  }
  
  if(sum(is.na(iconv(var.nms)))>0){
    na.inx <- is.na(iconv(var.nms));
    nms <- paste(var.nms[na.inx], collapse="; ");
    AddErrMsg(paste("No special letters (i.e. Latin, Greek) are allowed in feature names!", nms, collapse=" "));
    return(0);
  }
  
  # only keep alphabets, numbers, ",", "." "_", "-" "/"
  smpl.nms <- CleanNames(smpl.nms, "sample_name");
  
  # keep a copy of original names for saving tables 
  orig.var.nms <- var.nms;
  var.nms <- CleanNames(var.nms, "var_name"); # allow space, comma and period
  names(orig.var.nms) <- var.nms;
  
  current.msg <<- msg;
  # now create the dataSet
  dataSet$orig.var.nms <- orig.var.nms;
  data <- data.frame(apply(data, 2, function(x) as.numeric(as.character(x))))
  # now reassgin the dimension names
  colnames(data) <- smpl.nms;
  rownames(data) <- var.nms;
  dataSet$data.proc <- data
  dataSet$data.raw <- data
  dataSet$data.annotated <- ""
  dataSet$data.missed <- ""
  dataSet$data.filtered <- ""
  dataSet$name <- fileName;
  dataSet$de.method <- "NA"
  dataSet$type <- omics.type;
  if(omics.type == "rna_b"){
    readableType <- "Transcriptomics";
  }else if (omics.type == "met_t" || omics.type == "met_u"){
    readableType <- "Metabolomics";
  }else if (omics.type == "mic_m"){
    readableType <- "Microbiome";
  }else if (omics.type == "prot"){
    readableType <- "Proteomics";
  }else if (omics.type == "mirna"){
    readableType <- "miRNA";
  }else{
    readableType <-  omics.type;
  }
  dataSet$readableType <- readableType;
  dataSet$enrich_ids = rownames(dataSet$data.proc)
  names(dataSet$enrich_ids) = rownames(dataSet$data.proc)
  dataSet$meta <- meta.info[which(rownames(meta.info) %in% colnames(dataSet$data.proc)), ];
  # update current dataset
  RegisterData(dataSet);
  return(1)
}

SanityCheckMeta <- function(){

  rdtSet <- .get.rdt.set();
  sel.nms <- names(mdata.all)
  data.list = list();

  
  for(i in 1:length(sel.nms)){
    dataSet <- qs::qread(sel.nms[i])
    data.list[[i]] <- dataSet$meta;
    mdata.all[[i]] <- 1;
  }

  samples_intersect <- intersect_rownames(data.list);
  meta.info <- rdtSet$dataSet$meta.info[samples_intersect,];

  rdtSet$dataSet$meta.info <- meta.info;
  rdtSet$dataSet.origin <- rdtSet$dataSet;
  for(i in 1:length(sel.nms)){
    dataSet <- qs::qread(sel.nms[i])
    dataSet$meta <- rdtSet$dataSet$meta.info;
    dataSet$data.proc <- dataSet$data.proc.origin <- dataSet$data.proc[,samples_intersect];
    RegisterData(dataSet);
  }
  .set.rdt.set(rdtSet)

disc.vec <- paste(names(rdtSet$dataSet$disc.inx)[which(rdtSet$dataSet$disc.inx)],collapse=", ")  
cont.vec <- paste(names(rdtSet$dataSet$cont.inx)[which(rdtSet$dataSet$cont.inx)],collapse=", ")  
na.vec <- na.check(meta.info)
 return(c(ncol(meta.info),length(which(rdtSet$dataSet$disc.inx)),disc.vec,
         length(which(rdtSet$dataSet$cont.inx)),cont.vec,names(meta.info)[1],length(unique(meta.info[,1])),paste(unique(meta.info[,1]),collapse=", "),na.vec ));
}

intersect_rownames <- function(df_list) {
  # Find the intersection of row names across all data frames
  row_names <- Reduce(intersect, lapply(df_list, row.names))
  
  return(row_names)
}


NormalizingDataOmics <-function (data, norm.opt="NA", colNorm="NA", scaleNorm="NA"){
  msg <- ""
  rnms <- rownames(data)
  cnms <- colnames(data)

  # column(sample)-wise normalization
  if(colNorm=="SumNorm"){
    data<-t(apply(data, 2, SumNorm));
    rownm<-"Normalization to constant sum";
  }else if(colNorm=="MedianNorm"){
    data<-t(apply(data, 2, MedianNorm));
    rownm<-"Normalization to sample median";
  }else{
    # nothing to do
    rownm<-"N/A";
  }
  
  if(norm.opt=="log"){
    min.val <- min(data[data>0], na.rm=T)/10;
    numberOfNeg = sum(data<=0, na.rm = TRUE) + 1; 
    totalNumber = length(data)
    data[data<=0] <- min.val;
    data <- log2(data);
    msg <- paste(msg, "Log2 transformation.", collapse=" ");
  }else if(norm.opt=="vsn"){
    require(limma);
    data <- normalizeVSN(data);
    msg <- paste(msg, "VSN normalization.", collapse=" ");
  }else if(norm.opt=="quantile"){
    require('preprocessCore');
    data <- normalize.quantiles(as.matrix(data), copy=TRUE);
    msg <- paste(msg, "Quantile normalization.", collapse=" ");
  }else if(norm.opt=="combined"){
    require(limma);
    data <- normalizeVSN(data);
    require('preprocessCore');
    data <- normalize.quantiles(as.matrix(data), copy=TRUE);
    msg <- paste(msg, "VSN followed by quantile normalization.", collapse=" ");
  }else if(norm.opt=="logcount"){ # for count data, do it in DE analysis, as it is dependent on design matrix
    require(edgeR);
    nf <- calcNormFactors(data);
    y <- voom(data,plot=F,lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    msg <- paste(msg, "Limma based on log2-counts per million transformation.", collapse=" ");
  } else if(norm.opt=="rle"){
    data <- edgeRnorm(data,method="RLE");
    msg <- c(msg, paste("Performed RLE Normalization"));
  }else if(norm.opt=="TMM"){
    data <- edgeRnorm(data,method="TMM");
    msg <- c(msg, paste("Performed TMM Normalization"));
  }else if(norm.opt=="clr"){
    data <- apply(data, 2, clr_transform);
    msg <- "Performed centered-log-ratio normalization.";
  }else if(norm.opt=='LogNorm'){
    min.val <- min(abs(data[data!=0]))/10;
    data<-apply(data, 2, LogNorm, min.val);
  }else if(norm.opt=='CrNorm'){
    norm.data <- abs(data)^(1/3);
    norm.data[data<0] <- - norm.data[data<0];
    data <- norm.data;
  }
  
  
  # scaling
  if(scaleNorm=='MeanCenter'){
    data<-apply(data, 1, MeanCenter);
    scalenm<-"Mean Centering";
  }else if(scaleNorm=='AutoNorm'){
    data<-apply(data, 1, AutoNorm);
    scalenm<-"Autoscaling";
  }else if(scaleNorm=='ParetoNorm'){
    data<-apply(data, 1, ParetoNorm);
    scalenm<-"Pareto Scaling";
  }else if(scaleNorm=='RangeNorm'){
    data<-apply(data, 1, RangeNorm);
    scalenm<-"Range Scaling";
  }else if(scaleNorm=="colsum"){
    data <- sweep(data, 2, colSums(data), FUN="/")
    data <- data*10000000;
    #msg <- c(msg, paste("Performed total sum normalization."));
  }else if(scaleNorm=="upperquartile"){
    suppressMessages(library(edgeR))
    otuUQ <- edgeRnorm(data,method="upperquartile");
    data <- as.matrix(otuUQ$counts);
    #msg <- c(msg, paste("Performed upper quartile normalization"));
  }else if(scaleNorm=="CSS"){
    suppressMessages(library(metagenomeSeq))
    #biom and mothur data also has to be in class(matrix only not in phyloseq:otu_table)
    data1 <- as(data,"matrix");
    dataMR <- newMRexperiment(data1);
    data <- cumNorm(dataMR,p=cumNormStat(dataMR));
    data <- MRcounts(data,norm = T);
    #msg <- c(msg, paste("Performed cumulative sum scaling normalization"));
  }else{
    scalenm<-"N/A";
  }
  if(scaleNorm %in% c('MeanCenter', 'AutoNorm', 'ParetoNorm', 'RangeNorm')){
    data <- t(data)
  }
  
  data <- as.data.frame(data)
  rownames(data) <- rnms;
  colnames(data) <- cnms;
  msg.vec <<- msg;
  return(data)
}

