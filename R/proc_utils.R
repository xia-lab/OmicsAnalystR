
#'Mapping from different gene IDs
#'@param q.vec Input the genes to be mapped.
#'@param org Input the name of the organism.
#'@param q.type Input the query-type.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
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

queryGeneDB <- function(table.nm, data.org){
  require('RSQLite')
  
  conv.db <- dbConnect(SQLite(), paste(sqlite.path, data.org, "_genes.sqlite", sep="")); 
  db.map <- dbReadTable(conv.db, table.nm)
  dbDisconnect(conv.db); cleanMem();
  
  return(db.map)
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
#'License: GNU GPL (>= 2)
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
