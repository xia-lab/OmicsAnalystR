##################################################
## R scripts for NetworkAnalyst 
## Various utility methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################


# new range [a, b]
rescale2NewRange <- function(qvec, a, b){
  qvec = replace(qvec, qvec == 0, 1)
  q.min <- min(qvec);
  q.max <- max(qvec);
  if(length(qvec) < 50){
    a <- a*2;
  }
  if(q.max == q.min){
    new.vec <- rep(8, length(qvec));
  }else{
    coef.a <- (b-a)/(q.max-q.min);
    const.b <- b - coef.a*q.max;
    new.vec <- coef.a*qvec + const.b;
  }
  return(new.vec);
}

GetExtendRange<-function(vec, unit=10){
  var.max <- max(vec);
  var.min <- min(vec);
  exts <- (var.max - var.min)/unit;
  c(var.min-exts, var.max+exts);
}

# given a data with duplicates, dups is the one with duplicates
RemoveDuplicates <- function(data, lvlOpt, quiet=T){
  
  all.nms <- rownames(data);
  colnms <- colnames(data);
  dup.inx <- duplicated(all.nms);
  dim.orig  <- dim(data);
  data <- apply(data, 2, as.numeric); # force to be all numeric
  dim(data) <- dim.orig; # keep dimension (will lost when only one item) 
  rownames(data) <- all.nms;
  colnames(data) <- colnms;
  if(sum(dup.inx) > 0){
    uniq.nms <- all.nms[!dup.inx];
    uniq.data <- data[!dup.inx,,drop=F];
    
    dup.nms <- all.nms[dup.inx];
    uniq.dupnms <- unique(dup.nms);
    uniq.duplen <- length(uniq.dupnms);
    
    for(i in 1:uniq.duplen){
      nm <- uniq.dupnms[i];
      hit.inx.all <- which(all.nms == nm);
      hit.inx.uniq <- which(uniq.nms == nm);
      
      # average the whole sub matrix 
      if(lvlOpt == "mean"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, mean, na.rm=T);
      }else if(lvlOpt == "median"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, median, na.rm=T);
      }else if(lvlOpt == "max"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, max, na.rm=T);
      }else{ # sum
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, sum, na.rm=T);
      }
    }
    if(!quiet){
      if(numOfLists == 1){
        current.msg <<- paste(current.msg, paste("A total of ", sum(dup.inx), " of duplicates were replaced by their ", lvlOpt, ".", sep=""), collapse="\n");
      }else{
        current.msg <<- paste(current.msg, paste0("<b>", listInxU, "</b> : ", length(data), " genes;"), collapse="\n");
      }
    }
    return(uniq.data);
  }else{
    if(!quiet){
      if(numOfLists == 1){
        current.msg <<- paste(current.msg, "All IDs are unique.", collapse="\n");
      }else{
        current.msg <<- paste(current.msg, paste0("<b>", listInxU, "</b> : ", length(data), " genes;"), collapse="\n");
      }
    }
    return(data);
  }
} 


# utils to remove from
# within, leading and trailing spaces
# remove /
ClearFactorStrings<-function(query){
  chars <- substr(query, 0, 1);
  num.inx<- chars >= '0' & chars <= '9';
  if(all(num.inx[!(is.na(num.inx))])){
    query = as.numeric(query);
    query <- factor(query, levels=sort(unique(query)));
  }else{
   query<-factor(query, levels= unique(query))
  }
  return (query);
}

# borrowed from Hmisc
all.numeric <- function (x, what = c("test", "vector"), extras = c(".", "NA")){
  what <- match.arg(what)
  old <- options(warn = -1)
  on.exit(options(old));
  x <- sub("[[:space:]]+$", "", x);
  x <- sub("^[[:space:]]+", "", x);
  inx <- x %in% c("", extras);
  xs <- x[!inx];
  isnum <- !any(is.na(as.numeric(xs)))
  if (what == "test") 
    isnum
  else if (isnum) 
    as.numeric(x)
  else x
}

# Adds an error message
AddErrMsg <- function(msg){
  msg.vec <<- c(msg.vec, msg);
  print(msg);
}

UnzipUploadedFile<-function(zip_file){
    dataName <- unzip(zip_file, list = TRUE)$Name;
    if(length(dataName) > 1){
      # test if "__MACOSX" or ".DS_Store"
      osInx <- grep('MACOSX',dataName,perl=TRUE);
      if(length(osInx) > 0){
        dataName <- dataName[-osInx];
      }
      dsInx <- grep('DS_Store',dataName,perl=TRUE);
      if(length(dsInx) > 0){
        dataName <- dataName[-dsInx];
      }
      if(length(dataName) != 1){
        current.msg <<- "More than one data files found in the zip file.";
        return("NA");
      }
    }
    a <- try(unzip(zip_file));
    if(class(a) == "try-error" | length(a)==0){
      current.msg <<- "Failed to unzip the uploaded files!";
      return ("NA");
    }
    return(dataName);
}

# note, this may leads to duplicate names, use make.unque as last step
CleanNames <- function(query, type){

  if(type=="sample_name"){
    query <- gsub("[^[:alnum:]./_-]", "", query);
  }else{
    query <- gsub("[^[:alnum:][:space:],'./_-]", "", query)
  }
  return(make.unique(query));
}

generate_colors <- function(n_colors, coltype="default", filenm=NULL) {
    if(coltype == "colorblind"){
        palette = c("#9F2CB9", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        if(n_colors <= 0) {
            stop("Number of colors must be a positive integer")
        }
        n_palette <- length(palette)
        if(n_colors <= n_palette) {
            colors <- palette[1:n_colors]
        } else {
            n_repeats <- ceiling(n_colors/n_palette)
            colors <- rep(palette, n_repeats)[1:n_colors]
        }
    }else{
        pal18 <- c("#e6194B", "#3cb44b", "#4363d8", "#ffff00", "#f032e6", "#ffe119", "#911eb4", "#f58231", "#bfef45", "#fabebe", "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075");
        if(n_colors <= 18){ # update color and respect default
          colors <- pal18[1:n_colors];
        }else{
          colors <- colorRampPalette(pal18)(n_colors);
        }
    }
  if(is.null(filenm)){
    return(colors);
  }else{
    library(RJSONIO)
    sink(filenm);
    cat(toJSON(colors));
    sink();
    return(filenm);
  }
}

generate_continuous_colors <- function(n, primary_color="green", filenm=NULL) {
  colors <- colorRampPalette(c("white", primary_color))(n)
  if(is.null(filenm)){
    return(colors);
  }else{
    library(RJSONIO)
    sink(filenm);
    cat(toJSON(colors));
    sink();
    return(filenm);
  }
}

unitAutoScale <- function(df){
    df <- as.data.frame(df)
    row.nms <- rownames(df);
    col.nms <- colnames(df);
    df<-apply(df, 2, AutoNorm);
    rownames(df) <- row.nms;
    colnames(df) <- col.nms;
    maxVal <- max(abs(df))
    df<- df/maxVal
    return(df)
}


# loading mirfamily library accroding to the species. The names for set.ids are the same as set.ids.
LoadmiRFamLib <- function(){
  mirfamily.rda <- paste(lib.path, "mirfamily.rda", sep="");
  load(mirfamily.rda);    
  print(paste("adding library: ", mirfamily.rda));
  current.mset <- mirfam[[dataSet$org]];
  set.ids <- names(current.mset);
  names(set.ids) <- names(current.mset);
  current.setlink <<- "http://www.mirbase.org";
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

doSymbol2EntrezMapping <- function(entrez.vec){
    gene.map <-  queryGeneDB("entrez", data.org);
    gene.map[] <- lapply(gene.map, as.character)

    hit.inx <- match(entrez.vec, gene.map[, "symbol"]);
    symbols <- gene.map[hit.inx, "gene_id"];

    # if not gene symbol, use id by itself
    na.inx <- is.na(symbols);
    symbols[na.inx] <- entrez.vec[na.inx];
    return(symbols);
}

LoadCellLib<-function(){
  motif.path <- paste(lib.path, data.org, "/cell_set.rds", sep="");
  motif_set<-readRDS(motif.path);
  current.mset <- motif_set;
  current.mset = lapply(current.mset, function(x){as.character(x)})
  set.ids<- names(current.mset); 
  names(set.ids) <- names(current.mset)
  current.setlink <<- "";
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

LoadGOLib<-function(onto){
  go.path <- paste(lib.path, data.org, "/go_", tolower(onto), ".rds", sep="");
  if(tolower(onto) %in% c("panthbp","panthcc","panthmf","go_bp")){
    go_bp <- readRDS(go.path);
    
    if(is.null(names(go_bp))){ # new go lib does not give names
      names(go_bp) <- c("link", "term", "sets");
    }
    current.link <- go_bp$link;
    current.mset <- go_bp$sets;
    set.ids<- names(current.mset); 
    names(set.ids) <- names(current.mset) <- go_bp$term;
  }else if(tolower(onto) == "mf"){
    go_mf <- readRDS(go.path);
    if(is.null(names(go_mf))){
      names(go_mf) <- c("link", "term", "sets");
    }
    current.link <- go_mf$link;
    current.mset <- go_mf$sets;
    set.ids<- names(current.mset); 
    names(set.ids) <- names(current.mset) <- go_mf$term;
  }else{
    go_cc <- readRDS(go.path);
    if(is.null(names(go_cc))){
      names(go_cc) <- c("link", "term", "sets");
    }
    current.link <- go_cc$link;
    current.mset <- go_cc$sets;
    set.ids<- names(current.mset); 
    names(set.ids) <- names(current.mset) <- go_cc$term;
  }
  current.setlink <<- current.link;
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

checkEntrezMatches <- function(entrez.vec){
  gene.map <-  queryGeneDB("entrez", data.org);
  gene.map[] <- lapply(gene.map, as.character)
  
  hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
  
  return(length(hit.inx));
}


# loading mirfamily library accroding to the species. The names for set.ids are the same as set.ids.
LoadmiRFamLib <- function(){
  mirfamily.rda <- paste(lib.path, "mirfamily.rda", sep="");

  load(mirfamily.rda);    
  print(paste("adding library: ", mirfamily.rda));
  current.mset <- mirfam[[data.org]];
  set.ids <- names(current.mset);
  names(set.ids) <- names(current.mset);
  current.setlink <<- "http://www.mirbase.org";
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}


# loading miRNA functional annotation library Tam 2.0 (human only)
LoadFuncLib <- function(){
  func.rda <- paste(lib.path, data.org, "/tam_func.rda", sep="");
  load(func.rda);
  print(paste("adding library: ", func.rda));
  current.mset <- tam_func$sets;
  set.ids <- names(current.mset);
  names(set.ids) <- names(current.mset) <- tam_func$term;
  current.setlink <<- "http://www.lirmed.com/tam2/";
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}


# loading miRNA hmdd disease annotation library Tam 2.0 (human only)
LoadHMDDLib <- function(){
  hmdd.rda <- paste(lib.path, data.org, "/tam_hmdd.rda", sep="");
  load(hmdd.rda);
  print(paste("adding library: ", hmdd.rda));
  current.mset <- tam_hmdd$sets;
  set.ids <- names(current.mset);
  names(set.ids) <- names(current.mset) <- tam_hmdd$term;
  current.setlink <<- "http://www.lirmed.com/tam2/";
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

cleanMem <- function(n=8) { for (i in 1:n) gc() }

###########
# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    capture.output(format(utils::object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  print(lapply(dataSet, object.size));
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

# shorthand
ShowMemoryUse <- function(..., n=40) {
  library(pryr);
  sink(); # make sure print to screen
  print(mem_used());
  print(sessionInfo());
  print(.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n));
  print(warnings());
}

color_scale <- function(c1="grey", c2="red") {
  pal <- colorRampPalette(c(c1, c2))
  colors <- pal(100)
  return(colors)
}


# in public web, this is done by microservice
.perform.computing <- function(){
  dat.in <- qs::qread("dat.in.qs"); 
  dat.in$my.res <- dat.in$my.fun();
  qs::qsave(dat.in, file="dat.in.qs");    
}

fast.write <- function(dat, file, row.names=TRUE){
    tryCatch(
        {
           if(is.data.frame(dat)){
                # there is a rare bug in data.table (R 3.6) which kill the R process in some cases 
                data.table::fwrite(dat, file, row.names=row.names);
           }else{
                write.csv(dat, file, row.names=row.names);  
           }
        }, error=function(e){
            print(e);
            fast.write.csv(dat, file, row.names=row.names);   
        }, warning=function(w){
            print(w);
            fast.write.csv(dat, file, row.names=row.names); 
        });
}

rowcolFt =  function(x, fac, var.equal, which = 1L) {
  
  if(!(which %in% c(1L, 2L)))
    stop(sQuote("which"), " must be 1L or 2L.")
  
  if(which==2L)
    x = t(x)

  if (typeof(x) == "integer")
      x[] <- as.numeric(x)

  sqr = function(x) x*x
  
  stopifnot(length(fac)==ncol(x), is.factor(fac), is.matrix(x))
  x   <- x[,!is.na(fac), drop=FALSE]
  fac <- fac[!is.na(fac)]

  ## Number of levels (groups)
  k <- nlevels(fac)

  ## xm: a nrow(x) x nlevels(fac) matrix with the means of each factor
  ## level
  xm <- matrix(
     sapply(levels(fac), function(fl) rowMeans(x[,which(fac==fl), drop=FALSE])),
     nrow = nrow(x),
     ncol = nlevels(fac))

  ## x1: a matrix of group means, with as many rows as x, columns correspond to groups 
  x1 <- xm[,fac, drop=FALSE]

  ## degree of freedom 1
  dff    <- k - 1

  if(var.equal){
    ## x0: a matrix of same size as x with overall means
    x0 <- matrix(rowMeans(x), ncol=ncol(x), nrow=nrow(x))
  
    ## degree of freedom 2
    dfr    <- ncol(x) - dff - 1

    ## mean sum of squares
    mssf   <- rowSums(sqr(x1 - x0)) / dff
    mssr   <- rowSums(sqr( x - x1)) / dfr

    ## F statistic
    fstat  <- mssf/mssr

  } else{

    ## a nrow(x) x nlevels(fac) matrix with the group size  of each factor
    ## level
    ni <- t(matrix(tapply(fac,fac,length),ncol=nrow(x),nrow=k))

    ## wi: a nrow(x) x nlevels(fac) matrix with the variance * group size of each factor
    ## level
    sss <- sqr(x-x1)
    x5 <- matrix(
       sapply(levels(fac), function(fl) rowSums(sss[,which(fac==fl), drop=FALSE])),
       nrow = nrow(sss),
       ncol = nlevels(fac))          
    wi <- ni*(ni-1) /x5

    ## u : Sum of wi
    u  <- rowSums(wi)

    ## F statistic
    MR <- rowSums(sqr((1 - wi/u)) * 1/(ni-1))*1/(sqr(k)-1)
    fsno <- 1/dff * rowSums(sqr(xm - rowSums(wi*xm)/u) * wi)
    fsdeno <- 1+ 2* (k-2)*MR
    fstat <- fsno/fsdeno

    ## degree of freedom 2: Vector with length nrow(x)
    dfr <- 1/(3 * MR)
  
  }
  
  res = data.frame(statistic = fstat,
                   p.value   = pf(fstat, dff, dfr, lower.tail=FALSE),
                   row.names = rownames(x))

  attr(res, "df") = c(dff=dff, dfr=dfr)
  return(res)
}

rowcoltt =  function(x, fac, tstatOnly, which, na.rm) {
    
  dyn.load(.getDynLoadPath());
  
  if (!missing(tstatOnly) && (!is.logical(tstatOnly) || is.na(tstatOnly)))
      stop(sQuote("tstatOnly"), " must be TRUE or FALSE.")
  
  f = checkfac(fac)
  if ((f$nrgrp > 2) || (f$nrgrp <= 0))
    stop("Number of groups is ", f$nrgrp, ", but must be >0 and <=2 for 'rowttests'.")

  if (typeof(x) == "integer")
      x[] <- as.numeric(x)

  cc = .Call("rowcolttests", x, f$fac, f$nrgrp, which-1L, na.rm)
    
  res = data.frame(statistic = cc$statistic,
                   dm        = cc$dm,
                   row.names = dimnames(x)[[which]])

  if (!tstatOnly)
    res = cbind(res, p.value = 2*pt(abs(res$statistic), cc$df, lower.tail=FALSE))

  attr(res, "df") = cc$df    
  return(res)
}

checkfac = function(fac) {

  if(is.numeric(fac)) {
    nrgrp = as.integer(max(fac, na.rm=TRUE)+1)
    fac   = as.integer(fac)
  }
  ## this must precede the factor test
  if(is.character(fac))
    fac = factor(fac)

  if (is.factor(fac)) {
    nrgrp = nlevels(fac)
    fac   = as.integer(as.integer(fac)-1)
  } 
  if(!is.integer(fac))
    stop("'fac' must be factor, character, numeric, or integer.")
  
  if(any(fac<0, na.rm=TRUE))
    stop("'fac' must not be negative.")
    
  return(list(fac=fac, nrgrp=nrgrp))
}

.getDynLoadPath <- function() {

    path = "../../rscripts/OmicsAnalystR/src/OmicsAnalyst.so";
    path=normalizePath(path)
    
    return(path)
}


## fast T-tests/F-tests using genefilter
PerformFastUnivTests <- function(data, cls, var.equal=TRUE, nonpar=F){
    print("Performing fast univariate tests ....2");

    # note, feature in rows for gene expression
    data <- as.matrix(data);
    if(length(levels(cls)) > 2){
        res <- try(rowcolFt(data, cls, var.equal = var.equal));
    }else{
        res <- try(rowcoltt(data, cls, FALSE, 1L, FALSE));
    }  

    if(class(res) == "try-error") {
        res <- cbind(NA, NA);
    }else{
        # res <- cbind(res$statistic, res$p.value);
        # make sure row names are kept
        res <- as.matrix(res[, c("statistic", "p.value")]);
    }

    return(res);
}


fviz_silhouette <- function (sil.obj, titlenm="Silhouette plot",ord=1, label = FALSE, print.summary = FALSE, ...) 
{
    if (inherits(sil.obj, c("eclust", "hcut", "pam", "clara", 
        "fanny"))) {
        df <- as.data.frame(sil.obj$silinfo$widths, stringsAsFactors = TRUE)
    }
    else if (inherits(sil.obj, "silhouette")) 
        df <- as.data.frame(sil.obj[, 1:3], stringsAsFactors = TRUE)
    else stop("Don't support an oject of class ", class(sil.obj))
    df <- df[order(df$cluster, -df$sil_width), ]
    if (!is.null(rownames(df))) 
        df$name <- factor(rownames(df), levels = rownames(df))
    else df$name <- as.factor(1:nrow(df))
    df$cluster <- as.factor(df$cluster)
    mapping <- aes_string(x = "name", y = "sil_width", color = "cluster", 
        fill = "cluster")
    p <- ggplot(df, mapping) + geom_bar(stat = "identity") + 
        labs(y = "Silhouette width Si", x = "", title = paste0(titlenm, 
            "\n Average silhouette width: ", round(mean(df$sil_width), 
                2))) + ggplot2::ylim(c(NA, 1)) + geom_hline(yintercept = mean(df$sil_width), 
        linetype = "dashed", color = "red") + theme_bw()
    p <- ggpubr::ggpar(p, ...)
    if (!label) 
        p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    else if (label) 
        p <- p + theme(axis.text.x = element_text(angle = 45))
    ave <- tapply(df$sil_width, df$cluster, mean)
    n <- tapply(df$cluster, df$cluster, length)
    sil.sum <- data.frame(cluster = names(ave), size = n, ave.sil.width = round(ave, 
        2), stringsAsFactors = TRUE)
    if (print.summary) 
        print(sil.sum)
    p
}

SumNorm<-function(x){
  1000*x/sum(x, na.rm=T);
}

# normalize by median
MedianNorm<-function(x){
  x/median(x, na.rm=T);
}


# normalize to zero mean and unit variance
AutoNorm<-function(x){
  (x - mean(x))/sd(x, na.rm=T);
}

# normalize to zero mean but variance/SE
ParetoNorm<-function(x){
  (x - mean(x))/sqrt(sd(x, na.rm=T));
}

# normalize to zero mean but variance/SE
MeanCenter<-function(x){
  x - mean(x);
}

# normalize to zero mean but variance/SE
RangeNorm<-function(x){
  if(max(x) == min(x)){
    x;
  }else{
    (x - mean(x))/(max(x)-min(x));
  }
}


.discretisation <- function(eigenVectors) {
  
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  
  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
    
  }
  
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}

.discretisationEigenVectorData <- function(eigenVector) {
  
  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1
  
  return(Y)
  
}


# Calculate the mutual information between vectors x and y.
.mutualInformation <- function(x, y) {
  classx <- unique(x)
  classy <- unique(y)
  nx <- length(x)
  ncx <- length(classx)
  ncy <- length(classy)
  
  probxy <- matrix(NA, ncx, ncy)
  for (i in 1:ncx) {
    for (j in 1:ncy) {
      probxy[i, j] <- sum((x == classx[i]) & (y == classy[j])) / nx
    }
  }
  
  probx <- matrix(rowSums(probxy), ncx, ncy)
  proby <- matrix(colSums(probxy), ncx, ncy, byrow=TRUE)
  result <- sum(probxy * log(probxy / (probx * proby), 2), na.rm=TRUE)
  return(result)
}

# Calculate the entropy of vector x.
.entropy <- function(x) {
  class <- unique(x)
  nx <- length(x)
  nc <- length(class)
  
  prob <- rep.int(NA, nc)
  for (i in 1:nc) {
    prob[i] <- sum(x == class[i])/nx
  }
  
  result <- -sum(prob * log(prob, 2))
  return(result)
}

.repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  if (is.null(dim(X))) {
    mx = length(X)
    nx = 1
  } else {
    mx = dim(X)[1]
    nx = dim(X)[2]
  }
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}


# based on phyloseq post: https://github.com/joey711/shiny-phyloseq/blob/master/panels/paneldoc/Transform.md
clr_transform <- function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}


edgeRnorm = function(x,method){

  # Perform edgeR-encoded normalization, using the specified method (...)
  nf <- edgeR::calcNormFactors(x, method=method);
  y <- voom(x,plot=F,lib.size=colSums(x)*nf);
  z <- y$E; # copy per million

  return(z)
}

# generalize log, tolerant to 0 and negative values
LogNorm<-function(x, min.val){
  log10((x + sqrt(x^2 + min.val^2))/2)
}

fast.write.csv <- function(dat, file, row.names=TRUE){
    tryCatch(
        {
           if(is.data.frame(dat)){
                # there is a rare bug in data.table (R 3.6) which kill the R process in some cases 
                data.table::fwrite(dat, file, row.names=row.names);
           }else{
                write.csv(dat, file, row.names=row.names);  
           }
        }, error=function(e){
            print(e);
            write.csv(dat, file, row.names=row.names);   
        }, warning=function(w){
            print(w);
            write.csv(dat, file, row.names=row.names); 
        });
}



saveSet <- function(obj=NA, set="", output=1){
    
    #if(globalConfig$anal.mode == "api"){ 
      qs:::qsave(obj, paste0(set, ".qs"));
    #}else{
    #  if(set == ""){
    #    set <- obj$objName;
    #  }
    #  if(set == "dataSet"){
    #    dataSet <<- obj;
    #  }else if(set == "analSet"){
    #    analSet <<- obj;
    #  }else if(set == "imgSet"){
    #    imgSet <<- obj;
    #  }else if(set == "paramSet"){
    #    paramSet <<- obj;
    #  }else if(set == "msgSet"){
    #    msgSet <<- obj;
    #  }else if(set == "cmdSet"){
    #    cmdSet <<- obj;
    #  }else if(set == "infoSet"){
    #    infoSet <<- obj;
    #  }

    #}
      return(output);

}

readSet <- function(obj=NA, set=""){
    #if(globalConfig$anal.mode == "api"){
     # path <- "";
     # if(exists('user.path')){
     #   path <- user.path;
     # }

     # if(path != ""){
     # obj <- load_qs(paste0(path, set, ".qs"));
     # }else{
      obj <- qs:::qread(paste0(set, ".qs"));
     # }
    #}
    return(obj);
}


#'Replace infinite numbers
#'@description Replace -Inf, Inf to 99999 and -99999
#'@param bdata Input matrix to clean numbers
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'
CleanNumber <-function(bdata){
  if(sum(bdata==Inf)>0){
    inx <- bdata == Inf;
    bdata[inx] <- NA;
    bdata[inx] <- 999999;
  }
  if(sum(bdata==-Inf)>0){
    inx <- bdata == -Inf;
    bdata[inx] <- NA;
    bdata[inx] <- -999999;
  }
  bdata;
}

# get qualified inx with at least number of replicates
GetDiscreteInx <- function(my.dat, min.rep=2){
  good.inx <- apply(my.dat, 2, function(x){
     x = x[x!="NA"]
    good1.inx <- length(x) > length(unique(x));
    good2.inx <- min(table(x)) >= min.rep;
    return (good1.inx & good2.inx);
  });
  return(good.inx);
}

GetNumbericalInx <- function(my.dat){
  good.inx <- apply(my.dat, 2, function(x){
    isNum = as.numeric(as.character(x[x!="NA"]))
    return(all(!is.na(as.numeric(as.character(isNum)))));
  });
  return(good.inx);
}

na.check <- function(mydata){
  na.idx <- apply(mydata,2,function(x) "NA" %in% x)
  if(all(!na.idx)){
    return("None")
  }
  na.num <- apply(mydata,2,function(x) length(which(x=="NA")))
  naInfo <- data.frame(names(mydata)[na.idx],num = na.num[na.num>0])
  naInfo <- apply(naInfo, 1, function(x) paste0(x[1]," (",x[2],")"))
  naInfo <- paste(naInfo,collapse = ", ")
  return(naInfo)
}

# need to obtain the full path to convert (from imagemagik) for cropping images
GetBashFullPath<-function(){
  path <- system("which bash", intern=TRUE);
    
  if((length(path) == 0) && (typeof(path) == "character")){
    print("Could not find bash in the PATH!");
    return("NA");
  }
  return(path);
}


# a utility function to get pheatmap image size (before saving to PNG)
# https://stackoverflow.com/questions/61874876/get-size-of-plot-in-pixels-in-r
get_pheatmap_dims <- function(dat, annotation, view.type, width, cellheight = 15, cellwidth = 15){
  png("NUL"); # trick to avoid open device in server 
  heat_map <- pheatmap::pheatmap(dat, annotation=annotation, cellheight = cellheight, cellwidth = cellwidth);
  h <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"));
  w  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"));
  dev.off();

  # further refine 
  myW <- ncol(dat)*20 + 200;  
  if(myW < 650){
      myW <- 650;
  }   
  myW <- round(myW/72,2);
  if(w < myW){
    w <- myW;
  }

  if(view.type == "overview"){
    if(is.na(width)){
      if(w > 9){
        w <- 9;
      }
    }else if(width == 0){
      if(w > 7.2){
        w <- 7.2;
      }
      
    }else{
      w <- 7.2;
    }
    if(h > w){
      h <- w;
    }
  }

  return(list(height = h, width = w));
}

readDataset <- function(fileName=""){
    if(globalConfig$anal.mode == "api"){
      if(exists('user.path')){
        path <- user.path;
        obj <- load_qs(paste0(path, fileName));
      }else{
        obj <- qs:::qread(fileName);
      }
    }else{
        obj <- qs:::qread(fileName);
       #obj <- dataSets[[fileName]];
    }

    return(obj);
}

saveDataQs <-function(data, name, dataName){
    qs::qsave(data, file=paste0(dataName, "_data/", name));
}

readDataQs <-function(name, dataName){
    dat <- qs::qread(file=paste0(dataName, "_data/", name));
}


makeReadable <- function(str){
    result <- switch(str,
                 pct = "Percent",
                 abs = "Absolute",
                 log = "Log2",
                 rle = "RLE",
                 array = "Microarray",
                 count= "RNA-Seq",
                 hsa = "H. sapiens (human)",
                 mmu = "M. musculus (mouse)",
                 rno = "R. norvegicus (rat)",
                 cel = "C. elegans (roundworm)",
                 dme = "D. melanogaster (fruitfly)",
                 dre = "D. rerio (zebrafish)",
                 sce = "S. cerevisiae (yeast)",
                 eco = "E. coli",
                 ath = "A. thaliana (Arabidopsis)",
                 bta = "B. taurus (cow)",
                 gga = "G. gallus (chicken)",
                 mun = "M. unguiculatus (Mongolian gerbil)",
                 bsu = "B. subtilis",
                 pae = "P. aeruginosa",
                 mtb = "M. tuberculosis",
                 smm = "S. mansoni (schistosomiasis)",
                 tbr = "T. brucei (trypanosoma)",
                 pfa = "P. falciparum (malaria)",
                 cjo = "C. japonica (japanese quail)",
                 xla = "X. laevis (African clawed frog)",
                 ppr = "P. promelas (fathead minnow; custom)",
                 fhm = "P. promelas (fathead minnow; NCBI)",
                 nlf = "L. pipiens (northern leopard frog)",
                 omk = "O. mykiss (rainbow trout)",
                 ham = "H. americanus (American lobster)",
                 cdi = "C. dilutus",
                 dma = "D. magna",
                 rsu = "R. subcapitata",
                 haz = "H. azteca",
                 fcd = "F. candida",
                 "entrez" = "Entrez ID",
                 "refseq" = "RefSeq ID",
                   "gb" = "Genbank ID",
                   "symbol" = "Official Gene Symbol",
                   "embl_gene" = "Ensembl Gene ID",
                   "embl_transcript" = "Ensemble Transcript ID",
                   "embl_protein" = "Ensembl Protein ID",
                   "uniprot" = "Uniprot Accession ID",
                   "hgu95a" = "Affymetrix Human Genome U95 (chip hgu95a)",
                   "hgu95av2" = "Affymetrix Human Genome U95 (chip hgu95av2)",
                   "hgu95b" = "Affymetrix Human Genome U95 (chip hgu95b)",
                   "hgu95c" = "Affymetrix Human Genome U95 (chip hgu95c)",
                   "hgu95d" = "Affymetrix Human Genome U95 (chip hgu95d)",
                   "hgu95e" = "Affymetrix Human Genome U95 (chip hgu95e)",
                   "hgu133a" = "Affymetrix Human Genome U133 (chip hgu133a)",
                   "hgu133b" = "Affymetrix Human Genome U133 (chip hgu133b)",
                   "hgu133plus2" = "Affymetrix Human Genome U133plus2 (hgu133plus2)",
                   "hgu133plus2pm" = "Affymetrix Human Genome U133plus2_PM (hgu133plus2pm)",
                   "lumiht12v3" = "Illumina HumanHT-12 V3 BeadArray",
                   "lumiht12v4" = "Illumina HumanHT-12 V4 BeadArray",
                   "lumiref8v2" = "Illumina HumanRef-8 V2 BeadArray",
                   "lumiref8v3" = "Illumina HumanRef-8 V3 BeadArray",
                   "lumiwg6v2" = "Illumina HumanWG-6 V2 BeadArray",
                   "lumiwg6v3" = "Illumina HumanWG-6 V3 BeadArray",
                   "agi4100a" = "Agilent Human 1 cDNA Microarray (4100A)",
                   "agi4101a" = "Agilent Human 2 cDNA Microarray (4101A)",
                   "agi4110b" = "Agilent Human 1A cDNA Microarray (4110B)",
                   "agi4111a" = "Agilent Human 1B cDNA Microarray (4111A)",
                   "agi4112a" = "Agilent Human Genome Whole Microarray (4x44k/4112)",
                   "agi4845a" = "Agilent Human AMADID 026652 Microarray (4845A)",
                   "lumiwg6v1" = "Illumina MouseWG-6 v1.0 Bead Array",
                   "lumiwg6v11" = "Illumina MouseWG-6 v1.1 Bead Array",
                   "lumiwg6v2" = "Illumina MouseWG-6 v2.0 Bead Array",
                   "lumiref8v1" = "Illumina MouseRef-8 v1.0 Bead Array",
                   "lumiref8v2" = "Illumina MouseRef-8 v2.0 Bead Array",
                   "mgu74a" = "Affymetrix Murine Genome U74v2 (chip mgu74a)",
                   "mgu74av2" = "Affymetrix Murine Genome U74v2 (chip mgu74av2)",
                   "mgu74b" = "Affymetrix Murine Genome U74v2 (chip mgu74b)",
                   "mgu74bv2" = "Affymetrix Murine Genome U74v2 (chip mgu74bv2)",
                   "mgu74c" = "Affymetrix Murine Genome U74v2 (chip mgu74c)",
                   "mgu74cv2" = "Affymetrix Murine Genome U74v2 (chip mgu74cv2)",
                   "moe430a" = "Affymetrix Mouse Expression Set 430 (chip moe430a)",
                   "moe430b" = "Affymetrix Mouse Expression Set 430 (chip moe430b)",
                   "moe430_2" = "Affymetrix GeneChip Mouse Genome 430 2.0",
                   "mgi_st1" = "Affymetrix Mouse Gene 1.0 ST Array",
                   "mgu4101a" = "Agilent Mouse Array (chip mgug4104a)",
                   "mgu4120a" = "Agilent Mouse Array (chip mgug4120a)",
                   "mgu4121a" = "Agilent Mouse Array (chip mgug4121a)",
                   "mgu4122a" = "Agilent Mouse Array (chip mgug4122a)",
                   "kegg" = "KEGG",
                    "reactome" = "Reactome",
                    "go_bp" = "GO:BP",
                    "go_mf" = "GO:MF",
                    "go_cc" = "GO:CC",
                    "panth" = "PANTHER Slim",
                    "motif_set" = "Motif",
                 str)
}
