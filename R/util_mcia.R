perform_mcia <- function (df.list, cia.nf = 2, cia.scan = FALSE, nsc = T, svd=TRUE) 
{
  library(ade4)
  df.list <- lapply (df.list, function(x) {
    if (inherits(x, "ExpressionSet")) {
      r <- getdata(x)
    } else {
      r <- x
    }
    return(r)
  })

  for (i in names(df.list)) {
    df <- df.list[[i]]
    minn <- min(df)
    ind <- apply(df, 1, function(x) all(x == minn))
    if (any(ind)) 
      stop(paste("There are features in data.frame ", i, 
                 " do not\n        expressed in all observations, please remove these features"))
  }

  N <- sapply(df.list, ncol)
  df.list <- lapply(df.list, as.matrix)
  if (length(unique(N)) != 1) 
    stop("Nonequal number of individual across data.frames")
  infi <- sapply(df.list, function(x) any(is.infinite(x)))
  if (any(infi)) 
    stop("Infinite numeric in the data.frames")
  na <- sapply(df.list, function(x) any(is.na(x)))
  if (any(na)) 
    stop("NAs in the data.frames")
  if (is.null(names(df.list))) 
    names(df.list) <- paste("df", 1:length(df.list), sep = "")
  
  # ====================================
  # ===== lapack function which is called by svd fails to converge in some cases
  # ===== This function is used to replace svd when this happens
  # =================================
  
  mcoaEnv <- environment(mcoa)
  fakeEnv <- new.env(parent = mcoaEnv)
  mcoa2 <- ade4::mcoa
  environment(mcoa2) <- fakeEnv
  
  if (is.logical(svd)) {
    if (svd)
      assign("svd", base::svd, fakeEnv)
    else 
      assign("svd", function(df) {
                res <- list()
                m <- tcrossprod(df, df)
                em <- eigen(m)
                em$values[em$values < 0] <- 1e-30
                res$d <- sqrt(em$values)
                res$u <- em$vectors
                res$v <- t(apply(t(df) %*% em$vectors, 1, function(x) x/sqrt(em$values)))
                return(res)}, 
             fakeEnv)
  } else
    stop("logical value required for svd")
  # =========================================
  pairwise.rv <- function(data.list) {
    rv <- function(m1, m2) {
      nscm1 <- crossprod(as.matrix(m1))
      nscm2 <- crossprod(as.matrix(m2))
      rv <- sum(nscm1*nscm2)/(sum(nscm1*nscm1)*sum(nscm2*nscm2))^0.5
      return(rv)
    }
    n <- length(data.list)
    a <- combn(1:n, 2, FUN=function(x) 
      rv(data.list[[x[1]]], data.list[[x[2]]]), simplify=TRUE)
    m <- matrix(1, n, n)
    m[lower.tri(m)] <- a
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    colnames(m) <- rownames(m) <- names(data.list)
    return(m)
  }
  
  df.list <- lapply(df.list, as.data.frame)
  if (nsc) {
    df.list <- lapply(df.list, array2ade4, pos = TRUE)
    coa.list <- lapply(df.list, dudi.nsc, scannf = FALSE, nf = cia.nf)
    coa.list.t <- lapply(coa.list, ade4:::t.dudi)
    dfl <- lapply(coa.list, function(x) x$tab)
    RV <- pairwise.rv(dfl)
    ktcoa <- ktab.list.dudi(coa.list.t)
  }
  if (!nsc) {
    df.list <- lapply(df.list, t)
    df.list <- lapply(df.list, array2ade4, pos = TRUE)
    coa1 <- dudi.coa(df.list[[1]], scannf = FALSE, nf = cia.nf)
    coa.list <- lapply(df.list[-1], made4:::dudi.rwcoa, rowweights = coa1$lw, 
                       scannf = FALSE, nf = cia.nf)
    coa.list.name <- names(coa.list)
    coa.list$coa1 <- coa1
    coa.list <- coa.list[c("coa1", coa.list.name)]
    names(coa.list)[1] <- names(df.list)[1]
    for (i in 1:length(coa.list)) {
      coa.list[[i]]$lw <- round(coa.list[[i]]$lw, digits = 8)
    }
    dfl <- lapply(coa.list, function(x) x$tab)
    dfl <- lapply(dfl, t) 
    RV <- pairwise.rv(dfl)
    ktcoa <- ktab.list.dudi(coa.list)
  }
  mcoin <- try(mcoa2(X = ktcoa, nf = cia.nf, scannf = FALSE), silent=TRUE) # ...
  if (inherits (mcoin, "try-error")) {
    cat("'svd' fail to convergence, 'eigen' used to perform singular value decomposition")
    assign("svd", function(df) {
              res <- list()
              m <- tcrossprod(df, df)
              em <- eigen(m)
              em$values[em$values < 0] <- 1e-30
              res$d <- sqrt(em$values)
              res$u <- em$vectors
              res$v <- t(apply(t(df) %*% em$vectors, 1, function(x) x/sqrt(em$values)))
              return(res)}, 
            fakeEnv)
    mcoin <- mcoa2(X = ktcoa, nf = cia.nf, scannf = FALSE)
  }
  tab <- scalewt(mcoin$Tco, wt = ktcoa$cw, center = F, scale = T)
  colnames(tab) <- paste("Axis", 1:ncol(tab), sep = "")
  mcoin$Tlw <- ktcoa$lw
  mcoin$Tcw <- ktcoa$cw
  mcoin$blo <- ktcoa$blo
  mcoin$Tc1 <- tab
  mcoin$RV <- RV
  call <- match.call()
  mciares <- list(call = call, mcoa = mcoin, coa = coa.list)
  class(mciares) <- "mcia"
  return(mciares)
}


"getdata" <- function(arraydata) {
        # Edited from vsn function getIntensityMatrix() from W. Huber
        # To run ade4 the arraydata needs to be in a data.frame format.
        y = switch(class(arraydata),
                    matrix    = { if (!is.numeric(arraydata))
                                     stop("Arraydata was found to be a matrix, but is not numeric.")
                                   data.frame(arraydata)
                                 },
                   data.frame = {  if (!all(sapply(arraydata, is.numeric)))
                                     stop("Arraydata was found to be a data.frame, but contains non-numeric columns.")
                                   arraydata
                                 },
                  ExpressionSet = { data.frame(exprs(arraydata))
                                },
                  marrayRaw = {
                               # if (require(affy, quietly = TRUE)) {
                                  nrslides = as.integer(ncol(arraydata@maRf))
                                  nrspots = as.integer(nrow(arraydata@maRf))
                                  tmp = matrix(NA, nrow = nrspots, ncol = 2 * nrslides)
                                  tmp[, (1:nrslides) * 2 - 1] = arraydata@maGf - arraydata@maGb
                                  tmp[, (1:nrslides) * 2] = arraydata@maRf - arraydata@maRb
                                  tmp.names = vector(mode = "character", length = 2 * nrslides)
                                  tmp.names[(1:nrslides) * 2 - 1] = paste("G",colnames(arraydata@maGf),sep="_")
                                  tmp.names[(1:nrslides) * 2] = paste("R",colnames(arraydata@maRf),sep="_")
                                  colnames(tmp) = tmp.names
                                #  }
                                as.data.frame(tmp)
                               },
                   stop(paste("Arraydata has class ", class(arraydata), ". Permitted are: matrix, data.frame, ExpressionSet, marrayRaw", sep=""))
        )  ## end of switch statement
  
        return(y)
}

# This is copied form made4 package in bioconductor
# v1.61

"array2ade4" <-
function(dataset, pos=FALSE,  trans=FALSE){

        # Allows matrix, data.frame, ExpressionSet, marrayRaw to be read as data.frame
        if (!is.data.frame(dataset)) dataset<-getdata(dataset)
 
        if (any(is.na(dataset)))
             stop("Arraydata must not contain NA values. Use impute.knn in library(impute), KNNimpute from Troyanskaya et al., 2001 or LSimpute from Bo et al., 2004 to impute missing values\n")

   
	# COA needs table of positive data, will add real no to make +ve
	if(pos){
               if (any(dataset < 0)) {
                   num<-round(min(dataset)-1)
                   dataset<-dataset+abs(num)
               }
	}

        if(trans) {
               # Transpose matrix  (as BGA, CIA expects the samples to be in the rows)
               # dudi.nsc should not be transposed, use t.dudi instead to ensure row weight are equal
               # There is a horrible bug is dudi.pca/coa etc, if a dataset with vars>>cases is given
               # It can end abruptly crashing the session. This is a bug in sweep
               # There will now use t.dudi rather than transpose the data
              
               # using t convert data.frame to matrix and messes up affymetrix probe ID names
               # It changes all of the "-" to "." in probeids like AFFX-CreX-5_at
               # So save names and change the names of the transposed matrix

               colnam= colnames(dataset)
               rownam = rownames(dataset)               
               dataset<-t(dataset)		
               dimnames(dataset) = list(colnam, rownam) 
               
               # convert matrix to data.frame for ade4
               dataset <- as.data.frame(dataset)
               
               if (!is.data.frame(dataset)) stop("Problems checking dataset")
        }
        
        data.out<-dataset        
        return(data.out)
}
