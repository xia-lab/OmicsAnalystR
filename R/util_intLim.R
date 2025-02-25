 ##################################################
## R scripts for OmicsAnalyst
## Description: Pairwise linear modeling
## Author: Yao, yao.lu5@mail.mcgill.ca
###################################################

#' IntLim.Anal 
#' @param imgName image name
#' @param imgFormat image format
#' @param analysis.var variable of analysis
#' @param ref reference group
#' @param thresh threshold
#' @param contrast.cls contrast group
#' @export
IntLim.Anal <- function( imgName="NA", imgFormat="png",  
                         analysis.var, ref = NULL, thresh=0.05,pval_type="raw",
                         contrast.cls = "anova",dt1,dt2,outcome=1,topNum=1000){
 
  dataSet1 <- readDataset(dt1);
  dataSet2 <- readDataset(dt2);
  reductionSet <- .get.rdt.set();
 
  if(!exists('adj.vec')){
    adj.bool = F;
    vars <- analysis.var;
    covariates.vec <- c() #for report generation purpose only
  }else{
    if(length(adj.vec) > 0){
      adj.bool = T;
      vars <- c(analysis.var, adj.vec)
      covariates.vec <- adj.vec;
      if(length(adj.vec) == 1){
        if(adj.vec == ""){
          adj.bool = F;
        }
      }
    }else{
      adj.bool = F;
      vars <- analysis.var;
      covariates.vec <- c()
    }
  }
   dt1 = as.matrix(dataSet1$data.proc)[1:min(nrow(dataSet1$data.proc),topNum),]
   dt2 = as.matrix(dataSet2$data.proc)[1:min(nrow(dataSet2$data.proc),topNum),]

  meta.info=reductionSet$dataSet$meta.info
  cls.type <-  reductionSet[["dataSet"]][["meta.types"]]
  idx <- which(cls.type[colnames(meta.info)]=="cont")
  if(length(idx)>0){
    meta.info[,idx] <- apply(meta.info[,idx,drop=F],2,as.numeric)
  }
  if(cls.type[analysis.var]=="disc"){
    meta.info <- meta.info[meta.info[[analysis.var]] %in% c(ref,contrast.cls),,drop=F]
    meta.info[,analysis.var] = droplevels(meta.info[,analysis.var])
    dt1 = dt1[,match(rownames(meta.info),colnames(dt1))]
    dt2 = dt2[,match(rownames(meta.info),colnames(dt2))]
    continuous = F
  }else{
    continuous = T
  }
   if(outcome==1){
    inputData =list(outcome=dt1,independentArray = dt2,sampleMetaData = meta.info)  
    independent.var.type=2
  }else{
    inputData =list(outcome=dt2,independentArray = dt1,sampleMetaData = meta.info)  
    independent.var.type=1
  }
  
   
  myres <- RunIntLim(inputData, stype= analysis.var,  covar=covariates.vec, 
                     continuous=continuous,  save.covar.pvals = F, suppressWarnings = TRUE)
 if(!is.list(myres)){
   return(0)
 }
 
 reductionSet$intLim <- list(stype=analysis.var,continuous=continuous,outcomeArray=dt1,independentArray = dt2,sampleMetaData=meta.info)
  
  reductionSet <- ProcessResults(reductionSet,inputResults=myres,
                             pvalcutoff=thresh, pval.type=pval_type)
 
  .set.rdt.set(reductionSet)
  return(c(1,nrow(reductionSet$intLim_sigmat)))

}


RunIntLim <- function(inputData,stype="", covar=c(), 
                      continuous = FALSE, 
                      save.covar.pvals=FALSE, suppressWarnings = FALSE){
  
  if(!continuous & length(unique(stats::na.omit(inputData$sampleMetaData[,stype]))) != 2) {
    stop(paste("IntLim currently requires only two categories.  Make sure the column",
               stype,"only has two unique values. Did you mean to set",
               "continuous to TRUE?"))
  }

  myres <- NULL
 
  if(length(inputData$outcome)>0 && length(inputData$independentArray)>0) {
    myres <- RunLM(inputData,type =inputData$sampleMetaData[,stype],covar=covar, 
                   continuous = continuous, 
                   save.covar.pvals=save.covar.pvals,
                   suppressWarnings = suppressWarnings)
  }else{
    AddMsg("One type of analyte data is missing.")
    return(0)
  }
  
  return(myres)
}
  

RunLM <- function(inputData,type, covar=c(), 
                  continuous=FALSE, save.covar.pvals = FALSE, keep.highest.pval = FALSE,
                  suppressWarnings = FALSE) {
 
  # Initialize types 1 and 2 and warning list.
  type1 <- inputData$outcome
  type2 <- inputData$independentArray
  sampleMetaData <- inputData$sampleMetaData
  # Initialize covariate matrix, type, and messages.
  covarMatrix <- as.matrix(sampleMetaData[,covar])
  stype <- type 
  mymessage<-list()
  
  if(length(covar)>0){
    # If one of the covariate names contains a plus sign, change it because this
    # will cause problems when constructing the formula.
    adjNames <- RemovePlusInCovars(covar, colnames(sampleMetaData))
    covar <- adjNames$covar
    colnames(sampleMetaData) <- adjNames$sampleDataColnames
    
    # Convert covariates to matrix. Ensure that matrix is one-hot encoded.
    f <- paste('~ 0 + ', paste(covar, collapse = ' + '))
    dat <-  sampleMetaData[,covar,drop=F]
   
    covarMatrix <- stats::model.matrix(stats::as.formula(f), data = dat)
    covar <- colnames(covarMatrix)
    
    # Since names will be changed now, we need to remove plus signs again. For example,
    # if we have a variable 'treatment' that has values 'med1', 'med2', and 'med1+med2',
    # the column names will now include 'treatmentmed1' and 'treatmentmed1+med2'.
    adjNames <- RemovePlusInCovars(covar, colnames(covarMatrix))
    covar <- adjNames$covar
    colnames(covarMatrix) <- adjNames$sampleDataColnames
  }
  
  # Find all standard deviations.
  type1sd <- as.numeric(apply(type1,1,function(x){stats::sd(as.numeric(x),na.rm=TRUE)}))
  type2sd <- as.numeric(apply(type2,1,function(x){stats::sd(as.numeric(x),na.rm=TRUE)}))
  covarsd <- as.numeric(apply(covarMatrix,2,function(x){
    return(stats::sd(x,na.rm=TRUE))}))
  if(methods::is(stype, "character")){
    stype <- as.numeric(as.factor(stype))
  }else if(methods::is(stype, "factor")){
    stype <- as.numeric(stype)
  }
  stypesd <- stats::sd(stype,na.rm=TRUE)
  
  # If the standard deviation of the phenotype is zero, then stop.
  if(stypesd == 0){
    AddMsg("stype variable has a standard deviation of zero. Cannot run.")
    return(0)
  }
  # If the standard deviation of analyte type 1 is zero, then remove and add a warning.
  if (any(type1sd == 0)) {
    toremove <- type1sd == 0  # Logical vector of rows to remove
    namestoremove <- rownames(type1)[toremove]
    type1 <- type1[!toremove, , drop = FALSE]  # Drop removed rows
     
  }
  
  # Check and remove analytes with zero standard deviation for type 2
  if (any(type2sd == 0)) {
    toremove <- type2sd == 0
    namestoremove <- rownames(type2)[toremove]
    type2 <- type2[!toremove, , drop = FALSE]
    
  }
  
  # Check and remove covariates with zero standard deviation
  if (any(covarsd == 0)) {
    toremove <- covarsd == 0
    namestoremove <- colnames(covarMatrix)[toremove]
    covarMatrix <- covarMatrix[, !toremove, drop = FALSE]
    covar <- colnames(covarMatrix)  # Update covariate names
     
  }
  
  myres <- getStatsAllLM(outcome = type1, independentVariable = type2,
                           type = type, covar = covar, covarMatrix = covarMatrix,
                            continuous = continuous, save.covar.pvals = save.covar.pvals,
                             suppressWarnings = suppressWarnings)
 
  
  return(myres)
}


getStatsAllLM <- function(outcome, independentVariable, type, covar, covarMatrix, 
                          continuous, save.covar.pvals, 
                          suppressWarnings = FALSE) {
  outcomeArrayData <- data.frame(outcome)
  independentArrayData <- data.frame(independentVariable)
  num <-  nrow(independentArrayData) 
  # Set up formula and interaction term.
  form.add <- "Y ~ a + type + a:type"
  interactionTerm <- "a:type"
  
  # Set numprog.
  numprog <- round(num*0.1)
  
  # Add covariates to the formula.
  if (length(covar) > 0) {
    form.add <- paste(form.add, paste(covar, collapse = " + "), sep = " + ")
  }
  
  # Initialize stats to collect.
  list.pvals <- list()
  list.coefficients <- list()
  list.rsquared <- list()
  list.covariate.pvals <- list()
  list.covariate.coefficients <- list()
  
  # Run each model.
  warnings <- list()
  for (i in 1:num) {
    # Set up clinical data.
    a <- as.numeric(independentArrayData[i, ])
    if (is.null(covar)|length(covar)==0) {
      clindata <- data.frame(a,type)
    } else {
      clindata <- data.frame(a, type, covarMatrix)
      colnames(clindata)[3:ncol(clindata)] <- covar
    }
    
    # Change type for continuous data (factor to numeric)
    if(continuous){
      clindata[2] <- lapply(clindata[2], as.character)
      clindata[2] <- lapply(clindata[2], as.numeric)
    }
    
    # Run all models for this independent analyte.
    mlin <- getstatsOneLM(stats::as.formula(form.add), clindata = clindata,
                          arraydata = outcomeArrayData, 
                          analytename = rownames(independentArrayData)[i],
                          suppressWarnings = suppressWarnings)
  
    mlin <- mlin[["mlin"]]

    term.pvals <- rownames(mlin$p.value.coeff)
    
    # Return the primary p-values and coefficients.
    index.interac <- grep(interactionTerm, term.pvals)
    term.coefficient <- rownames(mlin$coefficients)
    index.coefficient <-  grep(interactionTerm, term.coefficient)
    p.val.vector <- as.vector(mlin$p.value.coeff[index.interac,])
    coefficient.vector <- as.vector(mlin$coefficients[index.coefficient,])
    
    term.rsquared <- rownames(mlin$r.squared.val)
    r.squared.vector <- as.vector(mlin$r.squared.val)
    
    if (numprog != 0){
      if (i %% numprog == 0) {
        progX <- round(i/num*100)
        message(paste(progX,"% complete"))
      }
    }
    list.pvals[[i]] <-  p.val.vector
    list.coefficients[[i]] <- coefficient.vector
    list.rsquared[[i]] <- r.squared.vector
    
    if(save.covar.pvals == TRUE){
      # Save covariate p-values.
      covariate.pvals <- lapply(term.pvals, function(covariate){
        return(mlin$p.value.coeff[covariate,])
      })
      covariate.pvals.df <- do.call("cbind", covariate.pvals)
      rownames(covariate.pvals.df) <- paste(rownames(independentArrayData)[i], 
                                            rownames(covariate.pvals.df),
                                            sep="__")
      colnames(covariate.pvals.df) <- term.pvals
      list.covariate.pvals[[i]] <- covariate.pvals.df
      
      # Save covariate coefficients.
      covariate.coefficients <- lapply(term.coefficient, function(covariate){
        return(mlin$coefficients[covariate,])
      })
      covariate.coefficients.df <- do.call("cbind", covariate.coefficients)
      rownames(covariate.coefficients.df) <- paste(rownames(independentArrayData)[i],
                                                   rownames(covariate.coefficients.df), 
                                                   sep="__")
      colnames(covariate.coefficients.df) <- term.pvals
      list.covariate.coefficients[[i]] <- covariate.coefficients.df
    }
  }
    failed.calc = which(unlist(lapply(list.coefficients,function(x) length(x)==0)))

list.coefficients[failed.calc] <- list(rep(0,nrow(outcomeArrayData)))
list.pvals[failed.calc] <- list(rep(1,nrow(outcomeArrayData)))
list.rsquared[failed.calc] <- list(rep(0,nrow(outcomeArrayData)))
  # Convert the stats into matrix form.
  mat.pvals <- do.call(rbind, list.pvals)
  mat.coefficients <- do.call(rbind, list.coefficients)
  mat.rsquared <- do.call(rbind, list.rsquared)
  covariate.pvals <- do.call(rbind, list.covariate.pvals)
  covariate.coefficients <- do.call(rbind, list.covariate.coefficients)
  covariate.pvalsadj <- covariate.pvals
  if(!is.null(covariate.pvalsadj)){
    covariate.pvalsadj <- do.call(cbind, lapply(1:ncol(covariate.pvals), function(i){
      return(stats::p.adjust(covariate.pvals[,i], method = 'fdr'))
    }))
    colnames(covariate.pvalsadj) <- colnames(covariate.pvals)
  }
  
  # Adjust p-values.
  row.pvt <- dim(mat.pvals)[1]
  col.pvt <- dim(mat.pvals)[2]
  myps <- as.vector(mat.pvals)
  mypsadj <- stats::p.adjust(myps, method = 'fdr')
  mat.pvalsadj <- matrix(mypsadj, row.pvt, col.pvt)
  
  # Assign names to results.
 
    colnames(mat.pvals) <- colnames(mat.pvalsadj) <-   
      colnames(mat.coefficients) <- colnames(mat.rsquared) <- rownames(outcomeArrayData)
 
    rownames(mat.pvals) <- rownames(mat.pvalsadj) <-  
    rownames(mat.coefficients) <-  
    rownames(mat.rsquared) <- rownames(independentArrayData)
  
  
  # Add the matrices to a final list.
  list.mat <- list()
  list.mat[["interaction.pvalues"]] <- as.matrix(mat.pvals)
  list.mat[["interaction.adj.pvalues"]] <- as.matrix(mat.pvalsadj)
  list.mat[["interaction.coefficients"]] <- as.matrix(mat.coefficients)
  list.mat[["model.rsquared"]] <- as.matrix(mat.rsquared)
  list.mat[["covariate.pvalues"]] <- as.data.frame(covariate.pvals)
  list.mat[["covariate.adj.pvalues"]] <- as.data.frame(covariate.pvalsadj)
  list.mat[["covariate.coefficients"]] <- as.data.frame(covariate.coefficients)
 
  return(list.mat)
}

getstatsOneLM <- function(form, clindata, arraydata, analytename, suppressWarnings = FALSE) {
  
  # Transpose data matrix
  YY <- t(arraydata)
  N <- nrow(YY)
  
  # Precompute means and sums
  EY <- colMeans(YY)
  SYY <- colSums(YY^2) - N * EY^2
  
  # Prepare design matrix
  clindata$Y <- YY[, 1]  # Add dummy Y column for design matrix
  X <- stats::model.matrix(form, clindata)
  XtX <- crossprod(X)
  
  # Solve XtX directly (faster)
  ixtx <- tryCatch({
    solve(XtX)
  }, error = function(e) {
    if (!suppressWarnings) warning(paste("Using pseudoinverse for", analytename))
    MASS::ginv(XtX)
  })
  
  # Coefficients and residuals
  bhat <- ixtx %*% t(X) %*% YY
  yhat <- X %*% bhat
  errors <- YY - yhat
  sse <- colSums(errors^2)
  
  # Degrees of freedom and statistics
  rdf <- ncol(X) - 1
  edf <- N - rdf - 1
  mse <- sse / edf
  var.y <- SYY
  r.squared <- 1 - (sse / var.y)
  
  # Standard errors and p-values
  stderror.coeff <- sqrt(diag(ixtx) %o% mse)
  t.coeff <- bhat / stderror.coeff
  p.val.coeff <- 2 * pt(-abs(t.coeff), df = edf)
  
  list(
    mlin = list(
      coefficients = bhat,
      p.value.coeff = p.val.coeff,
      r.squared.val = r.squared
    ),
    warnings = NULL
  )
}

RemovePlusInCovars <- function(covar=c(), sampleDataColnames){
  # Find which covariates have plus signs.
  which_plus <- which(grepl("+", covar, fixed = TRUE) == TRUE)
  oldCovars <- covar
  
  # Replace each plus sign with "plus".
  covar <- unlist(lapply(1:length(covar), function(i){
    retval <- covar[i]
    if(i %in% which_plus){
      splitOnPlus <- strsplit(covar[i], "+", fixed = TRUE)
      retval <- paste(splitOnPlus[[1]], collapse = "plus")
    }
    return(retval)
  }))
  
  # Replace the plus signs in the sampleMetaData column names as well.
  sampleDataColnames <- unlist(lapply(1:length(sampleDataColnames),
                                      function(i){
                                        retval <- sampleDataColnames[i]
                                        if(retval %in% oldCovars){
                                          whichMatch <- which(oldCovars == retval)
                                          retval <- covar[whichMatch]
                                        }
                                        return(retval)
                                      }))
  # Return new values.
  return(list(covar = covar, sampleDataColnames = sampleDataColnames))
}


ProcessResults <- function(reductionSet,inputResults, 
                           pvalcutoff=0.01,
                           pval.type="raw",
                           coeffPercentile=0,
                           rsquaredCutoff = 0,
                           coefficient = "interaction"){
   
  
    # Store all results in shorter variables.
    mydat <-inputResults$interaction.pvalues
    mydat.adjust <-inputResults$interaction.adj.pvalues
    mydat.interac <- inputResults$interaction.coefficients
    mydat.rsq <- inputResults$model.rsquared 

      # Call continuous function if applicable.
    if( reductionSet$intLim$continuous == 1){
      filtResults <- ProcessResultsContinuous(inputResults,
                                              coeffPercentile,
                                              pvalcutoff, rsquaredCutoff,
                                              coefficient)
      which_int <- which(grepl("a:", colnames(filtResults), fixed = TRUE) == TRUE)
      colnames(filtResults)[which_int] <- "a:type"
      which_type <- which(lapply(colnames(filtResults), function(name){
        retval <- FALSE
        if(substr(name, 1, 4) == "type"){
          retval <- TRUE
        }
        return(retval)
      }) == TRUE)
      colnames(filtResults)[which_type] <- "type"
      
    } else{
       # Create melted matrices to filter by inputs.
      finmydat <- reshape2::melt(mydat)
      finmydat.adj <- reshape2::melt(mydat.adjust)
      finmydat.coef <- reshape2::melt(mydat.interac)
      finmydat.rsq <- reshape2::melt(mydat.rsq)
       finmydat = data.frame(
                Analyte1=as.character(finmydat[,"Var1"]),
                Analyte2 = as.character(finmydat[,"Var2"]),
                interaction_coeff=finmydat.coef[,"value"],
                rsquared =finmydat.rsq[,"value"],
                Pval=finmydat[,"value"],
                FDRadjPval =finmydat.adj[,"value"] 
              )
      fast.write.csv(finmydat, file=paste0("IntLim_raw.csv"))
      
      if(nrow(finmydat)>20000){
        top_idx <- order(abs(finmydat$interaction_coeff), decreasing = TRUE)[1:20000]
        finmydat <- finmydat[top_idx, , drop = FALSE]
      }
    if(pval.type=="raw"){
      finmydat$pLog = -log(finmydat$Pval,10)
    }else{
      finmydat$pLog = -log(finmydat$FDRadjPval,10)
    }
      rownames(finmydat) <- paste(as.character(finmydat$Analyte1), as.character(finmydat$Analyte2), sep = "__")
      
      reductionSet$intLim_filtres <- finmydat 
      
      fast.write.csv(finmydat, file=paste0("IntLim_filtres.csv"))
   
        # Filter by p-value cutoff.
      if(pvalcutoff != 1) { #(no filtering)
        if(pval.type=="raw"){
          keepers2 <- which(finmydat$Pval <= pvalcutoff)
          finmydat <- finmydat[keepers2,] 
        }else{
          keepers2 <- which(finmydat$FDRadjPval <= pvalcutoff)
         finmydat <- finmydat[keepers2,] 
        }
      }
      # Filter by coefficient.
     
      if(coeffPercentile > 0){
        quantiles <- getQuantileForCoefficient(finmydat$Coef, coeffPercentile)
        
        keepers <- (finmydat$Coef > quantiles[2]) | (finmydat$Coef < quantiles[1])
        
        finmydat <- finmydat[keepers, , drop = FALSE] 
      }

     
      # Filter by r-squared.
      if(rsquaredCutoff>0){
        finmydat <- finmydat[finmydat$rsquared>=rsquaredCutoff,]
      }
       
      reductionSet$intLim_sigmat <- finmydat 
      
      fast.write.csv(finmydat, file=paste0("IntLim_sigmat.csv"))
      
      reductionSet$intLim_sigpos = length(which(finmydat$interaction_coeff>0))
      reductionSet$intLim_signeg = length(which(finmydat$interaction_coeff<0))
      reductionSet$intLim_nosig = length(nrow(reductionSet$intLim_filtres)-nrow(finmydat))
     }
     
    # Print and return the results.
    message(paste(nrow(finmydat), 'pairs found given cutoffs'))
    return(reductionSet)
  } 


GetVolcanoDnMat <- function(reductionSet=NA){
    reductionSet <- .get.rdt.set();
    vcn <- reductionSet$intLim_filtres; 
    blue.inx <- which(!rownames(vcn) %in% rownames(reductionSet$intLim_sigmat));

    if(sum(blue.inx)>0){
      xs <- vcn$interaction_coeff[blue.inx]
      ys <- vcn$pLog[blue.inx];   
        return(as.matrix(cbind(xs, ys)));
    }else{
      return(as.matrix(cbind(-1, -1)));
    }
  }
  
  GetVolcanoUpLftMat <- function(reductionSet=NA){
    reductionSet <- .get.rdt.set();
    vcn <- reductionSet$intLim_sigmat; 
    red.inx <- which(vcn$interaction_coeff<0);
    if(sum(red.inx)>0){
      xs <- vcn$interaction_coeff[red.inx]
      ys <- vcn$pLog[red.inx]; 
      return(as.matrix(cbind(xs, ys)));
    }else{
      return(as.matrix(cbind(-1, -1)));
    }
  }
  
  GetVolcanoUpRgtMat <- function(reductionSet=NA){
    reductionSet <- .get.rdt.set();
    vcn <- reductionSet$intLim_sigmat; 
    red.inx <-  which(vcn$interaction_coeff>0);
    if(sum(red.inx)>0){
      xs <- vcn$interaction_coeff[red.inx]
      ys <- vcn$pLog[red.inx]; 
      return(as.matrix(cbind(xs, ys)));
    }else{
      return(as.matrix(cbind(-1, -1)));
    }
  }
  
  GetVolcanoUpLftIDs <- function(mSetObj=NA){
    reductionSet <- .get.rdt.set();
    vcn <- reductionSet$intLim_sigmat; 
    red.inx <-  which(vcn$interaction_coeff<0);
    if(sum(red.inx)>0){
      return(rownames(vcn)[red.inx]);
    }else{
      return("NA");
    }
  }
  
  GetVolcanoUpRgtIDs <- function(reductionSet=NA){
    reductionSet <- .get.rdt.set();
    vcn <- reductionSet$intLim_sigmat; 
    red.inx <-  which(vcn$interaction_coeff>0);
    if(sum(red.inx)>0){
      return(rownames(vcn)[red.inx]);
    }else{
      return("NA");
    }
  }
  
  GetVolcanoDnIDs <- function(mSetObj=NA){
    reductionSet <- .get.rdt.set();
    vcn <- reductionSet$intLim_filtres; 
    blue.inx <- which(!rownames(vcn) %in% rownames(reductionSet$intLim_sigmat));
    if(sum(blue.inx)>0){
      return(rownames(vcn)[blue.inx]);
    }else{
      return("NA");
    }
  }
  
 

PlotPairCorr <- function(reductionSet=NA,imgName,corrID,dpi=72,format="png"){
 
    reductionSet <- .get.rdt.set();
   
  imgName <- paste(imgName, "dpi", dpi, ".", format, sep = "")
   library(ggplot2)
 
    # Set type.
    stype <- reductionSet$intLim$stype
    sampleMetaData <- reductionSet$intLim$sampleMetaData
    
    # Check whether continuous or discrete.
    unique_stypes <- unique(sampleMetaData[,stype])
    
    # For continuous data, plot the marginal effects graph.
    if(length(unique_stypes) > 2){
      MarginalEffectsGraph(
        dataframe = MarginalEffectsGraphDataframe(inputResults = inputResults,
                                                  inputData = inputData,
                                                  outcomeAnalyteOfInterest = outcomeAnalyteOfInterest,
                                                  independentAnalyteOfInterest = independentAnalyteOfInterest,
                                                  outcome = outcome,
                                                  independentVariable = independentVariable), 
        title = paste("Marginal Effects -", independentAnalyteOfInterest,
                      "and", outcomeAnalyteOfInterest), xlab = independentAnalyteOfInterest,
        ylab = outcomeAnalyteOfInterest)
    }else{
       
      independentAnalyteOfInterest = unlist(strsplit(corrID,split = "__"))[1]
      outcomeAnalyteOfInterest = unlist(strsplit(corrID,split = "__"))[2]
 
      outcomeData <-  reductionSet$intLim$outcomeArray
      independentData <- reductionSet$intLim$independentArray
      sOutcome<-as.numeric(outcomeData[rownames(outcomeData)==outcomeAnalyteOfInterest,])
      sIndependent<-as.numeric(independentData[rownames(independentData)==independentAnalyteOfInterest,])
      
        
      data<-data.frame(x=sIndependent,y=sOutcome,sample=colnames(independentData),
                       label=sampleMetaData[,stype])
   
      p<-ggplot(data, aes(x = x, y = y, color = label)) +
        geom_point(size = 3) + # Dots represent the samples
        geom_smooth(method = "lm", se = FALSE) + 
        ggsci::scale_color_aaas(alpha = 0.85) +            
        labs(title = "",
             x = independentAnalyteOfInterest,
             y = outcomeAnalyteOfInterest,
             color = stype) +
        theme_minimal()+
      theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), # Add panel border
             plot.margin = margin(10, 10, 10, 10),
      axis.title = element_text(size = 11, color = "black"),
    axis.text = element_text(size = 10, color = "black")) 
      
    }

   w <- 7.5
     h<-6
     Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, width = w, height = h, type = format, bg = "white")
     print(p)
     dev.off()

  }

#####feature detailed table
GetVolcMat <- function() {
  rdtSet <- .get.rdt.set()
 
  if (is.null(rdtSet$intLim_sigmat)) {
    stop("IntLIM reatult table not found.")
  }
 
  varPart_matrix <- as.matrix(subset(rdtSet$intLim_sigmat, select = -c(Analyte1,Analyte2,pLog))) # Removing the symbol column
  
  return(varPart_matrix)
}

 
GetVolcIds <- function() {
  rdtSet <- .get.rdt.set()
  
  if (is.null(rdtSet$intLim_sigmat)) {
    stop("IntLIM reatult table not found.")
  }
 
  varPart_ids <- rownames(rdtSet$intLim_sigmat)
  
  return(varPart_ids)
}

 
GetVolcSymbols1 <- function() {
  rdtSet <- .get.rdt.set()
  
   if (is.null(rdtSet$intLim_sigmat)) {
    stop("IntLIM reatult table not found.")
  }
   
  varPart_symbols <- rdtSet$intLim_sigmat[,"Analyte1"]
  
  return(varPart_symbols)
}

GetVolcSymbols2 <- function() {
  rdtSet <- .get.rdt.set()
  
   if (is.null(rdtSet$intLim_sigmat)) {
    stop("IntLIM reatult table not found.")
  }
   
  varPart_symbols <- rdtSet$intLim_sigmat[,"Analyte2"]
  
  return(varPart_symbols)
}
 
GetVolcColNames <- function() {
  rdtSet <- .get.rdt.set()
   if (is.null(rdtSet$intLim_sigmat)) {
    stop("IntLIM reatult table not found.")
  }
 
  varPart_colnames <- setdiff(colnames(rdtSet$intLim_sigmat),c("Analyte1","Analyte2","pLog")) # Exclude the symbol column
  
  return(varPart_colnames)
}

 
GetVolcFileName <- function() {
  rdtSet <- .get.rdt.set()
  
   if (is.null(rdtSet$intLim_sigmat)) {
    stop("IntLIM reatult table not found.")
  }
  
  return("IntLim_sigmat.csv")
}

