##################################################
## R script for OmicsAnalyst
## Description: Compute dimension reduction
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Ewald (jessica.ewald@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################

reduce.dimension <- function(reductionOpt, diabloMeta="", diabloPar=0.2){  
  infoSet <- readSet(infoSet, "infoSet");
  ncomps = 5;
  sel.nms <- names(mdata.all)[mdata.all==1];
  data.list = list();
  omics.type = vector();
  featureNms <- vector();
  uniqFeats <- vector();
  for(i in 1:length(sel.nms)){
  
    dataSet = readDataset(sel.nms[i])
    omics.type <- c(omics.type, dataSet$type)
    data.list[[dataSet$type]] <- dataSet$data.proc

    if(i == 1){       
      comp.res1 = dataSet$comp.res
      enrich.nms1 = dataSet$enrich_ids
      comp.res.inx1 = rep(1, nrow(comp.res1));
      featureNms <- rownames(dataSet$data.proc);
      omics.vec <- rep(dataSet$type, nrow(dataSet$data.proc));
      uniqFeats <- paste0(rownames(dataSet$data.proc),"_", dataSet$type)
    } else {
      comp.res1 = rbind(comp.res1, dataSet$comp.res)
      enrich.nms1 = c(enrich.nms1, dataSet$enrich_ids);
      comp.res.inx1 = c(comp.res.inx1, rep(i, nrow(dataSet$comp.res)));
      featureNms <- c(featureNms, rownames(dataSet$data.proc));
      omics.vec <- c(omics.vec,rep(dataSet$type, nrow(dataSet$data.proc)));
      uniqFeats <- c(uniqFeats, paste0(rownames(dataSet$data.proc),"_", dataSet$type))
    }
  }
  reductionSet <- .get.rdt.set();
  reductionSet$comp.res <- comp.res1;
  reductionSet$enrich_ids <- enrich.nms1;
  reductionSet$comp.res.inx <- comp.res.inx1;
  reductionSet$meta <- dataSet$meta;
  reductionSet$uniqFeats <- uniqFeats;
  reductionSet$reductionOpt <- reductionOpt;
  reductionSet$featureNms <- featureNms;
  reductionSet$omics.vec <- omics.vec;

  if(reductionOpt == "mcia") {
    
    mcoin <- mcia(data.list, cia.nf=ncomps)
    
    pos.xyz = mcoin$mcoa$SynVar;
    
    #setting rownames because mcia may modify the names (i.e "-")
    rownames(pos.xyz) <- rownames(reductionSet$meta);
    colnames(pos.xyz) <- c(paste0("Factor", 1:ncomps));
    
    loading.pos.xyz = mcoin$mcoa$Tco;
    loading.pos.xyz$ids = featureNms;
    loading.pos.xyz$type <- omics.vec;
    # get sample and weight names
    names = rownames(pos.xyz)
    
    var.exp <- t(mcoin$mcoa$cov2);
    var.exp <- round(var.exp, digits = 3);
    rownames(var.exp) <- colnames(pos.xyz);
  } else if (reductionOpt == "mofa") {
    # set up model
    data.list <- lapply(data.list, as.matrix)
    for(i in c(1:length(omics.type))){
      rownames(data.list[[i]]) <- paste0(rownames(data.list[[i]]), "_", omics.type[i])
    }

    if(.on.public.web){
        infoSet$paramSet$reductionOptGlobal <- reductionOpt;
        saveSet(infoSet);
        reductionSet$enrich.nms1 <- enrich.nms1;
        qs::qsave(reductionSet, "rdt.set.qs");
        saveRDS(data.list, file = "mofaInput.rds");
        return(2);
    } else {
        library(MOFA2)
        MOFAobject <- create_mofa(data.list);
        data_opts <- get_default_data_options(MOFAobject);
        model_opts <- get_default_model_options(MOFAobject);
        model_opts$num_factors <- 5;
        train_opts <- get_default_training_options(MOFAobject);

        MOFAobject <- prepare_mofa(
          object = MOFAobject,
          data_options = data_opts,
          model_options = model_opts,
          training_options = train_opts
        );

        model <- run_mofa(MOFAobject, save_data = FALSE, use_basilisk = TRUE, outfile="mofa_model.hdf5");

        factors <- get_factors(model, as.data.frame = T);
        pos.xyz <- reshape2::dcast(factors, sample ~ factor, value.var = "value")
        rownames(pos.xyz) <- pos.xyz$sample
        pos.xyz <- pos.xyz[,-1]

        weights <- get_weights(model, as.data.frame = T);
        loading.pos.xyz <- reshape2::dcast(weights, feature ~ factor, value.var = "value")
        loading.pos.xyz$ids <- as.character(loading.pos.xyz$feature)
        loading.pos.xyz <- loading.pos.xyz[,-1]
        loading.pos.xyz$ids <- gsub("_.*", "", loading.pos.xyz$ids)
        loading.pos.xyz$type <- omics.vec;

        var.exp <- model@cache[["variance_explained"]][["r2_per_factor"]][[1]]/100;
        var.exp <- round(var.exp, digits = 3);
    }
    
  } else if (reductionOpt == "diablo"){ # pos pars to tune: value from 0-1 inside matrix, which metadata to predict
    library(mixOmics)
    #print(diabloMeta);
    #print("===diablo");
    diablo.meta.type <- reductionSet$dataSet$meta.types[diabloMeta];
    
    if(diablo.meta.type == "disc"){
      Y <- reductionSet$meta[,diabloMeta];
      
      design = matrix(diabloPar, ncol = length(data.list), nrow = length(data.list), # default diabloPar was 0.2
                      dimnames = list(names(data.list), names(data.list)))
      diag(design) = 0;
      data.list <- lapply(data.list, t)
      model = block.splsda(X = data.list, Y = Y, ncomp = ncomps, design = design)
    } else {
      meta.var <- reductionSet$meta[,diabloMeta];
      Y <- matrix(c(as.numeric(as.character(meta.var))));
      rownames(Y) <- rownames(reductionSet$meta);
      
      design = matrix(diabloPar, ncol = length(data.list), nrow = length(data.list), # default diabloPar was 0.2
                      dimnames = list(names(data.list), names(data.list)));
      diag(design) = 0;
      data.list <- lapply(data.list, t)
      model = block.spls(X = data.list, Y = Y, ncomp = ncomps, design = design, mode = "regression")
    }
    
    # must calculate centroid factor scores
    variates <- model$variates
    variates$Y <- NULL
    variates <- lapply(variates, function(df){
      x_min <- min(df[,1])
      x_max <- max(df[,1])
      y_min <- min(df[,2])
      y_max <- max(df[,2])
      z_min <- min(df[,3])
      z_max <- max(df[,3])
      f4_min <- min(df[,4])
      f4_max <- max(df[,4])
      f5_min <- min(df[,5])
      f5_max <- max(df[,5])
      df[,1] <- (df[,1] - x_min)/ (x_max - x_min) - 0.5
      df[,2] <- (df[,2] - y_min)/ (y_max - y_min) - 0.5
      df[,3] <- (df[,3] - z_min)/ (z_max - z_min) - 0.5
      df[,4] <- (df[,4] - f4_min)/ (f4_max - f4_min) - 0.5
      df[,5] <- (df[,5] - f5_min)/ (f5_max - f5_min) - 0.5
      df
    })
    pos.xyz <- lapply(c(Factor1='comp1', Factor2='comp2', Factor3='comp3', Factor4='comp4', Factor5='comp5'), function(w){
      xORy <- lapply(variates, function(v) v[,w, drop=FALSE])
      xORy <- Reduce(x = xORy, f = cbind)
      xORy <- rowMeans(xORy)
    })
    pos.xyz <- as.data.frame(pos.xyz)
    
    # concatenate feature weights
    loading.pos.xyz <- data.frame()
    for(i in c(1:length(omics.type))){
      temp.mat <- as.data.frame(model[["loadings"]][[i]])
      rnms <- rownames(temp.mat);
      temp.mat <-as.data.frame(unitAutoScale(temp.mat));
      rownames(temp.mat) <- rnms;
      temp.mat$ids <- rownames(temp.mat);
      temp.mat$type <- omics.type[i]
      loading.pos.xyz <- rbind(loading.pos.xyz, temp.mat)
    }
    colnames(loading.pos.xyz) <- c(paste0("Factor", 1:ncomps), "ids", "type");
    var.exp <- model$prop_expl_var;
    var.exp$Y <- NULL;
    var.exp <- as.matrix(as.data.frame(var.exp));
    var.exp <- round(var.exp, digits = 3);
    rownames(var.exp) <- colnames(pos.xyz);
    loading.pos.xyz$type <- omics.vec;
  }
  
  # preserve original order
  loading.pos.xyz <- loading.pos.xyz[match(uniqFeats, paste0(loading.pos.xyz$ids, "_", loading.pos.xyz$type)), ]
  loading.pos.xyz$label <-  invert_named_vector(enrich.nms1)[as.character(loading.pos.xyz$ids)];
  pos.xyz <- pos.xyz[match(rownames(reductionSet$meta), rownames(pos.xyz)), ];

  #update colnames to "Loading"
  colnames(loading.pos.xyz)[c(1:ncomps)] <- c(paste0("Loading", 1:ncomps))

  reductionSet$pos.xyz <- pos.xyz;
  reductionSet$loading.pos.xyz <- loading.pos.xyz;
  reductionSet$var.exp <- var.exp;
  fileNm <- paste0("loading_result_", reductionOpt);
  reductionSet$loading.file.nm <- fileNm;
  fast.write.csv(loading.pos.xyz,file=fileNm);
  
  hit.inx <- match(featureNms, unname(enrich.nms1));
  loadingSymbols <- names(enrich.nms1[hit.inx]);
  reductionSet$loading.enrich <- loadingSymbols
  reductionSet$loading.names <- featureNms
  reductionSet$omicstype <- names(data.list)

  infoSet$paramSet$reductionOptGlobal <- reductionOpt;
  saveSet(infoSet);
  .set.rdt.set(reductionSet);
  
  return(1)
}


#used to get MOFA results
GetRdtQs <- function(){
    res <- qs::qread("rdt.set.qs");
    result.set <<- res;
    return(1);
}

mcia <- function (df.list, cia.nf = 2, cia.scan = FALSE, nsc = T, svd=TRUE) 
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


PlotDimredVarexp <- function(imgNm, dpi=72, format="png"){
  load_cairo();
  library(see)
  load_ggplot();
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx]
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  reductionSet <- .get.rdt.set();
  
  df <- reductionSet$var.exp;
  df <- reshape2::melt(df)
  colnames(df) <- c("Component", "Dataset", "value")
  df$Component <- gsub("Factor","", df$Component);
  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i]);
    df$Dataset <- gsub(dataSet$type,dataSet$readableType, df$Dataset);
  }
  min_r2 = 0
  max_r2 = max(df$value)
  
  p1 <- ggplot(df, aes_string(y="value", x="Component", group="Dataset")) + 
    geom_line(aes(color=Dataset),linewidth=2) +
    scale_fill_okabeito() +
    scale_color_okabeito() +
    labs(x="Component #", y="Var. (%)", title="") + theme_minimal() +
    theme(legend.text=element_text(size=11), legend.position = c(0.9, 0.95), legend.title=element_text(size=0));
  
  
  Cairo(file=imgNm, width=10, height=10, type=format, bg="white", unit="in", dpi=dpi);
  print(p1)
  dev.off();

  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$dimred_varexp <- imgNm;
  saveSet(infoSet);
}

PlotDimredFactors <- function(meta, pc.num = 5, imgNm, dpi=72, format="png"){
  
  load_cairo();
  load_ggplot();
  library(GGally)
  library(see)
  library(grid)
  
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  sel.nms <- names(mdata.all)
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readDataset(sel.nms[i])
    data.list[[i]] <- dat$data.proc
  }
  reductionSet <- .get.rdt.set();
  
  pclabels <- paste0("Component ", 1:pc.num);
  
  # draw plot
  Cairo::Cairo(file = imgNm, unit="in", dpi=dpi, width=10, height=10, type=format, bg="white");
  
  data <- as.data.frame(reductionSet$pos.xyz[,1:pc.num]);
  meta.info <- reductionSet$meta;
  meta.info <- meta.info[match(rownames(data), rownames(meta.info)),,drop=F]
  
  
  inx <- which(colnames(meta.info) == meta)
  cls <- meta.info[, inx];
  cls.type <- reductionSet$dataSet$meta.types[inx] ##### UPDATE THIS AFTER SUPPORT COMPLEX META
  
  if (cls.type == "disc"){ ## code to execute if metadata class is discrete
    
    # make plot
    p <- ggpairs(data, 
                 lower = list(continuous = wrap("points")), 
                 upper = list(continuous = wrap("density")),
                 diag = list(continuous = wrap("densityDiag", alpha = 0.5, color = NA)),
                 columnLabels = pclabels, mapping = aes(color = cls))
    
    auxplot <- ggplot(data.frame(cls = cls),aes(x=cls,y=cls,color=cls)) + 
      theme_bw() + geom_point(size = 6) + theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=11)) + 
      scale_fill_okabeito() + 
      scale_color_okabeito() +
      guides(col = guide_legend(nrow = 1))
    p <- p + theme_bw() + 
      scale_fill_okabeito() + 
      scale_color_okabeito() +
      theme(plot.margin = unit(c(0.25, 0.25, 0.6, 0.25), "in"))
    mylegend <- grab_legend(auxplot)
    
  } else { ## code to excute if metadata class is continuous
    
    colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(20));
    num.cls <- as.numeric(as.character(cls));
    cols <- colors[as.numeric(cut(num.cls,breaks = 20))];
    
    # make plot
    p <- ggpairs(data, lower = list(continuous = wrap("points", color = cols)), 
                 upper = list(continuous = wrap("density", color = "#505050")),
                 diag = list(continuous = wrap("densityDiag", fill = "#505050", color = NA)),
                 columnLabels = pclabels)
    
    auxplot <- ggplot(data.frame(cls = num.cls), aes(x=cls, y=cls, color=cls)) + 
      theme_bw() + geom_point(size = 6) + 
      theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=11)) + 
      guides(col = guide_legend(nrow = 1))
    
    p <- p + theme_bw() + theme(plot.margin = unit(c(0.25, 0.25, 0.8, 0.25), "in"))
    mylegend <- grab_legend(auxplot)
    
  }
  
  grid.newpage()
  grid.draw(p)
  vp = viewport(x=5, y=0.3, width=.35, height=.3, default.units = "in") ## control legend position
  pushViewport(vp)
  grid.draw(mylegend)
  upViewport()
  dev.off()

  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$dimred_factors <- imgNm;
  saveSet(infoSet);
}