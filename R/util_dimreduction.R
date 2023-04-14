##################################################
## R script for OmicsAnalyst
## Description: Compute dimension reduction
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Ewald (jessica.ewald@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################

reduce.dimension <- function(reductionOpt, diabloMeta="", diabloPar="0.1"){  
  ncomps = 5;
  sel.nms <- names(mdata.all)[mdata.all==1];
  data.list = list()
  omics.type = vector();
  featureNms <- vector();
  uniqFeats <- vector();
  
  for(i in 1:length(sel.nms)){
    
    dataSet = qs::qread(sel.nms[i])
    omics.type <- c(omics.type, dataSet$type)
    data.list[[dataSet$type]] <- dataSet$data.proc
    
    if(i == 1){       
      comp.res1 = dataSet$comp.res
      enrich.nms1 = dataSet$enrich_ids
      comp.res.inx1 = rep(1, nrow(comp.res1));
      featureNms <- rownames(dataSet$data.proc);
      omics.vec <- rep(dataSet$type, length(featureNms));
      uniqFeats <- paste0(rownames(dataSet$data.proc),"_", dataSet$type)
    } else {
      comp.res1 = rbind(comp.res1, dataSet$comp.res)
      enrich.nms1 = c(enrich.nms1, dataSet$enrich_ids);
      comp.res.inx1 = c(comp.res.inx1, rep(i, nrow(dataSet$comp.res)));
      featureNms <- c(featureNms, rownames(dataSet$data.proc));
      omics.vec <- c(omics.vec,rep(dataSet$type, length(featureNms)));
      uniqFeats <- c(uniqFeats, paste0(rownames(dataSet$data.proc),"_", dataSet$type))
      
    }
  }
  
  reductionSet <- .get.rdt.set();
  reductionSet$comp.res = comp.res1
  reductionSet$enrich_ids = enrich.nms1
  reductionSet$comp.res.inx = comp.res.inx1
  reductionSet$meta = dataSet$meta
  
  if(reductionOpt == "mcia") {
    
    library(omicade4)
    mcoin <- mcia(data.list, cia.nf=ncomps)
    
    pos.xyz = mcoin$mcoa$SynVar;
    colnames(pos.xyz) <- c(paste0("Factor", 1:ncomps))
    
    loading.pos.xyz = mcoin$mcoa$Tco;
    colnames(pos.xyz) <- c(paste0("Factor", 1:ncomps))
    loading.pos.xyz$ids = featureNms;
    
    # get sample and weight names
    names = rownames(pos.xyz)
    
    var.exp <- t(mcoin$mcoa$cov2);
    var.exp <- round(var.exp, digits = 3);
    rownames(var.exp) <- colnames(pos.xyz);
  } else if (reductionOpt == "mofa") {
    tmp_dir <- tempdir();
    do.call(file.remove, list(list.files(tmp_dir, full.names = TRUE, recursive = TRUE)));
    
    library(MOFA2)
    
    # set up model
    data.list <- lapply(data.list, as.matrix)
    for(i in c(1:length(omics.type))){
      rownames(data.list[[i]]) <- paste0(rownames(data.list[[i]]), "_", omics.type[i])
    }
    
    MOFAobject <- create_mofa(data.list);
    data_opts <- get_default_data_options(MOFAobject);
    model_opts <- get_default_model_options(MOFAobject);
    model_opts$num_factors <- ncomps;
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
    
    var.exp <- model@cache[["variance_explained"]][["r2_per_factor"]][[1]]/100;
    var.exp <- round(var.exp, digits = 3);
    
  } else if (reductionOpt == "diablo"){ # pos pars to tune: value from 0-1 inside matrix, which metadata to predict
    library(mixOmics)
    diablo.meta.type <- reductionSet$dataSet$meta.types[diabloMeta]
    diabloPar <- as.numeric(diabloPar);
    
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
    
  }
  
  # preserve original order
  loading.pos.xyz <- loading.pos.xyz[match(featureNms, loading.pos.xyz$ids), ]
  loading.pos.xyz$label <-  invert_named_vector(enrich.nms1)[as.character(loading.pos.xyz$ids)];
  pos.xyz <- pos.xyz[match(rownames(reductionSet$meta), rownames(pos.xyz)), ]
  
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
  
  reductionOptGlobal <<- reductionOpt
  .set.rdt.set(reductionSet);
  
  return(1)
}