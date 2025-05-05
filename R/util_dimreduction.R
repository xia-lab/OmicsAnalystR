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
    sanitized_names <- gsub("[[:cntrl:]]|[^[:ascii:]]", "_", rownames(dataSet$data.proc), perl = TRUE)
    rownames(dataSet$data.proc) <- sanitized_names;
    data.list[[dataSet$type]] <- dataSet$data.proc

    if(i == 1){       
      comp.res1 = dataSet$comp.res
      enrich.nms1 = dataSet$enrich_ids
      comp.res.inx1 = rep(1, nrow(comp.res1));
      featureNms <- rownames(dataSet$data.proc);
      omics.vec <- rep(dataSet$type, nrow(dataSet$data.proc));
      uniqFeats <- paste0(rownames(dataSet$data.proc),"_", dataSet$type)
      filenms <- sel.nms[i]
    } else {
      comp.res1 = rbind(comp.res1, dataSet$comp.res)
      enrich.nms1 = c(enrich.nms1, dataSet$enrich_ids);
      comp.res.inx1 = c(comp.res.inx1, rep(i, nrow(dataSet$comp.res)));
      featureNms <- c(featureNms, rownames(dataSet$data.proc));
      omics.vec <- c(omics.vec,rep(dataSet$type, nrow(dataSet$data.proc)));
      uniqFeats <- c(uniqFeats, paste0(rownames(dataSet$data.proc),"_", dataSet$type))
      filenms <- c(filenms,sel.nms[i])
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
  reductionSet$filenms <- filenms;

  if(reductionOpt == "mcia") {
    
    mcoin <- run.mcia(data.list, cia.nf=ncomps)
    
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
    if(.on.public.web){
        saveSet(infoSet);
        reductionSet$enrich.nms1 <- enrich.nms1;
        qs::qsave(reductionSet, "rdt.set.qs");
        saveRDS(data.list, file = "mofaInput.rds");
        return(2);
    } else {
        if(!exists("run_mofa")){ # public web on same user dir
            compiler::loadcmp("../../rscripts/OmicsAnalystR/R/mofa_core.Rc");   
            compiler::loadcmp("../../rscripts/OmicsAnalystR/R/util_mofa.Rc");    
        }

        #library(MOFA2);
        # set up model
        # Sanitize the row names of each matrix to remove non-ASCII characters and append the omics type.
        data.list <- lapply(data.list, function(matrix, omics) {
            # Replace non-ASCII characters with an underscore or another ASCII character
            sanitized_names <- gsub("[[:cntrl:]]|[^[:ascii:]]", "_", rownames(matrix), perl = TRUE)
            # Append the omics type to the sanitized row names
            rownames(matrix) <- paste0(sanitized_names, "_", omics)
            return(as.matrix(matrix))  # Ensure the sanitized data is still in matrix form
        }, omics.type);  # Use the corresponding omics.type for each matrix

        MOFAobject <- create_mofa_from_matrix(data.list);
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

        model <- run_mofa(MOFAobject, save_data = FALSE, outfile="mofa_model.hdf5");

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

        # flag to be consistent with public
        reductionSet[["mofa.complete"]] <- TRUE;
    }

  } else if (reductionOpt == "diablo"){ # pos pars to tune: value from 0-1 inside matrix, which metadata to predict
    library(mixOmics)
    diablo.meta.type <- reductionSet$dataSet$meta.types[diabloMeta];
    reductionSet$diabloMeta <- diabloMeta;
    reductionSet$diabloPar <- diabloPar;
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
  loading.pos.xyz$filenm <-   filenms
  #update colnames to "Loading"
  colnames(loading.pos.xyz)[c(1:ncomps)] <- c(paste0("Loading", 1:ncomps))

  reductionSet[[reductionOpt]]$pos.xyz <- pos.xyz;
  reductionSet[[reductionOpt]]$loading.pos.xyz <- loading.pos.xyz;
  reductionSet[[reductionOpt]]$var.exp <- var.exp;
  fileNm <- paste0("loading_result_", reductionOpt);
  reductionSet[[reductionSet$reductionOpt]]$loading.file.nm <- fileNm;
  infoSet$imgSet[[reductionOpt]]$loading.pos.xyz <- loading.pos.xyz;
  fast.write.csv(loading.pos.xyz,file=fileNm);
  
  hit.inx <- match(featureNms, unname(enrich.nms1));
  loadingSymbols <- names(enrich.nms1[hit.inx]);
  reductionSet[[reductionOpt]]$loading.enrich <- loadingSymbols
  reductionSet[[reductionOpt]]$loading.names <- featureNms
  reductionSet$omicstype <- names(data.list)

  reductionSet$reductionOpt <- reductionOpt;
  saveSet(infoSet);
  .set.rdt.set(reductionSet);

  return(1)
}


#used to get MOFA results
GetRdtQs <- function(){
    res <- qs::qread("rdt.set.qs");    
#    rdt.set <<- res;
    return(1);
}

run.mcia <- function(df.list, cia.nf = 2, cia.scan = FALSE, nsc = T, svd=TRUE){
    if(!exists("perform_mcia")){ # public web on same user dir
        compiler::loadcmp("../../rscripts/OmicsAnalystR/R/util_mcia.Rc");    
    }
    return(perform_mcia(df.list, cia.nf, cia.scan, nsc, svd));
}

PlotDimredVarexp <- function(imgNm, dpi=72, format="png"){
  infoSet <- readSet(infoSet, "infoSet");
  load_cairo();
  library(see)
  load_ggplot();
  sel.inx <- mdata.all==1;
 
  sel.nms <- names(mdata.all)[sel.inx]
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
 
  reductionSet <- .get.rdt.set();
  df <- reductionSet[[reductionSet$reductionOpt]]$var.exp;

  # reshape deprecated, use data.table
  #df <- reshape2::melt(df) 

  library(data.table);
  df <- as.data.frame(df)
  df$myID <- rownames(df);
  df <- as.data.frame(melt(as.data.table(df), "myID")); 

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
    labs(x="Component #", y="Var. (%)", title="") + theme_minimal(base_size=15) +
    theme(legend.text=element_text(size=16), legend.position = c(0.9, 0.95), legend.title=element_text(size=0));

  
  Cairo(file=imgNm, width=10, height=10, type=format, bg="white", unit="in", dpi=dpi);
  print(p1)
  dev.off();

  infoSet$imgSet[[paste0("dimred_varexp_", reductionSet$reductionOpt)]]<- imgNm;
  saveSet(infoSet);
}

PlotDimredFactors <- function(meta, pc.num = 5, imgNm, dpi=72, format="png"){
  #save.image("factorsDimRed.RData");
  infoSet <- readSet(infoSet, "infoSet");
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
  
  data <- as.data.frame(reductionSet[[reductionSet$reductionOpt]]$pos.xyz[,1:pc.num]);
  meta.info <- reductionSet$meta;
  meta.info <- meta.info[match(rownames(data), rownames(meta.info)),,drop=F]
  
  
  inx <- which(colnames(meta.info) == meta)
  cls <- meta.info[, inx];
  cls.type <- reductionSet$dataSet$meta.types[inx] ##### UPDATE THIS AFTER SUPPORT COMPLEX META
  base_size=15;
  
  if(is.null(cls.type)){
    cls.type <- "disc";
  }

  if (cls.type == "disc"){ ## code to execute if metadata class is discrete
    
    # make plot
    p <- ggpairs(data, 
                 lower = list(continuous = wrap("points")), 
                 upper = list(continuous = wrap("density")),
                 diag = list(continuous = wrap("densityDiag", alpha = 0.5, color = NA)),
                 columnLabels = pclabels, mapping = aes(color = cls))
    
    auxplot <- ggplot(data.frame(cls = cls),aes(x=cls,y=cls,color=cls)) + 
      theme_bw(base_size=base_size) + geom_point(size = 6) + theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=15)) + 
      scale_fill_okabeito() + 
      scale_color_okabeito() +
      guides(col = guide_legend(nrow = 1))
    p <- p + theme_bw(base_size=base_size) + 
      scale_fill_okabeito() + 
      scale_color_okabeito() +
      theme(plot.margin = unit(c(0.25, 0.25, 0.6, 0.25), "in"));

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
      theme_bw(base_size=base_size) + geom_point(size = 6) + 
      theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=15)) + 
      guides(col = guide_legend(nrow = 1))
    
    p <- p + theme_bw(base_size=base_size) + theme(plot.margin = unit(c(0.25, 0.25, 0.8, 0.25), "in"))
    mylegend <- grab_legend(auxplot)
    
  }
  
  grid.newpage()
  grid.draw(p)
  vp = viewport(x=5, y=0.3, width=.35, height=.3, default.units = "in") ## control legend position
  pushViewport(vp)
  grid.draw(mylegend)
  upViewport()
  dev.off()
  
  infoSet$imgSet[[paste0("dimred_factors_", reductionSet$reductionOpt)]]<- imgNm;
  saveSet(infoSet);
}
