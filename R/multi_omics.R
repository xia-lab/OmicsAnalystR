

DoIntegrativeAnalysis <- function(method, sign="both", threshold=0.6, nComp){
  load_igraph();
  intRes <- DoDimensionReductionIntegrative(method);
  threshold <- as.numeric(threshold)
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx]
  data <- list()
  labels <- vector();
  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i])
    dat <- dataSet$data.proc
    df <- data.frame(dat, stringsAsFactors = FALSE)
    df <- t(df)
    data[[sel.nms[i]]] <- df
    labels <- c(labels, dataSet$enrich_ids)
  }
  reductionSet<-.get.rdt.set();
  res <- reductionSet$dim.res

  # Isolate mixOmics::network + igraph in subprocess
  net_result <- tryCatch({
    rsclient_isolated_exec(
      func_body = function(input_data) {
        library(mixOmics)
        library(igraph)
        res <- input_data$dim.res
        threshold <- input_data$threshold
        net.res <- mixOmics::network(res, cutoff = threshold, save = "jpeg")
        cor_edge_list <- igraph::as_data_frame(net.res$gR, 'edges')
        gc(verbose = FALSE, full = TRUE)
        return(list(cor_edge_list = cor_edge_list))
      },
      input_data = list(dim.res = res, threshold = threshold),
      packages = c("mixOmics", "igraph", "qs"),
      timeout = 300,
      output_type = "qs"
    )
  }, error = function(e) {
    AddErrMsg(paste("Integrative analysis failed:", e$message))
    NULL
  })
  if (is.list(net_result) && isFALSE(net_result$success)) { AddErrMsg(net_result$message); return(0) }
  if (is.null(net_result)) return(0)

  cor_edge_list <- net_result$cor_edge_list
  if(sign == "both"){
    cor.inx <- abs(cor_edge_list$weight) > threshold
  }else if(sign == "positive"){
    cor.inx <- cor_edge_list$weight > threshold
  }else{
    cor.inx <- cor_edge_list$weight < -threshold
  }
  only_sig <- cor_edge_list[cor.inx, ];
  new_g <- igraph::graph_from_data_frame(only_sig, FALSE);

  type.list <- list()
  for(i in 1:length(sel.nms)){
    type.list[[sel.nms[[i]]]] <- unique(cor_edge_list[,i]);
  }
  ProcessGraphFile(new_g, labels, type.list, TRUE);
}

NormalizeDataWrapper <-function (nm, opt, colNorm){
  if(nm == "NA"){
    sel.nms <- names(mdata.all)
  }else{
    sel.nms <- c(nm);
  }
  for(i in 1:length(sel.nms)){
    dataSet = readDataset(sel.nms[i])
    data <- NormalizingDataOmics(dataSet$data.filtered, opt, colNorm, "NA")
    dataSet$data.proc <- data;
    if(exists("m2m",dataSet)){
      data.norm.taxa <- lapply(dataSet$dataSet$data.filt.taxa, function(x) {
        NormalizingDataOmics(x, opt, colNorm, "NA")
      }
      )
      dataSet$data.proc.taxa <- data.norm.taxa
    }
    fast.write(data, paste0("table_", dataSet$name));
    RegisterData(dataSet)
  }
  return(1)
}

ScaleDataWrapper <-function (nm, scaleNorm){
  if(nm == "NA"){
    sel.nms <- names(mdata.all)
  }else{
    sel.nms <- c(nm);
  }
  for(i in 1:length(sel.nms)){
    dataSet = readDataset(sel.nms[i])
    data <- NormalizingDataOmics(dataSet$data.proc, "NA", "NA", scaleNorm)
    dataSet$data.proc <- data;
    if(exists("m2m",dataSet)){
      data.norm.taxa <- lapply(dataSet$dataSet$data.proc.taxa, function(x) {
        NormalizingDataOmics(x, "NA", "NA", scaleNorm)
      })
      dataSet$data.proc.taxa <- data.norm.taxa
    }
    RegisterData(dataSet)
  }
  return(1);
}

# Per-dataset normalization for multi-omics harmonization. Reads data.proc (the
# harmonization-FILTERED matrix) and writes data.proc, so it composes with the
# variance filter — unlike NormalizeDataWrapper, which reads data.filtered and would
# undo the harmonization filter. `opt`: "auto" detects raw-vs-normalized per layer
# (raw sequencing counts -> log2-CPM via limma::voom; raw intensity/concentration ->
# log2; already-normalized layers left unchanged), or an explicit NormalizingDataOmics
# norm.opt ("logcount","log","NA",...). Records dataSet$normInfo for transparency.
NormalizeDataMultiOmics <- function(nm, opt = "auto"){
  if(nm == "NA"){
    sel.nms <- names(Filter(function(x) isTRUE(x == 1L), mdata.all))
  }else{
    sel.nms <- c(nm);
  }
  if(length(sel.nms) == 0L) return(0);
  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i])
    mat <- dataSet$data.proc
    if(is.null(mat)) next
    use.opt <- opt
    if(identical(opt, "auto")){
      # Per-layer declared data state (set in data-prep) drives the choice; the
      # auto-detector is only the fallback for an UNdeclared layer. "normalized" ->
      # scale only (Mode A); "raw" -> transform by omics type (Mode B); else detect.
      state <- tryCatch(tolower(as.character(dataSet$dataState)), error = function(e) "")
      type  <- tryCatch(dataSet$type, error = function(e) "")
      if(identical(state, "normalized")){
        use.opt <- "NA"
      } else if(identical(state, "raw")){
        use.opt <- .OmicsTypeNormOpt(mat, type)
      } else {
        use.opt <- .OmicsDetectNormOpt(mat, type)
      }
    }
    if(!identical(use.opt, "NA")){
      # log / voom on raw data: clip the few stray negatives (artifacts) to 0 first,
      # otherwise log2 of a negative value returns NaN.
      if(use.opt %in% c("log", "logcount", "clr")){
        mat <- as.matrix(mat); mat[!is.na(mat) & mat < 0] <- 0
      }
      dataSet$data.proc <- NormalizingDataOmics(mat, use.opt, "NA", "NA")
    }
    dataSet$normInfo <- use.opt
    message("[NormalizeDataMultiOmics] ", sel.nms[i], " -> ", use.opt)
    RegisterData(dataSet)
  }
  return(1);
}

# Omics-type -> normalization map for RAW data. Used directly when a layer is DECLARED
# raw (Mode B), and by the auto-detector once it concludes a layer is raw.
#   Microbiome (ASV/OTU) -> CLR (centered log-ratio): feature-level microbiome stays at
#     ASV/OTU here (no taxa collapse, unlike MMP) and is sparse, zero-inflated and
#     compositional; plain log2 leaves the zero-inflation + a mismatched scale that
#     dominates MCIA/DIABLO. CLR centres each sample by its geometric mean.
#   Sequencing (RNA-seq/miRNA) -> voom/log2-CPM, but ONLY for integer count data; data
#     that is not integer counts (TPM/FPKM/CPM, already library-normalized) -> log2.
#   Intensity / concentration / generic -> log2.
.OmicsTypeNormOpt <- function(mat, omicsType = ""){
  t <- tolower(as.character(omicsType))
  if(grepl("mic", t)) return("clr")
  if(grepl("rna|mirna|seq|count|gene", t)){
    m  <- suppressWarnings(matrix(as.numeric(as.matrix(mat)), nrow = nrow(mat)))
    nz <- m[is.finite(m) & m > 0]
    if(length(nz) > 2000L) nz <- nz[seq_len(2000L)]
    if(length(nz) > 0L && mean(nz == round(nz)) > 0.95) return("logcount")
    return("log")
  }
  "log"
}

# Detect a layer's data STATE and pick the normalization to bring it to log-space
# (the auto-scale step runs on every layer afterwards regardless). Returns the
# NormalizingDataOmics opt ("logcount" voom for counts, "clr" for microbiome ASV/OTU
# tables, "log" log2 for intensities) when the layer is RAW, or "NA" (skip) when it is
# already normalized or auto-scaled.
.OmicsDetectNormOpt <- function(mat, omicsType = ""){
  m <- suppressWarnings(matrix(as.numeric(as.matrix(mat)), nrow = nrow(mat)))
  if(all(is.na(m))) return("NA")

  # (1) Already AUTO-SCALED (z-scored): the DEFINING signature is per-feature mean ~ 0
  #     and sd ~ 1. This is unambiguous, so detect it directly and skip.
  rmean <- suppressWarnings(rowMeans(m, na.rm = TRUE))
  rsd   <- suppressWarnings(apply(m, 1, stats::sd, na.rm = TRUE))
  ok    <- is.finite(rmean) & is.finite(rsd) & rsd > 0
  if(sum(ok) >= 5L && stats::median(abs(rmean[ok])) < 0.15 &&
     abs(stats::median(rsd[ok]) - 1) < 0.25) return("NA")

  # (2) Already NORMALIZED + CENTERED (e.g. log-CPM): a SUBSTANTIAL negative fraction.
  #     A few stray negatives (imputation / sparse outliers) in raw data do NOT count --
  #     that wrongly skipped a raw microbiome layer (max ~15000, 0.01% negatives).
  if(suppressWarnings(mean(m < 0, na.rm = TRUE)) > 0.10) return("NA")

  # (3) Small dynamic range (max < 20) => already in log / normalized space => skip.
  #     log2/ln-transformed data never exceeds ~25 even for counts in the millions, so
  #     a large max is the reliable signal that a layer is still RAW.
  mx <- suppressWarnings(max(m, na.rm = TRUE))
  if(!is.finite(mx) || mx < 20) return("NA")

  # (4) max in (20, 50]: small RAW counts vs already-log -> distinguish by integer-ness
  #     (counts are whole numbers, log values are not). max > 50 is unambiguously RAW.
  is.raw <- TRUE
  if(mx <= 50){
    nz <- m[m > 0 & !is.na(m)]
    if(length(nz) > 2000L) nz <- nz[seq_len(2000L)]
    is.raw <- length(nz) > 0L && mean(nz == round(nz)) > 0.95
  }
  if(!is.raw) return("NA")
  .OmicsTypeNormOpt(mat, omicsType)
}

# Generic PREVALENCE filter for count-type data. `data` is features x samples: keep
# features detected (value >= min.count) in at least `min.prev` fraction of samples.
# Sparse features (e.g. microbiome ASVs nonzero in a handful of samples) survive
# abundance/rank filters but are the main cause of the density spike + within-fold zero
# variance in integration — prevalence is the right cut for counts. Never reduces below
# `min.keep` features (returns the input unchanged when the cut would leave too few).
# Defaults match MicrobiomeAnalyst's MMP filter (ApplyAbundanceFilter "prevalence", 4, 0.2):
# >= 4 reads in >= 20% of samples.
FilterByPrevalence <- function(data, min.prev = 0.2, min.count = 4, min.keep = 10L){
  data <- as.matrix(data)
  nS <- ncol(data)
  if(nS < 5L || nrow(data) <= min.keep) return(data)
  prev <- rowSums(data >= min.count, na.rm = TRUE) / nS
  keep <- is.finite(prev) & prev >= min.prev
  if(sum(keep) >= min.keep) data[keep, , drop = FALSE] else data
}

FilterDataMultiOmicsHarmonization <- function(dataName,filterMethod, filterPercent = 0){
  filterPercent <- suppressWarnings(as.numeric(filterPercent));
  if (length(filterPercent) == 0L || is.na(filterPercent)) filterPercent <- 0;
  # "all datasets" = the ACTIVE datasets (mdata.all==1), NOT every entry — a
  # deselected / stale dataset (e.g. manual->AI contamination) must not be
  # re-filtered, and iterating one with mismatched samples drops int.mat to a vector.
  if (length(dataName) == 0L || is.na(dataName) || tolower(as.character(dataName)) %in% c("na","all")) {
    sel.nms <- names(mdata.all)[vapply(mdata.all, function(x) isTRUE(x == 1), logical(1))]
  } else {
    sel.nms <- as.character(dataName);
  }
  if (length(sel.nms) == 0L) { AddErrMsg("No active dataset available to filter."); return(0); }

  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i])
    # STRUCTURAL FIX: filtering and scaling are sequential and share the SAME source —
    # the in-memory data.proc set by the preceding step. The previous code read the
    # on-disk annotated matrix (ov_qs_read(data.annotated.path)), which the AI
    # multi-omics load does not reliably write, so EVERY filter method/dataset failed
    # while scaling (which uses data.proc) worked. Use data.proc, exactly like scaling.
    int.mat <- dataSet$data.proc
    if (is.null(int.mat)) { AddErrMsg(paste0("No data matrix (data.proc) for ", dataSet$name, " — cannot filter.")); return(0); }
    int.mat <- as.matrix(int.mat);
    suppressWarnings(storage.mode(int.mat) <- "numeric");

    # Basic raw-data filtering (before normalization). For COUNT data (microbiome ASV/OTU,
    # sequencing) apply the generic prevalence filter; then drop constant features (zero
    # variance) for every layer — they carry no signal and break DIABLO CV folds.
    .ot <- tryCatch(tolower(as.character(dataSet$type)), error = function(e) "")
    if(grepl("mic|rna|mirna|seq|count|gene", .ot)){
      int.mat <- FilterByPrevalence(int.mat)   # MMP-aligned defaults: >=4 reads in >=20% samples
    }
    .fvar <- suppressWarnings(apply(int.mat, 1, stats::var, na.rm = TRUE))
    .nzv  <- is.finite(.fvar) & .fvar > 0
    if(sum(!.nzv) > 0L && sum(.nzv) >= 2L){
      message("[filter] ", sel.nms[i], " (", .ot, "): removed ", sum(!.nzv), " constant feature(s)")
      int.mat <- int.mat[.nzv, , drop = FALSE]
    }

    if(filterMethod == "variance"){
      # Variance can't rank already-scaled data: z-scoring sets every feature's variance to
      # ~1, so var-of-variances ~ 0. IQR, by contrast, reflects each feature's distribution
      # SHAPE and still varies after scaling — so fall back to IQR ranking there. Non-scaled
      # data is ranked by variance as requested.
      featVar <- apply(int.mat, 1, var, na.rm = TRUE);
      vfv <- suppressWarnings(var(featVar, na.rm = TRUE));
      if(is.na(vfv) || vfv < 0.001){
        message("[filter] near-uniform feature variance (already scaled) -> ranking by IQR instead of variance");
        res  <- tryCatch(PerformFeatureFilter(t(int.mat), "iqr", filterPercent, "", T), error = function(e) NULL);
        data <- if(is.null(res)) int.mat else t(res$data);
      } else {
        data <- FilterDataByVariance(int.mat, filterPercent);
      }
    }else{
      # iqr / other robust ranking — works on already-scaled data too.
      res  <- tryCatch(PerformFeatureFilter(t(int.mat), filterMethod, filterPercent, "", T), error = function(e) NULL);
      data <- if(is.null(res)) int.mat else t(res$data);
    }
    # Defensive: a helper returning a status string instead of a matrix -> keep the input
    # (a layer is passed through, never aborts the whole filter).
    if(any(class(data) == "character")){
      data <- int.mat;
    }
    dataSet$data.proc <- data;
    if(exists("m2m",dataSet)){
      data.norm.taxa <- lapply(dataSet$dataSet$data.annotated.taxa, function(x) {
        if(filterMethod == "variance"){
          data <- FilterDataByVariance(x, filterPercent);
        }else{
          filterPercent <- filterPercent/100;
          res <- PerformFeatureFilter(t(x), filterMethod, filterPercent, "", T);
          data <- t(res$data);
        }
        return(data);
      })
      dataSet$data.proc.taxa <- data.norm.taxa
    }

    RegisterData(dataSet)
  }
  return(1)
}

FilterDataByVariance <- function(data, filterPercent){
  data <- as.matrix(data);
  if (is.null(dim(data)) || nrow(data) < 2L || ncol(data) < 2L) return(data);  # nothing meaningful to filter
  suppressWarnings(storage.mode(data) <- "numeric");
  featVar <- apply(data, 1, var, na.rm = TRUE);
  # Always remove zero-variance features (essential for downstream methods like DIABLO)
  nonzero <- featVar > 0
  if (sum(nonzero) < nrow(data)) {
    message(paste0("Removed ", sum(!nonzero), " zero-variance features"))
    data <- data[nonzero, , drop = FALSE]
    featVar <- featVar[nonzero]
  }
  if(length(featVar) == 0 || suppressWarnings(var(featVar, na.rm = TRUE)) < 0.001 || is.na(suppressWarnings(var(featVar, na.rm = TRUE)))){
    # Near-uniform feature variance (e.g. already z-scored / normalized): variance ranking is
    # meaningless and the quantile cut would drop every feature — pass the data through
    # unchanged so the filter still yields a matrix for the scaling step.
    return(data);
  }
  varThresh <- quantile(featVar, (filterPercent/100), na.rm = TRUE);
  featKeep <- which(featVar > varThresh);
  data <- data[featKeep, , drop = FALSE];
  return(data);
}

#'Plot PCA plot for multi-omics samples
#'@description 
#'@param imgNm name of the image to output
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'

PlotMultiPCA <- function(imgNm, dpi=150, format="png",factor="1", interactive=F){
  try(RecordRCommand(paste0("PlotMultiPCA(\"", imgNm, "\")")), silent = TRUE)
  load_cairo();
  load_ggplot();
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  sel.nms <- names(mdata.all)
  fig.list <- list()
  
  pca.list<- list()
  pct <- list();
  all_data <- list()

  for(i in 1:length(sel.nms)){
    print(i)
    dataSet = readDataset(sel.nms[i])
    x <- dataSet$data.proc
   
    pca <- prcomp(t(na.omit(x)), center=T, scale=T);
    imp.pca<-summary(pca)$importance;
    xlabel <- paste0("PC1"," (", 100*round(imp.pca[2,][1], 3), "%)")
    ylabel <- paste0("PC2"," (", 100*round(imp.pca[2,][2], 3), "%)")
    names <- colnames(x);
    pca.res <- as.data.frame(pca$x);
    pca.res <- pca.res[,c(1,2)]
    
    xlim <- GetExtendRange(pca.res$PC1);
    ylim <- GetExtendRange(pca.res$PC2);
    
    Factor <- dataSet$meta[,1];
    pca.rest <- pca.res
    pca.rest$Conditions <- Factor
    pca.rest$names <- rownames(pca.res)
    pca.rest$dataset <- dataSet$name # This is the new column specifying the dataset/source
    all_data[[i]] <- pca.rest

    #pcafig <- ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions)) +
    #  geom_point(size=3, alpha=0.5) + 
    #  xlim(xlim) + 
    #  ylim(ylim) + 
    #  xlab(xlabel) + 
    #  ylab(ylabel) +
    #  ggtitle(dataSet$name)+
    #  theme_bw()+
    #  theme(text=element_text(size=13), plot.title = element_text(hjust = 0.5));
    
    #fig.list[[i]] <- pcafig
    
    #computing 3d pca
    #pos.xyz = data.frame(x=pca$x[,1], y=pca$x[,2], z=pca$x[,3]);
    #loading.pos.xyz = data.frame(pca$rotation[,c(1:3)]);
    #loadingNames = rownames(pca$rotation);
    
    #pos.xyz <- as.data.frame(pos.xyz);
    #pos.xyz <- unitAutoScale(pos.xyz);
    #rownames(pos.xyz) = rownames(dataSet$meta);
    #dataSet$pos.xyz = pos.xyz;
    
    #loadingNames <- rownames(loading.pos.xyz);
    #loading.pos.xyz <- as.data.frame(loading.pos.xyz);
    #loading <- unitAutoScale(loading.pos.xyz);
    #rownames(loading) <- loadingNames;
    #nm <- paste0("pca_", dataSet$type);
    #pca.list[[nm]][["score"]] <- pos.xyz * 1000;
    #pca.list[[nm]][["loading"]] <- loading* 1000;    
    
    #pct2Nm <- paste0("pca_", dataSet$type)
    #pct[[pct2Nm]] <- unname(round(imp.pca[2,],3))[c(1:3)]*100;
  }

  pca.list$pct2 <- pct;
  ov_qs_save(pca.list, file="pca.scatter.qs");
  
  #h<-6*round(length(sel.nms)/2)
  h<-6
  #library("ggpubr")
  #p1 <- ggarrange(plotlist=fig.list, ncol = 2, nrow = round(length(fig.list)/2))
  combined_data <- do.call(rbind, all_data)

  # Method-standard: persist the per-sample PCA coordinates behind this figure so
  # Refine can re-plot from data and the scores are portable to any tool.
  if (exists("WfSaveFigureData"))
    tryCatch(WfSaveFigureData("oa_multi_pca", combined_data), error = function(e) NULL)

  p1 <- ggplot(combined_data, aes(x=PC1, y=PC2, color=Conditions)) +
    geom_point(size=3, alpha=0.5) + 
    facet_wrap(~ dataset, scales = "free") + # Use facet_wrap or facet_grid
    theme_bw() +
    theme(text=element_text(size=13))


  if(interactive){
    library(plotly);
        m <- list(
                l = 50,
                r = 50,
                b = 20,
                t = 20,
                pad = 0.5
            )

    ggp_build <- layout(ggplotly(p1), autosize = FALSE, width = 1000, height = 600, margin = m)
    return(ggp_build);
  }else{
  Cairo(file=imgNm, width=14, height=h, type=format, bg="white", unit="in", dpi=dpi);
  print(p1)
  dev.off();
  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$qc_multi_pca <- imgNm;
  saveSet(infoSet);
  }
}

PlotMultiDensity <- function(imgNm, dpi=150, format="png",factor="1", interactive=F){
  try(RecordRCommand(paste0("PlotMultiDensity(\"", imgNm, "\")")), silent = TRUE)
  load_cairo();
  load_ggplot();
  dpi <- as.numeric(dpi)
  factor <- as.numeric(factor)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  sel.nms <- names(mdata.all)
  fig.list <- list()
  # Pre-allocate list to avoid sequential rbind (memory optimization)
  df.list <- vector("list", length(sel.nms))

  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i])
    dat <- dataSet$data.proc
    st <- stack(as.data.frame(dat))
    st$Dataset <- rep(sel.nms[i], nrow(dat))
    df.list[[i]] <- st  # Store in list instead of sequential rbind
  }

  # Combine all data frames at once (instead of sequential rbind)
  merged.df <- do.call(rbind, df.list)
  
  type<-merged.df$type
  merged.df$ind <- paste0(merged.df$ind, "_", merged.df$type)
  # Method-standard: persist the plotted density data so the AI "Refine" control can
  # re-plot from data and users can regenerate it in any tool (key on the figure stem,
  # since this fn draws several stable-named densities: raw / norm / scaled / qc).
  if (exists("WfSaveFigureData"))
    tryCatch(WfSaveFigureData(sub("_dpi.*$", "", basename(imgNm)), merged.df), error = function(e) NULL)
    g =ggplot(merged.df, aes(x=values)) +
    geom_line(aes(color=Dataset, group=ind), stat="density", alpha=0.1) + 
    geom_line(aes(color=Dataset), stat="density", alpha=0.7, size=3) +
    theme_bw() +
    theme(text=element_text(size=13))



  if(interactive){
    library(plotly);
        m <- list(
                l = 50,
                r = 50,
                b = 20,
                t = 20,
                pad = 0.5
            )

    ggp_build <- layout(ggplotly(g), autosize = FALSE, width = 1000, height = 600, margin = m)
    return(ggp_build);
  }else{
  Cairo(file=imgNm, width=10, height=6, type=format, bg="white", dpi=dpi, unit="in");
  print(g)
  dev.off();
  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$qc_multi_density <- imgNm;
  saveSet(infoSet);
  }
} 

CheckMetaIntegrity <- function(){
  try(RecordRCommand("CheckMetaIntegrity()"), silent = TRUE)
  sel.nms <- names(mdata.all)

  data.list <- list()
  cnms <- list()
  metas <- list();
  for(i in 1:length(sel.nms)){
    dat = readDataset(sel.nms[i])
    cnms[[i]] <- colnames(dat$data.proc);
    metas[[i]] <- dat$meta;
  }
  
  if(length(metas) == 0){
    AddErrMsg('No metadata found. Provide a metadata file, or rows starting with "#CLASS".');
    return(0)
  }

  for(i in 1:length(sel.nms)){
    for(j in 1:length(sel.nms)){
      bool <- SameElements(cnms[[i]], cnms[[j]])
      if(!bool){
        ov <- length(intersect(cnms[[i]], cnms[[j]]))
        AddErrMsg(paste0(
          "Sample IDs are not consistent across the omics tables — every table must use the SAME ",
          "sample IDs as its column headers, matching the metadata row names. '", sel.nms[i], "' has ",
          length(cnms[[i]]), " samples and '", sel.nms[j], "' has ", length(cnms[[j]]),
          "; they share only ", ov, ". '", sel.nms[i], "' e.g.: ",
          paste(utils::head(cnms[[i]], 4), collapse=", "), " | '", sel.nms[j], "' e.g.: ",
          paste(utils::head(cnms[[j]], 4), collapse=", "), "."));
        return(0)
      }
      boolMeta <- identical(metas[[i]],metas[[j]])
      if(!boolMeta){
        AddErrMsg("The metadata is not consistent across the uploaded omics tables (the same samples must map to the same metadata in every table).");
        return(0)
      }
    }
  }
  return(1)
  
}

SameElements <- function(a, b) return(identical(sort(a), sort(b)));

#'Plot t-sne plot for multi-omics samples
#'@description 
#'@param imgNm name of the image to output
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'

PlotMultiTsne <- function(imgNm, dpi=150, format="png",factor="1"){
  dpi<-as.numeric(dpi)

  sel.nms <- names(mdata.all)[mdata.all == 1]
  data_list <- list()
  meta_list <- list()
  for (i in 1:length(sel.nms)) {
    dataSet <- readDataset(sel.nms[i])
    data_list[[sel.nms[i]]] <- dataSet$data.proc
    meta_list[[sel.nms[i]]] <- dataSet$meta
  }

  # Isolate ggpubr + Rtsne in subprocess
  params <- list(
    imgNm = imgNm, dpi = dpi, format = format, factor = factor, sel.nms = sel.nms
  )

  tsne_result <- tryCatch({
    rsclient_isolated_exec(
      func_body = function(input_data) {
        library(ggpubr)
        library(ggplot2)
        library(Rtsne)
        library(Cairo)

        data_list <- input_data$data_obj$data_list
        meta_list <- input_data$data_obj$meta_list
        params <- input_data$params

        imgNm <- paste(params$imgNm, "dpi", params$dpi, ".", params$format, sep = "")
        sel.nms <- params$sel.nms
        fig.list <- list()

        for (i in 1:length(sel.nms)) {
          x <- t(data_list[[sel.nms[i]]])
          max.perx <- floor((nrow(x) - 1) / 3)
          if (max.perx > 30) max.perx <- 30
          tsne_out <- Rtsne::Rtsne(x, pca = FALSE, perplexity = max.perx, theta = 0.0, check_duplicates = FALSE)
          pca.res <- as.data.frame(tsne_out$Y)
          colnames(pca.res) <- c("tsne1", "tsne2")
          xlim <- range(pca.res[, 1])
          xlim <- xlim + c(-1, 1) * diff(xlim) * 0.1
          ylim <- range(pca.res[, 2])
          ylim <- ylim + c(-1, 1) * diff(ylim) * 0.1
          Factor <- meta_list[[sel.nms[i]]][, 1]
          pca.rest <- pca.res
          pca.rest$Conditions <- Factor
          pca.rest$names <- rownames(pca.res)
          pcafig <- ggplot2::ggplot(pca.rest, ggplot2::aes(x = tsne1, y = tsne2, color = Conditions)) +
            ggplot2::geom_point(size = 3, alpha = 0.5) +
            ggplot2::xlim(xlim) + ggplot2::ylim(ylim) +
            ggplot2::xlab("tsne1") + ggplot2::ylab("tsne2") +
            ggplot2::theme_bw()
          fig.list[[i]] <- pcafig
        }

        h <- 6 * round(length(fig.list) / 2)
        Cairo::Cairo(file = imgNm, width = 14, height = h, type = params$format, bg = "white", unit = "in", dpi = params$dpi)
        p1 <- ggpubr::ggarrange(plotlist = fig.list, ncol = 2, nrow = round(length(fig.list) / 2), labels = sel.nms)
        print(p1)
        dev.off()
        gc(verbose = FALSE, full = TRUE)
        return(1)
      },
      input_data = list(
        data_obj = list(data_list = data_list, meta_list = meta_list),
        params = params
      ),
      packages = c("ggpubr", "ggplot2", "Rtsne", "Cairo", "qs"),
      timeout = 300,
      output_type = "qs"
    )
  }, error = function(e) {
    AddErrMsg(paste("PlotMultiTsne failed:", e$message))
    NULL
  })
  if (is.list(tsne_result) && isFALSE(tsne_result$success)) { AddErrMsg(tsne_result$message); return(0) }
  if (is.null(tsne_result)) return(0)
  return(tsne_result)
}
