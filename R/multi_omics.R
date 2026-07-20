

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

# Per-dataset normalization for multi-omics harmonization. Reads data.proc and writes
# data.proc, so it composes with the surrounding filter steps. `opt`: "auto" resolves a
# per-layer plan from the layer's data STATE and omics TYPE, or an explicit
# NormalizingDataOmics option ("log","logcpm","clr","pqn","vst","none",...) chosen in the
# manual UI. The workflow accepts BOTH raw abundance tables (full type-specific transform)
# AND already-normalized matrices from an individual-omics workflow. A normalized layer is
# left unchanged here (it is only IQR-selected then auto-scaled downstream); a layer that
# arrives already auto-scaled is flagged, because IQR feature selection can no longer rank
# its features. Records dataSet$normInfo + dataSet$inputState for the transparency table.
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
    type <- tryCatch(dataSet$type, error = function(e) "")
    if(identical(opt, "auto")){
      declared <- tryCatch(tolower(as.character(dataSet$dataState)), error = function(e) "")
      state <- .OmicsLayerState(mat, type, declared)
      dataSet$inputState <- state
      if(identical(state, "raw")){
        tag <- .OmicsTypeNormPlan(mat, type)          # raw -> transform by omics type
      } else {
        tag <- "NA"                                    # normalized/scaled -> do not re-normalize
        if(identical(state, "scaled")){
          AddErrMsg(paste0("Layer '", sel.nms[i], "' looks already auto-scaled (per-feature mean~0, sd~1). ",
                           "Feed a NORMALIZED (not scaled) matrix so IQR feature selection can rank features."))
        }
      }
    } else {
      tag <- opt                                       # explicit method from the manual UI
      dataSet$inputState <- "manual"
    }
    if(!identical(tag, "NA") && !identical(tag, "none") && !identical(tag, "")){
      dataSet$data.proc <- .OmicsApplyNorm(mat, tag)
      # Safety net: no normalization path may leave residual NAs. The prefilter retains
      # below-LOD NAs for the QRILC paths (raw metabolomics / proteomics); the auto tags
      # impute them, but an explicit "vsn" override (not the auto "vsn_qrilc" plan) would
      # not — fill any survivors from the low tail so scaling stays defined.
      if(anyNA(as.matrix(dataSet$data.proc))){
        m <- qrilc_impute(as.matrix(dataSet$data.proc))
        dataSet$data.proc <- as.data.frame(m)
      }
    }
    dataSet$normInfo <- tag
    message("[NormalizeDataMultiOmics] ", sel.nms[i], " (", type, ") state=", dataSet$inputState, " -> ", tag)
    RegisterData(dataSet)
  }
  return(1);
}

# log2 an intensity matrix and impute the below-LOD entries with QRILC (MNAR). MS matrices
# encode non-detects as 0, and log2(0) = -Inf would crash imputeLCMD::impute.QRILC. So the
# order is strictly: (1) convert every 0 / non-positive / non-finite value to NA; (2) log2
# the remaining positive values (never sees a 0, so no -Inf); (3) QRILC-impute the NA from
# the lower tail. Used by the metabolomics path; proteomics reuses qrilc_impute directly on
# the VSN-transformed (glog2) matrix — see the "vsn_qrilc" branch in .OmicsApplyNorm.
.OmicsLogQrilc <- function(m){
  m[!is.finite(m) | m <= 0] <- NA                      # (1) 0 / below-LOD -> NA (avoid -Inf)
  m <- log2(m)                                         # (2) log2 on positive values only
  qrilc_impute(m)                                      # (3) QRILC impute the NA (MNAR)
}

# Apply one normalization TAG to a raw/normalized feature x sample matrix.
#   log/logcount/logcpm/clr/vst -> delegate to NormalizingDataOmics (negatives clipped
#     to 0 first so the log is defined); pqn -> PQN dilution correction then log2;
#     pqn_log_qrilc -> PQN -> log2 -> QRILC below-LOD impute (the auto metabolomics path);
#     vsn_qrilc -> VSN (per-sample calibration + glog2, normalize+transform in one) -> QRILC
#     below-LOD impute (the auto proteomics path); any other string (log, vsn, quantile, ...)
#     is passed straight through to NormalizingDataOmics.
.OmicsApplyNorm <- function(mat, tag){
  mat <- as.matrix(mat)
  if(tag %in% c("log", "logcount", "logcpm", "clr", "vst")){
    mat[!is.na(mat) & mat < 0] <- 0
    return(NormalizingDataOmics(mat, tag, "NA", "NA"))
  }
  if(identical(tag, "pqn")){
    mat[!is.na(mat) & mat < 0] <- 0
    return(NormalizingDataOmics(mat, "log", "pqn", "NA"))
  }
  if(identical(tag, "pqn_log_qrilc")){                 # metabolomics: PQN -> log2 -> QRILC
    m <- mat; m[!is.na(m) & m < 0] <- 0
    m <- PQNNorm(m)
    m <- .OmicsLogQrilc(m)
    d <- as.data.frame(m); rownames(d) <- rownames(mat); colnames(d) <- colnames(mat)
    return(d)
  }
  if(identical(tag, "vsn_qrilc")){                     # proteomics: VSN (steps 2-3) -> QRILC (step 4)
    d <- NormalizingDataOmics(mat, "vsn", "NA", "NA")  # per-sample calibration + glog2 in one op
    m <- as.matrix(d)
    if(anyNA(m)) m <- qrilc_impute(m)                  # impute residual below-LOD entries (MNAR)
    d <- as.data.frame(m); rownames(d) <- rownames(mat); colnames(d) <- colnames(mat)
    return(d)
  }
  NormalizingDataOmics(mat, tag, "NA", "NA")
}

# Omics-type -> normalization PLAN tag for RAW data. Used when a layer is raw (declared,
# or concluded by the state classifier). .OmicsApplyNorm consumes the tag.
#   Microbiome (ASV/OTU) -> "clr": total-sum scaling + additive pseudo-count + centered
#     log-ratio. Sparse, zero-inflated and compositional; plain log2 would leave the
#     zero-inflation and a scale that dominates the joint projection.
#   Sequencing (RNA-seq/miRNA) -> "logcpm" (TMM-normalized log2-CPM, prior.count=2) for
#     integer count data; already library-normalized values (TPM/FPKM/CPM) -> "log".
#   Metabolomics (incl. LC-MS / GC-MS labelled as metabolomics) -> "pqn_log_qrilc": PQN
#     dilution correction, log2, QRILC below-LOD impute (MS has genuine below-LOD zeros).
#   Proteomics -> "vsn_qrilc": VSN variance-stabilizing normalization (limma) normalizes each
#     sample AND applies a glog2 transform in one step (tolerates zeros, no -Inf); the residual
#     below-detection NAs it leaves are then QRILC-imputed (MNAR) — the dedicated missing-value
#     step for proteomics.
#   Generic / type-agnostic layer (e.g. "generic_layer1", an arbitrary table or a platform we
#     do not special-case such as NMR) -> "vsn": variance-stabilizing normalization is the safe
#     assumption-free default — it normalizes each sample AND applies a glog2 transform in one
#     step and tolerates zeros. The data STATE still gates it (an already-normalized generic
#     layer is skipped and only scaled).
.OmicsTypeNormPlan <- function(mat, omicsType = ""){
  t <- tolower(as.character(omicsType))
  # Generic / type-agnostic layer FIRST — checked before the biological-type patterns because
  # "generic" contains the substring "gene" (RNA). (The slot code is "generic_layer", not
  # "generic_omics", so it also avoids the "mic" in "omics" — this is belt-and-suspenders.)
  if(grepl("generic", t)) return("vsn")
  if(grepl("mic", t)) return("clr")
  if(grepl("rna|mirna|seq|count|gene", t)){
    m  <- suppressWarnings(matrix(as.numeric(as.matrix(mat)), nrow = nrow(mat)))
    nz <- m[is.finite(m) & m > 0]
    if(length(nz) > 2000L) nz <- nz[seq_len(2000L)]
    if(length(nz) > 0L && mean(nz == round(nz)) > 0.95) return("logcpm")
    return("log")
  }
  if(grepl("prot", t)) return("vsn_qrilc")                             # proteomics: VSN -> QRILC
  if(grepl("met", t)) return("pqn_log_qrilc")                          # metabolomics (incl. LC-MS/GC-MS)
  "log"                                                                # generic / unknown -> data-state log
}

# Classify a layer's INPUT STATE: "raw" (abundance/counts/intensities -> full transform),
# "normalized" (already normalized by an upstream single-omics workflow -> leave as-is,
# only IQR-select + auto-scale downstream), or "scaled" (already z-scored -> unusable for
# variance-based feature selection; the caller warns). A declared state ("raw"/"normalized")
# is honoured; "auto"/unset falls back to the statistical signature.
.OmicsLayerState <- function(mat, omicsType = "", declared = ""){
  d <- tolower(as.character(declared))
  if(length(d) != 1L) d <- ""
  if(d %in% c("raw", "normalized", "scaled")) return(d)

  m <- suppressWarnings(matrix(as.numeric(as.matrix(mat)), nrow = nrow(mat)))
  if(all(is.na(m))) return("normalized")

  # (1) Already AUTO-SCALED (z-scored): per-feature mean ~ 0 and sd ~ 1 — unambiguous.
  rmean <- suppressWarnings(rowMeans(m, na.rm = TRUE))
  rsd   <- suppressWarnings(apply(m, 1, stats::sd, na.rm = TRUE))
  ok    <- is.finite(rmean) & is.finite(rsd) & rsd > 0
  if(sum(ok) >= 5L && stats::median(abs(rmean[ok])) < 0.15 &&
     abs(stats::median(rsd[ok]) - 1) < 0.25) return("scaled")

  # (2) Already NORMALIZED + CENTERED (e.g. log-CPM): a SUBSTANTIAL negative fraction.
  #     A few stray negatives (imputation / sparse outliers) in raw data do NOT count.
  if(suppressWarnings(mean(m < 0, na.rm = TRUE)) > 0.10) return("normalized")

  # Integer-count signature: log/normalized values are essentially never all whole
  # numbers, so integer non-negative data is RAW counts regardless of magnitude — this
  # rescues a low-depth count layer (e.g. sparse microbiome, max < 20) that the
  # dynamic-range heuristic alone would misread as already-log.
  nz <- m[m > 0 & !is.na(m)]
  if(length(nz) > 2000L) nz <- nz[seq_len(2000L)]
  is.int <- length(nz) > 0L && mean(nz == round(nz)) > 0.95 &&
            suppressWarnings(min(m, na.rm = TRUE)) >= 0

  # (3) Small dynamic range (max < 20) => already in log / normalized space, unless the
  #     data is integer counts (handled above).
  mx <- suppressWarnings(max(m, na.rm = TRUE))
  if(!is.finite(mx) || mx < 20) return(if(is.int) "raw" else "normalized")

  # (4) max in (20, 50]: small RAW counts vs already-log -> integer-ness decides.
  #     max > 50 is unambiguously RAW.
  if(mx <= 50 && !is.int) return("normalized")
  "raw"
}

# Back-compat shim: return the normalization tag for a layer concluded RAW, else "NA".
.OmicsDetectNormOpt <- function(mat, omicsType = ""){
  if(identical(.OmicsLayerState(mat, omicsType), "raw")) .OmicsTypeNormPlan(mat, omicsType) else "NA"
}

# Generic PREVALENCE filter for count-type data. `data` is features x samples: keep
# features detected (value >= min.count) in at least `min.prev` fraction of samples.
# Sparse features (e.g. microbiome ASVs nonzero in a handful of samples) survive
# abundance/rank filters but are the main cause of the density spike + within-fold zero
# variance in integration — prevalence is the right cut for counts. Never reduces below
# `min.keep` features (returns the input unchanged when the cut would leave too few).
# Defaults: keep features with >= min.count reads in >= min.prev of samples. The 20%
# prevalence matches MicrobiomeAnalyst's MMP filter (ApplyAbundanceFilter "prevalence",
# smpl.perc = 0.2); the count threshold is 2 (slightly more permissive than MMP's 4).
FilterByPrevalence <- function(data, min.prev = 0.2, min.count = 2, min.keep = 10L){
  data <- as.matrix(data)
  nS <- ncol(data)
  if(nS < 5L || nrow(data) <= min.keep) return(data)
  prev <- rowSums(data >= min.count, na.rm = TRUE) / nS
  keep <- is.finite(prev) & prev >= min.prev
  if(sum(keep) >= min.keep) data[keep, , drop = FALSE] else data
}

# Microbiome prevalence filter on RELATIVE abundance: keep taxa present above
# `min.relabund` (fraction of a sample's total, default 0.01%) in at least `min.prev`
# of samples (default 10%). Robust to sequencing-depth differences (unlike a raw-count
# threshold). Never reduces below `min.keep`.
.PrevalenceFilterMic <- function(data, min.prev = 0.10, min.relabund = 1e-4, min.keep = 10L){
  data <- as.matrix(data); nS <- ncol(data)
  if(nS < 5L || nrow(data) <= min.keep) return(data)
  cs <- colSums(data, na.rm = TRUE); cs[cs <= 0] <- NA
  rel <- sweep(data, 2, cs, FUN = "/")
  prev <- rowSums(rel > min.relabund, na.rm = TRUE) / nS
  keep <- is.finite(prev) & prev >= min.prev
  if(sum(keep) >= min.keep) data[keep, , drop = FALSE] else data
}

# edgeR::filterByExpr for sequencing counts. Defaults keep genes with >= 10 reads in
# >= 70% of samples (min.count = 10, min.prop = 0.7). Never reduces below `min.keep`.
.FilterByExpr <- function(data, min.keep = 10L){
  data <- as.matrix(data)
  if(nrow(data) <= min.keep || !requireNamespace("edgeR", quietly = TRUE)) return(data)
  keep <- tryCatch(edgeR::filterByExpr(data), error = function(e) rep(TRUE, nrow(data)))
  if(sum(keep) >= min.keep) data[keep, , drop = FALSE] else data
}

# Metabolomics "80% rule": drop features present (finite, > 0) in fewer than `min.present`
# of samples. Applied per group when a grouping vector is supplied (kept if it clears the
# threshold in ANY group), else globally. Never reduces below `min.keep`.
.Rule80Filter <- function(data, group = NULL, min.present = 0.80, min.keep = 10L){
  data <- as.matrix(data); nS <- ncol(data)
  if(nS < 5L || nrow(data) <= min.keep) return(data)
  present <- is.finite(data) & data > 0
  if(!is.null(group) && length(group) == nS && length(unique(group)) > 1L){
    keep <- rep(FALSE, nrow(data))
    for(g in unique(group)){
      cols <- which(group == g)
      keep <- keep | (rowSums(present[, cols, drop = FALSE], na.rm = TRUE) / length(cols) >= min.present)
    }
  } else {
    keep <- rowSums(present, na.rm = TRUE) / nS >= min.present
  }
  if(sum(keep) >= min.keep) data[keep, , drop = FALSE] else data
}

# STEP 1 — low-abundance / quality pre-filter, per omics type, applied to RAW layers
# only (a normalized layer coming from an upstream workflow is passed through untouched
# except for constant-feature removal). Reads/writes data.proc so it composes with the
# rest of the harmonization chain. Records dataSet$filterInfo$n_before.
PrefilterOmicsByType <- function(dataName){
  if (length(dataName) == 0L || is.na(dataName) || tolower(as.character(dataName)) %in% c("na","all")) {
    sel.nms <- names(mdata.all)[vapply(mdata.all, function(x) isTRUE(x == 1), logical(1))]
  } else {
    sel.nms <- as.character(dataName);
  }
  if (length(sel.nms) == 0L) { AddErrMsg("No active dataset available to pre-filter."); return(0); }
  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i])
    int.mat <- dataSet$data.proc
    if(is.null(int.mat)) next
    int.mat <- as.matrix(int.mat); suppressWarnings(storage.mode(int.mat) <- "numeric")
    n0.feat <- nrow(int.mat)
    .ot <- tryCatch(tolower(as.character(dataSet$type)), error = function(e) "")
    .dec <- tryCatch(tolower(as.character(dataSet$dataState)), error = function(e) "")
    state <- .OmicsLayerState(int.mat, .ot, .dec)
    if(identical(state, "raw")){
      # Generic / type-agnostic layer FIRST (before the "mic"/"gene" substring patterns) —
      # a plain missing-fraction drop, no platform-specific abundance filter.
      if(grepl("generic", .ot)){
        keep <- rowMeans(!is.finite(int.mat) | int.mat == 0, na.rm = TRUE) <= 0.5
        if(sum(keep) >= 10L) int.mat <- int.mat[keep, , drop = FALSE]
      } else if(grepl("mic", .ot)){
        int.mat <- .PrevalenceFilterMic(int.mat)                 # >0.01% in >=10% of samples
      } else if(grepl("rna|mirna|seq|count|gene", .ot)){
        int.mat <- .FilterByExpr(int.mat)                        # edgeR filterByExpr
      } else if(grepl("met", .ot)){
        grp <- tryCatch(.OmicsPrimaryGroup(dataSet), error = function(e) NULL)
        int.mat <- .Rule80Filter(int.mat, group = grp)           # 80% rule (metabolomics)
      } else if(grepl("prot", .ot)){
        grp <- tryCatch(.OmicsPrimaryGroup(dataSet), error = function(e) NULL)
        int.mat <- .Rule80Filter(int.mat, group = grp, min.present = 0.70)  # valid-value filter
      } else {
        keep <- rowMeans(!is.finite(int.mat) | int.mat == 0, na.rm = TRUE) <= 0.5
        if(sum(keep) >= 10L) int.mat <- int.mat[keep, , drop = FALSE]
      }
    }
    # constant / zero-variance features carry no signal and break downstream CV folds
    .fvar <- suppressWarnings(apply(int.mat, 1, stats::var, na.rm = TRUE))
    .nzv  <- is.finite(.fvar) & .fvar > 0
    if(sum(!.nzv) > 0L && sum(.nzv) >= 2L) int.mat <- int.mat[.nzv, , drop = FALSE]
    # impute residual missing values (per-feature half-minimum) so normalization/scaling
    # are defined and NA features are not silently dropped. Keep zeros/NAs for a layer whose
    # normalization plan runs the MNAR QRILC path downstream — raw metabolomics ("pqn_log_qrilc")
    # AND raw proteomics ("vsn_qrilc"), which impute below-detection values in the dedicated
    # missing-value step; every other layer (generic -> VSN, all normalized layers) is filled
    # here, since none of them impute NA downstream (VSN tolerates zeros but leaves NA as-is).
    qrilc.path <- identical(state, "raw") &&
                  .OmicsTypeNormPlan(int.mat, .ot) %in% c("pqn_log_qrilc", "vsn_qrilc")
    if(!qrilc.path && anyNA(int.mat)){
      for(r in which(rowSums(is.na(int.mat)) > 0L)){
        v <- int.mat[r, ]; mn <- suppressWarnings(min(v[v > 0 & !is.na(v)], na.rm = TRUE))
        int.mat[r, is.na(v)] <- if(is.finite(mn)) mn / 2 else 0
      }
    }
    dataSet$data.proc <- int.mat
    dataSet$filterInfo <- list(n_before = n0.feat, n_after = nrow(int.mat), state = state)
    message("[PrefilterOmicsByType] ", sel.nms[i], " (", .ot, ", ", state, "): ", n0.feat, " -> ", nrow(int.mat))
    RegisterData(dataSet)
  }
  return(1)
}

# Best-effort primary grouping vector for a dataset (metabolomics 80% rule). Tries the
# per-dataset analysis variable, then a shared primary metadata factor; NULL -> global rule.
.OmicsPrimaryGroup <- function(dataSet){
  g <- tryCatch(dataSet$analysis.var, error = function(e) NULL)
  if(!is.null(g) && length(g) == ncol(dataSet$data.proc)) return(as.character(g))
  meta <- tryCatch(dataSet$meta.info, error = function(e) NULL)
  if(is.null(meta)) meta <- tryCatch(dataSet$meta, error = function(e) NULL)
  if(!is.null(meta) && is.data.frame(meta) && ncol(meta) >= 1L &&
     nrow(meta) == ncol(dataSet$data.proc)) return(as.character(meta[[1]]))
  NULL
}

# Per-feature IQR (Q3 - Q1) for a features x samples matrix, using the vectorized
# matrixStats::rowIQRs (C, much faster than apply on a large layer) when it is available,
# else a base-R fallback so matrixStats is not a hard dependency.
.OmicsRowIQR <- function(m){
  m <- as.matrix(m)
  if(requireNamespace("matrixStats", quietly = TRUE)){
    matrixStats::rowIQRs(m, na.rm = TRUE)
  } else {
    suppressWarnings(apply(m, 1, stats::IQR, na.rm = TRUE))
  }
}

FilterDataMultiOmicsHarmonization <- function(dataName,filterMethod, filterPercent = 0, topN = 0){
  filterPercent <- suppressWarnings(as.numeric(filterPercent));
  if (length(filterPercent) == 0L || is.na(filterPercent)) filterPercent <- 0;
  topN <- suppressWarnings(as.integer(topN));
  if (length(topN) == 0L || is.na(topN)) topN <- 0L;
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
    n0.feat <- nrow(int.mat);   # features entering this filter pass (for the transparency table)

    # High-variance feature SELECTION (spec Step 4) runs AFTER normalization; the
    # type-specific abundance/prevalence pre-filter lives in PrefilterOmicsByType (Step 1)
    # and must NOT be repeated here (a count threshold on log/CLR values is meaningless).
    # Drop any residual constant features for every layer — they break downstream CV folds.
    .ot <- tryCatch(tolower(as.character(dataSet$type)), error = function(e) "")
    .fvar <- suppressWarnings(apply(int.mat, 1, stats::var, na.rm = TRUE))
    .nzv  <- is.finite(.fvar) & .fvar > 0
    if(sum(!.nzv) > 0L && sum(.nzv) >= 2L){
      message("[filter] ", sel.nms[i], " (", .ot, "): removed ", sum(!.nzv), " constant feature(s)")
      int.mat <- int.mat[.nzv, , drop = FALSE]
    }

    if(topN > 0L){
      # High-variance selection by IQR (spec Step 4). Rank features by interquartile range
      # (Q3-Q1, robust to single-sample outliers, ranks meaningfully on normalized data).
      # First discard features with IQR == 0 (or non-finite): they have no middle-50% spread,
      # add nothing but flat noise dimensions, and shouldn't survive even the keep-all case.
      # (The z-scoring divisor is SD, already guarded by the zero-variance filter above; this
      # keeps the metric and the valid-feature count consistent.) Then keep exactly the top N
      # valid features, or all of them when a layer has <= N.
      iqrv  <- .OmicsRowIQR(int.mat)
      valid <- which(is.finite(iqrv) & iqrv > 0)
      if(length(valid) == 0L){
        data <- int.mat                                  # degenerate: nothing has spread
      } else if(length(valid) <= topN){
        data <- int.mat[valid, , drop = FALSE]           # <= N valid features -> keep all valid
      } else {
        ord  <- order(iqrv[valid], decreasing = TRUE)
        data <- int.mat[valid[ord[seq_len(topN)]], , drop = FALSE]
      }
    }else if(filterMethod == "variance"){
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

    # Per-layer feature counts for the Harmonization transparency table. PrefilterOmicsByType
    # sets n_before (the true pre-filter count) and the input state; this step updates n_after
    # to the retained count after high-variance selection while preserving both.
    .fi <- tryCatch(dataSet$filterInfo, error = function(e) NULL)
    n_before <- if (!is.null(.fi) && !is.null(.fi$n_before)) .fi$n_before else n0.feat
    .st <- if (!is.null(.fi) && !is.null(.fi$state)) .fi$state else NA
    dataSet$filterInfo <- list(n_before = n_before, n_after = nrow(data),
                               state = .st, topN = topN)

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

# Harmonization sanity check: per-layer boxplot of the (final, scaled) data.proc. After
# harmonization every layer should sit near zero with a comparable spread — a layer that
# stays shifted or skewed flags an upstream filtering / normalization problem.
PlotMultiBoxplot <- function(imgNm, dpi=150, format="png", factor="1", interactive=F){
  try(RecordRCommand(paste0("PlotMultiBoxplot(\"", imgNm, "\")")), silent = TRUE)
  load_cairo();
  load_ggplot();
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  sel.nms <- names(Filter(function(x) isTRUE(x == 1L), mdata.all))
  if(length(sel.nms) == 0L) sel.nms <- names(mdata.all)
  df.list <- vector("list", length(sel.nms))
  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i])
    dat <- tryCatch(as.matrix(dataSet$data.proc), error = function(e) NULL)
    if(is.null(dat)) next
    st <- utils::stack(as.data.frame(dat))
    st$Dataset <- sel.nms[i]
    df.list[[i]] <- st
  }
  merged.df <- do.call(rbind, df.list)
  if (exists("WfSaveFigureData"))
    tryCatch(WfSaveFigureData(sub("_dpi.*$", "", basename(imgNm)), merged.df), error = function(e) NULL)
  g <- ggplot(merged.df, aes(x=Dataset, y=values, fill=Dataset)) +
       geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
       geom_boxplot(outlier.size=0.3, alpha=0.85) +
       labs(x=NULL, y="Value (after scaling)") +
       theme_bw() +
       theme(text=element_text(size=13), legend.position="none",
             axis.text.x=element_text(angle=30, hjust=1))
  if(interactive){
    library(plotly);
    return(ggplotly(g));
  }else{
    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    print(g)
    dev.off();
    infoSet <- readSet(infoSet, "infoSet");
    infoSet$imgSet$qc_multi_boxplot <- imgNm;
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
