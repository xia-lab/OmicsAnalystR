##################################################
## R script for OmicsAnalyst
## Description: GSEA functions
## Author: G. Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#'Perform Gene Set Enrichment Analysis test on single or multiple gene expression matrices
#'@param dataSetObj file name of the data, .txt format
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
PerformGSEA<- function(dataName, file.nm, fun.type,omics.type="", input.type="loading",loading.comp=1, mode = "multi"){
  rdtSet <- .get.rdt.set();
  setres <- .loadEnrichLib(fun.type, data.org);
  current.geneset <- setres$current.geneset;
  require("fgsea");

  normalize_ids <- function(x) {
    x <- gsub("^[A-Za-z]+:", "", x);
    x <- gsub("\\s+", "", x);
    x
  }
  
  map_any_to_kegg <- function(qvec) {
    db.path.qs <- paste0(lib.path, "compound_db.qs");
    db.path.rds <- paste0(lib.path, "compound_db.rds");
    cmpd.db <- NULL;
    if (file.exists(db.path.qs)) {
      cmpd.db <- qs::qread(db.path.qs);
    } else if (file.exists(db.path.rds)) {
      cmpd.db <- readRDS(db.path.rds);
    } else {
      return(rep(NA_character_, length(qvec)));
    }

    out <- rep(NA_character_, length(qvec));
    qlow <- tolower(qvec);
    fields <- list(
      kegg_id = cmpd.db$kegg_id,
      kegg = cmpd.db$kegg,
      hmdb = cmpd.db$hmdb,
      pubchem = cmpd.db$pubchem,
      chebi = cmpd.db$chebi,
      metlin = cmpd.db$metlin,
      name = cmpd.db$name
    );
    for (fld in names(fields)) {
      vals <- fields[[fld]];
      inx <- match(qlow, tolower(vals));
      hit <- !is.na(inx);
      out[hit & is.na(out)] <- cmpd.db$kegg_id[inx[hit]];
    }

    # Fuzzy name mapping (OmicsNet approach) using synonyms
    syn.db.path <- paste0(lib.path, "syn_nms.qs");
    if (exists("MatchCompoundNames", mode = "function") && file.exists(syn.db.path)) {
      syn.db <- qs::qread(syn.db.path);
      todo <- is.na(out) & !is.na(qvec) & qvec != "";
      if (any(todo)) {
        match.res <- MatchCompoundNames(qvec[todo], cmpd.db, syn.db);
        hit.inx <- match.res$hit.inx;
        ok <- hit.inx > 0;
        out[todo][ok] <- cmpd.db$kegg_id[hit.inx[ok]];
      }
    }
    out
  }
  
  if(input.type == "loading"){
    loading.pos.xyz <- rdtSet[[reductionSet$reductionOpt]]$loading.pos.xyz.orig;
    loading.pos.xyz <- loading.pos.xyz[loading.pos.xyz$omicstype == omics.type,]
    rankedVec <- loading.pos.xyz[,loading.comp];
    names(rankedVec) <- loading.pos.xyz$ids;
  }else{
    dataSet <- readDataset(dataName);
    comp.res <- dataSet$comp.res;
    rankedVec <- NULL;

    if (!is.null(comp.res)) {
      if ("coefficient" %in% colnames(comp.res)) {
        rankedVec <- comp.res[, "coefficient"];
      } else if ("logFC" %in% colnames(comp.res)) {
        rankedVec <- comp.res[, "logFC"];
      } else if ("stat" %in% colnames(comp.res)) {
        rankedVec <- comp.res[, "stat"];
      } else if ("t" %in% colnames(comp.res)) {
        rankedVec <- comp.res[, "t"];
      } else {
        ave.idx <- match("AveExpr", colnames(comp.res));
        if (!is.na(ave.idx) && ave.idx > 1) {
          comp.mat <- as.matrix(comp.res[, 1:(ave.idx - 1), drop = FALSE]);
          mode(comp.mat) <- "numeric";
          rankedVec <- apply(comp.mat, 1, function(x) max(abs(x), na.rm = TRUE));
        } else {
          num.cols <- which(sapply(comp.res, is.numeric));
          if (length(num.cols) > 0) {
            rankedVec <- comp.res[, num.cols[1]];
          }
        }
      }
    }

    if (is.null(rankedVec)) {
      stop("GSEA failed: unable to locate a ranking column in comparison results.");
    }


    if ("ids" %in% colnames(comp.res)) {
      names(rankedVec) <- comp.res$ids;
    } else {
      names(rankedVec) <- rownames(comp.res);
    }

    # Map feature IDs to enrichment IDs (e.g., Entrez/KEGG) when available
    # NOTE: comp.res$ids may already contain post-mapping names (from data.proc rownames)
    # Strategy: Use the VALUES of rawToEntrez/enrich_ids since comp.res may have mapped names
    prefer_kegg <- FALSE;
    if (!is.null(dataSet$rawToEntrez)) {
      # Get all mapped KEGG IDs from rawToEntrez values
      all_kegg_ids <- unname(dataSet$rawToEntrez);
      kegg_ids_only <- all_kegg_ids[grepl("^C\\d{5}$", all_kegg_ids)];

      # For each feature in rankedVec, check if it's already a KEGG ID or needs mapping
      mapped_raw <- character(length(names(rankedVec)));
      for (i in seq_along(names(rankedVec))) {
        feat_id <- names(rankedVec)[i];

        # Strategy 1: Feature is already a KEGG ID
        if (grepl("^C\\d{5}$", feat_id)) {
          mapped_raw[i] <- feat_id;
        } else {
          # Strategy 2: Try matching feat_id against rawToEntrez NAMES (original mapping table)
          idx <- match(feat_id, names(dataSet$rawToEntrez));
          if (!is.na(idx)) {
            mapped_raw[i] <- dataSet$rawToEntrez[idx];
          } else {
            # Strategy 3: Check if feat_id is itself a VALUE in rawToEntrez (already mapped)
            idx2 <- match(feat_id, dataSet$rawToEntrez);
            if (!is.na(idx2)) {
              mapped_raw[i] <- feat_id;
            } else {
              mapped_raw[i] <- NA_character_;
            }
          }
        }
      }

      valid_mapped <- sum(!is.na(mapped_raw) & mapped_raw != "");
      if (valid_mapped > 0) {
        names(rankedVec) <- mapped_raw;
        prefer_kegg <- sum(grepl("^C\\d{5}$", mapped_raw), na.rm = TRUE) > 0;
      }
    }
    if (!prefer_kegg && !is.null(dataSet$enrich_ids)) {
      # Same three-strategy approach for enrich_ids
      mapped <- character(length(names(rankedVec)));
      for (i in seq_along(names(rankedVec))) {
        feat_id <- names(rankedVec)[i];

        if (grepl("^C\\d{5}$", feat_id)) {
          mapped[i] <- feat_id;
        } else {
          idx <- match(feat_id, names(dataSet$enrich_ids));
          if (!is.na(idx)) {
            mapped[i] <- dataSet$enrich_ids[idx];
          } else {
            idx2 <- match(feat_id, dataSet$enrich_ids);
            if (!is.na(idx2)) {
              mapped[i] <- feat_id;
            } else {
              mapped[i] <- NA_character_;
            }
          }
        }
      }

      valid_mapped <- sum(!is.na(mapped) & mapped != "");
      if (valid_mapped > 0) {
        names(rankedVec) <- mapped;
      }
    }

    # Filter out NA, empty strings, and "NA" string before dedup
    valid_idx <- !is.na(names(rankedVec)) & names(rankedVec) != "" & names(rankedVec) != "NA";
    rankedVec <- rankedVec[valid_idx];

    tmp <- tapply(rankedVec, names(rankedVec), function(x) max(x, na.rm = TRUE));
    rankedVec <- as.numeric(tmp);
    names(rankedVec) <- names(tmp);
  }

  # Normalize IDs (handles KEGG-style prefixes like "hsa:")
  names(rankedVec) <- normalize_ids(names(rankedVec))
  current.geneset <- lapply(current.geneset, function(x) normalize_ids(unlist(x)))

  scoreType <- "std";
  if (all(rankedVec > 0, na.rm = TRUE)) {
    scoreType <- "pos";
  } else if (all(rankedVec < 0, na.rm = TRUE)) {
    scoreType <- "neg";
  }

  # Filter pathways by overlap with rankedVec
  overlap_sizes <- sapply(current.geneset, function(gs) sum(gs %in% names(rankedVec)))
  if (is.null(names(rankedVec)) || length(names(rankedVec)) == 0) {
    return(0);
  }

  # Fallback mapping using enrich_ids VALUES (useful for metabolites where IDs are KEGG/HMDB)
  if (sum(overlap_sizes >= 1) == 0 && !is.null(dataSet$enrich_ids)) {
    alt <- names(rankedVec)
    alt <- unname(dataSet$enrich_ids[match(alt, unname(dataSet$enrich_ids))])
    if (sum(!is.na(alt)) > 0) {
      names(rankedVec) <- alt
      rankedVec <- rankedVec[!is.na(names(rankedVec)) & names(rankedVec) != ""]
      tmp <- tapply(rankedVec, names(rankedVec), function(x) max(x, na.rm = TRUE))
      rankedVec <- as.numeric(tmp)
      names(rankedVec) <- names(tmp)
      names(rankedVec) <- normalize_ids(names(rankedVec))
      overlap_sizes <- sapply(current.geneset, function(gs) sum(gs %in% names(rankedVec)))
    }
  }

  # Metabolite fallback: map common names to KEGG compound IDs when geneset is KEGG compounds
  if (sum(overlap_sizes >= 1) == 0) {
    gs_ids <- unique(unlist(current.geneset))
    kegg_gs <- sum(grepl("^C\\d{5}$", gs_ids)) > 0
    if (kegg_gs) {
      kegg_map <- map_any_to_kegg(names(rankedVec))
      names(rankedVec) <- kegg_map
      keep_inx <- grepl("^C\\d{5}$", names(rankedVec))
      rankedVec <- rankedVec[keep_inx]
      if (length(rankedVec) > 0) {
        tmp <- tapply(rankedVec, names(rankedVec), function(x) max(x, na.rm = TRUE))
        rankedVec <- as.numeric(tmp)
        names(rankedVec) <- names(tmp)
      }
      overlap_sizes <- sapply(current.geneset, function(gs) sum(gs %in% names(rankedVec)))
    }
  }

  # Final fallback: map KEGG compound genesets to common names and match by name
  if (sum(overlap_sizes >= 1) == 0) {
    gs_ids <- unique(unlist(current.geneset))
    kegg_gs <- sum(grepl("^C\\d{5}$", gs_ids)) > 0
    if (kegg_gs) {
      cmpd.db <- readRDS(paste0(lib.path, "compound_db.rds"));
      kegg_to_name <- setNames(cmpd.db$name, cmpd.db$kegg_id)
      current.geneset <- lapply(current.geneset, function(gs) {
        nm <- unname(kegg_to_name[gs])
        nm <- nm[!is.na(nm) & nm != ""]
        normalize_ids(nm)
      })
      names(rankedVec) <- normalize_ids(names(rankedVec))
      overlap_sizes <- sapply(current.geneset, function(gs) sum(gs %in% names(rankedVec)))
    }
  }

  minSize <- 15
  maxSize <- 500
  if (sum(overlap_sizes >= minSize & overlap_sizes <= maxSize) == 0) {
    minSize <- 5
  }
  if (sum(overlap_sizes >= minSize & overlap_sizes <= maxSize) == 0) {
    minSize <- 1
  }

  filt.geneset <- current.geneset[which(overlap_sizes >= minSize & overlap_sizes <= maxSize)]

  if (length(rankedVec) == 0 || is.null(names(rankedVec)) || all(is.na(names(rankedVec))) || all(names(rankedVec) == "")) {
    return(0);
  }

  if(mode == "simple"){
    fgseaRes <- fgsea(pathways = filt.geneset, 
                      stats = rankedVec,
                      minSize=minSize,
                      maxSize=maxSize,
                      scoreType=scoreType,
                      nperm=10000);
  } else {
    fgseaRes <- fgsea(pathways = filt.geneset, 
                      stats = rankedVec,
                      minSize=minSize,
                      maxSize=maxSize,
                      scoreType=scoreType);
  }
  
  fgseaRes <- fgseaRes[!duplicated(fgseaRes$pathway),];
  rownames(fgseaRes) <- make.names(fgseaRes$pathway, unique=TRUE)
  fgseaRes <- fgseaRes[,c("size","ES", "pval", "pathway", "padj")]
  fgseaRes <- fgseaRes[order(fgseaRes$pval),]
  if(nrow(fgseaRes[which(fgseaRes$pval < 0.05),])<20 ){
    if(nrow(fgseaRes)>20){
      fgseaRes <- fgseaRes[c(1:20),]
    }
  }else{
    fgseaRes <- fgseaRes[which(fgseaRes$pval < 0.05),]
  } 
  
  current.mset <- current.geneset[fgseaRes$pathway]
  current.mset <- current.mset[!duplicated(names(current.mset))]
  
  ora.vec <- names(rankedVec)
  ora.nms <- ora.vec #doEntrez2SymbolMapping(ora.vec)
  hits.query <- lapply(current.mset, 
                       function(x) {
                         unique(ora.nms[ora.vec%in%unlist(x)]);
                       });
  
  set.num <- unlist(lapply(current.mset, function(x){length(unique(x))}), use.names=TRUE);
  names(hits.query) <- names(current.mset);
  hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=TRUE);
  
  qs::qsave(hits.query, "hits_query.qs");
  if (nrow(fgseaRes) > 0) {
    hit.idx <- match(fgseaRes$pathway, names(hit.num));
    total.idx <- match(fgseaRes$pathway, names(set.num));
    fgseaRes$hits <- ifelse(is.na(hit.idx), 0, hit.num[hit.idx]);
    fgseaRes$total <- ifelse(is.na(total.idx), 0, set.num[total.idx]);
  } else {
    fgseaRes$hits <- numeric(0);
    fgseaRes$total <- numeric(0);
  }

  fgseaRes.pre <- fgseaRes
  
  fgseaRes <- fgseaRes[which(fgseaRes$hits>1),];
  fgseaRes <- fgseaRes[which(fgseaRes$hits<500),];
  fgseaRes <- fgseaRes[which(fgseaRes$total<2000),];
  
  fgseaRes=fgseaRes[order(fgseaRes$pval),];
  if(nrow(fgseaRes[which(fgseaRes$pval < 0.05),])<20 ){
    if(nrow(fgseaRes)>20){
      fgseaRes <- fgseaRes[c(1:20),];
    }
  }else{
    fgseaRes <- fgseaRes[which(fgseaRes$pval < 0.05),];
  } 
  
  if (nrow(fgseaRes) == 0 && nrow(fgseaRes.pre) > 0) {
    fgseaRes <- fgseaRes.pre[order(fgseaRes.pre$pval), ];
    if (nrow(fgseaRes) > 20) {
      fgseaRes <- fgseaRes[1:20, ];
    }
  }
  
  fgseaRes <- data.frame(fgseaRes, stringsAsFactors=FALSE);
  
  #get gene symbols
  current.msg <<- "Functional analysis was completed";
  
  # write json
  fun.anot <- hits.query; 
  fun.pval <- fgseaRes[,3]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  fun.padj <- fgseaRes[,5]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  es.num <- fgseaRes[,2]; if(length(es.num) ==1) { es.num <- matrix(es.num) };
  setids <- setres$current.setids;
  if (is.null(setids) || length(setids) == 0) {
    setids <- names(current.geneset);
  }
  fun.ids <- as.vector(setids[names(fun.anot)]);
  if (length(fun.ids) == 1) { fun.ids <- matrix(fun.ids) };
  if (all(is.na(fun.ids))) {
    fun.ids <- names(fun.anot);
  }
  
  json.res <- list(
    fun.link = setres$current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    pathname = fgseaRes[,"pathway"],
    es.num = es.num,
    hit.num = fgseaRes[,"hits"],
    total = fgseaRes[,"total"]
  );
  

  json.res$org <- data.org
  json.res$analType <- anal.type
  json.res$naviString <- "GSEA";
  
  json.mat <- RJSONIO::toJSON(json.res);
  json.nm <- paste(file.nm, ".json", sep="");
  partialToBeSaved <- c(partialToBeSaved, c(json.nm, "current_geneset.qs"))
  sink(json.nm)
  cat(json.mat);
  sink();
  
  # write csv
  
  ftype <- fun.type
  if(fun.type %in% c("bp", "mf", "cc")){
    ftype <- paste0("go_", fun.type);
  }
  
  csvDf <- data.frame(Name=fgseaRes$pathway, Total=fgseaRes$total, Hits=fgseaRes$hits, EnrichmentScore=fgseaRes$ES, Pval=fgseaRes$pval, Padj=fgseaRes$padj);
  fast.write(csvDf, file=paste0(file.nm, ".csv"));
  
  return(1);
}