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

  # Remove near-zero-variance features from each omics block.
  # data.proc is features(rows) x samples(cols). Compute per-row variance.
  # Features constant across samples cause failures in DIABLO BER CV folds
  # and in MOFA/MCIA factor extraction.
  for (nm in names(data.list)) {
    mat <- data.list[[nm]]  # features x samples
    feat.var <- apply(mat, 1, var, na.rm = TRUE)
    keep <- !is.na(feat.var) & feat.var > 1e-10
    if (sum(!keep) > 0) {
      removed <- rownames(mat)[!keep]
      data.list[[nm]] <- mat[keep, , drop = FALSE]
      # Update feature tracking vectors
      rm.inx <- featureNms %in% removed & omics.vec == nm
      featureNms <- featureNms[!rm.inx]
      omics.vec <- omics.vec[!rm.inx]
      uniqFeats <- uniqFeats[!rm.inx]
      comp.res1 <- comp.res1[!rm.inx, , drop = FALSE]
      enrich.nms1 <- enrich.nms1[!rm.inx]
      comp.res.inx1 <- comp.res.inx1[!rm.inx]
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
      print("s7")
    loading.pos.xyz = mcoin$mcoa$Tco;
    loading.pos.xyz$ids = featureNms;
    loading.pos.xyz$type <- omics.vec;
    # get sample and weight names
    names = rownames(pos.xyz)
 
    var.exp <- t(mcoin$mcoa$cov2);
    var.exp <- round(var.exp, digits = 3);
    rownames(var.exp) <- colnames(pos.xyz);
  } else if (reductionOpt == "mofa") {

    # Isolate MOFA2 in subprocess
    data.list.mofa <- lapply(seq_along(data.list), function(i) {
        mat <- data.list[[i]]
        sanitized_names <- gsub("[[:cntrl:]]|[^[:ascii:]]", "_", rownames(mat), perl = TRUE)
        rownames(mat) <- paste0(sanitized_names, "_", omics.type[i])
        as.matrix(mat)
      })
      names(data.list.mofa) <- names(data.list)

      mofa_result <- tryCatch({
        rsclient_isolated_exec(
          func_body = function(input_data) {
            library(MOFA2)
            library(reshape2)
            data.list <- input_data$data_obj
            num_factors <- input_data$num_factors
            MOFAobject <- MOFA2::create_mofa_from_matrix(data.list)
            data_opts <- MOFA2::get_default_data_options(MOFAobject)
            model_opts <- MOFA2::get_default_model_options(MOFAobject)
            model_opts$num_factors <- num_factors
            train_opts <- MOFA2::get_default_training_options(MOFAobject)
            MOFAobject <- MOFA2::prepare_mofa(
              object = MOFAobject, data_options = data_opts,
              model_options = model_opts, training_options = train_opts
            )
            model <- tryCatch({
              MOFA2::run_mofa(MOFAobject, save_data = FALSE, use_basilisk = FALSE)
            }, error = function(e) {
              if (grepl("mofapy2", e$message, ignore.case = TRUE)) {
                MOFA2::run_mofa(MOFAobject, save_data = FALSE, use_basilisk = TRUE)
              } else stop(e)
            })
            factors <- MOFA2::get_factors(model, as.data.frame = TRUE)
            pos.xyz <- reshape2::dcast(factors, sample ~ factor, value.var = "value")
            rownames(pos.xyz) <- pos.xyz$sample
            pos.xyz <- pos.xyz[, -1]
            weights <- MOFA2::get_weights(model, as.data.frame = TRUE)
            loading.pos.xyz <- reshape2::dcast(weights, feature ~ factor, value.var = "value")
            loading.pos.xyz$ids <- as.character(loading.pos.xyz$feature)
            loading.pos.xyz <- loading.pos.xyz[, -1]
            loading.pos.xyz$ids <- gsub("_.*", "", loading.pos.xyz$ids)
            var.exp <- model@cache[["variance_explained"]][["r2_per_factor"]][[1]] / 100
            var.exp <- round(var.exp, digits = 3)
            gc(verbose = FALSE, full = TRUE)
            return(list(pos.xyz = pos.xyz, loading.pos.xyz = loading.pos.xyz, var.exp = var.exp))
          },
          input_data = list(data_obj = data.list.mofa, num_factors = ncomps),
          packages = c("MOFA2", "reshape2", "basilisk", "qs"),
          timeout = 600,
          output_type = "qs"
        )
      }, error = function(e) {
        AddErrMsg(paste("MOFA2 reduction failed:", e$message))
        NULL
      })
      if (is.list(mofa_result) && isFALSE(mofa_result$success)) { AddErrMsg(mofa_result$message); return(0) }
      if (is.null(mofa_result)) return(0)

      pos.xyz <- mofa_result$pos.xyz
      loading.pos.xyz <- mofa_result$loading.pos.xyz
      loading.pos.xyz$type <- omics.vec
      var.exp <- mofa_result$var.exp
      reductionSet[["mofa.complete"]] <- TRUE

  } else if (reductionOpt == "diablo"){ # pos pars to tune: value from 0-1 inside matrix, which metadata to predict
    diablo.meta.type <- reductionSet$dataSet$meta.types[diabloMeta];
    reductionSet$diabloMeta <- diabloMeta;
    reductionSet$diabloPar <- diabloPar;

    # Isolate mixOmics in subprocess
    rsclient_input <- list(
        data.list = data.list,
        meta = reductionSet$meta,
        diabloMeta = diabloMeta,
        diabloPar = diabloPar,
        diablo.meta.type = diablo.meta.type,
        ncomps = ncomps,
        omics.type = omics.type,
        omics.vec = omics.vec
      )

      diablo_result <- tryCatch({
        rsclient_isolated_exec(
          func_body = function(input_data) {
            library(mixOmics)
            library(BiocParallel)
            register(SerialParam())
            data.list <- input_data$data.list
            meta <- input_data$meta
            diabloMeta <- input_data$diabloMeta
            diabloPar <- input_data$diabloPar
            diablo.meta.type <- input_data$diablo.meta.type
            ncomps <- input_data$ncomps
            omics.type <- input_data$omics.type
            omics.vec <- input_data$omics.vec

            if (diablo.meta.type == "disc") {
              Y <- meta[, diabloMeta]
              design <- matrix(diabloPar, ncol = length(data.list), nrow = length(data.list),
                              dimnames = list(names(data.list), names(data.list)))
              diag(design) <- 0
              data.list <- lapply(data.list, function(x) {
                xt <- t(x)
                feat.var <- apply(xt, 2, var, na.rm = TRUE)
                xt[, !is.na(feat.var) & feat.var > 0, drop = FALSE]
              })
              model <- mixOmics::block.splsda(X = data.list, Y = Y, ncomp = ncomps, design = design)
            } else {
              meta.var <- meta[, diabloMeta]
              Y <- matrix(c(as.numeric(as.character(meta.var))))
              rownames(Y) <- rownames(meta)
              design <- matrix(diabloPar, ncol = length(data.list), nrow = length(data.list),
                              dimnames = list(names(data.list), names(data.list)))
              diag(design) <- 0
              data.list <- lapply(data.list, function(x) {
                xt <- t(x)
                feat.var <- apply(xt, 2, var, na.rm = TRUE)
                xt[, !is.na(feat.var) & feat.var > 0, drop = FALSE]
              })
              model <- mixOmics::block.spls(X = data.list, Y = Y, ncomp = ncomps, design = design, mode = "regression")
            }

            # Calculate centroid factor scores
            variates <- model$variates
            variates$Y <- NULL
            variates <- lapply(variates, function(df) {
              for (col in 1:ncol(df)) {
                col_min <- min(df[, col])
                col_max <- max(df[, col])
                df[, col] <- (df[, col] - col_min) / (col_max - col_min) - 0.5
              }
              df
            })
            pos.xyz <- lapply(c(Factor1 = 'comp1', Factor2 = 'comp2', Factor3 = 'comp3', Factor4 = 'comp4', Factor5 = 'comp5'), function(w) {
              xORy <- lapply(variates, function(v) v[, w, drop = FALSE])
              xORy <- Reduce(x = xORy, f = cbind)
              rowMeans(xORy)
            })
            pos.xyz <- as.data.frame(pos.xyz)

            # Get loadings
            loading.list <- vector("list", length(omics.type))
            for (i in 1:length(omics.type)) {
              temp.mat <- as.data.frame(model[["loadings"]][[i]])
              rnms <- rownames(temp.mat)
              temp.mat <- as.data.frame(apply(temp.mat, 2, function(x) {
                x_range <- max(x) - min(x)
                if (x_range == 0) return(rep(0, length(x)))
                (x - min(x)) / x_range
              }))
              rownames(temp.mat) <- rnms
              temp.mat$ids <- rownames(temp.mat)
              temp.mat$type <- omics.type[i]
              loading.list[[i]] <- temp.mat
            }
            loading.pos.xyz <- do.call(rbind, loading.list)
            colnames(loading.pos.xyz) <- c(paste0("Factor", 1:ncomps), "ids", "type")

            var.exp <- model$prop_expl_var
            var.exp$Y <- NULL
            var.exp <- as.matrix(as.data.frame(var.exp))
            var.exp <- round(var.exp, digits = 3)
            rownames(var.exp) <- colnames(pos.xyz)
            loading.pos.xyz$type <- omics.vec

            # Save model for fallback PlotDiabloBER/Circos/Loading
            qs::qsave(model, "diablo_model.qs", preset = "fast")

            # Generate BER diagnostic plot
            ber_img <- NULL
            opt.comp <- NULL
            ber_table <- NULL
            tryCatch({
              if (diablo.meta.type == "disc") {
                # Limit ncomp for perf to avoid zero-variance in folds (same as MicrobiomeAnalyst)
                ncomp_perf <- min(model$ncomp[1], 8)
                res_perf <- model
                res_perf$ncomp <- rep(ncomp_perf, length(model$ncomp))
                if (!is.null(res_perf$keepX)) {
                  res_perf$keepX <- lapply(res_perf$keepX, function(k) k[seq_len(min(length(k), ncomp_perf))])
                }
                if (!is.null(res_perf$keepA)) {
                  res_perf$keepA <- lapply(res_perf$keepA, function(k) k[seq_len(min(length(k), ncomp_perf))])
                }
                perf.res <- mixOmics:::perf(res_perf, validation = 'Mfold', folds = 5, nrepeat = 1, dist = 'max.dist', near.zero.var = TRUE)
                tmp_ber <- "diablo_berdpi150.png"
                if (!is.null(perf.res$choice.ncomp)) {
                  opt.comp <- median(perf.res$choice.ncomp$WeightedVote)
                }
                # Extract BER table
                tryCatch({
                  wv <- perf.res$WeightedVote.error.rate
                  if (!is.null(wv) && is.list(wv)) {
                    first_mat <- wv[[1]]
                    if (is.matrix(first_mat)) {
                      ber_table <- data.frame(Component = colnames(first_mat), stringsAsFactors = FALSE)
                      n_dist <- length(wv)
                      for (dist_nm in names(wv)) {
                        mat <- wv[[dist_nm]]
                        for (rn in rownames(mat)) {
                          col_name <- if (n_dist == 1) rn else paste0(dist_nm, ".", rn)
                          ber_table[[col_name]] <- round(as.numeric(mat[rn, ]), 4)
                        }
                      }
                    }
                  }
                  if (is.null(ber_table) && !is.null(perf.res$error.rate)) {
                    er <- perf.res$error.rate
                    if (is.list(er)) {
                      first_el <- er[[1]]
                      if (is.numeric(first_el) && !is.matrix(first_el)) {
                        comp_names <- names(first_el)
                        if (is.null(comp_names)) comp_names <- paste0("comp", seq_along(first_el))
                        ber_table <- data.frame(Component = comp_names, stringsAsFactors = FALSE)
                        for (dist_nm in names(er)) ber_table[[dist_nm]] <- round(as.numeric(er[[dist_nm]]), 4)
                      } else if (is.matrix(first_el)) {
                        ber_table <- data.frame(Component = colnames(first_el), stringsAsFactors = FALSE)
                        n_dist <- length(er)
                        for (dist_nm in names(er)) {
                          mat <- er[[dist_nm]]
                          for (rn in rownames(mat)) {
                            col_name <- if (n_dist == 1) rn else paste0(dist_nm, ".", rn)
                            ber_table[[col_name]] <- round(as.numeric(mat[rn, ]), 4)
                          }
                        }
                      }
                    }
                  }
                }, error = function(e) {
                  message("BER table extraction error: ", e$message)
                })

                # Generate ggplot BER line plot
                if (!is.null(ber_table)) {
                  library(see)
                  dt <- data.table::as.data.table(ber_table)
                  dt_long <- data.table::melt(dt, id.vars = "Component", variable.name = "Metric", value.name = "Error.Rate")
                  dt_long <- as.data.frame(dt_long)
                  p1 <- ggplot2::ggplot(dt_long, ggplot2::aes(x = Component, y = Error.Rate, group = Metric)) +
                    ggplot2::geom_line(ggplot2::aes(color = Metric), linewidth = 2) +
                    see::scale_fill_okabeito() + see::scale_color_okabeito() +
                    ggplot2::labs(x = "Component #", y = "Error Rate", title = "") +
                    ggplot2::theme_minimal(base_size = 15) +
                    ggplot2::theme(legend.text = ggplot2::element_text(size = 16),
                                   legend.position = c(0.9, 0.95),
                                   legend.title = ggplot2::element_text(size = 0))
                  Cairo::Cairo(file = tmp_ber, width = 8, height = 7, type = "png", bg = "white", unit = "in", dpi = 150)
                  print(p1)
                  dev.off()
                  ber_img <- tmp_ber
                } else {
                  Cairo::Cairo(file = tmp_ber, width = 8, height = 7, type = "png", bg = "white", unit = "in", dpi = 150)
                  plot(perf.res)
                  dev.off()
                  ber_img <- tmp_ber
                }
              }
            }, error = function(e) {
              message("BER diagnostic failed (non-fatal): ", e$message,
                      "\nNote: Cross-validation performance (BER) may fail due to small sample size ",
                      "or features with zero variance within individual folds. ",
                      "This does not affect the DIABLO model results.")
              tryCatch(dev.off(), error = function(e2) {})
            })

            # Generate Loading plot
            loading_img <- NULL
            tryCatch({
              library(grid)
              library(gridExtra)
              library(cowplot)
              fig.list <- list()
              ncomp_plot <- min(ncomps, 3)
              for (cc in 1:ncomp_plot) {
                local({
                  comp_idx <- cc
                  fig.list[[comp_idx]] <<- cowplot::as_grob(function() {
                    par(mar = c(4, 12, 2, 2))
                    mixOmics::plotLoadings(model, ndisplay = 10, comp = comp_idx, contrib = "max",
                                           method = "median", size.name = 1.1, legend = TRUE)
                  })
                })
              }
              tmp_loading <- "diablo_loadingdpi150.png"
              h <- 8 * length(fig.list)
              Cairo::Cairo(file = tmp_loading, width = 13, height = h, type = "png", bg = "white", unit = "in", dpi = 150)
              gridExtra::grid.arrange(grobs = fig.list, nrow = length(fig.list))
              dev.off()
              loading_img <- tmp_loading
            }, error = function(e) {
              message("Loading plot failed (non-fatal): ", e$message)
              tryCatch(dev.off(), error = function(e2) {})
            })

            # Generate circos JSON for interactive chord diagram (same logic as GenerateDiabloCircosJson)
            tryCatch({
              block_names <- names(model$X)
              X_proj <- lapply(model$X, function(x) x[, which(apply(x, 2, var) > 0), drop = FALSE])
              cutoff <- 0.5
              maxEdges <- 100
              all_edges <- list()
              # Compute all pairwise cross-block correlations
              for (bi in 1:(length(X_proj)-1)) {
                for (bj in (bi+1):length(X_proj)) {
                  cor_cross <- cor(X_proj[[bi]], X_proj[[bj]])
                  sig_idx <- which(abs(cor_cross) > cutoff, arr.ind = TRUE)
                  if (nrow(sig_idx) > 0) {
                    for (k in 1:nrow(sig_idx)) {
                      all_edges[[length(all_edges)+1]] <- list(
                        source = rownames(cor_cross)[sig_idx[k,1]],
                        target = colnames(cor_cross)[sig_idx[k,2]],
                        corr = round(cor_cross[sig_idx[k,1], sig_idx[k,2]], 4),
                        type1 = block_names[bi], type2 = block_names[bj],
                        label1 = rownames(cor_cross)[sig_idx[k,1]],
                        label2 = colnames(cor_cross)[sig_idx[k,2]])
                    }
                  }
                }
              }
              # Cap at maxEdges by top absolute correlation
              if (length(all_edges) > maxEdges) {
                corrs <- sapply(all_edges, function(e) abs(e$corr))
                all_edges <- all_edges[order(corrs, decreasing=TRUE)[1:maxEdges]]
              }
              jsonlite::write_json(list(DIABLO = all_edges), "diablo_circos.json", auto_unbox = TRUE, pretty = FALSE)
            }, error = function(e) message("[DIABLO circosJSON] ", e$message))

            gc(verbose = FALSE, full = TRUE)
            return(list(
              pos.xyz = pos.xyz, loading.pos.xyz = loading.pos.xyz, var.exp = var.exp,
              ber_img = ber_img, loading_img = loading_img, opt.comp = opt.comp, ber_table = ber_table
            ))
          },
          input_data = rsclient_input,
          packages = c("mixOmics", "qs", "Cairo", "grid", "gridExtra", "cowplot", "ggplot2", "see", "data.table"),
          timeout = 300,
          output_type = "qs"
        )
      }, error = function(e) {
        AddErrMsg(paste("DIABLO reduction failed:", e$message))
        NULL
      })
      if (is.list(diablo_result) && isFALSE(diablo_result$success)) { AddErrMsg(diablo_result$message); return(0) }
      if (is.null(diablo_result)) return(0)

      pos.xyz <- diablo_result$pos.xyz
      loading.pos.xyz <- diablo_result$loading.pos.xyz
      var.exp <- diablo_result$var.exp

      # Register pre-generated BER/Loading images
      if (!is.null(diablo_result$ber_img)) {
        infoSet$imgSet[["diablo_ber"]] <- diablo_result$ber_img
        reductionSet[["diablo"]]$ber_done <- TRUE
      }
      if (!is.null(diablo_result$opt.comp)) {
        reductionSet[["diablo"]]$opt.ncomp <- diablo_result$opt.comp
      }
      if (!is.null(diablo_result$loading_img)) {
        infoSet$imgSet[["diablo_loading"]] <- diablo_result$loading_img
        reductionSet[["diablo"]]$loading_done <- TRUE
      }
      if (!is.null(diablo_result$ber_table)) {
        reductionSet[["diablo"]]$ber_table <- diablo_result$ber_table
        tryCatch({
          arrow_path <- "diablo_ber_table.arrow"
          df_ber <- diablo_result$ber_table
          if (file.exists(arrow_path)) { unlink(arrow_path); Sys.sleep(0.01) }
          arrow::write_feather(df_ber, arrow_path, compression = "uncompressed")
          Sys.sleep(0.02)
        }, error = function(e) {
          warning(paste("BER Arrow export failed:", e$message))
        })
      }

  }

  # preserve original order
    
  loading.pos.xyz <- loading.pos.xyz[match(uniqFeats, paste0(loading.pos.xyz$ids, "_", loading.pos.xyz$type)), ]
  loading.pos.xyz$label <-  invert_named_vector(enrich.nms1)[as.character(loading.pos.xyz$ids)];
  pos.xyz <- pos.xyz[match(rownames(reductionSet$meta), rownames(pos.xyz)), ];
  #loading.pos.xyz$filenm <-   filenms
  #update colnames to "Loading"
  colnames(loading.pos.xyz)[c(1:ncomps)] <- c(paste0("Loading", 1:ncomps))
  # Initialize the list element before assigning properties
  reductionSet[[reductionOpt]] <- list()
  reductionSet[[reductionOpt]]$pos.xyz <- pos.xyz;
  reductionSet[[reductionOpt]]$loading.pos.xyz <- loading.pos.xyz;
  reductionSet[[reductionOpt]]$var.exp <- var.exp;
  fileNm <- paste0("loading_result_", reductionOpt);
  reductionSet[[reductionOpt]]$loading.file.nm <- fileNm;
  infoSet$imgSet[[reductionOpt]]$loading.pos.xyz <- loading.pos.xyz;
  fast.write.csv(loading.pos.xyz,file=fileNm);

  # Export to Arrow for Java DataTable zero-copy access
  tryCatch({
    arrow_path <- paste0("loading_", reductionOpt, ".arrow")
    df <- loading.pos.xyz
    # Add row_names_id column
    df <- cbind(row_names_id = rownames(df), df)
    # Convert factors to character
    for (col in names(df)) {
      if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])
    }
    if (file.exists(arrow_path)) {
      unlink(arrow_path)
      Sys.sleep(0.01)
    }
    arrow::write_feather(df, arrow_path, compression = "uncompressed")
    Sys.sleep(0.02)
  }, error = function(e) {
    warning(paste("Loading Arrow export failed:", e$message))
  })

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

  # Isolate ade4 in subprocess
  mcia_script_path <- NULL
    candidates <- c(
      "../../rscripts/OmicsAnalystR/R/util_mcia.R",
      "../OmicsAnalystR/R/util_mcia.R",
      file.path(getwd(), "util_mcia.R")
    )
    for (candidate in candidates) {
      if (file.exists(candidate)) {
        mcia_script_path <- normalizePath(candidate)
        break
      }
    }
    mcia_script_rc <- NULL
    if (!is.null(mcia_script_path)) {
      mcia_script_rc <- sub("\\.R$", ".Rc", mcia_script_path)
      if (!file.exists(mcia_script_rc)) mcia_script_rc <- NULL
    }

    mcia_result <- tryCatch({
      rsclient_isolated_exec(
        func_body = function(input_data) {
          library(ade4)
          df.list <- input_data$data_obj
          if (!is.null(input_data$mcia_script_rc) && file.exists(input_data$mcia_script_rc)) {
            compiler::loadcmp(input_data$mcia_script_rc, envir = globalenv())
          } else if (!is.null(input_data$mcia_script_path) && file.exists(input_data$mcia_script_path)) {
            source(input_data$mcia_script_path, local = FALSE)
          } else {
            AddErrMsg("Cannot find util_mcia.R or util_mcia.Rc script"); return(0);
          }
          result <- perform_mcia(df.list, cia.nf = input_data$cia.nf,
                                 cia.scan = input_data$cia.scan,
                                 nsc = input_data$nsc, svd = input_data$svd)
          gc(verbose = FALSE, full = TRUE)
          return(result)
        },
        input_data = list(
          data_obj = df.list,
          cia.nf = cia.nf,
          cia.scan = cia.scan,
          nsc = nsc,
          svd = svd,
          mcia_script_path = mcia_script_path,
          mcia_script_rc = mcia_script_rc
        ),
        packages = c("ade4", "qs"),
        timeout = 180,
        output_type = "qs"
      )
    }, error = function(e) {
      AddErrMsg(paste("MCIA analysis failed:", e$message))
      NULL
    })
    if (is.list(mcia_result) && isFALSE(mcia_result$success)) { AddErrMsg(mcia_result$message); return(0) }
    if (is.null(mcia_result)) return(0)
    return(mcia_result)
}

PlotDimredVarexp <- function(imgNm, dpi=150, format="png"){
  infoSet <- readSet(infoSet, "infoSet");
  load_cairo();
  library(see)
  load_ggplot();
  sel.inx <- mdata.all==1;
 
  sel.nms <- names(mdata.all)[sel.inx]
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
 
  reductionSet <- .get.rdt.set();
  if (is.null(reductionSet$reductionOpt) || is.null(reductionSet[[reductionSet$reductionOpt]])) {
    AddErrMsg("Dimension reduction results not available. Please run the analysis first.");
    return(0);
  }
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

  
  Cairo(file=imgNm, width=8, height=7, type=format, bg="white", unit="in", dpi=dpi);
  print(p1)
  dev.off();

  infoSet$imgSet[[paste0("dimred_varexp_", reductionSet$reductionOpt)]]<- imgNm;
  saveSet(infoSet);
}

# Cumulative R-squared plot for MOFA/MCIA — stacked bar chart with cumulative line
PlotCumR2 <- function(imgNm, dpi=150, format="png") {
  infoSet <- readSet(infoSet, "infoSet");
  load_cairo(); load_ggplot();
  library(data.table)
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="")

  reductionSet <- .get.rdt.set()
  if (is.null(reductionSet$reductionOpt) || is.null(reductionSet[[reductionSet$reductionOpt]])) {
    AddErrMsg("Dimension reduction results not available."); return(0)
  }
  var.exp <- reductionSet[[reductionSet$reductionOpt]]$var.exp
  if (is.null(var.exp)) { AddErrMsg("Variance data not available."); return(0) }

  sel.nms <- names(mdata.all)[mdata.all == 1]
  df <- as.data.frame(var.exp * 100)  # convert proportions to percentages
  for (i in seq_along(sel.nms)) {
    dataSet <- readDataset(sel.nms[i])
    colnames(df) <- gsub(dataSet$type, dataSet$readableType, colnames(df))
  }
  view.cols <- colnames(df)
  df$Total <- rowSums(df[, view.cols, drop=FALSE])
  df$Factor <- if (!is.null(rownames(var.exp))) rownames(var.exp) else paste0("Factor", 1:nrow(df))
  df$Factor <- gsub("Factor", "Factor ", df$Factor)
  df$Cumulative <- cumsum(df$Total)

  df_long <- melt(as.data.table(df[, c("Factor", view.cols), drop=FALSE]),
                  id.vars="Factor", variable.name="View", value.name="Variance")
  df_long$Factor <- factor(df_long$Factor, levels=df$Factor)

  p <- ggplot(df_long, aes(x=Factor, y=Variance, fill=View)) +
    geom_bar(stat="identity", width=0.6) +
    geom_line(data=df, aes(x=Factor, y=Cumulative, group=1),
              inherit.aes=FALSE, linewidth=1.2, color="black", linetype="dashed") +
    geom_point(data=df, aes(x=Factor, y=Cumulative),
               inherit.aes=FALSE, size=3, color="black") +
    geom_text(data=df, aes(x=Factor, y=Cumulative, label=sprintf("%.1f%%", Cumulative)),
              inherit.aes=FALSE, vjust=-0.8, size=3.5) +
    labs(x="", y="Variance Explained (%)", title="Cumulative Variance Explained") +
    theme_minimal(base_size=12) +
    theme(plot.title=element_text(hjust=0.5, size=13), panel.grid.major.x=element_blank())

  Cairo(file=imgNm, unit="in", dpi=dpi, width=8, height=6, type=format, bg="white")
  print(p)
  dev.off()
  infoSet$imgSet[[paste0("cumr2_", reductionSet$reductionOpt)]] <- imgNm
  saveSet(infoSet)
  return(1)
}

PlotDimredFactors <- function(meta, pc.num = 5, imgNm, dpi=150, format="png"){
  infoSet <- readSet(infoSet, "infoSet");
  load_cairo();
  load_ggplot();
  library(see)

  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");

  reductionSet <- .get.rdt.set();

  if (is.null(reductionSet$reductionOpt) || is.null(reductionSet[[reductionSet$reductionOpt]])) {
    AddErrMsg("Dimension reduction results not available. Please run the analysis first.");
    return(0);
  }

  # For MOFA/MCIA/DIABLO: plot variance explained heatmap instead of GGally ggpairs
  if (reductionSet$reductionOpt %in% c("mofa", "mcia", "diablo")) {
    sel.inx <- mdata.all == 1
    sel.nms <- names(mdata.all)[sel.inx]

    library(data.table)
    df <- as.data.frame(reductionSet[[reductionSet$reductionOpt]]$var.exp)

    # Replace internal omics type names (columns) with readable names
    for (i in 1:length(sel.nms)) {
      dataSet <- readDataset(sel.nms[i])
      colnames(df) <- gsub(dataSet$type, dataSet$readableType, colnames(df))
    }

    df$Factor <- rownames(df)
    df_long <- as.data.frame(melt(as.data.table(df), id.vars = "Factor", variable.name = "View", value.name = "Variance"))
    df_long$Variance <- df_long$Variance * 100
    df_long$Factor <- gsub("Factor", "Factor ", df_long$Factor)

    p1 <- ggplot(df_long, aes(x = Factor, y = View, fill = Variance)) +
      geom_tile(color = "grey30", linewidth = 0.8) +
      geom_text(aes(label = sprintf("%.2f%%", Variance)), size = 4, color = "black") +
      scale_fill_gradient(low = "white", high = "#C0392B", name = "Var. (%)") +
      labs(x = "", y = "", title = "Variance Explained per Factor") +
      theme_minimal(base_size = 15) +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13),
        axis.text.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5, size = 16),
        panel.grid = element_blank()
      )

    Cairo::Cairo(file = imgNm, width = 8, height = 7, type = format, bg = "white", unit = "in", dpi = dpi)
    print(p1)
    dev.off()

    infoSet$imgSet[[paste0("dimred_factors_", reductionSet$reductionOpt)]] <- imgNm
    saveSet(infoSet)
    return(1)
  }

  # For non-MOFA methods: GGally ggpairs scatter/density plot
  sel.nms <- names(mdata.all)
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readDataset(sel.nms[i])
    data.list[[i]] <- dat$data.proc
  }

  pclabels <- paste0("Component ", 1:pc.num);

  data <- as.data.frame(reductionSet[[reductionSet$reductionOpt]]$pos.xyz[,1:pc.num]);
  meta.info <- reductionSet$meta;
  meta.info <- meta.info[match(rownames(data), rownames(meta.info)),,drop=F]

  inx <- which(colnames(meta.info) == meta)
  cls <- meta.info[, inx];
  cls.type <- reductionSet$dataSet$meta.types[inx] ##### UPDATE THIS AFTER SUPPORT COMPLEX META

  if(is.null(cls.type)){
    cls.type <- "disc";
  }

  # Isolate GGally in subprocess
  ggally_result <- tryCatch({
      rsclient_isolated_exec(
        func_body = function(input_data) {
          library(GGally)
          library(ggplot2)
          library(Cairo)
          library(grid)
          library(see)
          library(RColorBrewer)

          data <- input_data$data
          cls <- input_data$cls
          cls.type <- input_data$cls.type
          pclabels <- input_data$pclabels
          imgNm <- input_data$imgNm
          dpi <- input_data$dpi
          format <- input_data$format
          pc.num <- input_data$pc.num
          base_size <- 15

          Cairo::Cairo(file = imgNm, unit = "in", dpi = dpi, width = 10, height = 10, type = format, bg = "white")

          if (cls.type == "disc") {
            p <- GGally::ggpairs(data,
                     lower = list(continuous = GGally::wrap("points")),
                     upper = list(continuous = GGally::wrap("density")),
                     diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.5, color = NA)),
                     columnLabels = pclabels, mapping = ggplot2::aes(color = cls))
            auxplot <- ggplot2::ggplot(data.frame(cls = cls), ggplot2::aes(x = cls, y = cls, color = cls)) +
              ggplot2::theme_bw(base_size = base_size) + ggplot2::geom_point(size = 6) +
              ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank(),
                             legend.text = ggplot2::element_text(size = 15)) +
              see::scale_fill_okabeito() + see::scale_color_okabeito() +
              ggplot2::guides(col = ggplot2::guide_legend(nrow = 1))
            p <- p + ggplot2::theme_bw(base_size = base_size) +
              see::scale_fill_okabeito() + see::scale_color_okabeito() +
              ggplot2::theme(plot.margin = ggplot2::unit(c(0.25, 0.25, 0.6, 0.25), "in"))
            mylegend <- GGally::grab_legend(auxplot)
          } else {
            colors <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(20))
            num.cls <- as.numeric(as.character(cls))
            cols <- colors[as.numeric(cut(num.cls, breaks = 20))]
            p <- GGally::ggpairs(data,
                     lower = list(continuous = GGally::wrap("points", color = cols)),
                     upper = list(continuous = GGally::wrap("density", color = "#505050")),
                     diag = list(continuous = GGally::wrap("densityDiag", fill = "#505050", color = NA)),
                     columnLabels = pclabels)
            auxplot <- ggplot2::ggplot(data.frame(cls = num.cls), ggplot2::aes(x = cls, y = cls, color = cls)) +
              ggplot2::theme_bw(base_size = base_size) + ggplot2::geom_point(size = 6) +
              ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank(),
                             legend.text = ggplot2::element_text(size = 15)) +
              ggplot2::guides(col = ggplot2::guide_legend(nrow = 1))
            p <- p + ggplot2::theme_bw(base_size = base_size) +
              ggplot2::theme(plot.margin = ggplot2::unit(c(0.25, 0.25, 0.8, 0.25), "in"))
            mylegend <- GGally::grab_legend(auxplot)
          }

          grid::grid.newpage()
          grid::grid.draw(p)
          vp <- grid::viewport(x = 5, y = 0.3, width = .35, height = .3, default.units = "in")
          grid::pushViewport(vp)
          grid::grid.draw(mylegend)
          grid::upViewport()
          dev.off()

          gc(verbose = FALSE)
          return(1)
        },
        input_data = list(
          data = data, cls = cls, cls.type = cls.type, pclabels = pclabels,
          imgNm = imgNm, dpi = dpi, format = format, pc.num = pc.num
        ),
        packages = c("GGally", "ggplot2", "Cairo", "grid", "see", "RColorBrewer", "qs"),
        timeout = 300,
        output_type = "qs"
      )
    }, error = function(e) {
      AddErrMsg(paste("PlotDimredFactors failed:", e$message))
      NULL
    })
    if (is.list(ggally_result) && isFALSE(ggally_result$success)) { AddErrMsg(ggally_result$message); return(0) }
    if (is.null(ggally_result)) return(0)

  infoSet$imgSet[[paste0("dimred_factors_", reductionSet$reductionOpt)]]<- imgNm;
  saveSet(infoSet);
}

# Extract BER table from perf() result - handles multiple mixOmics output formats
.extract_ber_table <- function(perf.res) {
  ber_table <- NULL
  tryCatch({
    # Try WeightedVote.error.rate first (list of matrices per distance)
    wv <- perf.res$WeightedVote.error.rate
    if (!is.null(wv) && is.list(wv)) {
      first_mat <- wv[[1]]
      if (is.matrix(first_mat)) {
        ber_table <- data.frame(Component = colnames(first_mat), stringsAsFactors = FALSE)
        n_dist <- length(wv)
        for (nm in names(wv)) {
          mat <- wv[[nm]]
          for (rn in rownames(mat)) {
            # If single distance type, use row name directly; otherwise prefix with distance name
            col_name <- if (n_dist == 1) rn else paste0(nm, ".", rn)
            ber_table[[col_name]] <- round(as.numeric(mat[rn, ]), 4)
          }
        }
      }
    }
    # Fallback: error.rate
    if (is.null(ber_table) && !is.null(perf.res$error.rate)) {
      er <- perf.res$error.rate
      if (is.list(er)) {
        first_el <- er[[1]]
        if (is.numeric(first_el) && !is.matrix(first_el)) {
          comp_names <- names(first_el)
          if (is.null(comp_names)) comp_names <- paste0("comp", seq_along(first_el))
          ber_table <- data.frame(Component = comp_names, stringsAsFactors = FALSE)
          for (nm in names(er)) ber_table[[nm]] <- round(as.numeric(er[[nm]]), 4)
        } else if (is.matrix(first_el)) {
          ber_table <- data.frame(Component = colnames(first_el), stringsAsFactors = FALSE)
          n_dist <- length(er)
          for (nm in names(er)) {
            mat <- er[[nm]]
            for (rn in rownames(mat)) {
              col_name <- if (n_dist == 1) rn else paste0(nm, ".", rn)
              ber_table[[col_name]] <- round(as.numeric(mat[rn, ]), 4)
            }
          }
        }
      }
    }
  }, error = function(e) {
    message("BER table extraction error: ", e$message)
  })
  return(ber_table)
}

# Regenerate DIABLO circos JSON with custom parameters
GenerateDiabloCircosJson <- function(cutoff=0.5, maxEdges=100) {
  cutoff <- as.numeric(cutoff)
  maxEdges <- as.integer(maxEdges)
  tryCatch(
    rsclient_isolated_exec(
      func_body = function(input_data) {
        suppressPackageStartupMessages(library(mixOmics))
        model <- qs::qread("diablo_model.qs")
        block_names <- names(model$X)
        X_proj <- lapply(model$X, function(x) x[, which(apply(x, 2, var) > 0), drop = FALSE])
        all_edges <- list()
        for (bi in 1:(length(X_proj)-1)) {
          for (bj in (bi+1):length(X_proj)) {
            cor_cross <- cor(X_proj[[bi]], X_proj[[bj]])
            sig_idx <- which(abs(cor_cross) > input_data$cutoff, arr.ind = TRUE)
            if (nrow(sig_idx) > 0) {
              for (k in 1:nrow(sig_idx)) {
                all_edges[[length(all_edges)+1]] <- list(
                  source = rownames(cor_cross)[sig_idx[k,1]],
                  target = colnames(cor_cross)[sig_idx[k,2]],
                  corr = round(cor_cross[sig_idx[k,1], sig_idx[k,2]], 4),
                  type1 = block_names[bi], type2 = block_names[bj],
                  label1 = rownames(cor_cross)[sig_idx[k,1]],
                  label2 = colnames(cor_cross)[sig_idx[k,2]])
              }
            }
          }
        }
        if (length(all_edges) > input_data$maxEdges) {
          corrs <- sapply(all_edges, function(e) abs(e$corr))
          all_edges <- all_edges[order(corrs, decreasing=TRUE)[1:input_data$maxEdges]]
        }
        jsonlite::write_json(list(DIABLO = all_edges), "diablo_circos.json", auto_unbox = TRUE, pretty = FALSE)
        TRUE
      },
      input_data = list(cutoff = cutoff, maxEdges = maxEdges),
      packages = c("mixOmics", "jsonlite", "qs"),
      timeout = 120,
      output_type = "qs"
    ),
  error = function(e) AddErrMsg(paste("[GenerateDiabloCircosJson]", e$message)))
  return(1)
}

# Plot DIABLO BER (Balanced Error Rate) diagnostic - performance vs number of components
PlotDiabloBER <- function(imgNm, dpi=150, format="png") {
  infoSet <- readSet(infoSet, "infoSet");
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");

  # Use pre-generated BER from reduce.dimension subprocess
  pregenerated <- paste0("diablo_ber", "dpi", dpi, ".", format)
  if (file.exists(pregenerated) && pregenerated != imgNm) {
    file.rename(pregenerated, imgNm)
  }

  infoSet$imgSet[["diablo_ber"]] <- imgNm
  saveSet(infoSet)
  return(1)
}

# Plot DIABLO Circos plot showing correlations between omics layers
GetBerTableRows <- function() {
  reductionSet <- .get.rdt.set()
  bt <- reductionSet[["diablo"]]$ber_table
  if (is.null(bt)) return(0)
  return(nrow(bt))
}

GetBerTableColNames <- function() {
  reductionSet <- .get.rdt.set()
  bt <- reductionSet[["diablo"]]$ber_table
  if (is.null(bt)) return("")
  # Return non-Component column names as semicolon-separated string
  cols <- colnames(bt)[colnames(bt) != "Component"]
  return(paste(cols, collapse=";"))
}

GetBerTableComp <- function(row) {
  reductionSet <- .get.rdt.set()
  bt <- reductionSet[["diablo"]]$ber_table
  return(bt$Component[as.integer(row)])
}

GetBerTableValues <- function(row) {
  reductionSet <- .get.rdt.set()
  bt <- reductionSet[["diablo"]]$ber_table
  r <- as.integer(row)
  # Return all numeric columns (everything except Component)
  cols <- colnames(bt)[colnames(bt) != "Component"]
  return(as.numeric(bt[r, cols]))
}

PlotDiabloLoading <- function(imgNm, dpi=150, format="png") {
  infoSet <- readSet(infoSet, "infoSet");
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");

  # Use pre-generated loading from reduce.dimension subprocess
  pregenerated <- paste0("diablo_loading", "dpi", dpi, ".", format)
  if (file.exists(pregenerated) && pregenerated != imgNm) {
    file.rename(pregenerated, imgNm)
  }

  infoSet$imgSet[["diablo_loading"]] <- imgNm
  saveSet(infoSet)
  return(1)
}

PlotDiabloCircos <- function(imgNm, dpi=150, format="png", cutoff=0.7) {
  infoSet <- readSet(infoSet, "infoSet")
  dpi <- as.numeric(dpi)
  cutoff <- as.numeric(cutoff)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="")

  # Read from diablo_circos.json (same data as web chord diagram)
  json_path <- "diablo_circos.json"
  if (!file.exists(json_path)) {
    GenerateDiabloCircosJson(cutoff)
  }
  edges <- tryCatch(
    jsonlite::read_json(json_path, simplifyVector = TRUE)$DIABLO,
    error = function(e) { AddErrMsg(paste("[PlotDiabloCircos] JSON read failed:", e$message)); NULL }
  )

  tryCatch({
    require(circlize)
    if (is.null(edges) || length(edges) == 0) {
      AddErrMsg("[PlotDiabloCircos] No edges in circos JSON")
    } else {
      sources <- data.frame(name = edges$source, type = edges$type1, stringsAsFactors = FALSE)
      targets <- data.frame(name = edges$target, type = edges$type2, stringsAsFactors = FALSE)
      nodes <- unique(rbind(sources, targets))
      nodes <- nodes[order(nodes$type, nodes$name), ]

      type_colors <- c("#ED7D31", "#70AD47", "#FFC000", "#9B59B6")
      type_names <- unique(nodes$type)
      sector_colors <- setNames(
        type_colors[match(nodes$type, type_names)],
        nodes$name
      )

      link_colors <- ifelse(as.numeric(edges$corr) > 0,
                            rgb(255, 80, 80, alpha = 153, maxColorValue = 255),
                            rgb(80, 80, 255, alpha = 153, maxColorValue = 255))

      Cairo::Cairo(file = imgNm, width = 10, height = 10,
                   type = format, bg = "white", unit = "in", dpi = dpi)

      circos.clear()
      n_nodes <- nrow(nodes)
      gap_size <- if (n_nodes > 80) 0.2 else if (n_nodes > 40) 0.5 else 1
      type_gap <- if (n_nodes > 80) 2 else if (n_nodes > 40) 4 else 8
      # Add larger gap between data types
      type_boundaries <- which(diff(match(nodes$type, type_names)) != 0)
      gap_vec <- rep(gap_size, n_nodes)
      if (length(type_boundaries) > 0) gap_vec[type_boundaries] <- type_gap
      gap_vec[n_nodes] <- type_gap
      circos.par(gap.after = gap_vec, start.degree = 90, cell.padding = c(0, 0, 0, 0))
      circos.initialize(factors = factor(nodes$name, levels = nodes$name),
                        xlim = c(0, 1))

      label_cex <- if (n_nodes > 80) 0.35 else if (n_nodes > 40) 0.45 else 0.6
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
        sector.name <- CELL_META$sector.index
        circos.text(CELL_META$xcenter, CELL_META$ycenter,
                    sector.name, facing = "clockwise", niceFacing = TRUE,
                    adj = c(0, 0.5), cex = label_cex)
      }, bg.col = sector_colors, bg.border = "grey50",
      track.height = 0.05)

      for (i in seq_len(nrow(edges))) {
        circos.link(edges$source[i], c(0.2, 0.8),
                    edges$target[i], c(0.2, 0.8),
                    col = link_colors[i], border = NA)
      }

      legend("bottomleft", legend = type_names, fill = type_colors[1:length(type_names)],
             title = "Data Type", bty = "n", cex = 0.9)
      legend("bottomright", legend = c("Positive", "Negative"),
             fill = c(rgb(255, 80, 80, alpha = 153, maxColorValue = 255),
                      rgb(80, 80, 255, alpha = 153, maxColorValue = 255)),
             title = "Correlation", bty = "n", cex = 0.9)

      circos.clear()
      dev.off()
    }
  }, error = function(e) { AddErrMsg(paste("[PlotDiabloCircos]", e$message)) })

  infoSet$imgSet[["diablo_circos"]] <- imgNm
  saveSet(infoSet)
  return(1)
}
