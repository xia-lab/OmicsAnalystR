##################################################
## R scripts for OmicsAnalyst
## Description: Related to correlation analysis
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#default feature selection based on sig genes

DoFeatSelectionForCorr <- function(type="default", retainedNumber=20, retainedComp=3) {
 sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
#if(!exists("mem.featSelectionForCorr")){
 #   require("memoise");
  #  mem.featSelectionForCorr <<- memoise(.do.feat.selectionForCorr);
 # }
  return(.do.feat.selectionForCorr(type,retainedNumber,retainedComp,sel.nms));
 
}


.do.feat.selectionForCorr <- function(type="default", retainedNumber=20, retainedComp=3,sel.nms){
  #print(c(type,"DoFeatSelectionForCorr"))
  sel.dats <- list();
  labels <- vector();
  reductionSet <- .get.rdt.set()
 
  if(type %in% c("default","custom","var")){
    # Pre-allocate lists to avoid sequential rbind (memory optimization)
    data.list <- vector("list", length(sel.nms))
    sig.list <- vector("list", length(sel.nms))
    dataset.list <- vector("list", length(sel.nms))
    type.vec <- vector("character", length(sel.nms))

    # First pass: collect data and sig matrices
    for(i in 1:length(sel.nms)){
      dataName = sel.nms[i];
      dataSet <- readDataset(dataName);
      data.list[[i]] <- dataSet$data.proc
      dataset.list[[i]] <- dataSet  # Store entire dataset for later use
      type.vec[i] <- dataSet$type

      if(type == "default"){
        sig.mat <- dataSet$sig.mat
      }else if(type == "var"){
        sig.mat <- dataSet$varPar$varPart.df[1:dataSet$varPar$topNum,];
        rownames(sig.mat) <-  sig.mat$ID
      }else{
        sig.mat <- dataSet$custom.sig.mat
      }

      #if more than 1000 sig features, then limit to top 1000 only
      if(nrow(sig.mat) > 1000){
        sig.mat <- sig.mat[c(1:1000),];
      }

      sig.list[[i]] <- sig.mat
    }

    # Combine all data matrices once (instead of sequential rbind)
    all.mat <- do.call(rbind, data.list)

    # Second pass: apply filtering logic
    for(i in 1:length(sel.nms)){
      dataSet <- dataset.list[[i]]  # Retrieve stored dataset
      sig.mat <- sig.list[[i]]

      inxAll = which(rownames(all.mat) %in% rownames(sig.mat));
      inx = which(rownames(data.list[[i]]) %in% rownames(sig.mat));

      all.mat <- all.mat[inxAll, ];
      dat <- data.list[[i]][inx, ];
      rownames(dat) <- paste0(rownames(dat), "_", type.vec[i]);
      sel.dats[[i]] <- dat
      labels <- c(labels, rownames(dat));
            
      if(exists("m2m",dataSet)){ 
        
        all.mat.taxa <- dataSet$data.proc.taxa

        if(type == "default"){
          sig.mat.taxa <- dataSet$sig.mat.tax
        }else{
          sig.mat <- dataSet$custom.sig.mat
        }
        inxAll <- list()
        inx <- list()
        sel.dats.taxa <- list()
        labels.taxa <- list()
        for(j in 1:length(all.mat.taxa)){
          inxAll[[j]] = which(rownames(all.mat.taxa[[j]]) %in% rownames(sig.mat.taxa[[j]]));
          inx[[j]] = which(rownames(dataSet$data.proc.taxa[[j]]) %in% rownames(sig.mat.taxa[[j]]));
          all.mat.taxa[[j]] <- all.mat.taxa[[i]][inxAll[[j]], ]
          sel.dats.taxa[[j]] <- dataSet$data.proc.taxa[[j]][inx[[j]], ]
          
        }
        names(sel.dats.taxa)<- colnames(dataSet$taxa_table)
        reductionSet$selDatsCorr.taxa <- sel.dats.taxa
        reductionSet$micidx<-i
        reductionSet$residx<-3-i
      }
      
    }
  }else{
    sel.dats <- list();
    reductionSet$corr.axis.nms <- list();
    print(sel.nms)
    for(i in 1:length(sel.nms)){
      dataName = sel.nms[i]
      dataSet <- readDataset(dataName);
      inx = which(reductionSet[[type]]$loading.pos.xyz$ids %in% rownames(dataSet$data.proc));
      loading.df <- reductionSet[[type]]$loading.pos.xyz[inx, ]
      
      if(retainedNumber > nrow(loading.df)){
        numToKeep <- nrow(loading.df);
      }else{
        numToKeep <- retainedNumber
      }
      
      for(j in 1:retainedComp){
        if(j == 1){
          loading <- loading.df[,1]
          names(loading) <- rownames(loading.df)
          loading <- loading[order(-abs(loading))]
          reductionSet$corr.axis.nms[[j]] <-names(loading)[c(1:numToKeep)]
          toKeep <- names(loading)[c(1:numToKeep)]
        }else{
          loading <- loading.df[,j]
          names(loading) <- rownames(loading.df)
          loading <- loading[order(-abs(loading))]
          reductionSet$corr.axis.nms[[j]] <-names(loading)[c(1:numToKeep)]
          toKeep <- c(toKeep, names(loading)[c(1:numToKeep)])
        }
      }
       library(stringr)
      
      # Check if all elements start with "X"
      all_start_with_x <- all(str_detect(toKeep, "^X"))
      if(all_start_with_x){
        toKeep <- substring(toKeep, 2, nchar(toKeep))
      }
      print(type)
      dat <- dataSet$data.proc
    if(type=="mofa"){
     rownames(dat) <- paste0(rownames(dat),"_",dataSet$type)
       toKeep <- paste0(toKeep, "_", dataSet$type)
        dat <- dat[rownames(dat) %in% toKeep, ]
      }else if(type=="mcia"){

      toKeep <- gsub(paste0("\\.", dataSet$type, "$"), "", toKeep)
      dat <- dat[rownames(dat) %in% toKeep, ]
      rownames(dat) <- paste0(rownames(dat), "_", dataSet$type);
    }else{
      dat <- dat[rownames(dat) %in% toKeep, ]
      rownames(dat) <- paste0(rownames(dat), "_", dataSet$type);
    }
     
      sel.dats[[i]] <- dat
    }
   
  }
 
  names(sel.dats) = sel.nms
  reductionSet$selDatsCorr <- sel.dats;
  reductionSet$feat.sel.type <- type;
   .set.rdt.set(reductionSet);
  return(1)
}


DoCorrelationFilter <- function(corSign="both", crossOmicsOnly="false", networkInfer="NA", threshold.inter=0.5,
                                threshold.intra=0.9, numToKeep=2000, updateRes="false", taxlvl="genus", datagem="agora"){

  # igraph/dplyr/reshape2 in subprocess
  reductionSet <- .get.rdt.set();
    sel.inx <- mdata.all == 1;
    sel.nms <- names(mdata.all)[sel.inx];

    # Validate correlation matrix exists
    if (is.null(reductionSet$corr.mat.path) || !file.exists(reductionSet$corr.mat.path)) {
      msg.vec <<- "Correlation matrix not found. Please run DoOmicsCorrelation first."
      return(0)
    }

    # Collect datasets
    dataSetList <- lapply(sel.nms, readDataset)
    names(dataSetList) <- sel.nms

    has_taxa <- exists("selDatsCorr.taxa", reductionSet)

    data_for_sub <- list(
      sel.nms = sel.nms,
      dataSetList = dataSetList,
      selDatsCorr = reductionSet$selDatsCorr,
      labels = reductionSet$labels,
      corr.mat.path = reductionSet$corr.mat.path,
      has_taxa = has_taxa
    )

    if (has_taxa) {
      data_for_sub$selDatsCorr.taxa <- reductionSet$selDatsCorr.taxa
      data_for_sub$corr.mat.taxa <- reductionSet$corr.mat.taxa
      data_for_sub$micidx <- reductionSet$micidx
      data_for_sub$residx <- reductionSet$residx
    }

    filter_result <- tryCatch({
      rsclient_isolated_exec(
        func_body = function(input_data) {
          require(igraph)
          require(dplyr)
          require(reshape2)
          require(qs)

          data_obj <- input_data$data_obj
          params <- input_data$params

          sel.nms <- data_obj$sel.nms
          dataSetList <- data_obj$dataSetList
          selDatsCorr <- data_obj$selDatsCorr
          labels <- data_obj$labels

          tryCatch({
            if (!data_obj$has_taxa) {
              # =========== Standard correlation filter (non-taxa path) ===========
              labels_vec <- unlist(lapply(dataSetList, function(x) x$enrich_ids))
              types <- unlist(lapply(dataSetList, function(x) rep(x$type, length(x$enrich_ids))))
              type_df <- data.frame(name = labels_vec, type = types)
              type_lookup <- setNames(type_df$type, type_df$name)

              corr.mat <- ov_qs_read(data_obj$corr.mat.path)
              corr.p.mat <- ov_qs_read("corr.p.mat.qs")
              corr.p.mat <- reshape2::melt(corr.p.mat)

              g <- igraph::graph_from_adjacency_matrix(corr.mat, mode = "undirected",
                                                        diag = FALSE, weighted = 'correlation')

              toMatch <- unlist(lapply(dataSetList, function(x) x$type))
              pattern <- paste0("^(.*)((", paste(toMatch, collapse = "|"), "))$")
              igraph::V(g)$type <- gsub(pattern, "\\2", igraph::V(g)$name)
              igraph::V(g)$label <- gsub(pattern, "", igraph::V(g)$name)

              edge_list <- igraph::as_data_frame(igraph::simplify(g, remove.loops = TRUE, edge.attr.comb = "max"), "edges")

              v1 <- igraph::V(g)$name[igraph::V(g)$type == unique(types)[1]]
              v2 <- igraph::V(g)$name[igraph::V(g)$type == unique(types)[2]]

              inter_inx <- igraph::V(g)[igraph::ends(g, igraph::E(g))[, 1]]$type != igraph::V(g)[igraph::ends(g, igraph::E(g))[, 2]]$type
              intra_inx <- igraph::V(g)[igraph::ends(g, igraph::E(g))[, 1]]$type == igraph::V(g)[igraph::ends(g, igraph::E(g))[, 2]]$type
              inter_g <- igraph::delete_edges(g, igraph::E(g)[intra_inx])
              intra_g <- igraph::delete_edges(g, igraph::E(g)[inter_inx])

              ov_qs_save(list(corr.graph.inter = inter_g, corr.graph.intra = intra_g), "corr.graph.qs")

              cor.list <- list(all = NULL, inter = NULL, intra = NULL)

              if (params$corSign == "both") {
                toRm.inter <- igraph::E(inter_g)[!abs(correlation) > params$threshold.inter]
                toRm.intra <- igraph::E(intra_g)[!abs(correlation) > params$threshold.intra]
              } else if (params$corSign == "positive") {
                toRm.inter <- igraph::E(inter_g)[!correlation > params$threshold.inter]
                toRm.intra <- igraph::E(intra_g)[!correlation > params$threshold.intra]
              } else {
                toRm.inter <- igraph::E(inter_g)[!correlation < -params$threshold.inter]
                toRm.intra <- igraph::E(intra_g)[!correlation < -params$threshold.intra]
              }

              inter_g_sub <- igraph::delete_edges(inter_g, igraph::E(inter_g)[toRm.inter])
              intra_g_sub <- igraph::delete_edges(intra_g, igraph::E(intra_g)[toRm.intra])

              edge_list_inter <- igraph::as_edgelist(inter_g_sub)
              edge_list_intra <- igraph::as_edgelist(intra_g_sub)

              cor.list$inter <- data.frame(edge_list_inter, as.numeric(igraph::E(inter_g_sub)$correlation))
              cor.list$intra <- data.frame(edge_list_intra, as.numeric(igraph::E(intra_g_sub)$correlation))
              colnames(cor.list$inter) <- c("source", "target", "correlation")
              colnames(cor.list$intra) <- c("source", "target", "correlation")
              cor.list$all <- rbind(cor.list$inter, cor.list$intra)

              ov_qs_save(cor.list, file = "cor.list.qs")

              if (params$crossOmicsOnly == "true") {
                cor_edge_list <- cor.list$inter
              } else {
                cor_edge_list <- dplyr::bind_rows(cor.list$inter, cor.list$intra)
              }

              numToKeep <- params$numToKeep
              if (numToKeep > length(unique(cor_edge_list$correlation))) {
                numToKeep <- length(unique(cor_edge_list$correlation))
              }

              if (nrow(cor_edge_list) >= 3) {
                top.edge <- sort(abs(unique(cor_edge_list$correlation)))[1:numToKeep]
                top.inx <- match(abs(cor_edge_list$correlation), top.edge)
                cor_edge_list <- cor_edge_list[!is.na(top.inx), , drop = FALSE]

                new_g <- igraph::graph_from_data_frame(cor_edge_list, directed = FALSE)
                new_g <- igraph::simplify(new_g, edge.attr.comb = "mean")

                igraph::V(new_g)$type <- gsub(pattern, "\\2", igraph::V(new_g)$name)
                toMatch2 <- unlist(lapply(dataSetList, function(x) paste0("_", x$type)))
                pattern2 <- paste0("(", paste(toMatch2, collapse = "|"), ")$")
                igraph::V(new_g)$featureId <- gsub(pattern2, "", igraph::V(new_g)$name)

                type.list <- list()
                for (i in 1:length(sel.nms)) {
                  type.list[[sel.nms[[i]]]] <- unique(cor_edge_list[, i])
                }

                cor_edge_list <- cor_edge_list %>%
                  dplyr::left_join(corr.p.mat, by = c("source" = "Var1", "target" = "Var2")) %>%
                  dplyr::mutate(pval = value)
                cor_edge_list$value <- NULL

                cor_edge_list$label1 <- gsub(paste(paste0("_", toMatch, "$"), collapse = "|"), "", cor_edge_list$source)
                cor_edge_list$label2 <- gsub(paste(paste0("_", toMatch, "$"), collapse = "|"), "", cor_edge_list$target)
                cor_edge_list$label1 <- names(labels)[match(cor_edge_list$label1, labels)]
                cor_edge_list$label2 <- names(labels)[match(cor_edge_list$label2, labels)]

                write.csv(corr.mat, "corNet.csv", row.names = FALSE)

                # ProcessGraphFile logic
                overall.graph <- new_g
                nms <- igraph::V(new_g)$name
                if (length(nms) < 1) {
                  nms <- igraph::V(new_g)$id
                  new_g <- igraph::set_vertex_attr(new_g, "name", value = nms)
                }

                lblsNm <- names(labels)
                names(lblsNm) <- unname(labels)
                lbls <- unname(lblsNm[igraph::V(new_g)$featureId])
                node.data <- data.frame(nms, lbls)
                new_g <- igraph::set_vertex_attr(new_g, "label", value = lbls)
                seed.proteins <- nms

                if (!is.null(type.list) && is.null(igraph::V(new_g)$type)) {
                  typeVec <- rep("NA", length(nms))
                  for (i in 1:length(type.list)) {
                    inx <- nms %in% type.list[[i]]
                    typeVec[inx] <- names(type.list)[i]
                  }
                  new_g <- igraph::set_vertex_attr(new_g, "type", value = typeVec)
                }

                e <- igraph::as_edgelist(new_g)
                edge.data <- data.frame(Source = e[, 1], Target = e[, 2])

                comps <- igraph::decompose(new_g, min.vertices = 3)
                if (length(comps) == 0) {
                  return(list(success = 0, msg = "No subnetworks containing at least 3 edges are identified"))
                }

                comp_sizes <- sapply(comps, function(x) igraph::vcount(x))
                comps <- comps[order(comp_sizes, decreasing = TRUE)]
                names(comps) <- paste0("subnetwork", 1:length(comps))

                net.stats <- data.frame(Node = character(length(comps)),
                                        Edge = integer(length(comps)),
                                        Query = integer(length(comps)),
                                        stringsAsFactors = FALSE)
                rownames(net.stats) <- names(comps)

                for (j in 1:length(comps)) {
                  g <- comps[[j]]
                  nd.res <- ""
                  for (i in 1:length(sel.nms)) {
                    dataSet <- dataSetList[[sel.nms[i]]]
                    lbl <- dataSet$readableType
                    if (sum(igraph::V(g)$type == dataSet$type) > 0 && !grepl(lbl, nd.res)) {
                      nd.res <- paste0(lbl, ": ", sum(igraph::V(g)$type == dataSet$type), "; ", nd.res)
                    }
                  }
                  net.stats[j, ] <- c(nd.res, igraph::ecount(g), 0)
                }

                ov_qs_save(overall.graph, "overall.graph.qs")
                ov_qs_save(comps, "ppi.comps.qs")
                ov_qs_save(node.data, "node.data.qs")
                ov_qs_save(edge.data, "edge.data.qs")
                ov_qs_save(net.stats, "net.stats.qs")

                gc(verbose = FALSE, full = TRUE)

                return(list(
                  success = 1,
                  corNet = cor_edge_list,
                  threshold.inter = params$threshold.inter,
                  threshold.intra = params$threshold.intra,
                  crossOmicsOnly = params$crossOmicsOnly,
                  taxlvl = "Feature",
                  seed.proteins = seed.proteins,
                  current.net.nm = names(comps)[1],
                  ppi.net = list(db.type = "abc", order = 1, seeds = nms, table.nm = " ",
                                 node.data = node.data, edge.data = edge.data)
                ))
              } else {
                return(list(success = 0, msg = "Less than 3 correlations have been identified using current parameters"))
              }

            } else {
              # =========== Taxa-based correlation path ===========
              selDatsCorr.taxa <- data_obj$selDatsCorr.taxa
              corr.mat.taxa <- data_obj$corr.mat.taxa
              micidx <- data_obj$micidx
              residx <- data_obj$residx

              taxlvl <- gsub("(^[[:alpha:]])", "\\U\\1", params$taxlvl, perl = TRUE)

              if (!is.null(corr.mat.taxa) && taxlvl %in% names(corr.mat.taxa)) {
                corr.mat <- corr.mat.taxa[[taxlvl]]
              } else {
                return(list(success = 0, msg = paste("Taxa correlation matrix not found for level:", taxlvl)))
              }

              sel.nms <- data_obj$sel.nms
              dataSet_mic <- dataSetList[[sel.nms[micidx]]]
              dataSet_res <- dataSetList[[sel.nms[residx]]]

              taxa_names <- unique(dataSet_mic$taxa_table[, taxlvl])
              metab_names <- dataSet_res$enrich_ids
              labels <- c(setNames(taxa_names, taxa_names), metab_names)

              mic_rows <- rownames(corr.mat) %in% names(taxa_names)
              res_rows <- !mic_rows

              corr.mat.inter <- corr.mat[mic_rows, res_rows, drop = FALSE]
              corr.mat.intra1 <- corr.mat[mic_rows, mic_rows, drop = FALSE]
              corr.mat.intra2 <- corr.mat[res_rows, res_rows, drop = FALSE]

              cor_g_inter <- igraph::graph_from_incidence_matrix(corr.mat.inter, directed = FALSE, weighted = 'correlation')
              cor_g_intra1 <- igraph::graph_from_incidence_matrix(corr.mat.intra1, directed = FALSE, weighted = 'correlation')
              cor_g_intra2 <- igraph::graph_from_incidence_matrix(corr.mat.intra2, directed = FALSE, weighted = 'correlation')

              cor_edge_list_inter1 <- igraph::as_data_frame(cor_g_inter, 'edges')
              cor_edge_list_inter1 <- cor_edge_list_inter1[!is.na(cor_edge_list_inter1$correlation), ]
              cor_edge_list_intra1 <- igraph::as_data_frame(cor_g_intra1, 'edges')
              cor_edge_list_intra1 <- cor_edge_list_intra1[!is.na(cor_edge_list_intra1$correlation), ]
              cor_edge_list_intra2 <- igraph::as_data_frame(cor_g_intra2, 'edges')
              cor_edge_list_intra2 <- cor_edge_list_intra2[!is.na(cor_edge_list_intra2$correlation), ]

              cor_edge_list_inter <- rbind(cor_edge_list_inter1, cor_edge_list_intra2)
              cor_edge_list_intra <- cor_edge_list_intra1

              cor_edge_list_inter <- cor_edge_list_inter[cor_edge_list_inter$correlation != 1, ]
              cor_edge_list_intra <- cor_edge_list_intra[cor_edge_list_intra$correlation != 1, ]

              cor.list <- list(
                all = rbind(cor_edge_list_inter, cor_edge_list_intra),
                inter = cor_edge_list_inter,
                intra = cor_edge_list_intra
              )
              ov_qs_save(cor.list, file = paste0("cor.list", taxlvl, ".qs"))

              if (params$corSign == "both") {
                cor.inx.inter <- abs(cor_edge_list_inter$correlation) > params$threshold.inter
                cor.inx.intra <- abs(cor_edge_list_intra$correlation) > params$threshold.intra
              } else if (params$corSign == "positive") {
                cor.inx.inter <- cor_edge_list_inter$correlation > params$threshold.inter
                cor.inx.intra <- cor_edge_list_intra$correlation > params$threshold.intra
              } else {
                cor.inx.inter <- cor_edge_list_inter$correlation < -params$threshold.inter
                cor.inx.intra <- cor_edge_list_intra$correlation < -params$threshold.intra
              }

              cor_edge_list_inter <- cor_edge_list_inter[cor.inx.inter, ]
              cor_edge_list_intra <- cor_edge_list_intra[cor.inx.intra, ]

              if (params$crossOmicsOnly == "true") {
                cor_edge_list <- cor_edge_list_inter
              } else {
                cor_edge_list <- dplyr::bind_rows(cor_edge_list_inter, cor_edge_list_intra)
              }

              numToKeep <- params$numToKeep
              if (numToKeep > length(unique(cor_edge_list$correlation))) {
                numToKeep <- length(unique(cor_edge_list$correlation))
              }

              top.edge <- sort(abs(unique(cor_edge_list$correlation)))[1:numToKeep]
              top.inx <- match(abs(cor_edge_list$correlation), top.edge)
              cor_edge_list <- cor_edge_list[!is.na(top.inx), , drop = FALSE]

              if (nrow(cor_edge_list) < 3) {
                return(list(success = 0, msg = paste0("Less than 3 correlations identified (inter threshold: ",
                                                       params$threshold.inter, ", intra threshold: ", params$threshold.intra, ")")))
              }

              new_g <- igraph::graph_from_data_frame(cor_edge_list, directed = FALSE)
              new_g <- igraph::simplify(new_g, edge.attr.comb = "mean")

              cor_g_inter_all <- igraph::graph_from_data_frame(cor_edge_list_inter, directed = FALSE)
              cor_g_intra_all <- igraph::graph_from_data_frame(cor_edge_list_intra, directed = FALSE)
              ov_qs_save(list(corr.graph.inter = cor_g_inter_all, corr.graph.intra = cor_g_intra_all), "corr.graph.qs")

              # ProcessGraphFile logic
              overall.graph <- new_g
              nms <- igraph::V(new_g)$name
              if (length(nms) < 1) {
                nms <- igraph::V(new_g)$id
                new_g <- igraph::set_vertex_attr(new_g, "name", value = nms)
              }

              lblsNm <- names(labels)
              names(lblsNm) <- unname(labels)
              lbls <- lblsNm[igraph::V(new_g)$name]
              node.data <- data.frame(nms, lbls)
              new_g <- igraph::set_vertex_attr(new_g, "label", value = lbls)
              seed.proteins <- nms

              e <- igraph::get.edgelist(new_g)
              edge.data <- data.frame(Source = e[, 1], Target = e[, 2])

              comps <- igraph::decompose(new_g, min.vertices = 3)
              if (length(comps) == 0) {
                return(list(success = 0, msg = "No subnetworks containing at least 3 edges are identified"))
              }

              comp_sizes <- sapply(comps, function(x) igraph::vcount(x))
              comps <- comps[order(comp_sizes, decreasing = TRUE)]
              names(comps) <- paste0("subnetwork", 1:length(comps))

              net.stats <- data.frame(Node = character(length(comps)),
                                      Edge = integer(length(comps)),
                                      Query = integer(length(comps)),
                                      stringsAsFactors = FALSE)
              rownames(net.stats) <- names(comps)

              for (j in 1:length(comps)) {
                g <- comps[[j]]
                nd.res <- ""
                lbl_taxa <- taxlvl
                taxa_count <- sum(grepl(paste0("^", dataSet_mic$type), igraph::V(g)$name))
                if (taxa_count > 0 && !grepl(lbl_taxa, nd.res)) {
                  nd.res <- paste0(lbl_taxa, ": ", taxa_count, "; ", nd.res)
                }
                lbl_metab <- "Metabolite"
                metab_count <- sum(grepl(paste0("^", dataSet_res$type), igraph::V(g)$name))
                if (metab_count > 0 && !grepl(lbl_metab, nd.res)) {
                  nd.res <- paste0(lbl_metab, ": ", metab_count, "; ", nd.res)
                }
                net.stats[j, ] <- c(nd.res, igraph::ecount(g), 0)
              }

              ov_qs_save(overall.graph, "overall.graph.qs")
              ov_qs_save(comps, "ppi.comps.qs")
              ov_qs_save(node.data, "node.data.qs")
              ov_qs_save(edge.data, "edge.data.qs")
              ov_qs_save(net.stats, "net.stats.qs")

              gc(verbose = FALSE, full = TRUE)

              return(list(
                success = 1,
                taxlvl = taxlvl,
                datagem = params$datagem,
                threshold.inter = params$threshold.inter,
                threshold.intra = params$threshold.intra,
                crossOmicsOnly = params$crossOmicsOnly,
                seed.proteins = seed.proteins,
                current.net.nm = names(comps)[1],
                ppi.net = list(db.type = "abc", order = 1, seeds = nms, table.nm = " ",
                               node.data = node.data, edge.data = edge.data)
              ))
            }

          }, error = function(e) {
            return(list(success = -1, msg = paste("Correlation filter failed:", e$message)))
          })
        },
        input_data = list(
          data_obj = data_for_sub,
          params = list(
            corSign = corSign,
            crossOmicsOnly = crossOmicsOnly,
            networkInfer = networkInfer,
            threshold.inter = as.numeric(threshold.inter),
            threshold.intra = as.numeric(threshold.intra),
            numToKeep = as.numeric(numToKeep),
            updateRes = updateRes,
            taxlvl = taxlvl,
            datagem = datagem
          )
        ),
        packages = c("igraph", "dplyr", "reshape2", "qs"),
        timeout = 600
      )
    }, error = function(e) {
      msg.vec <<- paste("Correlation filter failed:", e$message)
      NULL
    })
    if (is.list(filter_result) && isFALSE(filter_result$success)) { AddErrMsg(filter_result$message); return(0) }

    if (is.null(filter_result)) return(0)

    if (filter_result$success == 1) {
      # Restore state to Master session
      reductionSet$corr.graph.path <- "corr.graph.qs"
      reductionSet$cor.list.path <- "cor.list.qs"
      reductionSet$threshold.inter <- filter_result$threshold.inter
      reductionSet$threshold.intra <- filter_result$threshold.intra
      reductionSet$crossOmicsOnly <- filter_result$crossOmicsOnly
      reductionSet$taxlvl <- filter_result$taxlvl

      if (!is.null(filter_result$datagem)) {
        reductionSet$datagem <- filter_result$datagem
      }
      if (!is.null(filter_result$corNet)) {
        reductionSet$corNet <- filter_result$corNet
      }

      .set.rdt.set(reductionSet)

      # Restore globals from qs files
      overall.graph <<- ov_qs_read("overall.graph.qs")
      ppi.comps <<- ov_qs_read("ppi.comps.qs")
      net.stats <<- ov_qs_read("net.stats.qs")
      seed.proteins <<- filter_result$seed.proteins
      seed.genes <<- filter_result$seed.proteins
      seed.expr <<- rep(0, nrow(filter_result$ppi.net$node.data))
      current.net.nm <<- filter_result$current.net.nm
      net.nmu <<- filter_result$current.net.nm
      ppi.net <<- filter_result$ppi.net
      data.idType <<- "NA"

      return(1)

    } else if (filter_result$success == 0) {
      msg.vec <<- filter_result$msg
      return(0)

    } else {
      msg.vec <<- filter_result$msg
      return(0)
    }

}

GenerateNetworkJson <- function(fileName="omicsanalyst_net_0.json"){
    intres <- convertIgraph2JSON(current.net.nm , fileName);
    return(intres);
}

ExportOmicsPairs <- function(fileName, type){
  
  cor.list <- ov_qs_read("cor.list.qs");
  cor.obj <- cor.list[[type]];
  colnames(cor.obj) = c("Id1", "Id2", "Correlation");
  write.table(cor.obj,row.name=F, file=fileName);
  
}


DoOmicsCorrelation <- function(cor.method="univariate",cor.stat="pearson",ifAll="true",metaSel,group){
  labels <- vector();
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];

  m2midx<-0
  for(i in 1:length(sel.nms)){
    dataName = sel.nms[i]
    dataSet <- readDataset(dataName);
    labels <- c(labels, dataSet$enrich_ids);
    if(exists("m2m",dataSet)){
      labels.taxa <- lapply(dataSet$data.taxa, function(x) rownames(x))
      labels.taxa <-  lapply(labels.taxa, function(x) setNames(x,x))
    }else if(exists("labels.taxa")){
      labels.taxa <- lapply(labels.taxa, function(x) c(x,dataSet$enrich_ids))
    }
  }

  reductionSet <- .get.rdt.set();
  reductionSet$cor.stat <- cor.stat;
  reductionSet$cor.method <- cor.method;
  sel.dats <- reductionSet$selDatsCorr;
  load_igraph();

  if(ifAll != "true"){
    meta.info = reductionSet$dataSet$meta.info
    sampleInclude = row.names(meta.info[which(meta.info[[metaSel]]==group),])

    if(length(sampleInclude)<4){
         return(-1)
    }else{
     sel.dats <- lapply(sel.dats,function(x){
         return(x[,sampleInclude])
     })
   }
  }

  residx <- reductionSet$residx
  selDatsCorr.taxa <- if (exists("selDatsCorr.taxa", reductionSet)) reductionSet$selDatsCorr.taxa else NULL

  # Hmisc/parmigene/ppcor in subprocess
  corr_result <- tryCatch({
      rsclient_isolated_exec(
        func_body = function(input_data) {
          sel.dats <- input_data$sel.dats
          selDatsCorr.taxa <- input_data$selDatsCorr.taxa
          cor.method <- input_data$cor.method
          cor.stat <- input_data$cor.stat
          residx <- input_data$residx

          corr.mat.taxa <- NULL
          corr.p.mat <- NULL

          if (cor.method == "univariate") {
            require(Hmisc)

            if (!is.null(selDatsCorr.taxa)) {
              corr.mat.taxa <- lapply(selDatsCorr.taxa, function(x) {
                transposed_taxa <- cbind(t(x), t(sel.dats[[residx]]))
                result <- cor(transposed_taxa, method = cor.stat)
                rm(transposed_taxa)
                return(result)
              })
            }

            transposed_data <- cbind(t(sel.dats[[1]]), t(sel.dats[[2]]))

            if (tolower(cor.stat) == "kendall") {
              corr.mat <- cor(transposed_data, method = "kendall", use = "pairwise.complete.obs")
              n_vars <- ncol(transposed_data)
              corr.p.mat <- matrix(NA, nrow = n_vars, ncol = n_vars)
              rownames(corr.p.mat) <- colnames(corr.p.mat) <- colnames(corr.mat)
              for (i in 1:(n_vars - 1)) {
                for (j in (i + 1):n_vars) {
                  valid_idx <- complete.cases(transposed_data[, i], transposed_data[, j])
                  n_valid <- sum(valid_idx)
                  tau <- corr.mat[i, j]
                  z <- tau / sqrt((2 * (2 * n_valid + 5)) / (9 * n_valid * (n_valid - 1)))
                  corr.p.mat[i, j] <- corr.p.mat[j, i] <- 2 * pnorm(-abs(z))
                }
              }
              diag(corr.p.mat) <- 0
            } else {
              res <- Hmisc::rcorr(as.matrix(transposed_data), type = cor.stat)
              corr.mat <- res$r
              corr.p.mat <- res$P
              rm(res)
            }
            rm(transposed_data)

          } else if (cor.method == "MI") {
            require(parmigene)
            combined_data <- rbind(sel.dats[[1]], sel.dats[[2]])
            res <- parmigene::knnmi.all(combined_data, k = 5)
            rm(combined_data)
            scale <- 1 / max(res)
            corr.mat <- res * scale
            rm(res, scale)

            if (!is.null(selDatsCorr.taxa)) {
              corr.mat.taxa <- list()
              for (j in 1:length(selDatsCorr.taxa)) {
                combined_taxa <- rbind(selDatsCorr.taxa[[j]], sel.dats[[residx]])
                res <- parmigene::knnmi.all(combined_taxa, k = 5)
                rm(combined_taxa)
                scale <- 1 / max(res)
                corr.mat.taxa[[j]] <- res * scale
                rm(res, scale)
              }
            }

          } else {
            require(ppcor)
            sel.res <- cbind(t(sel.dats[[1]]), t(sel.dats[[2]]))
            res <- ppcor::pcor(sel.res, method = cor.stat)
            corr.mat <- res$estimate
            corr.p.mat <- res$p.value
            rm(res)
            rownames(corr.mat) <- colnames(corr.mat) <- colnames(sel.res)
            rownames(corr.p.mat) <- colnames(corr.p.mat) <- colnames(sel.res)
            rm(sel.res)

            if (!is.null(selDatsCorr.taxa)) {
              sel.res <- lapply(selDatsCorr.taxa, function(x) {
                cbind(t(x), t(sel.dats[[residx]]))
              })
              res <- lapply(sel.res, function(x) { ppcor::pcor(x, method = cor.stat) })
              corr.mat.taxa <- lapply(res, function(x) x$estimate)
              for (j in 1:length(corr.mat.taxa)) {
                rownames(corr.mat.taxa[[j]]) <- colnames(sel.res[[j]])
                colnames(corr.mat.taxa[[j]]) <- colnames(sel.res[[j]])
              }
              rm(sel.res, res)
            }
          }

          gc(verbose = FALSE, full = TRUE)
          list(corr.mat = corr.mat, corr.p.mat = corr.p.mat, corr.mat.taxa = corr.mat.taxa)
        },
        input_data = list(
          sel.dats = sel.dats,
          selDatsCorr.taxa = selDatsCorr.taxa,
          cor.method = cor.method,
          cor.stat = cor.stat,
          residx = residx
        ),
        packages = if (cor.method == "univariate") c("Hmisc")
                   else if (cor.method == "MI") c("parmigene")
                   else c("ppcor"),
        timeout = 300
      )
    }, error = function(e) {
      msg.vec <<- paste("Correlation analysis failed:", e$message)
      NULL
    })
    if (is.list(corr_result) && isFALSE(corr_result$success)) { AddErrMsg(corr_result$message); return(0) }

    if (is.null(corr_result)) return(0)

    # Save results to disk
    ov_qs_save(corr_result$corr.mat, "corr.mat.qs")
    if (!is.null(corr_result$corr.p.mat)) {
      ov_qs_save(corr_result$corr.p.mat, "corr.p.mat.qs")
    }
    if (!is.null(corr_result$corr.mat.taxa)) {
      reductionSet$corr.mat.taxa <- corr_result$corr.mat.taxa
    }
    rm(corr_result)
  reductionSet$labels <- labels
  reductionSet$corr.mat.path <- "corr.mat.qs"
  rm(labels)
  gc(verbose = FALSE)

  .set.rdt.set(reductionSet);

  return(1);
}

PlotCorrViolin <- function(imgNm, dpi=150, format="png", corNetOpt="default"){
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");

  reductionSet <- .get.rdt.set();

  # ggpubr/ggplot2/Cairo/igraph in subprocess
  # Prepare data for subprocess
    if (corNetOpt == "intLim") {
      intLim_filtres <- reductionSet$intLim_filtres
      intLim_sigmat <- reductionSet$intLim_sigmat
      intLim_pvalcutoff <- reductionSet$intLim$pvalcutoff
      graphs <- NULL
      threshold_inter <- NULL
      threshold_intra <- NULL
    } else {
      intLim_filtres <- NULL
      intLim_sigmat <- NULL
      intLim_pvalcutoff <- NULL
      if (!is.null(reductionSet$corr.graph.path) && file.exists(reductionSet$corr.graph.path)) {
        graphs <- ov_qs_read(reductionSet$corr.graph.path)
      } else {
        msg.vec <<- "No correlation graph available for violin plot"
        return(0)
      }
      threshold_inter <- reductionSet$threshold.inter
      threshold_intra <- reductionSet$threshold.intra
    }

    plot_result <- tryCatch({
      rsclient_isolated_exec(
        func_body = function(input_data) {
          require(ggpubr)
          require(ggplot2)
          require(Cairo)
          require(scales)
          require(igraph)

          imgNm <- input_data$imgNm
          dpi <- input_data$dpi
          format <- input_data$format
          corNetOpt <- input_data$corNetOpt
          intLim_filtres <- input_data$intLim_filtres
          intLim_sigmat <- input_data$intLim_sigmat
          intLim_pvalcutoff <- input_data$intLim_pvalcutoff
          graphs <- input_data$graphs
          threshold_inter <- input_data$threshold_inter
          threshold_intra <- input_data$threshold_intra

          fig.list <- list()

          if (corNetOpt == "intLim") {
            df_res <- intLim_filtres[, c(1, 2, 3, 5)]
            colnames(df_res) <- c("source", "target", "correlation", "pval")
            df_res <- df_res[rownames(df_res) %in% rownames(intLim_sigmat), ]
            df_res$type <- "Between-omics coefficient"
            threshold1 <- min(df_res$correlation[df_res$correlation > 0])
            threshold2 <- max(df_res$correlation[df_res$correlation < 0])
            sig.pos <- scales::oob_censor(df_res$correlation, c(threshold1, max(df_res$correlation)))
            sig.neg <- scales::oob_censor(df_res$correlation, c(min(df_res$correlation), threshold2))
            sig.pos.num <- length(na.omit(sig.pos))
            sig.neg.num <- length(na.omit(sig.neg))

            fig.list[[1]] <- ggplot2::ggplot(df_res, ggplot2::aes(x = type, y = correlation, fill = type)) +
              ggplot2::geom_violin(trim = FALSE, fill = "#d3d3d3", show.legend = FALSE) +
              ggplot2::labs(x = "Between-omics coefficient") +
              ggplot2::labs(y = paste0("coefficient (p-value <", intLim_pvalcutoff, ")")) +
              ggplot2::geom_jitter(height = 0, width = 0.05, alpha = 0, show.legend = FALSE) +
              ggplot2::theme(legend.position = "none") +
              ggplot2::scale_x_discrete(labels = NULL) +
              ggplot2::scale_y_continuous(limits = c(min(df_res$correlation), max(df_res$correlation))) +
              ggplot2::geom_hline(yintercept = threshold1, linetype = "dashed", color = "red", size = 0.5) +
              ggplot2::annotate("text", x = 0.5, y = threshold1, label = sig.pos.num, vjust = -1) +
              ggplot2::geom_hline(yintercept = threshold2, linetype = "dashed", color = "red", size = 0.5) +
              ggplot2::annotate("text", x = 0.5, y = threshold2, label = sig.neg.num, vjust = 1.5) +
              ggplot2::theme_bw() +
              ggplot2::theme(text = ggplot2::element_text(size = 13),
                             plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
                             panel.grid.minor = ggplot2::element_blank(),
                             panel.grid.major = ggplot2::element_blank())
            fig.list[[2]] <- ggplot2::ggplot() + ggplot2::theme_void()

          } else {
            for (i in 1:2) {
              if (i == 1) {
                g <- graphs$corr.graph.inter
                titleText <- "Between-omics correlation"
                threshold <- threshold_inter
              } else {
                g <- graphs$corr.graph.intra
                titleText <- "Intra-omics correlation"
                threshold <- threshold_intra
              }

              df_res <- data.frame(igraph::as_edgelist(g), as.numeric(igraph::E(g)$correlation))
              df_res <- df_res[!duplicated(df_res), ]
              colnames(df_res) <- c("source", "target", "correlation")
              df_res$type <- titleText

              sig.pos <- scales::oob_censor(df_res$correlation, c(threshold, 1))
              sig.neg <- scales::oob_censor(df_res$correlation, c(-1, -threshold))
              sig.pos.num <- length(na.omit(sig.pos))
              sig.neg.num <- length(na.omit(sig.neg))

              fig.list[[i]] <- ggplot2::ggplot(df_res, ggplot2::aes(x = type, y = correlation, fill = type)) +
                ggplot2::geom_violin(trim = FALSE, fill = "#d3d3d3", show.legend = FALSE) +
                ggplot2::labs(x = titleText) +
                ggplot2::geom_jitter(height = 0, width = 0.05, alpha = 0, show.legend = FALSE) +
                ggplot2::theme(legend.position = "none") +
                ggplot2::scale_x_discrete(labels = NULL) +
                ggplot2::scale_y_continuous(limits = c(-1, 1)) +
                ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "red", size = 0.5) +
                ggplot2::annotate("text", x = 0.5, y = threshold, label = sig.pos.num, vjust = -1) +
                ggplot2::geom_hline(yintercept = -threshold, linetype = "dashed", color = "red", size = 0.5) +
                ggplot2::annotate("text", x = 0.5, y = -threshold, label = sig.neg.num, vjust = 1.5) +
                ggplot2::theme_bw() +
                ggplot2::theme(text = ggplot2::element_text(size = 13),
                               plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
                               panel.grid.minor = ggplot2::element_blank(),
                               panel.grid.major = ggplot2::element_blank())
            }
          }

          Cairo::Cairo(file = imgNm, width = 10, height = 8, unit = "in", type = "png", bg = "white", dpi = dpi)
          p1 <- ggpubr::ggarrange(plotlist = fig.list, ncol = 2)
          print(p1)
          dev.off()

          gc(verbose = FALSE, full = TRUE)
          list(imgNm = imgNm, success = 1)
        },
        input_data = list(
          imgNm = imgNm,
          dpi = dpi,
          format = format,
          corNetOpt = corNetOpt,
          intLim_filtres = intLim_filtres,
          intLim_sigmat = intLim_sigmat,
          intLim_pvalcutoff = intLim_pvalcutoff,
          graphs = graphs,
          threshold_inter = threshold_inter,
          threshold_intra = threshold_intra
        ),
        packages = c("ggpubr", "ggplot2", "Cairo", "scales", "igraph"),
        timeout = 300
      )
    }, error = function(e) {
      msg.vec <<- paste("PlotCorrViolin failed:", e$message)
      NULL
    })
    if (is.list(plot_result) && isFALSE(plot_result$success)) { AddErrMsg(plot_result$message); return(0) }

    if (is.null(plot_result)) return(0)

    infoSet <- readSet(infoSet, "infoSet");
    infoSet$imgSet$correlation_distribution <- plot_result$imgNm;
    saveSet(infoSet);
    return(plot_result$success);
}


PlotDegreeHistogram <- function(imgNm, netNm = "NA", dpi=150, format="png"){
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");

  # igraph degree + ggplot in subprocess
  graph_to_plot <- overall.graph;
    input_data <- list(
      graph = graph_to_plot,
      imgNm = imgNm
    )
    isolated_func <- function(input_data) {
      library(igraph)
      library(ggplot2)
      library(Cairo)
      graph <- input_data$graph
      imgNm <- input_data$imgNm
      Cairo(file=imgNm, width=400, height=400, type="png", bg="white")
      G.degrees <- igraph::degree(graph)
      G.degree.histogram <- as.data.frame(table(G.degrees))
      G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])
      p <- ggplot(G.degree.histogram, aes(x = G.degrees, y = Freq)) +
        geom_point() +
        scale_x_continuous("Degree\n(nodes containing that amount of connections)",
                           breaks = c(1, 3, 10, 30, 100, 300),
                           trans = "log10") +
        scale_y_continuous("Frequency\n(number of nodes)",
                           breaks = c(1, 3, 10, 30, 100, 300, 1000),
                           trans = "log10") +
        ggtitle("Degree Distribution (log-log)") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      print(p)
      dev.off()
      return(list(success = 1))
    }
    tryCatch({
      rsclient_isolated_exec(
        func_body = isolated_func,
        input_data = input_data,
        packages = c("igraph", "ggplot2", "Cairo"),
        timeout = 300
      )
    }, error = function(e) {
      AddErrMsg(paste("PlotDegreeHistogram failed:", e$message))
      return(0)
    })
    # plot/write failure is non-fatal

  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$degree_distribution <- imgNm;
  saveSet(infoSet);
}

PlotBetweennessHistogram <- function(imgNm, netNm = "NA", dpi=150, format="png"){
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");

  if(netNm != "NA"){
    graph_to_plot <- ppi.comps[[netNm]];
  } else {
    graph_to_plot <- overall.graph;
  }

  # igraph betweenness + ggplot in subprocess
  input_data <- list(
      graph = graph_to_plot,
      imgNm = imgNm
    )
    isolated_func <- function(input_data) {
      library(igraph)
      library(ggplot2)
      library(Cairo)
      graph <- input_data$graph
      imgNm <- input_data$imgNm
      Cairo(file=imgNm, width=400, height=400, type="png", bg="white")
      G.degrees <- igraph::betweenness(graph)
      G.degree.histogram <- as.data.frame(table(G.degrees))
      G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])
      p <- ggplot(G.degree.histogram, aes(x = G.degrees, y = Freq)) +
        geom_point() +
        scale_x_continuous("Betweenness\n(nodes with that amount of betweenness)",
                           breaks = c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000),
                           trans = "log10") +
        scale_y_continuous("Frequency\n(number of nodes)",
                           breaks = c(1, 3, 10, 30, 100, 300, 1000),
                           trans = "log10") +
        ggtitle("Betweenness Distribution (log-log)") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      print(p)
      dev.off()
      return(list(success = 1))
    }
    tryCatch({
      rsclient_isolated_exec(
        func_body = isolated_func,
        input_data = input_data,
        packages = c("igraph", "ggplot2", "Cairo"),
        timeout = 300
      )
    }, error = function(e) {
      AddErrMsg(paste("PlotBetweennessHistogram failed:", e$message))
      return(0)
    })
    # plot/write failure is non-fatal

  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$betweenness_distribution <- imgNm;
  saveSet(infoSet);
}

GetNetworkTopology <- function(netnm){
  g <- ppi.comps[[netnm]];
  globalProperties <-list();
  globalProperties[["Diameter"]] <-diameter(g);
  globalProperties[["Radius"]] <-radius(g);
  globalProperties[["Average path length"]] <-signif(mean_distance(g), 3);
  globalProperties[["Clustering coefficient"]] <- transitivity(g, type="global");
  propertiesVector <- c(globalProperties[[1]], globalProperties[[2]], globalProperties[[3]], globalProperties[[4]]);
  return(propertiesVector)
}

#expand non-square adjancy matrix to square
expand.matrix <- function(A){
  m <- nrow(A)
  n <- ncol(A)
  B <- matrix(0,nrow = m, ncol = m)
  C <- matrix(0,nrow = n, ncol = n)
  cbind(rbind(B,t(A)),rbind(A,C))
}



###########generate chordgram
GenerateChordGram <- function(thresh=0.5,maxN,pval,imgName = "chordgram", format = "png", dpi = 300){
  #print(c(maxN,pval))
  plotjs <- paste0(imgName, ".json");
  reductionSet <- .get.rdt.set();

  # Load correlation matrices from disk
  corr.mat <- ov_qs_read("corr.mat.qs");
  corr.p.mat<- ov_qs_read("corr.p.mat.qs");

  # Convert to long format for filtering
  corr.mat.melted <- reshape2::melt(corr.mat)
  corr.p.mat.melted <- reshape2::melt(corr.p.mat)

  # CRITICAL: Free original large matrices immediately after melting
  # Original matrices: ~90 MB per user
  # Melted format: ~same size but will be filtered down
  rm(corr.mat, corr.p.mat)
  gc(verbose = FALSE)

  # Add p-values to correlation data
  corr.mat.melted$pval <- corr.p.mat.melted$value
  rm(corr.p.mat.melted)  # Clean up p-value melted data

  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];

  sel.dats <- reductionSet$selDatsCorr[sel.nms];

  # Filter to significant correlations only (reduces from ~90 MB to ~1-2 MB)
  corr.mat <- corr.mat.melted[abs(corr.mat.melted$value)>thresh & corr.mat.melted$pval < pval,]
  rm(corr.mat.melted)  # Clean up unfiltered data
  gc(verbose = FALSE)

  corr.mat <- corr.mat[corr.mat$Var1 %in% rownames(sel.dats[[1]]) & corr.mat$Var2 %in% rownames(sel.dats[[2]]),]
  corr.mat <- corr.mat[order(abs(corr.mat$value),decreasing = T),]

  if(nrow(corr.mat)>maxN){
    corr.mat <- corr.mat[1:maxN,]
  }

  dataSetList <- lapply(sel.nms, readDataset);
 
corr.mat$Var1 <- gsub(paste0("_", dataSetList[[1]]$type), "", corr.mat$Var1);
corr.mat$Var1 <- names(dataSetList[[1]]$enrich_ids)[match(corr.mat$Var1,dataSetList[[1]]$enrich_ids)]

corr.mat$Var2 <- gsub(paste0("_", dataSetList[[2]]$type), "", corr.mat$Var2);
corr.mat$Var2 <- names(dataSetList[[2]]$enrich_ids)[match(corr.mat$Var2,dataSetList[[2]]$enrich_ids)]


  reductionSet$chordGram <- corr.mat
  write.csv(corr.mat,"chord_diagram.csv",row.names=F)
  library(jsonlite)
  write_json(corr.mat,plotjs, pretty = TRUE)
  .set.rdt.set(reductionSet);
  return(nrow(corr.mat))
}


###############Functions for generating differential chord diagram
###############
DoOmicsDiffCorrelation <- function(cor.method="univariate",cor.stat="pearson",comp.meta="Diagnosis",selnm1,selnm2){
  return(mem.diffCorr(cor.method,cor.stat,comp.meta, selnm1,selnm2));
}


.do.omics.diffcorr <- function(cor.method="univariate",cor.stat="pearson",comp.meta="Diagnosis",selnm1,selnm2){
  labels <- vector();
  sel.nms <- c(selnm1,selnm2);
   
  m2midx<-0
  for(i in 1:length(sel.nms)){
    dataName = sel.nms[i]
    dataSet <- readDataset(dataName);
    labels <- c(labels, dataSet$enrich_ids);
    if(exists("m2m",dataSet)){
      labels.taxa <- lapply(dataSet$data.taxa, function(x) rownames(x))
      labels.taxa <-  lapply(labels.taxa, function(x) setNames(x,x))
    }else if(exists("labels.taxa")){
      labels.taxa <- lapply(labels.taxa, function(x) c(x,dataSet$enrich_ids))
    }
  }
  
  reductionSet <- .get.rdt.set();
  reductionSet$cor.stat <- cor.stat;
  reductionSet$cor.method <- cor.method;
  sel.dats <- reductionSet$selDatsCorr;
  sel.dats <- sel.dats[names(sel.dats) %in% sel.nms]
  sel.dats <- lapply(sel.dats, function(x){
    if(nrow(x)>500){
      return(x[1:500,])
    }else{
      return(x)
    }
  })
  load_igraph();
 
  residx <- reductionSet$residx
  
  meta.info =  reductionSet$dataSet$meta.info
  corr.ls <- split(rownames(meta.info),meta.info[[comp.meta]])  

  corr.ls <- corr.ls[unlist(lapply(corr.ls, length))>4]
  
  
 
  if(cor.method == "univariate"){
    library(Hmisc)
    # Optimize: transpose once per sample group instead of repeatedly
    corr.mat.ls <- lapply(corr.ls, function(samp){
      transposed_samp <- cbind(t(sel.dats[[1]][,samp]), t(sel.dats[[2]][,samp]))

      # rcorr only supports pearson and spearman
      if(tolower(cor.stat) == "kendall") {
        # For Kendall, compute correlation and p-values separately
        r_mat <- cor(transposed_samp, method="kendall", use="pairwise.complete.obs")

        # Compute p-values
        n_samples <- nrow(transposed_samp)
        n_vars <- ncol(transposed_samp)
        p_mat <- matrix(NA, nrow=n_vars, ncol=n_vars)
        rownames(p_mat) <- colnames(p_mat) <- colnames(r_mat)

        # Compute p-values for each pair
        for(i in 1:(n_vars-1)) {
          for(j in (i+1):n_vars) {
            valid_idx <- complete.cases(transposed_samp[,i], transposed_samp[,j])
            n_valid <- sum(valid_idx)
            tau <- r_mat[i,j]
            z <- tau / sqrt((2*(2*n_valid+5))/(9*n_valid*(n_valid-1)))
            p_mat[i,j] <- p_mat[j,i] <- 2 * pnorm(-abs(z))
          }
        }
        diag(p_mat) <- 0

        # Create result object compatible with rcorr structure
        res <- list(r = r_mat, n = n_samples, P = p_mat)
        class(res) <- "rcorr"
      } else {
        res = rcorr(as.matrix(transposed_samp), type = cor.stat);
      }

      rm(transposed_samp)  # Clean up transpose immediately
      return(res)
    })

    if(exists("selDatsCorr.taxa",reductionSet)){
      corr.mat.taxa <- lapply(reductionSet$selDatsCorr.taxa, function(x){
        # Transpose once and reuse
        transposed_taxa <- cbind(t(x), t(sel.dats[[residx]]))
        result <- cor(transposed_taxa, method=cor.stat)
        rm(transposed_taxa)
        return(result)
      })

      reductionSet$corr.mat.taxa <- corr.mat.taxa
    }

  }

  reductionSet$diffnet.mat.path <- "diffnet.mat.qs"

  # Save differential network matrices to disk
  ov_qs_save(corr.mat.ls, "diffnet.mat.qs");

  # CRITICAL: Free correlation list from memory after saving
  rm(corr.mat.ls)
  gc(verbose = FALSE)

  .set.rdt.set(reductionSet);
  return(1);
}

require("memoise");
mem.diffCorr <<- memoise(.do.omics.diffcorr);

GenerateDiffNet <- function(corr_thresh=0.7,p_thresh=0.05,imgName = "diffnet", format = "png", dpi = 300,dt1,dt2,topN = 100,layout="kk"){

  plotjs <- paste0(imgName, ".json");
  reductionSet <- .get.rdt.set();
  sel.dats <- reductionSet$selDatsCorr[c(dt1,dt2)];
  dataSetList <- lapply(c(dt1,dt2), readDataset);
  type1 <- dataSetList[[1]][["type"]]
  type2 <- dataSetList[[2]][["type"]]
  corr.mat.ls <- ov_qs_read("diffnet.mat.qs");
  corr.mat.ls <- lapply(corr.mat.ls, function(x){
    corr.mat<-reshape2::melt(x$r)
    p.mat<-reshape2::melt(x$P)
    corr.mat$pval <- p.mat$value
    corr.mat <- corr.mat[corr.mat$Var1 %in% rownames(sel.dats[[1]]) & corr.mat$Var2 %in% rownames(sel.dats[[2]]),]
    #corr.mat <- corr.mat[abs(corr.mat$value)>thresh,]
    return( corr.mat)
  })
  
  sig.idx <- unique(unlist(lapply(corr.mat.ls, function(x) return(which(abs(x$value)>corr_thresh&x$pval< p_thresh)) )))
  corr.mat.ls <- lapply(corr.mat.ls, function(x){
    x <- x[sig.idx,]
   # x$value[x$value<thresh] <-0
    names(x) <-c("source","target","corr","pval")
    x$weight <- abs(x$corr)
    x$source = as.character(x$source)
    x$target = as.character(x$target)
    return(x)
  })

corr.mat.ls <- lapply(corr.mat.ls,function(x){
   df=x
   df$type1 = type1
   df$type2 = type2
   df$label1 = gsub(paste0("_",type1,"$"),"",df$source)
   df$label2 = gsub(paste0("_",type2,"$"),"",df$target)
   df$label1 = names(dataSetList[[1]]$enrich_ids)[match(df$label1,dataSetList[[1]]$enrich_ids)]
   df$label2 = names(dataSetList[[2]]$enrich_ids)[match(df$label2,dataSetList[[2]]$enrich_ids)]
   df <- df[order(-df$weight,df$pval),]
   df <- df[which(df$weight>corr_thresh &df$pval < p_thresh), ] 
   if(nrow(df)>topN){
   df <- df[1:topN,]
   }

   return(df)
})
   
  rm = which(lapply(corr.mat.ls,nrow)==0)
  if(length(rm)>0){
corr.mat.ls <- corr.mat.ls[-rm]
  }
  reductionSet$diffList <- corr.mat.ls 
  .set.rdt.set(reductionSet);
  library(jsonlite)
  write_json(corr.mat.ls,plotjs, pretty = TRUE)
  return(1);


  library(igraph)
  library(ggplot2)
  library(ggraph)
  nodes <- data.frame(
    id=1:length( unique(c(corr.mat.ls[[1]]$from, corr.mat.ls[[1]]$to))),
   label = unique(c(corr.mat.ls[[1]]$from, corr.mat.ls[[1]]$to))
   )
  nodes$type <- ifelse(nodes$label %in% rownames(sel.dats[[1]]),type1,type2)
  
  
  edges.ls <-lapply(corr.mat.ls, function(x){
    x$from <- nodes$id[match(x$from,nodes$label)]
    x$to <- nodes$id[match(x$to,nodes$label)]
    return(x)
  })
  edges = edges.ls[[1]] 
  nodes$label <- gsub(type1,"",  nodes$label)
  nodes$label <- gsub(type2,"",  nodes$label)
  nodes$label <- gsub("_$","",  nodes$label)
  graphall <- graph_from_data_frame(edges[edges$weight>0,], vertices = nodes, directed = FALSE)

  min_edge_length = 1
   layout_fixed <- layout_with_kk(graphall)
  #layout_fixed <- layout.norm(layout_fixed, xmin = -min_edge_length, xmax = min_edge_length,  ymin = -min_edge_length, ymax = min_edge_length)
  unique_types <- unique(nodes$type)
  
  node_shapes <- setNames(c(21,22),unique_types)
  
   node_colors <-setNames(c("#00b300","#ffaa00"),unique_types)   # Node color by type
 
  layout_df <- data.frame(
    x = layout_fixed[, 1],
    y = layout_fixed[, 2],
    name = V(graphall)$name  # Node names for alignment
  )
  edges.ls <- lapply(edges.ls, function(x){
    x <- x[x$weight>corr_thresh & x$pval<p_thresh,]
    return(x)
  })
   
  edges.ls <- edges.ls[unlist(lapply(edges.ls,nrow))>0]
  subgraph_ls <- lapply(edges.ls, function(x){
    subg<-induced_subgraph(graphall, vids = V(graphall)[V(graphall) %in% x$from |V(graphall) %in% x$to])
    return(subg)
    })

   
  graph_res <-lapply(names(subgraph_ls), function(title){
    
    g<-plot_graph_with_fixed_layout(subgraph_ls[[title]],layout_df,title,node_colors,node_shapes)
    return(g)
  })

  combined_plot <-  patchwork::wrap_plots(graph_res, ncol = 2, guides = "collect") +
    patchwork::plot_layout(guides = "collect", heights = unit(c(1, 1), "null"), widths = unit(1, "null")) &
    theme(plot.margin = margin(15, 15, 15, 15))
 
  
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  
  if(length(graph_res)>2){
   h=8
  }else{

   h=6
  }
  Cairo::Cairo(file = imgName, type = format, dpi = dpi, width = 9, height = h, units = "in", bg = "white")
  print(combined_plot)
  dev.off()
  
  reductionSet$analSet$diffgraph.ls <- graph_res
  .set.rdt.set(reductionSet);
  return(1)
}

plot_graph_with_fixed_layout <- function(graph, layout_df, title,node_colors,node_shapes) {
  # Merge layout_df to preserve positions for all nodes
  layout_sub <- merge(data.frame(name = V(graph)$name), layout_df, by = "name", sort = FALSE)
  edge_colors <- ifelse(E(graph)$corr > 0, "red", "blue")  
  # Set 'name' explicitly as vertex labels
  E(graph)$color <- edge_colors  # Scale edge width
  V(graph)$color <- node_colors[V(graph)$type]  # Node fill color
  V(graph)$border.color <- V(graph)$color 
  V(graph)$degree <- degree(graph) 
  E(graph)$width <- E(graph)$weight 
 
  # Plot the graph
  ggraph(graph, layout = "manual", x = layout_sub$x, y = layout_sub$y) +
    # Edges: color and width based on correlation and weight
    geom_edge_link0(aes(edge_color = color, edge_width = width), alpha = 0.8) +
    
    # Nodes: Shape and fill dynamically mapped to 'type'
    geom_node_point(aes(shape = type, fill = type,size=4 ),  color = "black", stroke = 0.5) +
    
    # Node labels
     geom_node_text(aes(label = label), size = 2.5, color = "black", repel = TRUE) +
    
    # Scales for edge and node properties
    scale_edge_color_identity(guide = "none") +  # Use edge colors directly
    scale_edge_width(range = c(0.5, 1.2)) +  # Scale edge width
    scale_shape_manual(values = node_shapes) +  # Shapes by type
    scale_fill_manual(values = node_colors) +   # Fill colors by type
    scale_size_continuous(range = c(2, 6)) +
    ggtitle(title)+
    
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 13,margin = margin(t = 10, b = 10)),
          legend.position = "none" )
}


GetChordSymbols1 <- function() {
  rdtSet <- .get.rdt.set()
  symbols <- rdtSet$chordGram[,"Var1"]
   return(symbols)
}

GetChordSymbols2 <- function() {
  rdtSet <- .get.rdt.set()
  symbols <- rdtSet$chordGram[,"Var2"]
   return(symbols)
}

GetChordColNames <- function() {
  rdtSet <- .get.rdt.set()
   if (is.null(rdtSet$chordGram)) {
    AddErrMsg("correlation result table not found."); return(0);
  }
 
 chord_colnames <- setdiff(colnames(rdtSet$chordGram),c("Var1","Var2")) # Exclude the symbol column
  
  return(chord_colnames)
}

 
GetChordFileName <- function() {
  rdtSet <- .get.rdt.set()
  
  if (is.null(rdtSet$chordGram)) {
    AddErrMsg("correlation result table not found."); return(0);
  }

  return("chord_diagram.csv")
}

GetChordMat <- function() {
  rdtSet <- .get.rdt.set()
 
   if (is.null(rdtSet$chordGram)) {
    AddErrMsg("correlation result table not found."); return(0);
  }
 
  chord_matrix <- as.matrix(subset(rdtSet$chordGram, select = -c(Var1,Var2))) # Removing the symbol column
 
  return(chord_matrix)
}

GetchordIds <- function() {
  rdtSet <- .get.rdt.set()
  
  if (is.null(rdtSet$chordGram)) {
    AddErrMsg("correlation result table not found."); return(0);
  }
  ids <- rownames(rdtSet$chordGram)
  
  return(ids)
}

GetCorrNetSymbols1 <- function() {
  rdtSet <- .get.rdt.set()
  symbols <- rdtSet$corNet[,"source"]
  return(symbols)
}

GetCorrNetSymbols2 <- function() {
  rdtSet <- .get.rdt.set() 

  symbols <- rdtSet$corNet[,"target"]
  return(symbols)
}

GetCorrNetLabel1 <- function() {
  rdtSet <- .get.rdt.set() 
  symbols <- rdtSet$corNet[,"label1"]
  return(symbols)
}

GetCorrNetLabel2 <- function() {
  rdtSet <- .get.rdt.set()
  symbols <- rdtSet$corNet[,"label2"]
  return(symbols)
}

GetCorrNetColNames <- function() {
  rdtSet <- .get.rdt.set()
  if (is.null(rdtSet$corNet)) {
    AddErrMsg("correlation result table not found."); return(0);
  }
  
  corNet_colnames <- setdiff(colnames(rdtSet$corNet),c("source","target","label1","label2")) # Exclude the symbol column
  
  return(corNet_colnames)
}


GetCorrNetFileName <- function() {
  rdtSet <- .get.rdt.set()
  
  if (is.null(rdtSet$corNet)) {
    AddErrMsg("correlation result table not found."); return(0);
  }
  
  return("corNet.csv")
}

GetCorrNetMat <- function() {
  rdtSet <- .get.rdt.set()
  
  if (is.null(rdtSet$corNet)) {
    AddErrMsg("correlation result table not found."); return(0);
  }
  
  corNet_matrix <- as.matrix(subset(rdtSet$corNet, select = -c(source,target,label1,label2))) # Removing the symbol column
  #print(head( corNet_matrix ))
  return(corNet_matrix)
}

GetCorrNetIds <- function() {
  rdtSet <- .get.rdt.set()
  
  if (is.null(rdtSet$corNet)) {
    AddErrMsg("correlation result table not found."); return(0);
  }
  ids <- rownames(rdtSet$corNet)
  
  return(ids)
}

 
plotFeatCorr <- function(reductionSet=NA,imgName,feat1="PC.O.18.2.0.16.0.0",feat2="DNMT3A",selMeta,group,dpi=150,format="png"){
  
  reductionSet <- .get.rdt.set();
  
  imgName <- paste(imgName, "dpi", dpi, ".", format, sep = "")
  library(ggplot2)
  sel.inx <- mdata.all==1; 
  sel.nms <- names(mdata.all)[sel.inx];
  dataSetList <- lapply(sel.nms, readDataset);
  id1 = paste0(dataSetList[[1]]$enrich_ids[feat1],"_",dataSetList[[1]]$type)
  id2 = paste0(dataSetList[[2]]$enrich_ids[feat2],"_",dataSetList[[2]]$type)
  selDatsCorr <- reductionSet$selDatsCorr[sel.nms]
  data<-data.frame(x= t(selDatsCorr[[1]])[,id1] ,y=t(selDatsCorr[[2]])[,id2])
  if(group!="NA"){
  meta.info = reductionSet[["dataSet"]][["meta.info"]]
  selSamps = rownames(meta.info[meta.info[[selMeta]]==group,])
   data <- data[rownames(data) %in%selSamps, ]
  }
 
  
    
    p<-ggplot(data, aes(x = x, y = y)) +
      geom_point(size = 3) + # Dots represent the samples
      geom_smooth(method = "lm", se = FALSE) +             
      labs(title = "",
           x = feat1,
           y = feat2 ) +
      theme_minimal()+
      theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), # Add panel border
            plot.margin = margin(10, 10, 10, 10),
            axis.title = element_text(size = 11, color = "black"),
            axis.text = element_text(size = 10, color = "black")) 
    
  w <- 7.5
  h<-6
  Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, width = w, height = h, type = format, bg = "white")
  print(p)
  dev.off()
  
}



GetDiffNetSymbols1 <- function(group) {
   rdtSet <- .get.rdt.set()
  diffnet <- rdtSet$diffList[[group]]
  return(diffnet$source)
}

GetDiffNetSymbols2 <- function(group) {
  rdtSet <- .get.rdt.set() 
  diffnet <- rdtSet$diffList[[group]]
  return(diffnet$target)
}



GetDiffNetLabel1 <- function(group) {
   rdtSet <- .get.rdt.set()
  diffnet <- rdtSet$diffList[[group]]
  return(diffnet$label1)
}

GetDiffNetLabel2 <- function(group) {
  rdtSet <- .get.rdt.set() 
  diffnet <- rdtSet$diffList[[group]]
  return(diffnet$label2)
}



GetDiffNetColNames <- function(group) {
  rdtSet <- .get.rdt.set() 
  if (is.null(rdtSet$diffList[[group]])) {
    AddErrMsg("correlation result table not found."); return(0);
  }
  diffNet_colnames <- c("corr" , "pval")# Exclude the symbol column
  return(diffNet_colnames)
}


GetDiffNetFileName <- function(group) {
  rdtSet <- .get.rdt.set() 
  if (is.null(rdtSet$diffList[[group]])) {
    AddErrMsg("correlation result table not found."); return(0);
  }
  write.csv(rdtSet$diffList[[group]][,c("source", "target","corr" , "pval")],"diffNet.csv",rownames=F)
  return("diffNet.csv")
}

GetDiffNetMat <- function(group) {
  rdtSet <- .get.rdt.set()
diffNet <- rdtSet$diffList[[group]]
  if (is.null(diffNet)) {
    AddErrMsg("correlation result table not found."); return(0);
  }
  diffNet_matrix <- as.matrix(diffNet[,c("corr" , "pval")])  
  return(diffNet_matrix)
}

GetDiffNetIds <- function(group) {
  rdtSet <- .get.rdt.set()
  diffNet <- rdtSet$diffList[[group]]
  if (is.null(diffNet)) {
    AddErrMsg("correlation result table not found."); return(0);
  }
  ids <- 1:nrow(diffNet)
  
  return(ids)
}
