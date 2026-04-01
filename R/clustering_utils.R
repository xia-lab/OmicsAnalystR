##################################################
## R scripts for OmicsAnalyst
## Description: Related to clustering analysis
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

ComputeHeatmap <- function(fileNm, type){
  infoSet <- readSet(infoSet, "infoSet");
  reductionSet <- .get.rdt.set();
  if(type == "NA"){
    reductionSet$clustVec <- "NA";
  }
  .set.rdt.set(reductionSet);
  res.list <- list()
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i])
    res <- ComputePathHeatmapTable(dataSet);
    res.list[[i]] <- res;
  }
  cleanMem();
  require(rjson);
  res.list
  json.mat <- rjson::toJSON(res.list);
  sink(fileNm);
  cat(json.mat);
  sink();
  current.msg <<- "Data is now ready for heatmap visualization!";
  infoSet$paramSet$jsonNms$heatmap <- fileNm
  saveSet(infoSet);
  return(1)
}

ComputePathHeatmapTable <- function(dataSet){
  library(fastcluster);
  data <- dataSet$data.proc;
  rdtSet <- .get.rdt.set();
  
  sig.ids <- rownames(dataSet$data.proc);
  enrich.nms1 <- dataSet$enrich_ids;
  
  metadf <- rdtSet$dataSet$meta.info;
  
  res <- dataSet$comp.res[rownames(data),c(1:2)];
  colnames(res) <- c("statistic", "p.value");
  stat.pvals <- unname(as.vector(res[,2]));
  
  # scale each gene 
  dat <- t(scale(t(data)));
  
  rankPval <- order(as.vector(stat.pvals));
  stat.pvals <- stat.pvals[rankPval];
  dat <- dat[rankPval,];
  reductionSet <- .get.rdt.set();
  if (length(reductionSet$clustVec) > 1) {
    vec <- reductionSet$clustVec;
    names(vec) <- colnames(dat);
    veco <- vec[order(vec)];
    smpl.int.rk <- match(names(veco), colnames(dat));
    dat <- dat[, smpl.int.rk];
  }
  
  # now pearson and euclidean will be the same after scaleing
  dat.dist <- dist(dat); 
  
  orig.smpl.nms <- colnames(dat);
  orig.gene.nms <- rownames(dat);
  
  rest <- t(apply(dat, 1, function(x){as.numeric(cut(x, breaks=30))}));
  # do clustering and save cluster info
  # convert order to rank (score that can used to sort)
  if(nrow(dat)> 1){
    methods <- c("ward.D", "ave", "single", "complete")
    ranks <- list()
    
    t.dat.dist <- dist(t(dat))
    
    for (method in methods) {
      gene.ord <- fastcluster::hclust(dat.dist, method)$order
      ranks[[paste0("gene.", method, ".rk")]] <- match(orig.gene.nms, orig.gene.nms[gene.ord])
      
      smpl.ord <-  fastcluster::hclust(t.dat.dist, method)$order
      ranks[[paste0("smpl.", method, ".rk")]] <- match(orig.smpl.nms, orig.smpl.nms[smpl.ord])
    }
    
    gene.ward.rk <- ranks$gene.ward.D.rk
    gene.ave.rk <- ranks$gene.ave.rk
    gene.single.rk <- ranks$gene.single.rk
    gene.complete.rk <- ranks$gene.complete.rk
    
    smpl.ward.rk <- ranks$smpl.ward.D.rk
    smpl.ave.rk <- ranks$smpl.ave.rk
    smpl.single.rk <- ranks$smpl.single.rk
    smpl.complete.rk <- ranks$smpl.complete.rk
  }else{
    # force not to be single element vector which will be scaler
    #stat.pvals <- matrix(stat.pvals);
    gene.ward.rk <- gene.ave.rk <- gene.single.rk <- gene.complete.rk <- matrix(1);
    smpl.ward.rk <- smpl.ave.rk <- smpl.single.rk <- smpl.complete.rk <- 1:ncol(dat);
  }
  
  gene.cluster <- list(
    ward = gene.ward.rk,
    average = gene.ave.rk,
    single = gene.single.rk,
    complete = gene.complete.rk,
    pval = seq.int(c(1:length(gene.ward.rk)))
  );
  
  sample.cluster <- list(
    ward = smpl.ward.rk,
    average = smpl.ave.rk,
    single = smpl.single.rk,
    complete = smpl.complete.rk
  );
  
  meta.types <- reductionSet$dataSet$meta.types;
  
  # prepare meta info    
  # 1) convert meta.data info numbers
  # 2) match number to string (factor level)
  if(length(reductionSet$clustVec) > 1){
    metadf$Cluster <- reductionSet$clustVec;
    sample.cluster[[reductionSet$clustType]] <- smpl.int.rk;
    metadf <- metadf[smpl.int.rk,];
    meta.types <- c( meta.types,"disc");
    names(meta.types)[length(meta.types)] <- "Cluster";
  }
  meta <- metadf;
  meta <- meta[which(rownames(meta) %in% orig.smpl.nms), ,drop=F];
  grps <- colnames(metadf)
  nmeta <- meta.vec <- NULL;
  uniq.num <- 0;
  meta.grps <- vector();
  disc.inx <- rep(F, ncol(meta)*nrow(meta));  
  
  for (i in 1:ncol(meta)){
    cls <- meta[,i];
    grp.nm <- grps[i];
    meta.vec <- c(meta.vec, as.character(cls))
    # make sure each label are unqiue across multiple meta data
    #print(grp.nm);
    if( meta.types[grp.nm] == "disc"){
      ncls <- paste(grp.nm, as.numeric(cls)+99); # note, here to retain ordered factor
      disc.inx[c((nrow(meta)*(i-1)+1): (nrow(meta)*i))] <- T;
      sample.cluster[[grps[i]]] <- order(cls)
      
    }else{
      ncls <- as.numeric(cut(rank(as.numeric(as.character(cls))), breaks=30)); # note, here to retain ordered factor
      ord <- match(orig.smpl.nms, orig.smpl.nms[order(cls)]);
      sample.cluster[[grps[i]]] <- ord
      
    }
    meta.grps <- c(meta.grps, paste(grp.nm, rownames(meta))); 
    nmeta <- c(nmeta, ncls);
  }
  
  # convert back to numeric
  nmeta[disc.inx] <- as.numeric(as.factor(nmeta[disc.inx]))+99;
  unik.inx <- !duplicated(nmeta)   
  
  # get corresponding names
  #meta_anot <- meta.vec[unik.inx]; 
  #names(meta_anot) <- nmeta[unik.inx]; # name annotatation by their numbers
  
  meta_anot <- meta.vec; 
  names(meta_anot) <- meta.grps;
  
  nmeta <- matrix(nmeta, ncol=ncol(meta), byrow=F);
  colnames(nmeta) <- grps;
  
  # for each gene/row, first normalize and then tranform real values to 30 breaks 
  
  # note, use {} will lose order; use [[],[]] to retain the order
  
  gene.id = orig.gene.nms; if(length(gene.id) ==1) { gene.id <- matrix(gene.id) };
  hit.inx <- match(gene.id , unname(enrich.nms1));
  lbls <- names(enrich.nms1[hit.inx]);
  
  # re-order by p-value
  gene.id <- gene.id[rankPval];
  lbls <- lbls[rankPval];
  rownames(rest) <- NULL;
  list_of_lists <- apply(rest, 1, function(x) unname(as.list(x)))
  
  
  heatmap_file <- saveTopRowsHeatmapImage(dat, dataSet, lbls, topN=50, meta.info=meta, meta.types=meta.types);

  json.res <- list(
    data.name = dataSet$name,
    data.type = dataSet$type,
    gene.id = gene.id,
    gene.name = lbls,
    gene.cluster = gene.cluster,
    sample.cluster = sample.cluster,
    sample.names = orig.smpl.nms,
    meta = data.frame(nmeta),
    meta.anot = meta_anot,
    data = list_of_lists,
    org = data.org,
    pval = stat.pvals
  );
  if(!is.null(heatmap_file)){
    json.res$heatmap.image <- heatmap_file;
  }

  return(json.res);
}

saveTopRowsHeatmapImage <- function(dat, dataSet, lbls=NULL, topN = 50, meta.info=NULL, meta.types=NULL){
  if(is.null(dat) || nrow(dat) == 0 || ncol(dat) == 0 || is.null(dataSet)){
    return(NULL);
  }
  topN <- min(topN, nrow(dat));
  if(topN <= 0){
    return(NULL);
  }

  # Isolate pheatmap in subprocess
  ds_name <- tryCatch({
    if(is.null(dataSet$name)) "dataset"
    else if(is.function(dataSet$name)) "dataset"
    else as.character(dataSet$name)[1]
  }, error = function(e) "dataset")
  safeName <- gsub("[^A-Za-z0-9_]+", "_", ds_name);
  if(is.na(safeName) || safeName == ""){
    safeName <- "dataset";
  }
  heatmap_file <- paste0("heatmap_top50_", safeName, ".png");

  data_for_callr <- list(
    dat = dat,
    lbls = lbls,
    topN = topN,
    meta.info = meta.info,
    meta.types = meta.types,
    heatmap_file = heatmap_file,
    lib_paths = .libPaths()
  )

  isolated_func <- function(input_data) {
    if (!is.null(input_data$lib_paths)) {
      .libPaths(input_data$lib_paths)
    }

    dat <- input_data$dat
    lbls <- input_data$lbls
    topN <- input_data$topN
    meta.info <- input_data$meta.info
    meta.types <- input_data$meta.types
    heatmap_file <- input_data$heatmap_file

    mat <- dat[1:topN, , drop=FALSE];
    if(!is.null(lbls) && length(lbls) >= topN){
      rownames(mat) <- lbls[1:topN];
    }

    anno <- NULL;
    annotation_colors <- list();

    if(!is.null(meta.info) && nrow(meta.info) == ncol(mat)){
      primary <- names(meta.info)[1];
      if(!is.null(primary) && primary %in% names(meta.info)){
        col.order <- order(meta.info[[primary]]);
        mat <- mat[, col.order, drop=FALSE];
        if(!is.null(lbls) && length(lbls) >= ncol(mat)){
          lbls <- lbls[col.order];
        }
        meta.info <- meta.info[col.order, , drop=FALSE];
      }
      num.cols <- min(4, ncol(meta.info));
      if(num.cols > 0){
        anno <- as.data.frame(meta.info[, seq_len(num.cols), drop=FALSE]);
        rownames(anno) <- rownames(meta.info);
        for(col in colnames(anno)){
          type <- if(!is.null(meta.types)) meta.types[col] else NULL;
          is_cont <- !is.null(type) && grepl("cont", type, ignore.case=TRUE);

          if(is_cont){
            values <- suppressWarnings(as.numeric(as.character(anno[[col]])));
            if(!all(is.na(values))){
              anno[[col]] <- values;
              n_colors <- 100
              color_func <- grDevices::colorRampPalette(c("blue", "white", "red"))
              annotation_colors[[col]] <- color_func(n_colors)
            } else {
              anno[[col]] <- NULL;
            }
          } else {
            anno[[col]] <- as.character(anno[[col]]);
            values <- unique(anno[[col]]);
            values <- values[!is.na(values)];
            values <- values[order(values)];
            if(length(values) > 0){
              pal <- if(length(values) <= 8) RColorBrewer::brewer.pal(max(3, length(values)), "Set2")
                     else grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(values));
              annotation_colors[[col]] <- setNames(pal[seq_along(values)], as.character(values));
            } else {
              anno[[col]] <- NULL;
            }
          }
        }
        if(is.null(anno) || ncol(anno) == 0) anno <- NULL;
      }
    }

    width <- max(6, ncol(mat) * 0.15 + 2)
    height <- max(4, topN * 0.12 + 1.5)

    Cairo::Cairo(file = heatmap_file, unit="in", dpi=150, width=width, height=height, type="png", bg="white");
    pheatmap::pheatmap(mat,
                       color = rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256)),
                       fontsize=10, fontsize_row=8,
                       border_color = NA,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       scale = "column",
                       annotation_col = anno,
                       annotation_colors = if(length(annotation_colors)>0) annotation_colors else NULL,
                       show_rownames = TRUE,
                       show_colnames = TRUE,
                       display_numbers = FALSE);
    dev.off();

    return(heatmap_file)
  }

  result <- tryCatch({
    rsclient_isolated_exec(
      func_body = isolated_func,
      input_data = data_for_callr,
      packages = c("pheatmap", "Cairo", "RColorBrewer"),
      timeout = 120
    )
  }, error = function(e) {
    message(sprintf("saveTopRowsHeatmapImage failed: %s", e$message))
    return(NULL)
  })
  # plot/write failure is non-fatal

  return(result)
}

ComputeKmeans <- function(clusterNum="-1"){
   
  print(clusterNum)
 
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readDataset(sel.nms[i])
    data.list[[i]] <- dat$data.proc
  }
  reductionSet <- .get.rdt.set();
  reductionSet$omicstype <- sel.nms
  
  clusterNum <- as.numeric(clusterNum)
  if(clusterNum==-1){
    clusterNum =3
   }
  combined_data <-  do.call(rbind, data.list)
    
  res <- kmeans(t(combined_data), centers = clusterNum, nstart = 25)
   
  clust <- res$cluster
  
  kmNMI = .calNMI(clust, as.numeric(dat$cls));


  global_mean <- colMeans(t(combined_data))
  TSS <- sum( rowSums( (t(combined_data) - global_mean)^2 ) )
  BSS_i <- numeric(clusterNum)
  for (i in 1:clusterNum) {
    cluster_size <- res$size[i]
    cluster_centroid <- res$centers[i, ]
    dist2 <- sum((cluster_centroid - global_mean)^2)  # squared distance
    BSS_i[i] <- cluster_size * dist2
  }
    
  reductionSet$clustType <- "K-means"
  reductionSet$clustVec <- clust
  reductionSet$clustRes <- res 
  reductionSet$clustNmi <- kmNMI
  reductionSet$clustRes$frac_explained <-  BSS_i / TSS
  
  
  reductionSet$dataSet$clust.meta.types <- c("disc",reductionSet$dataSet$meta.types);
  names(reductionSet$dataSet$clust.meta.types)[1] <- "Cluster";
  #save results
  res.table <- data.frame(Cluster=clust,reductionSet$dataSet$meta.info);
   
  reductionSet$clustResTable <- res.table;
  
  write.csv(res.table, "kmeans_clustering_results.csv");
  
  .set.rdt.set(reductionSet);
  return(1)
}

# do k-means++
ComputeKmeansPP <- function(clusterNum="-1"){
   
  clusterNum <- as.numeric(clusterNum)
  if(clusterNum==-1){
    clusterNum <- 3
  }

  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat <- readDataset(sel.nms[i])
    data.list[[i]] <- dat$data.proc
  }
  combined_data <-  do.call(rbind, data.list); 
 
  # use kmeans++
  clust <- maotai::kmeanspp(t(combined_data), k = clusterNum); 

  kmNMI = .calNMI(clust, as.numeric(dat$cls));

  global_mean <- colMeans(t(combined_data))
  TSS <- sum( rowSums( (t(combined_data) - global_mean)^2 ) )
  BSS_i <- numeric(clusterNum)

   
  for (i in 1:clusterNum) {
    my.hits <- clust == i;
    cluster_size <- sum(my.hits);
    cluster_centroid <- colMeans(t(combined_data)[my.hits])
    dist2 <- sum((cluster_centroid - global_mean)^2)  # squared distance
    BSS_i[i] <- cluster_size * dist2
  }
 
  reductionSet <- .get.rdt.set();
  reductionSet$omicstype <- sel.nms;

  reductionSet$clustType <- "K-means++"
  reductionSet$clustVec <- clust
  reductionSet$clustNmi <- kmNMI
  reductionSet$clustRes$frac_explained <-  BSS_i / TSS
  
  
  reductionSet$dataSet$clust.meta.types <- c("disc",reductionSet$dataSet$meta.types);
  names(reductionSet$dataSet$clust.meta.types)[1] <- "Cluster";

   #save results
  res.table <- data.frame( Cluster=clust,reductionSet$dataSet$meta.info);
  
  reductionSet$clustResTable <- res.table;
  
  #write.csv(res.table, "kmeans_clustering_results.csv");
  write.csv(res.table, "kmeanspp_clustering_results.csv");
  .set.rdt.set(reductionSet);
  return(1)
}

 
ComputeSpectrum <- function(method="1", clusterNum="-1"){
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readDataset(sel.nms[i])
    data.list[[i]] <- dat$data.proc
  }
  reductionSet <- .get.rdt.set();
  reductionSet$omicstype <- sel.nms
  
  clusterNum <- as.numeric(clusterNum)
  if(clusterNum == -1){
    clusterNum = NULL;
    if(method == "eigengap"){
      method=1;
    }else if (method == "multimodality"){
      method=2;
    }
  }else{
    method=3;
  }

  spec_result <- tryCatch({
    rsclient_isolated_exec(
      func_body = function(input_data) {
        require(Spectrum)
        res <- Spectrum::Spectrum(input_data$data.list, method=input_data$method,
                                  fixk=input_data$clusterNum, maxk=10, missing=TRUE,
                                  fontsize=8, dotsize=2, showres=FALSE, silent=TRUE)
        list(assignments = res$assignments, res = res, similarity_matrix = res$similarity_matrix)
      },
      input_data = list(data.list = data.list, method = method, clusterNum = clusterNum),
      packages = c("Spectrum"),
      timeout = 300
    )
  }, error = function(e) {
    AddErrMsg(paste("ComputeSpectrum failed:", e$message))
    NULL
  })
  if (is.list(spec_result) && isFALSE(spec_result$success)) { AddErrMsg(spec_result$message); return(0) }
  if (is.null(spec_result)) return(0)
  res <- spec_result$res
  clust <- spec_result$assignments

  SNFNMI = .calNMI(clust, as.numeric(dat$cls));
 
  reductionSet$clustType <- "Spectrum"
  reductionSet$clustVec <- clust
  reductionSet$clustRes <- res
  reductionSet$clustDistMat <- res$similarity_matrix
  reductionSet$clustNmi <- SNFNMI

 reductionSet$dataSet$clust.meta.types <- c("disc",reductionSet$dataSet$meta.types);
  names(reductionSet$dataSet$clust.meta.types)[1] <- "Cluster";

  #save results
  res.table <- data.frame( Cluster=clust,reductionSet$dataSet$meta.info);
  reductionSet$clustResTable <- res.table;

  write.csv(res.table, "spectrum_clustering_results.csv");

  .set.rdt.set(reductionSet);
  return(1)
}

ComputePins <- function(method="kmeans", clusterNum="auto"){
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readDataset(sel.nms[i])
    data.list[[i]] <- t(dat$data.proc)
  }
  reductionSet <- .get.rdt.set();
  reductionSet$omicstype <- sel.nms
  clusterNum <- as.numeric(clusterNum)
  if(clusterNum == -1){
    clusterNum = NULL;
  }

  pins_result <- tryCatch({
    rsclient_isolated_exec(
      func_body = function(input_data) {
        require(PINSPlus)
        result <- PINSPlus::SubtypingOmicsData(dataList = input_data$data.list,
                                                 clusteringMethod = "hclust",
                                                 k = input_data$clusterNum,
                                                 verbose = FALSE, agreementCutoff = 0.8)
        list(cluster1 = result$cluster1, result = result)
      },
      input_data = list(data.list = data.list, clusterNum = clusterNum),
      packages = c("PINSPlus"),
      timeout = 300
    )
  }, error = function(e) {
    AddErrMsg(paste("ComputePins failed:", e$message))
    NULL
  })
  if (is.list(pins_result) && isFALSE(pins_result$success)) { AddErrMsg(pins_result$message); return(0) }
  if (is.null(pins_result)) return(0)
  result <- pins_result$result
  clust <- pins_result$cluster1

  SNFNMI = .calNMI(clust, as.numeric(dat$cls))
  
  reductionSet$clustType <- "Perturbation"
  reductionSet$clustVec <- clust
  reductionSet$clustRes <- result
  reductionSet$clustNmi <- SNFNMI
  
 reductionSet$dataSet$clust.meta.types <- c("disc",reductionSet$dataSet$meta.types);
  names(reductionSet$dataSet$clust.meta.types)[1] <- "Cluster";


  res.table <- data.frame(Cluster=clust,reductionSet$dataSet$meta.info);
  write.csv(res.table, "perturbation_clustering_results.csv");
  reductionSet$clustResTable <- res.table;

  .set.rdt.set(reductionSet);
  return(1)
}

# return gene and compounds highlighted in the pathway
GetClusterMembers<-function(clust){

    reductionSet <- .get.rdt.set();
    clustVec <- reductionSet$clustVec;
    sampleNames <- rownames(reductionSet$dataSet$meta.info);
    #print(sampleNames);
    #print(colnames(dataSet$proc));
    match.inx <- clustVec == clust;
    members <- sampleNames[match.inx];
    #print(members);
    return(cbind(paste0("Cluster ", clust), paste(unique(members), collapse="; ")));
}

GetDiagnosticSummary<- function(type){
  if(type %in% c("perturbation", "spectrum", "snf","kmeans")){
    reductionSet <- .get.rdt.set();
    clustNum <- length(unique(reductionSet$clustVec))
    nmi <- if (!is.null(reductionSet$clustNmi) && is.numeric(reductionSet$clustNmi)) {
      signif(reductionSet$clustNmi)
    } else {
      "NA"
    }
    return(c(clustNum, nmi))
  }else{
    return(c("","") )
  }
}

plotEig <- function(nclust=3, eig.vec){
  nkeig <- as.numeric(nclust)
  eig <- eig.vec
  proe <- eig/sum(eig)
  par(mar=c(3, 4, 1, 4))
  neig <- length(eig)
  po <- barplot(eig, plot=FALSE)
  barplot(eig,names.arg= 1:neig, col=c(rep("cyan", nkeig), rep("gray", neig-nkeig)), xlab="", ylab="Eigen Value")
  
  par(new=T) 
  plot(cbind(po, proe), frame.plot=FALSE, pch=20, axes=FALSE, xlab="Eigen Vector", ylab="")
  points(x=po, y=proe, pch=20)
  lines(x=po, y=proe)
  axis(side=4)
  mtext(side=4, text="Percentage", line=2.5, cex=0.8)
}

widget_file_size <- function(p, filenm) {
  d <- tempdir()
  filenm <- paste0(filenm, ".html")
  withr::with_dir(d, htmlwidgets::saveWidget(p,filenm))
  f <- file.path(d, filenm)
}

estimateClusters <- function (W, NUMC = 2:10) 
{
  if (min(NUMC) == 1) {
    warning("Note that we always assume there are more than one cluster.")
    NUMC <- NUMC[NUMC > 1]
  }
  W <- (W + t(W))/2
  diag(W) <- 0
  if (length(NUMC) <= 0) {
    warning(paste("Invalid NUMC provided, must be an integer vector", 
                  "with atleast one other number than 1.", "Using default NUMC=c(2,3,4,5)", 
                  sep = ""))
    NUMC <- 2:5
  }
  degs <- rowSums(W)
  degs[degs == 0] <- .Machine$double.eps
  D <- diag(degs)
  L <- D - W
  Di <- diag(1/sqrt(degs))
  L <- Di %*% L %*% Di
  eigs <- eigen(L)
  eigs_order <- sort(eigs$values, index.return = T)$ix
  eigs$values <- eigs$values[eigs_order]
  eigs$vectors <- eigs$vectors[, eigs_order]
  eigengap <- abs(diff(eigs$values))
  quality <- list()
  for (c_index in 1:length(NUMC)) {
    ck <- NUMC[c_index]
    UU <- eigs$vectors[, 1:ck]
    EigenvectorsDiscrete <- .discretisation(UU)[[1]]
    EigenVectors <- EigenvectorsDiscrete^2
    temp1 <- EigenVectors[do.call(order, lapply(1:ncol(EigenVectors), 
                                                function(i) EigenVectors[, i])), ]
    temp1 <- t(apply(temp1, 1, sort, TRUE))
    quality[[c_index]] <- (1 - eigs$values[ck + 1])/(1 - 
                                                       eigs$values[ck]) * sum(sum(diag(1/(temp1[, 1] + .Machine$double.eps)) %*% 
                                                                                    temp1[, 1:max(2, ck - 1)]))
  }
  t1 <- sort(eigengap[NUMC], decreasing = TRUE, index.return = T)$ix
  K1 <- NUMC[t1[1]]
  K12 <- NUMC[t1[2]]
  t2 <- sort(unlist(quality), index.return = TRUE)$ix
  K2 <- NUMC[t2[1]]
  K22 <- NUMC[t2[2]]
  output <- list(`Eigen-gap best` = K1, `Eigen-gap 2nd best` = K12, 
                 `Rotation cost best` = K2, `Rotation cost 2nd best` = K22, "eigs" =  (1-eigs$values)[NUMC])
  return(output)
}


.discretisation <- function(eigenVectors) {
  
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  
  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
    
  }
  
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}


ComputeSNF <- function(method="1", clusterNum="auto"){
  K = 20;		# number of neighbors, usually (10~30)
  alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
  T = 10; 	# Number of Iterations, usually (10~20)
  
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readDataset(sel.nms[i])
    data.list[[i]] <- dat$data.proc
  }
  reductionSet <- .get.rdt.set();
  reductionSet$omicstype <- sel.nms
  
  truelabel = as.numeric(dat$cls)
  Data1 <- t(data.list[[1]])
  Data2 <- t(data.list[[2]])
  clusterNum <- as.numeric(clusterNum)

  snf_result <- tryCatch({
    rsclient_isolated_exec(
      func_body = function(input_data) {
        require(SNFtool)
        K <- input_data$K; alpha <- input_data$alpha; T_iter <- input_data$T_iter
        Data1 <- input_data$Data1; Data2 <- input_data$Data2
        clusterNum <- input_data$clusterNum

        Dist1 <- SNFtool::dist2(as.matrix(Data1), as.matrix(Data1))
        Dist2 <- SNFtool::dist2(as.matrix(Data2), as.matrix(Data2))
        W1 <- SNFtool::affinityMatrix(Dist1, K, alpha)
        W2 <- SNFtool::affinityMatrix(Dist2, K, alpha)
        W <- SNFtool::SNF(list(W1, W2), K, T_iter)

        # estimateClusters is defined in clustering_utils.R, not SNFtool
        # Inline the eigen-gap estimation here
        W2m <- (W + t(W)) / 2; diag(W2m) <- 0
        degs <- rowSums(W2m); degs[degs == 0] <- .Machine$double.eps
        D <- diag(degs); L <- D - W2m; Di <- diag(1/sqrt(degs))
        L <- Di %*% L %*% Di
        eigs <- eigen(L)
        eigs_order <- sort(eigs$values, index.return = TRUE)$ix
        eigs$values <- eigs$values[eigs_order]
        eigengap <- abs(diff(eigs$values))
        NUMC <- 2:10
        t1 <- sort(eigengap[NUMC], decreasing = TRUE, index.return = TRUE)$ix
        auto_k <- NUMC[t1[1]]
        eigs_vals <- (1 - eigs$values)[NUMC]

        C <- if (clusterNum == -1) auto_k else clusterNum
        group <- SNFtool::spectralClustering(W, C)

        list(group = group, W = W, auto_k = auto_k, eigs_vals = eigs_vals)
      },
      input_data = list(Data1 = Data1, Data2 = Data2, K = K, alpha = alpha,
                        T_iter = T, clusterNum = clusterNum),
      packages = c("SNFtool"),
      timeout = 300
    )
  }, error = function(e) {
    AddErrMsg(paste("ComputeSNF failed:", e$message))
    NULL
  })
  if (is.list(snf_result) && isFALSE(snf_result$success)) { AddErrMsg(snf_result$message); return(0) }
  if (is.null(snf_result)) return(0)
  group <- snf_result$group
  W <- snf_result$W
  res <- estimateClusters(W, 2:10)
  SNFNMI = .calNMI(group, truelabel)
  reductionSet$clustType <- "SNF"
  reductionSet$clustVec <- group
  reductionSet$clustRes <- res
  reductionSet$clustDistMat <- W
  reductionSet$clustNmi <- SNFNMI 
   reductionSet$dataSet$clust.meta.types <- c("disc",reductionSet$dataSet$meta.types);
  names(reductionSet$dataSet$clust.meta.types)[1] <- "Cluster";
  res.table <- data.frame(Cluster=group,reductionSet$dataSet$meta.info);
  write.csv(res.table, "snf_clustering_results.csv");
  reductionSet$clustResTable <- res.table;

  .set.rdt.set(reductionSet);
  return(1)
}

.calNMI<-function (x, y) 
{
  x <- as.vector(x)
  y <- as.vector(y)
  mutual.info <- (.mutualInformation(x, y)/sqrt(.entropy(x) * 
                                                  .entropy(y)))
  return(max(0, mutual.info, na.rm = TRUE))
}

ComputeSilhouette <-function(type){
  
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readDataset(sel.nms[i])
    data.list[[i]] <- dat$data.proc
  }
  
  reductionSet <- .get.rdt.set();
  clusterVec <- reductionSet$clustVec
  # function to compute average silhouette for k clusters
  if(type == "spectrum"){
    clustMatrix <- reductionSet$clustRes$similarity_matrix
  }else if(type == "snf"){
    clustMatrix <- reductionSet$clustDistMat
  }

  sil_result <- tryCatch({
    rsclient_isolated_exec(
      func_body = function(input_data) {
        require(cluster)
        ss1 <- cluster::silhouette(input_data$clusterVec, cluster::daisy(t(input_data$data1)))
        ss2 <- cluster::silhouette(input_data$clusterVec, cluster::daisy(t(input_data$data2)))
        list(ss1 = ss1, ss2 = ss2)
      },
      input_data = list(clusterVec = clusterVec,
                        data1 = data.list[[1]], data2 = data.list[[2]]),
      packages = c("cluster"),
      timeout = 120
    )
  }, error = function(e) {
    AddErrMsg(paste("ComputeSilhouette failed:", e$message))
    NULL
  })
  if (is.list(sil_result) && isFALSE(sil_result$success)) { AddErrMsg(sil_result$message); return(0) }
  if (is.null(sil_result)) return(0)
  return(list(sil_result$ss1, sil_result$ss2))
}


PlotDiagnostic <- function(alg, imgName, dpi=150, format="png"){
  dpi <- as.numeric(dpi);
  imgNm <- paste(imgName, "dpi", dpi, ".", format, sep="");
  h <- 8;
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  reductionSet<-.get.rdt.set();

  if (alg == "perturbation") {
    # Perturbation branch uses ggpubr/tidyverse - isolate
    res <- reductionSet$clustRes
    auc1 <- res$dataTypeResult[[1]]$Discrepancy$AUC[-1]
    auc2 <- res$dataTypeResult[[2]]$Discrepancy$AUC[-1]
    n_clust <- length(unique(res$cluster1))

    tryCatch({
      rsclient_isolated_exec(
        func_body = function(input_data) {
          require(ggplot2)
          require(Cairo)
          require(tidyr)

          auc.df1 <- data.frame(K = seq.int(length(input_data$auc1)) + 1, evals = input_data$auc1)
          auc.df2 <- data.frame(K = seq.int(length(input_data$auc2)) + 1, evals = input_data$auc2)
          auc.df1$evals2 <- auc.df2$evals
          colnames(auc.df1) <- c("K", input_data$sel_nm1, input_data$sel_nm2)

          df <- tidyr::gather(auc.df1, key = "variable", value = "value", -K)

          p1 <- ggplot2::ggplot(df, ggplot2::aes(x = K, y = value)) +
            ggplot2::geom_point(ggplot2::aes(color = variable), size = 2) +
            ggplot2::geom_line(ggplot2::aes(color = variable)) +
            ggplot2::geom_vline(xintercept = input_data$n_clust, linetype = "dashed") +
            ggplot2::xlab("Number of clusters") + ggplot2::ylab("AUC") +
            ggplot2::theme_bw()

          Cairo::Cairo(file = input_data$imgNm, width = 10, height = 8,
                       type = "png", unit = "in", bg = "white", dpi = input_data$dpi)
          print(p1)
          dev.off()
          return(1)
        },
        input_data = list(auc1 = auc1, auc2 = auc2, sel_nm1 = sel.nms[[1]],
                          sel_nm2 = sel.nms[[2]], n_clust = n_clust,
                          imgNm = imgNm, dpi = dpi),
        packages = c("ggplot2", "Cairo", "tidyr", "qs"),
        timeout = 300
      )
    }, error = function(e) {
      AddErrMsg(paste("PlotDiagnostic failed:", e$message))
      return(0)
    })
    # plot/write failure is non-fatal
  } else {
    # Standard path or non-perturbation algorithms (base R plotting)
    load_cairo();
    if(alg %in% c("snf", "spectrum","kmeans") ){
      fig.list <- list()
      library(ggpubr);
    }

    Cairo(file=imgNm, width=10, height=h, type="png",unit="in", bg="white", dpi=dpi);
    if(alg == "spectrum"){
      if(!is.null(reductionSet$clustRes$eigenvector_analysis)){
        plotEig(length(unique(reductionSet$clustVec)), reductionSet$clustRes$eigenvector_analysis[,2]);
      }else{
        plotEig(length(unique(reductionSet$clustVec)), reductionSet$clustRes$eigensystem$values[1:10]);
      }
    }else if(alg == "perturbation"){
      res <- reductionSet$clustRes
      library(ggpubr)
      load_ggplot();
      xlabel="Number of clusters"
      ylabel="AUC"
      auc1 <- res$dataTypeResult[[1]]$Discrepancy$AUC[-1]
      auc.df1 <- data.frame(K=seq.int(length(auc1))+1, evals=auc1)
      auc2 <- res$dataTypeResult[[2]]$Discrepancy$AUC[-1]
      auc.df2 <- data.frame(K=seq.int(length(auc2))+1, evals=auc2)
      auc.df1$evals2 = auc.df2$evals
      colnames(auc.df1) = c("K", sel.nms[[1]], sel.nms[[2]])
      library("tidyverse")
      df <- auc.df1 %>%
        select(K,  sel.nms[[1]], sel.nms[[2]]) %>%
        gather(key = "variable", value = "value", -K)
      p1 <- ggplot(df, aes(x = K, y = value)) +
        geom_point(aes(color = variable), size=2) +
        geom_line(aes(color = variable)) +
        geom_vline(xintercept = length(unique(res$cluster1)),linetype="dashed")+ xlab(xlabel) +
        ylab(ylabel) +
        theme_bw()
      print(p1)
    }else if(alg == "snf"){
      plotEig(length(unique(reductionSet$clustVec)), reductionSet$clustRes[[5]])
    }else if(alg == "kmeans"){
      plotEig(length(unique(reductionSet$clustVec)), reductionSet$clustRes$frac_explained)
    }
    dev.off();
  }

  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$diagnostic_components <- imgNm;
  saveSet(infoSet);
  return(1);
}


PlotHeatmapDiagnosticPca <- function(imgNm, dpi=150, format="png",type="spectrum"){
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");

  sel.nms <- names(mdata.all)
  data.list <- list()
  cls.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readDataset(sel.nms[i])
    data.list[[sel.nms[i]]] <- dat$data.proc
    cls.list[[sel.nms[i]]] <- dat$cls
  }
  reductionSet <- .get.rdt.set();
  clust <- reductionSet$clustVec

  tryCatch({
    rsclient_isolated_exec(
      func_body = function(input_data) {
        require(ggpubr)
        require(ggplot2)
        require(Cairo)

        sel.nms <- input_data$sel.nms
        data.list <- input_data$data.list
        cls.list <- input_data$cls.list
        clust <- input_data$clust

        fig.list <- list()
        for (i in 1:length(sel.nms)) {
          x <- data.list[[sel.nms[i]]]
          pca <- prcomp(t(na.omit(x)))
          imp.pca <- summary(pca)$importance
          xlabel <- paste0("PC1", " (", 100 * round(imp.pca[2, ][1], 3), "%)")
          ylabel <- paste0("PC2", " (", 100 * round(imp.pca[2, ][2], 3), "%)")
          pca.res <- as.data.frame(pca$x)
          pca.res <- pca.res[, c(1, 2)]

          xlim <- range(pca.res$PC1)
          xlim <- xlim + c(-1, 1) * diff(xlim) * 0.1
          ylim <- range(pca.res$PC2)
          ylim <- ylim + c(-1, 1) * diff(ylim) * 0.1

          Factor <- as.factor(clust)
          pca.rest <- pca.res
          pca.rest$Cluster <- Factor
          pca.rest$names <- rownames(pca.res)
          Conditions <- as.factor(cls.list[[sel.nms[i]]])
          pca.rest$Conditions <- Conditions

          pcafig <- ggplot2::ggplot(pca.rest, ggplot2::aes(x = PC1, y = PC2, color = Cluster, shape = Conditions)) +
            ggplot2::geom_point(size = 3, alpha = 0.5) +
            ggplot2::xlim(xlim) + ggplot2::ylim(ylim) + ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) +
            ggplot2::theme_bw()
          fig.list[[i]] <- pcafig
        }

        h <- 6 * round(length(fig.list) / 2)
        Cairo::Cairo(file = input_data$imgNm, width = 14, height = h,
                     type = input_data$format, bg = "white", unit = "in", dpi = input_data$dpi)
        p1 <- ggpubr::ggarrange(plotlist = fig.list, ncol = 2,
                                 nrow = round(length(fig.list) / 2), labels = sel.nms)
        print(p1)
        dev.off()
        return(1)
      },
      input_data = list(sel.nms = sel.nms, data.list = data.list, cls.list = cls.list,
                        clust = clust, imgNm = imgNm, dpi = dpi, format = format),
      packages = c("ggpubr", "ggplot2", "Cairo", "qs"),
      timeout = 300
    )
  }, error = function(e) {
    AddErrMsg(paste("PlotHeatmapDiagnosticPca failed:", e$message))
    return(0)
  })
  # plot/write failure is non-fatal
}


PlotClusterHeatmap <- function(viewOpt="detailed", clustSelOpt="both", smplDist="pearson", clstDist="average", colorGradient="bwm",drawBorder=F, includeRowNames=T,imgName, format="png", dpi=150,width=NA){

  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  rdtSet <- .get.rdt.set();

  metaData <- rdtSet$clustResTable;
  if(is.null(metaData) || nrow(metaData) == 0){
    return(.set.rdt.set(rdtSet));
  }

  rdtSet$imgSet$clusterHeat <- imgName;

  # Also save to infoSet for report generation
  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$clusterHeat <- imgName;
  saveSet(infoSet);

  # Isolate pheatmap in subprocess
  metaData = metaData[order( metaData$Cluster),]
  smp.nms <- rownames(metaData);
  met <- sapply(metaData, function(x) as.integer(x))
  rownames(met) <- smp.nms;

  data_for_callr <- list(
    metaData = metaData,
    met = met,
    viewOpt = viewOpt,
    clustSelOpt = clustSelOpt,
    smplDist = smplDist,
    clstDist = clstDist,
    colorGradient = colorGradient,
    includeRowNames = includeRowNames,
    imgName = imgName,
    format = format,
    dpi = dpi,
    width = width,
    lib_paths = .libPaths()
  )

  isolated_func <- function(input_data) {
    if (!is.null(input_data$lib_paths)) {
      .libPaths(input_data$lib_paths)
    }

    met <- input_data$met
    metaData <- input_data$metaData

    colorGradient <- input_data$colorGradient
    if(colorGradient == "gbr"){
      colors <- colorRampPalette(c("green", "black", "red"), space="rgb")(256)
    } else if(colorGradient == "heat") {
      colors <- rev(heat.colors(10))
    } else if(colorGradient == "topo") {
      colors <- rev(topo.colors(10))
    } else if(colorGradient == "bwm") {
      colors <- colorRampPalette(c("#0000FF","white","#FF0000"))(256)
    } else {
      colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256));
    }

    clustSelOpt <- input_data$clustSelOpt
    if(clustSelOpt == "both"){
      rowBool <- TRUE; colBool <- TRUE;
    } else if(clustSelOpt == "row"){
      rowBool <- TRUE; colBool <- FALSE;
    } else if(clustSelOpt == "col"){
      rowBool <- FALSE; colBool <- TRUE;
    } else {
      rowBool <- FALSE; colBool <- FALSE;
    }

    met <- scale(met)

    w <- max(6, ncol(met) * 0.3 + 2)
    h <- max(4, nrow(met) * 0.15 + 1.5)

    Cairo::Cairo(file = input_data$imgName, unit="in", dpi=input_data$dpi,
                 width=w, height=h, type=input_data$format, bg="white");

    pheatmap::pheatmap(met,
                       fontsize=12, fontsize_row=8,
                       clustering_distance_rows = input_data$smplDist,
                       clustering_distance_cols = input_data$smplDist,
                       clustering_method = input_data$clstDist,
                       border_color = NA,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       scale = "column",
                       show_rownames = input_data$includeRowNames,
                       color = colors,
                       display_numbers = FALSE);
    dev.off();

    return(input_data$imgName)
  }

  tryCatch({
    rsclient_isolated_exec(
      func_body = isolated_func,
      input_data = data_for_callr,
      packages = c("pheatmap", "Cairo", "RColorBrewer"),
      timeout = 120
    )
  }, error = function(e) {
    message(sprintf("PlotClusterHeatmap failed: %s", e$message))
  })
  # plot/write failure is non-fatal

  return(.set.rdt.set(rdtSet));
}
