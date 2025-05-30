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

  return(json.res);
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
  library("Spectrum")
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
  
  res <- Spectrum::Spectrum(data.list, method=method, fixk=clusterNum, maxk=10, missing=TRUE, fontsize=8, dotsize=2, showres=F, silent=T)

  clust <- res$assignments
  
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
  
  library("PINSPlus")
  result <- SubtypingOmicsData(dataList = data.list, clusteringMethod = "hclust", k = clusterNum, verbose=F, agreementCutoff=0.8)
  clust <- result$cluster1
  
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
    return(c(clustNum, signif(reductionSet$clustNmi)))
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
  
  library("SNFtool")
  truelabel = as.numeric(dat$cls)
  Data1 <- t(data.list[[1]])
  Data2 <- t(data.list[[2]])
  
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  
  W = SNF(list(W1,W2), K, T);
  res <- estimateClusters(W, 2:10);
  clusterNum <- as.numeric(clusterNum)
  if(clusterNum == -1){
    C = res[[1]]
  }else{
    C = clusterNum
  }
  group = spectralClustering(W, C); 
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
  
  library(cluster)
  reductionSet <- .get.rdt.set();
  clusterVec <- reductionSet$clustVec
  # function to compute average silhouette for k clusters
  if(type == "spectrum"){
    clustMatrix <- reductionSet$clustRes$similarity_matrix
  }else if(type == "snf"){
    clustMatrix <- reductionSet$clustDistMat
  }
  
  ss1 <- cluster::silhouette(clusterVec, daisy(t(data.list[[1]])))
  ss2 <- cluster::silhouette(clusterVec, daisy(t(data.list[[2]])))
  
  return(list(ss1,ss2))
}


PlotDiagnostic <- function(alg, imgName, dpi=72, format="png"){
  dpi <- as.numeric(dpi);
  imgNm <- paste(imgName, "dpi", dpi, ".", format, sep="");
  load_cairo();
  if(alg %in% c("snf", "spectrum","kmeans") ){
    h=8
    fig.list <- list()
    library(ggpubr);
  }else{
    h=8
  }
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  reductionSet<-.get.rdt.set();
  
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
    #res <- ComputeSilhouette("snf");
    
    #sp2 <- fviz_silhouette(res[[2]], titlenm=sel.nms[[2]])+
    #theme(axis.title.y=element_blank())
    #lims <- layer_scales(sp2)$y$get_limits()
    
    #sp1 <- fviz_silhouette(res[[1]], titlenm=sel.nms[[1]])+ ylim(lims)+
    #theme(legend.position="none")
    
    #fig.list.sub<-list();
    #fig.list.sub[[1]] <- sp1
    #fig.list.sub[[2]] <- sp2
    
    #p11 <- ggarrange(plotlist=fig.list.sub, ncol = 2, nrow = 1)
    #fig.list[[1]] <- p11
    #fig.list[[1]] <-  function(){plotEig(length(unique(reductionSet$clustVec)), reductionSet$clustRes[[5]])}
    #p1 <- ggarrange(plotlist=fig.list, ncol = 1, nrow = 1)
    plotEig(length(unique(reductionSet$clustVec)), reductionSet$clustRes[[5]])
  }else if(alg == "kmeans"){
    plotEig(length(unique(reductionSet$clustVec)), reductionSet$clustRes$frac_explained)

 }
  dev.off();
  
  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$diagnostic_components <- imgNm;
  saveSet(infoSet);
  return(1);
}


PlotHeatmapDiagnosticPca <- function(imgNm, dpi=72, format="png",type="spectrum"){
  load_cairo();
  load_ggplot();
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  sel.nms <- names(mdata.all)
  data.list <- list()
  for(i in 1:length(sel.nms)){
    dat = readDataset(sel.nms[i])
    data.list[[i]] <- dat$data.proc
  }
  reductionSet <- .get.rdt.set();
  
  fig.list <- list();
  result <- reductionSet$clustRes
  clust <- reductionSet$clustVec
  
  fig.list <- list()
  for(i in 1:length(sel.nms)){
   print(sel.nms[i])
    dataSet = readDataset(sel.nms[i])
    x <- dataSet$data.proc
    pca <- prcomp(t(na.omit(x)));
    imp.pca<-summary(pca)$importance;
    xlabel <- paste0("PC1"," (", 100*round(imp.pca[2,][1], 3), "%)")
    ylabel <- paste0("PC2"," (", 100*round(imp.pca[2,][2], 3), "%)")
    names <- colnames(x);
    pca.res <- as.data.frame(pca$x);
    pca.res <- pca.res[,c(1,2)]
    
    xlim <- GetExtendRange(pca.res$PC1);
    ylim <- GetExtendRange(pca.res$PC2);
    
    Factor <- as.factor(clust)
    pca.rest <- pca.res
    pca.rest$Cluster <- Factor
    pca.rest$names <- rownames(pca.res)
    Conditions = as.factor(dataSet$cls)
    pca.rest$Conditions <- Conditions
    
    pcafig <- ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Cluster, shape=Conditions)) +
      geom_point(size=3, alpha=0.5) + 
      xlim(xlim) + ylim(ylim) + xlab(xlabel) + ylab(ylabel) +
      theme_bw()
    fig.list[[i]] <- pcafig
    
    pcafigcls <- ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions)) +
      geom_point(size=3, alpha=0.5) + 
      xlim(xlim) + ylim(ylim) + xlab(xlabel) + ylab(ylabel) +
      theme_bw()
    if(i == 1){
      
      #fig.list[[3]] <- pcafigcls
    }else if( i == 2){
      #fig.list[[4]] <- pcafigcls
    }
  }
  h <- 6*round(length(fig.list)/2)
  Cairo(file=imgNm, width=14, height=h, type=format, bg="white", unit="in", dpi=dpi);
  library("ggpubr")
  p1 <- ggarrange(plotlist=fig.list, ncol = 2, nrow = round(length(fig.list)/2), labels=sel.nms)
  print(p1)
  dev.off();
}


PlotClusterHeatmap <- function(viewOpt="detailed", clustSelOpt="both", smplDist="pearson", clstDist="average", colorGradient="bwm",drawBorder=F, includeRowNames=T,imgName, format="png", dpi=96,width=NA){
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  rdtSet <- .get.rdt.set();
   
  metaData <- rdtSet$clustResTable;
  metaData = metaData[order( metaData$Cluster),]
  smp.nms <- rownames(metaData);
  meta.num <- ncol(metaData)

  var.nms <- rownames(metaData);
  rdtSet$imgSet$clusterHeat <- imgName;

  met <- sapply(metaData, function(x) as.integer(x))
  rownames(met) <- smp.nms;



  # set up parameter for heatmap
  if(colorGradient=="gbr"){
    colors <- colorRampPalette(c("green", "black", "red"), space="rgb")(256);
    }else if(colorGradient == "heat"){
    colors <- heat.colors(256);
    }else if(colorGradient == "topo"){
    colors <- topo.colors(256);
    }else if(colorGradient == "gray"){
        colors <- grDevices::colorRampPalette(c("grey90", "grey10"), space="rgb")(256);
    }else if(colorGradient == "byr"){
        colors <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256));
    }else if(colorGradient == "viridis") {
        colors <- rev(viridis::viridis(10))
    }else if(colorGradient == "plasma") {
        colors <- rev(viridis::plasma(10))
    }else if(colorGradient == "npj"){
        colors <- c("#00A087FF","white","#E64B35FF")
    }else if(colorGradient == "aaas"){
        colors <- c("#4DBBD5FF","white","#E64B35FF");
    }else if(colorGradient == "d3"){
        colors <- c("#2CA02CFF","white","#FF7F0EFF");
    }else {
         colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256));
    }
   
  if(clustSelOpt == "both"){
    rowBool = T;
    colBool = T;
  }else if(clustSelOpt == "row"){
    rowBool = T;
    colBool = F;
  }else if(clustSelOpt == "col"){
    rowBool = F;
    colBool = T;
  }else{
    rowBool = F;
    colBool = F;
  }

    w = min(500,ncol(met)*100+50)
    h = min(2000,nrow(met)*14+50);
 
   met <- scale_mat(met,  "column")
 

  ##plot static
  plot_dims <- get_pheatmap_dims(met, NA, "overview", width);
  h <- plot_dims$height;
  w <- plot_dims$width;
  viewOpt <- "overview";
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
       displayText = metaData;
    if(viewOpt == "overview"){
       displayText = F;
    }
 
    
    pheatmap::pheatmap(met, 
                       fontsize=12, fontsize_row=8, 
                       clustering_distance_rows = smplDist,
                       clustering_distance_cols = smplDist,
                       clustering_method = clstDist, 
                        border_color = NA,#border.col,
                       cluster_rows = F, 
                       cluster_cols = F,
                       scale = "column",
                       show_rownames= includeRowNames,
                       color = colors,
                       display_numbers=displayText);
  dev.off();

  return(.set.rdt.set(rdtSet));
}

