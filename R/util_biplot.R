

GenerateBiplot <- function(method="RDA",topN=10,predictor,dataName, fileName = "biplot", format = "png",dpi = 300,colorGradient="d3") {
   library(vegan)
  library(ggplot2)
  library(ggrepel)
print(method)
  dataSet <- readDataset(dataName);
  rdtSet <- .get.rdt.set()
  meta <- rdtSet$dataSet$meta.info 
  feature_table <- t(dataSet$data.proc);
     if (!(exists("includeMeta") && length(includeMeta) == 0  )) {
      meta = meta[,c(predictor,includeMeta)]
    }else{
      meta = meta[,predictor,drop=F]
    }

idx <- which(rdtSet$dataSet$meta.types[colnames(meta)]=="cont")
  if(length(idx)>0){
    meta[,idx] <- apply(meta[,idx,drop=F],2,as.numeric)
  }
  
color_scale <- if (rdtSet$dataSet$meta.types[predictor] == "cont") {
   if(colorGradient == "gray"){
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
        colors <- scale_color_gradientn(colors = c("#2CA02CFF", "white", "#FF7F0EFF"));
    }else {
         colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256));
    }
} else {
   if(colorGradient == "gray"){
            colors <- c("grey90", "grey70", "grey50", "grey30", "grey10") 
    }else if(colorGradient == "byr"){
       colors <- RColorBrewer::brewer.pal(10, "RdYlBu") 
    }else if(colorGradient == "viridis") {
           colors <- viridis::viridis(10) 
    }else if(colorGradient == "plasma") {
        colors <-  viridis::plasma(10) 
    }else if(colorGradient == "npj"){
        colors <- ggsci::scale_color_npg("nrc")
    }else if(colorGradient == "aaas"){
        colors <- ggsci::scale_color_aaas("default")
    }else if(colorGradient == "d3"){
        colors <- ggsci::scale_color_d3("category10")
    }else {
         colors <- RColorBrewer::brewer.pal(10, "RdBu") 
    }
}
 
  if(method=="RDA"){
    cls.type <-  rdtSet[["dataSet"]][["meta.types"]]
    idx <- which(cls.type[colnames(meta)]=="cont")
    if(length(idx)>0){
      meta[,idx] <- apply(meta[,idx,drop=F],2,as.numeric)
    }
 
    rda_result=vegan::rda(feature_table~.,meta)#
    ii=summary(rda_result)  
    sp=as.data.frame(ii$species[,1:2])
    st=as.data.frame(ii$sites)
    yz=as.data.frame(ii$biplot[,1:2])
    importance <- apply(abs(sp), 1, sum)  # Sum of absolute scores for Axes 1 & 2
    top_features <- names(sort(importance, decreasing = TRUE)[1:topN])
    overall_test <- anova.cca(rda_result, step = 10000, by = "term")
    if(is.na(overall_test$F[1])){
      stats <- "The model is overfitted with no unconstrained (residual) component."
    }else{
      stats <- paste0('[Constrained Permutational ANOVA] F-value: ',round(overall_test$F[1],3),";p-value: ",overall_test$`Pr(>F)`[1])
    }
    
    scaling_factor <- max(abs(st$RDA1), abs(st$RDA2)) / max(abs(sp$RDA1), abs(sp$RDA2))
    sp=sp[rownames(sp) %in% top_features,]*scaling_factor
    scaling_factor <- max(abs(st$RDA1), abs(st$RDA2)) / max(abs(yz$RDA1),abs(yz$RDA2))
      yz=yz*scaling_factor
 
      p <- ggplot() +
      geom_point(data = st, aes(x = RDA1, y = RDA2, color = meta[, predictor]), size = 3, alpha = 0.8) +
      color_scale + 
      labs(color = predictor) + 
      theme_minimal() +
      geom_segment(data = sp,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                     arrow = arrow(angle=10,length = unit(0.3,"cm"),
            type = "closed"),linetype=1, size=0.3,colour = "#000078ff",alpha=1) +
      geom_text_repel(data = sp,aes(RDA1,RDA2,label=row.names(sp)),size=3,
                      color="#000078ff",segment.color = 'grey',fontface ='italic',alpha=1) +
     geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                     arrow = arrow(angle=10,length = unit(0.3,"cm"),
                                   type = "closed"),linetype=1, size=0.3,colour = "#006000ff",alpha=1) + 
     geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)),size=3,color='#006000ff',segment.color = 'grey',fontface ='plain')+ #bold
      labs(x=paste("RDA1 (", format(100 *ii$cont[[1]][2,1], digits=4), "%)", sep=""),
           y=paste("RDA2 (", format(100 *ii$cont[[1]][2,2], digits=4), "%)", sep=""))+
      geom_hline(yintercept=0,linetype=3,size=0.5,color='gray') + 
      geom_vline(xintercept=0,linetype=3,size=0.5,color='gray')+
      guides(shape=guide_legend(title=NULL,color="black"), fill = F)+
      theme(
axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank()   
      )
    
    
  }else{
    library(factoextra) 
    library(dplyr)
    pca <- prcomp(feature_table, scale=F);
    cls <- rdtSet$dataSet$meta.info[,predictor];
    cls.type <-  rdtSet[["dataSet"]][["meta.types"]][[predictor]];
    if(class(cls)=="integer"){
      cls <- as.factor(as.numeric(levels(cls))[cls]);
     meta[,predictor] <- as.numeric(meta[,predictor])
    }else{
      cls <- cls;
    }
    
    contrib <- get_pca_var(pca)$contrib
    contrib_pc1_pc2 <- contrib[, 1] + contrib[, 2]
    top_features <- names(sort(contrib_pc1_pc2, decreasing = TRUE))[1:topN]
    top_features <- top_features[!is.na(top_features)]
    
    scores <-  pca$x;
    choices = c(1, 2);
    lam <- pca$sdev[choices]
    n <- NROW(scores)
    lam <- lam * sqrt(n);
    ind_data <- data.frame(
      PC1 = scores[, choices[1]] / lam[1],
      PC2 = scores[, choices[2]] / lam[2] 
    )
    st=scores
    loadings <- pca$rotation
    var_data <- data.frame(
      PC1 = loadings[, choices[1]] * lam[1],
      PC2 = loadings[, choices[2]] * lam[2] 
    )%>% 
      filter(rownames(.) %in% top_features)
    
    scaling_factor <- max(abs(ind_data$PC1), abs(ind_data$PC2)) / max(abs(var_data$PC1), abs(var_data$PC2))
    
    var_data$PC1 <- var_data$PC1 * scaling_factor
    var_data$PC2 <- var_data$PC2 * scaling_factor

    data.dist <- dist(as.matrix(feature_table), method = 'euclidean');
    overall_test <- adonis2(data.dist ~ ., data = meta, permutations = 999)

    if(is.na(overall_test$F[1])){
      stats <- "The model is overfitted with no unconstrained (residual) component."
    }else{
      stats <- paste0('[PERMANOVA] F-value: ',round(overall_test$F[1],5),"; R-squared: ", signif(overall_test$R2, 5),";p-value: ",overall_test$`Pr(>F)`[1])
    }

    p <- ggplot() +
      geom_point(data = ind_data, aes(x = PC1, y = PC2, color = meta[, predictor]), size = 3, alpha = 0.8) +
      color_scale + 
      labs(color = predictor) + 
      theme_minimal() +
      geom_segment(data = var_data,aes(x = 0, y = 0, xend = PC1, yend = PC2), 
                   arrow = arrow(angle=10,length = unit(0.2,"cm"),
                                 type = "closed"),linetype=1, size=0.4,colour = "#000078ff",alpha=1) +
      geom_text_repel(data = var_data,aes(PC1,PC2,label=row.names(var_data)),size=3,
                      color="#000078ff",segment.color = 'grey',fontface ='italic',alpha=1) +
  labs( x = paste("PC1 (", round(summary(pca)$importance[2, 1] * 100, 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca)$importance[2, 2] * 100, 2), "%)", sep = ""))+
      geom_hline(yintercept=0,linetype=3,size=0.5,color='gray') + 
      geom_vline(xintercept=0,linetype=3,size=0.5,color='gray')+
      guides(shape=guide_legend(title=NULL,color="black"), fill = F)+
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    axis.line = element_blank()   
      )
           
  }
  
  # Plot the variance partitioning results and save the image
    imgName = paste(fileName, "dpi", dpi, ".", format, sep="");

  Cairo::Cairo(file = imgName, type = format, dpi = dpi, width = 7, height = 5, units = "in", bg = "white")
  print(p)
  dev.off()
  rdtSet$analSet$biplot_method <- method;
  rdtSet$analSet$biplot_fileName <- "biplot_ordination.csv";
  rdtSet$analSet$biplot_ordination <- st;
  
  fast.write.csv(st, file="biplot_ordination.csv")
  .set.rdt.set(rdtSet)
  return(stats)
}
