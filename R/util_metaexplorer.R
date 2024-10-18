
#'Generate heatmaps for metadata table
#'@description Plot a heatmap showing clustering patterns among the metadata
#'@param rdtSet Input the name of the created rdtSet (see InitDataObjects)
#'@param viewOpt high-level summary or plotting the names inside cell
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotMetaHeatmap <- function(viewOpt="detailed", clustSelOpt="both", smplDist="pearson", clstDist="average", colorGradient="bwm",includeRowNames=T,drawBorder=T, imgName, format="png", dpi=96,width=NA){
  plotjs <- paste0(imgName, ".json");
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  rdtSet <- .get.rdt.set();
  require(iheatmapr);

  metaData <- rdtSet$dataSet$meta.info;
  smp.nms <- rownames(metaData);
  meta.num <- ncol(metaData)

  var.nms <- rownames(metaData);
  rdtSet$imgSet$metahtmaptwo <- imgName;

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

    w = min(1200,ncol(met)*200+50)
    h = min(2000,nrow(met)*14+50);
 
   met <- scale_mat(met,  "column")
    p <- iheatmap(met,  name = " ", 
                  colors = colors, 
          colorbar_grid = setup_colorbar_grid(y_start = 0.85)) %>%
      add_col_labels(size = 0.1, font = list(size = 14)) 
 
    if (includeRowNames) {
      p <- p%>%
        add_row_labels(size = 0.2, font = list(size = 10), side = "right") 
    }
 
     if (rowBool) {
     if(smplDist=="correlation"){
      my.dist <- cor(t(met), method="pearson")
      my.dist <- 1-my.dist 
      my.dist <- as.dist(my.dist, diag = FALSE, upper = F)
    }else{
      my.dist = dist(met, method = smplDist)
    }
 
    dend_row <- hclust(my.dist, method = clstDist)
    p <- p %>% add_row_dendro(dend_row,  size = 0.15)  
  } 
     if (colBool) {
    if(smplDist=="correlation"){
      my.dist2 <- cor(met, method="pearson")
      my.dist2 <- 1-my.dist2 
      my.dist2 <- as.dist(my.dist2, diag = FALSE, upper = F)
    }else{
      my.dist2 = dist(t(met), method = smplDist)
    }
    dend_col <- hclust( my.dist2, method = clstDist)
    p <- p %>% add_col_dendro(dend_col,  size = ifelse(h>1300,0.02,0.05))
    } 
    saveRDS(p,"metadata_heatmap.rds");

     # Adjust the height and width (in pixels)
 
    as_list <- to_plotly_list(p)
    as_list[["layout"]][["width"]] <- w
    as_list[["layout"]][["height"]] <- max(h,500)
  
    as_json <- attr(as_list, "TOJSON_FUNC")(as_list)
    as_json <- paste0("{ \"x\":", as_json, ",\"evals\": [],\"jsHooks\": []}")
 
    print(plotjs)
    write(as_json, plotjs)

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

   p <- pheatmap::pheatmap(met, 
                       fontsize=12, fontsize_row=8, 
                       clustering_distance_rows = smplDist,
                       clustering_distance_cols = smplDist,
                       clustering_method = clstDist, 
                      # border_color = border.col,
                       cluster_rows = rowBool, 
                       cluster_cols = colBool,
                       scale = "column",
                       show_rownames= includeRowNames,
                       color = colors,
                       display_numbers=displayText,silent = TRUE);
       if(rowBool){
         p$tree_row$order <- rev(p$tree_row$order)
         rowBool <-  p$tree_row
      }else{
         met <- met[,ncol(met):1]
      }
    
    pheatmap::pheatmap(met, 
                       fontsize=12, fontsize_row=8, 
                       clustering_distance_rows = smplDist,
                       clustering_distance_cols = smplDist,
                       clustering_method = clstDist, 
                      # border_color = border.col,
                       cluster_rows = rowBool, 
                       cluster_cols = colBool,
                       scale = "column",
                       show_rownames= includeRowNames,
                       color = colors,
                       display_numbers=displayText);
  dev.off();

  return(.set.rdt.set(rdtSet));
}

#'Create high resolution static HeatMap for download only
#'@description Plot a heatmap showing clustering patterns among the metadata
#'@param rdtSet Input the name of the created rdtSet (see InitDataObjects)
#'@param viewOpt high-level summary or plotting the names inside cell
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotStaticMetaHeatmap <- function(rdtSet=NA, viewOpt="detailed", clustSelOpt="both", smplDist="pearson", clstDist="average", colorGradient="bwm",includeRowNames=T, imgName, format="png", dpi=96, width=NA){
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  rdtSet <- .get.mSet(rdtSet);

  metaData <- rdtSet$dataSet$meta.info;
  smp.nms <- rownames(metaData);
  meta.num <- ncol(metaData)

  var.nms <- rownames(metaData);
  rdtSet$imgSet$metahtmaptwo <- imgName;

  met <- sapply(metaData, function(x) as.integer(x))
  rownames(met) <- smp.nms;

  plot_dims <- get_pheatmap_dims(met, NA, viewOpt, width);
  h <- plot_dims$height;
  w <- plot_dims$width;

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

  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
       displayText = metaData;
    if(viewOpt == "overview"){
       displayText = F;
    }

   p <- pheatmap::pheatmap(met, 
                       fontsize=12, fontsize_row=8, 
                       clustering_distance_rows = smplDist,
                       clustering_distance_cols = smplDist,
                       clustering_method = clstDist, 
                      # border_color = border.col,
                       cluster_rows = rowBool, 
                       cluster_cols = colBool,
                       scale = "column",
                       show_rownames= includeRowNames,
                       color = colors,
                       display_numbers=displayText,silent = TRUE);
       if(rowBool){
         p$tree_row$order <- rev(p$tree_row$order)
         rowBool <-  p$tree_row
      }else{
         met <- met[,ncol(met):1]
      }
    
    pheatmap::pheatmap(met, 
                       fontsize=12, fontsize_row=8, 
                       clustering_distance_rows = smplDist,
                       clustering_distance_cols = smplDist,
                       clustering_method = clstDist, 
                      # border_color = border.col,
                       cluster_rows = rowBool, 
                       cluster_cols = colBool,
                       scale = "column",
                       show_rownames= includeRowNames,
                       color = colors,
                       display_numbers=displayText);
  dev.off();

  infoSet <- readSet(infoSet, "infoSet");
  infoSet$imgSet$metaHeatmap <- imgName;
  saveSet(infoSet);

  return(.set.rdt.set(rdtSet));
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
} 

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}
