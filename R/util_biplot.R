

PerformOrdination <- function(method,predictor,dataName) {

 
if(!exists("mem.ordination")){
    require("memoise");
    mem.ordination <<- memoise(.perform.ordination);
  }
  return(mem.ordination(method,predictor,dataName,includeMeta));
 
}


.perform.ordination <- function(method = "RDA",
                              predictor,
                              dataName,includeMeta) {
  print(includeMeta)
  # Load or prepare your data
  dataSet <- readDataset(dataName)
  rdtSet  <- .get.rdt.set()
  
  # Feature table, transposed, etc.
  feature_table <- t(dataSet$data.proc)
  colnames(feature_table) <- names(dataSet$enrich_ids)[
    match(colnames(feature_table), dataSet$enrich_ids)
  ]
  
  meta <- rdtSet$dataSet$meta.info
  
  # If you use some includeMeta logic
  if (!(exists("includeMeta") && length(includeMeta) == 0)) {
    meta <- meta[, c(predictor, includeMeta)]
  } else {
    meta <- meta[, predictor, drop = FALSE]
  }
 
  # Convert continuous metadata columns to numeric
  idx <- which(rdtSet$dataSet$meta.types[colnames(meta)] == "cont")
  if (length(idx) > 0) {
    meta[, idx] <- apply(meta[, idx, drop = FALSE], 2, as.numeric)
  }
  
  # Perform the ordination + significance tests
  if (method == "RDA") {
    # RDA
    library(vegan)
    rda_result <- vegan::rda(feature_table ~ ., data = meta)
    anova_res  <- anova.cca(rda_result, step = 10000, by = "term")
    
    if (is.na(anova_res$F[1])) {
      stats_msg <- "The model is overfitted with no unconstrained (residual) component."
    } else {
      stats_msg <- paste0("[Constrained Permutational ANOVA] F-value: ",
                          round(anova_res$F[1], 3),
                          "; p-value: ", anova_res$`Pr(>F)`[1])
    }
    
    # Store results
    rdtSet$analSet$RDA_res <- list(
      method     = "RDA",
      model      = rda_result,
      meta       = meta,
      featureTbl = feature_table,
      stats      = stats_msg
    )
    
  } else {
    # PCA
    library(factoextra)
    library(vegan)
    
    pca_result <- prcomp(feature_table, scale = FALSE)
    data_dist  <- dist(feature_table, method = "euclidean")
    adonis_res <- adonis2(data_dist ~ ., data = meta, permutations = 999)

    if (is.na(adonis_res$F[1])) {
      stats_msg <- "The model is overfitted with no unconstrained (residual) component."
    } else {
      stats_msg <- paste0("[PERMANOVA] F-value: ",
                          round(adonis_res$F[1], 5),
                          "; R-squared: ", signif(adonis_res$R2, 5),
                          "; p-value: ", adonis_res$`Pr(>F)`[1])
    }

    # Store results
    rdtSet$analSet$PCA_res <- list(
      method     = "PCA",
      model      = pca_result,
      meta       = meta,
      featureTbl = feature_table,
      stats      = stats_msg
    )
  }
  
  # Save it back
  .set.rdt.set(rdtSet)
  # Return some summary stats if you like:
  return(rdtSet$analSet$ordination_res$stats)
}




 
PlotBiplot <- function( method,topN=10,fileName = "biplot", format = "png",dpi = 300,colorGradient="d3") {
  
 rdtSet  <- .get.rdt.set()
  resnm <- paste0(method,"_res")

  if (is.null(rdtSet$analSet[[resnm]])) {
    stop("No ordination results found! Please run PerformOrdination() first.")
  }

  ordRes  <- rdtSet$analSet[[resnm]]

 
  meta          <- ordRes$meta
  feature_table <- ordRes$featureTbl
  stats         <- ordRes$stats 

color_scale <- if (rdtSet$dataSet$meta.types[colnames(meta)[1]] == "cont") {
    if (colorGradient == "gray") {
      scale_color_gradientn(colors = c("grey90", "grey10"))
    } else if (colorGradient == "byr") {
      scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(10, "RdYlBu")))
    } else if (colorGradient == "viridis") {
      # Using the built-in scale_color_viridis_c
      scale_color_viridis_c(option = "viridis")
    } else if (colorGradient == "plasma") {
     scale_color_viridis_c(option = "plasma")
    } else if (colorGradient == "npj") {
      scale_color_gradientn(colors = c("#00A087FF", "white", "#E64B35FF"))
    } else if (colorGradient == "aaas") {
      scale_color_gradientn(colors = c("#4DBBD5FF", "white", "#E64B35FF"))
    } else if (colorGradient == "d3") {
      scale_color_gradientn(colors = c("#2CA02CFF", "white", "#FF7F0EFF"))
    } else {
      # default
      scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(10, "RdBu")))
    }
} else {
   if (colorGradient == "gray") {
        scale_color_manual(values = c("grey90","grey70","grey50","grey30","grey10"))
    } else if (colorGradient == "byr") {
       scale_color_brewer(palette = "RdYlBu")
    } else if (colorGradient == "viridis") {
      # Discrete viridis
       scale_color_viridis_d(option = "viridis")
    } else if (colorGradient == "plasma") {
      scale_color_viridis_d(option = "plasma")
    } else if (colorGradient == "npj") {
       ggsci::scale_color_npg()
    } else if (colorGradient == "aaas") {
      ggsci::scale_color_aaas()
    } else if (colorGradient == "d3") {
      ggsci::scale_color_d3("category10")
    } else {
      # default
     scale_color_brewer(palette="Set1")
    }
}
 
  library(ggplot2)
  library(ggrepel)
  print(method)

  if (method == "RDA") {
    # Extract RDA model, summary, etc.
    rda_obj <- ordRes$model
    ii      <- summary(rda_obj)

     # Coordinates
    sp <- as.data.frame(ii$species[, 1:2, drop = FALSE])
    st <- as.data.frame(ii$sites)
    yz <- as.data.frame(ii$biplot[, 1:2, drop = FALSE])


   if(names(st)[2]=="PC1"){
     ylab="PC1"
     names(st)[2]<- names(sp)[2] <- names(yz)[2] <- "RDA2"
    }else{
     ylab="RDA2"  
    }
    # Decide topN features by their loadings
    importance   <- apply(abs(sp), 1, sum)
    top_features <- names(sort(importance, decreasing = TRUE))[1:topN]

    # Scale arrows
    scaling_factor <- max(abs(st$RDA1), abs(st$RDA2)) / max(abs(sp$RDA1), abs(sp$RDA2))
    sp <- sp[rownames(sp) %in% top_features, ] * scaling_factor
    
    scaling_factor2 <- max(abs(st$RDA1), abs(st$RDA2)) / max(abs(yz$RDA1), abs(yz$RDA2))
    yz <- yz * scaling_factor2
    
    # If your predictor is the 1st (or only) column in meta
    predictor <- colnames(meta)[1]

    p <- ggplot() +
      geom_point(data = st,
                 aes(x = RDA1, y = RDA2, color = meta[, predictor]),
                 size = 3, alpha = 0.8) +
      color_scale +
      labs(color = predictor) +
      theme_minimal() +
      # Feature arrows
      geom_segment(data = sp,
                   aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                   arrow = arrow(angle=10, length = unit(0.3, "cm"), type = "closed"),
                   linetype=1, size=0.3, colour="#000078ff", alpha=1) +
      geom_text_repel(data = sp,
                      aes(RDA1, RDA2, label = row.names(sp)),
                      size=3, color="#000078ff", segment.color='grey',
                      fontface='italic') +
      # Biplot arrows
      geom_segment(data = yz,
                   aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                   arrow = arrow(angle=10, length = unit(0.3, "cm"), type = "closed"),
                   linetype=1, size=0.3, colour="#006000ff", alpha=1) +
      geom_text_repel(data = yz,
                      aes(RDA1, RDA2, label=row.names(yz)),
                      size=3, color="#006000ff", segment.color='grey') +
      labs(x = paste0("RDA1 (", format(100*ii$cont[[1]][2,1], digits=4), "%)"),
           y = paste0("RDA2 (", format(100*ii$cont[[1]][2,2], digits=4), "%)"),
           title ="") +
      geom_hline(yintercept=0, linetype=3, size=0.5, color='gray') +
      geom_vline(xintercept=0, linetype=3, size=0.5, color='gray') +
      theme(
        panel.border = element_rect(color="black", fill=NA),
        axis.line = element_blank()
      )

  } else {
 
    # PCA route
    library(dplyr)
    pca_obj <- ordRes$model
    predictor <- colnames(meta)[1]
    
    # Identify topN via PC1+PC2 loadings
    contrib <- factoextra::get_pca_var(pca_obj)$contrib
    contrib_pc1_pc2 <- contrib[,1] + contrib[,2]
    top_features <- names(sort(contrib_pc1_pc2, decreasing=TRUE))[1:topN]

    # Scores for individuals
    scores  <- pca_obj$x
    choices <- c(1, 2)
    lam     <- pca_obj$sdev[choices]
    n       <- NROW(scores)
    lam     <- lam * sqrt(n)
    ind_data <- data.frame(
      PC1 = scores[, choices[1]] / lam[1],
      PC2 = scores[, choices[2]] / lam[2]
    )

    # Loadings for variables
    loadings <- pca_obj$rotation
    var_data <- data.frame(
      PC1 = loadings[, choices[1]] * lam[1],
      PC2 = loadings[, choices[2]] * lam[2]
    ) %>%
      filter(rownames(.) %in% top_features)
    
    # Scale the arrows so they fit
    scaling_factor <- max(abs(ind_data$PC1), abs(ind_data$PC2)) /
                      max(abs(var_data$PC1), abs(var_data$PC2))
    var_data$PC1 <- var_data$PC1 * scaling_factor
    var_data$PC2 <- var_data$PC2 * scaling_factor
    
    # Build plot
    p <- ggplot() +
      geom_point(data=ind_data,
                 aes(x=PC1, y=PC2, color=meta[, predictor]),
                 size=3, alpha=0.8) +
      color_scale +
      labs(color = predictor) +
      theme_minimal() +
      geom_segment(data = var_data,
                   aes(x=0, y=0, xend=PC1, yend=PC2),
                   arrow = arrow(angle=10, length = unit(0.2,"cm"), type="closed"),
                   linetype=1, size=0.4, colour="#000078ff", alpha=1) +
      geom_text_repel(data = var_data,
                      aes(PC1, PC2, label=rownames(var_data)),
                      size=3, color="#000078ff", segment.color='grey',
                      fontface='italic') +
      labs(
        x = paste0("PC1 (",
                   round(summary(pca_obj)$importance[2,1] * 100, 2), "%)"),
        y = paste0("PC2 (",
                   round(summary(pca_obj)$importance[2,2] * 100, 2), "%)"),
        title = ""
      ) +
      geom_hline(yintercept=0, linetype=3, size=0.5, color='gray') +
      geom_vline(xintercept=0, linetype=3, size=0.5, color='gray') +
      theme(
        panel.border = element_rect(color="black", fill=NA),
        axis.line = element_blank()
      )
  }

  # --- Save the plot ---
  imgName <- paste0(fileName, "dpi", dpi, ".", format)
  Cairo::Cairo(file = imgName,
               type = format,
               dpi = dpi,
               width = 7.5,
               height = 6,
               units = "in",
               bg = "white")
  print(p)
  dev.off()

  # Save the final coordinates if you want
  # For RDA: st is summary(rda_obj)$sites
  # For PCA: st is pca_obj$x
  # We'll store them generically as ordination coords
  if (method == "RDA") {
    st <- summary(ordRes$model)$sites
  } else {
    st <- ordRes$model$x
  }

  rdtSet$analSet$biplot_method     <- method
  rdtSet$analSet$biplot_fileName   <- "biplot_ordination.csv"
  rdtSet$analSet$biplot_ordination <- st

  fast.write.csv(st, file = "biplot_ordination.csv")
  .set.rdt.set(rdtSet)

  # Return the overall stats from the ordination
  return(stats)
}
