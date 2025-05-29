 ##################################################
## R scripts for OmicsAnalyst
## Description: Mediation analysis
## Author: Yao, yao.lu5@mail.mcgill.ca
###################################################

#' Mediation.Anal 
#' @param imgName image name   
#' @export


Mediation.Anal <- function(method="causal",predictor="7805",mediator="51606",outcome="",lvl1,lvl2,type1,type2){
  
  rdtSet <- .get.rdt.set();
  sel.nms <- names(mdata.all);
  dataSetList <- lapply(sel.nms, readDataset);
  types <- unlist(lapply(dataSetList,function(x) x[["type"]]))
  dataSetList <- dataSetList[match(c(type1,type2),types)]

  outcome <- unlist(strsplit(outcome,";"))
  predictor <- unlist(strsplit(predictor,";"))
  mediator <- unlist(strsplit(mediator,";"))
  
  predictor_matrix <- dataSetList[[1]]$data.proc[rownames(dataSetList[[1]]$data.proc) %in% predictor, ]
  mediator_matrix <- dataSetList[[2]]$data.proc[rownames(dataSetList[[2]]$data.proc) %in% mediator, ]
  
  input_data <- data.frame(t(rbind(predictor_matrix,mediator_matrix)))
  input_data$meta <- rdtSet$dataSet$meta.info[[outcome]][match(rownames(input_data),rownames(rdtSet$dataSet$meta.info))]
  input_data <- input_data[input_data$meta %in% c(lvl1,lvl2),]
  if(rdtSet$dataSet$meta.types[outcome]=="disc"){
    input_data$meta <- as.character(input_data$meta)
    input_data$meta[input_data$meta==lvl1] <- 0
    input_data$meta[input_data$meta==lvl2] <- 1
    input_data$meta<-as.factor(as.numeric(input_data$meta))
  }else{
    input_data$meta <- as.numeric(input_data$meta)
  }
 
 
  if(method=="causal"){
    result <- PerformCausalMediation(input_data,rdtSet$dataSet$meta.types[outcome])
    
  }else{
    
    
  }
  rdtSet$mediation<-list(method=method,,result<-result)
  
  return(1)
}

PerformCausalMediation <- function(input,outcome_type  ,  sims = 1000, boot = TRUE){
   
  
  names(input) <-c( "X", "M","Y")
  if (outcome_type == "disc" && !all(input$Y %in% c(0, 1))) {
    stop("For a binary outcome, Y must be coded 0/1.")
  }
  
  ## ---- 1. Mediator model (always linear here) ----
  model.M <- lm(M ~ X, data = input)
  
  ## ---- 2. Outcome model, conditional on type ----
  if (outcome_type == "disc") {
    model.Y <- glm(Y ~ X + M, data = input, family = binomial())
  } else {
    model.Y <- lm(Y ~ X + M, data = input)
  }
  
  ## ---- 3. Mediation analysis ----
  library(mediation)
  med.out <- mediate(model.M, model.Y,
                     treat     = "X",
                     mediator  = "M",
                     boot      = boot,
                     sims      = sims)
 
  ## ---- 4. Tidy results table ----
  if(outcome_type == "disc"){
    out_tbl <- data.frame(
      Effect = c("ACME (control)", "ACME (treated)", 
                 "ADE (control)", "ADE (treated)", 
                 "Total Effect", 
                 "Prop. Mediated (control)", "Prop. Mediated (treated)", 
                 "ACME (average)", "ADE (average)", "Prop. Mediated (average)"),
      
      Estimate = c(med.out$d0, med.out$d1,
                   med.out$z0, med.out$z1,
                   med.out$tau.coef,
                   med.out$n0, med.out$n1,
                   med.out$d.avg, med.out$z.avg, med.out$n.avg),
      
      CI.lower = c(med.out$d0.ci[1], med.out$d1.ci[1],
                   med.out$z0.ci[1], med.out$z1.ci[1],
                   med.out$tau.ci[1],
                   med.out$n0.ci[1], med.out$n1.ci[1],
                   med.out$d.avg.ci[1], med.out$z.avg.ci[1], med.out$n.avg.ci[1]),
      
      CI.upper = c(med.out$d0.ci[2], med.out$d1.ci[2],
                   med.out$z0.ci[2], med.out$z1.ci[2],
                   med.out$tau.ci[2],
                   med.out$n0.ci[2], med.out$n1.ci[2],
                   med.out$d.avg.ci[2], med.out$z.avg.ci[2],  med.out$n.avg.ci[2]),
      
      p.value = c(med.out$d0.p, med.out$d1.p,
                  med.out$z0.p, med.out$z1.p,
                  med.out$tau.p,
                  med.out$n0.p, med.out$n1.p,
                  med.out$d.avg.p, med.out$z.avg.p,  med.out$n.avg.p),
      
      stringsAsFactors = FALSE
    )
  }else{
    
  }
  
  
  return(out_tbl)
}