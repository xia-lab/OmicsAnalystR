##################################################
## Compound Name Approximate Matching Utilities
## Adapted from OmicsNet util_compound_match.R
## Performs fuzzy matching using synonym database
##################################################

#' Perform approximate compound name matching using synonym database
#' @description Given a query compound name, performs fuzzy matching against
#' a synonym database to find the best matches
#' @param q Input query compound name
#' @param cmpd.db Compound database (from compound_db.rds)
#' @param syn.db Synonym database (from syn_nms.qs)
#' @return A data frame with matched candidates (index, value, score) or NULL if no matches
#' @export
PerformCompoundApproxMatch <- function(q, cmpd.db, syn.db) {

  if(length(syn.db$syns.vec) != nrow(cmpd.db)) {
    min.size <- min(length(syn.db$syns.vec), nrow(cmpd.db));
    com.nms <- cmpd.db$name[1:min.size];
    syns.vec <- syn.db$syns.vec[1:min.size];
    syns.list <- syn.db$syns.list[1:min.size];
  } else {
    if("lipid" %in% colnames(cmpd.db)) {
      nonLipidInx <- cmpd.db$lipid == 0;
      com.nms <- cmpd.db$name[nonLipidInx];
      syns.vec <- syn.db$syns.vec[nonLipidInx];
      syns.list <- syn.db$syns.list[nonLipidInx];
    } else {
      com.nms <- cmpd.db$name;
      syns.vec <- syn.db$syns.vec;
      syns.list <- syn.db$syns.list;
    }
  }

  q.length <- nchar(q);
  s <- c(0, 0.1, 0.2);
  candidates <- NULL;

  for (j in s) {
    new.q <- q;
    if(q.length > 32){
      new.q <- substr(q, 1, 32);
    }

    matched.inx <- agrep(new.q, syns.vec, ignore.case=TRUE, max.distance=j, useBytes=TRUE);

    if(length(matched.inx) > 0) {
      candidates <- data.frame(
        index = vector(mode = "numeric", length=length(matched.inx)),
        value = vector(mode = "character", length=length(matched.inx)),
        score = vector(mode = "numeric", length=length(matched.inx)),
        stringsAsFactors = FALSE
      );

      for(n in 1:length(matched.inx)){
        nm.vec <- syns.list[[matched.inx[n]]];
        hit3.inx <- agrep(q, nm.vec, ignore.case=TRUE, max.distance=j, useBytes=TRUE);

        if(length(hit3.inx) > 0){
          hit3.nm <- vector(mode = "character", length=length(hit3.inx));
          hit3.score <- vector(mode = "numeric", length=length(hit3.inx));

          for(k in 1:length(hit3.inx)){
            idx <- hit3.inx[k];
            hit3.nm[k] <- nm.vec[idx];
            hit3.score[k] <- j + abs(nchar(nm.vec[idx]) - nchar(q)) / (10 * nchar(q));
          }

          matches2 <- c();
          if(length(grep("^[1-9a-z]{2}", q, ignore.case=TRUE)) > 0){
            matches2 <- grep(paste("^", substr(q, 1, 2), sep=""), hit3.nm, ignore.case=TRUE);
          } else if (length(grep("^[1-9a-z]", q, ignore.case=TRUE)) > 0){
            matches2 <- grep(paste("^", substr(q, 1, 1), sep=""), hit3.nm, ignore.case=TRUE);
          }

          if(length(matches2) > 0){
            hit3.score[matches2] <- hit3.score[matches2] - 0.05;
          }

          best.inx <- which(hit3.score == min(hit3.score))[1];
          candidates[n, 1] <- matched.inx[n];
          candidates[n, 2] <- com.nms[matched.inx[n]];
          candidates[n, 3] <- hit3.score[best.inx];
        }
      }

      rm.inx <- is.na(candidates[,2]) | candidates[,2] == "NA" | candidates[,2] == "";
      candidates <- candidates[!rm.inx, , drop=FALSE];
      candidates <- candidates[order(candidates[,3], decreasing=FALSE), , drop=FALSE];
      if(nrow(candidates) > 10){
        candidates <- candidates[1:10, , drop=FALSE];
      }

      if(nrow(candidates) > 0) {
        return(candidates);
      }
    }
  }

  return(NULL);
}

#' Match compound name with fuzzy matching
#' @description Wrapper function that performs compound name matching
#' @param query.vec Vector of query compound names
#' @param cmpd.db Compound database
#' @param syn.db Synonym database
#' @return List with hit.inx (indices), hit.values (matched names), and match.state (1=match, 0=no match)
#' @export
MatchCompoundNames <- function(query.vec, cmpd.db, syn.db) {

  n <- length(query.vec);
  hit.inx <- rep(0, n);
  hit.values <- rep("", n);
  match.state <- rep(0, n);

  for(i in 1:n) {
    q <- query.vec[i];
    candidates <- PerformCompoundApproxMatch(q, cmpd.db, syn.db);
    if(!is.null(candidates) && nrow(candidates) > 0) {
      best.idx <- candidates[1, 1];
      if(best.idx > 0 && best.idx <= nrow(cmpd.db)) {
        hit.inx[i] <- best.idx;
        hit.values[i] <- cmpd.db$name[best.idx];
        match.state[i] <- 1;
      }
    }
    # No per-item progress output here to avoid excessive stdout
  }

  return(list(
    hit.inx = hit.inx,
    hit.values = hit.values,
    match.state = match.state
  ));
}
