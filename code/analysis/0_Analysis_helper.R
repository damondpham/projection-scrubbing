get_ICC31_decomp <- function(x) {
  # x is measurements by subjects by variables
  # Get data dimensions.
  d <- dim(x)
  if (length(d) == 2) { x <- array(x, dim=c(dim(x), 1)); d <- dim(x) }
  stopifnot(length(d) == 3)
  nM <- d[1]
  nS <- d[2]

  sub_mean <- colMeans(x, na.rm=TRUE)
  visit_mean <- apply(x, c(1,3), mean, na.rm=TRUE)
  grand_mean <- colMeans(sub_mean, na.rm=TRUE)
  grand_mean2 <- array(rep(grand_mean, each=nM*nS), dim=dim(x))

  # SST <- apply((x - grand_mean2)^2, 3, sum, na.rm=TRUE)
  SSW <- apply((x - rep(sub_mean, each=nM))^2, 3, sum, na.rm=TRUE)
  SSB <- nM * colSums((sub_mean - rep(grand_mean, each=nS))^2, na.rm=TRUE)
  SSM <- nS * colSums((visit_mean - rep(grand_mean, each=nM))^2, na.rm=TRUE)
  SSR <- SSW - SSM

  # MST <- SST / (nM*nS - 1)
  # MSW <- SSW / ((nM-1)*(nS))
  MSB <- SSB / (nS-1)
  # MSM <- SSM / (nM-1)
  MSR <- SSR / ((nM-1)*(nS-1))

  ICC31 <- (MSB - MSR)/(MSB + (nM - 1) * MSR)
  cbind(MSB=MSB, MSR=MSR, ICC31=ICC31)
}

scrub_flag_overlap <- function(base) {
  pairs <- list(
    c("LevICA", "FD"),
    c("LevICA", "modFD"),
    c("LevICA", "DVARS"),
    c("LevICA", "LevFusedPCA"),
    c("FD", "DVARS"),
    c("modFD", "DVARS")
  )

  overlap <- array(NA, dim=c(8, length(subjects), length(pairs), 3))
  dimnames(overlap) <- list(
    c(paste0("v", seq(4), "_LR"), paste0("v", seq(4), "_RL")),
    subjects,
    vapply(pairs, function(x){paste(x, collapse="__")}, ""),
    c("first", "second", "both")
  )
  iters <- expand.grid(visit=seq(2), test=c(TRUE, FALSE), acquisition=c("LR", "RL"), subject=subjects)
  for (ii in seq(nrow(iters))) {
    # Get iteration info
    subject <- iters[ii, "subject"]
    acquisition <- iters[ii, "acquisition"]
    test <- iters[ii, "test"]
    visit <- iters[ii, "visit"]

    suffix <- paste0(subject, "_v", visit + (!test)*2, "_", acquisition)
    cat(suffix, "\n")

    scan <- list(
      Lev = file.path(dir_scrubMeas, "Lev", base, paste0("LEV_", suffix, ".rds")),
      DVARS = file.path(dir_scrubMeas, "DVARS", base, paste0("DVARS_", suffix, ".rds")),
      FD = file.path(dir_scrubMeas, "FD", paste0("FD_", suffix, ".rds"))
    )
    scan <- lapply(scan, readRDS)
    scan <- cbind(
      LevICA = scan$Lev$measure$ICA_kurt > median(scan$Lev$measure$ICA_kurt) * 5,
      LevFusedPCA = scan$Lev$measure$fusedPCA_kurt > median(scan$Lev$measure$fusedPCA_kurt) * 5,
      LevPCA = scan$Lev$measure$PCA_kurt > median(scan$Lev$measure$PCA_kurt) * 5,
      DVARS = (scan$DVARS$dv$DPD > 5) & (scan$DVARS$dv$ZD > qnorm(1-.05/(hcp_T-nDrop))),
      FD = scan$FD$og > .5,
      modFD = scan$FD$og_nfc_l4 > .5
    )
    colnames(scan) <- c("LevICA", "LevFusedPCA", "LevPCA", "DVARS", "FD", "modFD")
    scan <- lapply(pairs, function(p){scan[,p[1]] + scan[,p[2]] * 2})
    overlap[visit + test*2 + (acquisition=="RL")*4, which(subjects==subject),,] <- c(
      do.call(rbind, lapply(scan, function(x){c(sum(x==1), sum(x==2), sum(x==3))}))
    )
  }
  overlap
}

scrub_flag_overlap2 <- function(base) {
  pairs <- list(
    c("LevICA5", "DVARS"),
    c("LevICA4", "DVARS"),
    c("LevICA3", "DVARS"),
    c("LevICA2", "DVARS")
  )
  
  overlap <- array(NA, dim=c(8, length(subjects), length(pairs), 3))
  dimnames(overlap) <- list(
    c(paste0("v", seq(4), "_LR"), paste0("v", seq(4), "_RL")),
    subjects,
    vapply(pairs, function(x){paste(x, collapse="__")}, ""),
    c("first", "second", "both")
  )
  iters <- expand.grid(visit=seq(2), test=c(TRUE, FALSE), acquisition=c("LR", "RL"), subject=subjects)
  for (ii in seq(nrow(iters))) {
    # Get iteration info
    subject <- iters[ii, "subject"]
    acquisition <- iters[ii, "acquisition"]
    test <- iters[ii, "test"]
    visit <- iters[ii, "visit"]
    
    suffix <- paste0(subject, "_v", visit + (!test)*2, "_", acquisition)
    cat(suffix, "\n")
    
    scan <- list(
      Lev = file.path(dir_scrubMeas, "Lev", base, paste0("LEV_", suffix, ".rds")),
      DVARS = file.path(dir_scrubMeas, "DVARS", base, paste0("DVARS_", suffix, ".rds"))
    )
    scan <- lapply(scan, readRDS)
    scan <- cbind(
      LevICA5 = scan$Lev$measure$ICA_kurt > median(scan$Lev$measure$ICA_kurt) * 5,
      LevICA4 = scan$Lev$measure$ICA_kurt > median(scan$Lev$measure$ICA_kurt) * 4,
      LevICA3 = scan$Lev$measure$ICA_kurt > median(scan$Lev$measure$ICA_kurt) * 3,
      LevICA2 = scan$Lev$measure$ICA_kurt > median(scan$Lev$measure$ICA_kurt) * 2,
      DVARS = (scan$DVARS$dv$DPD > 5) & (scan$DVARS$dv$ZD > qnorm(1-.05/(hcp_T-nDrop)))
    )
    colnames(scan) <- c("LevICA5", "LevICA4", "LevICA3", "LevICA2", "DVARS")
    scan <- lapply(pairs, function(p){scan[,p[1]] + scan[,p[2]] * 2})
    overlap[visit + test*2 + (acquisition=="RL")*4, which(subjects==subject),,] <- c(
      do.call(rbind, lapply(scan, function(x){c(sum(x==1), sum(x==2), sum(x==3))}))
    )
  }
  overlap
}

fingerprint <- function(FC) {
  visits <- dimnames(FC)[[1]]
  visit_pairs <- list(
    LRt = c("test_LR1", "test_LR2"),
    RLt = c("test_RL1", "test_RL2"),
    LRr = c("retest_LR1", "retest_LR2"),
    RLr = c("retest_RL1", "retest_RL2")
  )
  fp_sim_meas <- "pearson"
  c_diags <- c("on", "off")

  fp <- array(NA, dim=c(
    length(labs)*length(c_diags),
    length(visit_pairs),
    length(subjects)*2
  ))
  dimnames(fp) <- list(
    outer(names(labs), c_diags, paste, sep="_"),
    names(visit_pairs),
    c(paste0("1s", subjects), paste0("2s", subjects))
  )

  # For each connection type (which network/subcortex, and on/off-diagonal)
  for (ll in seq(length(labs))) {
    for (dd in seq(length(c_diags))) {
      lab <- names(labs)[ll]
      c_diag <- c_diags[dd]
      cmask <- labs[[lab]][[c_diag]]
      cmask_name <- paste(lab, c_diag, sep="_")

      # For each visit pair
      for (vv in seq(length(visit_pairs))) {
        visit_pair <- visit_pairs[[vv]]

        # Get the scan-scan similarities.
        sim <- matrix(NA, nSubjects, nSubjects)
        for (ss in seq(nSubjects)) {
          sim[ss,] <- apply(
            FC[which(visits==visit_pair[2]),,cmask],
            1, cor,
            y=FC[which(visits==visit_pair[1]),ss,cmask], method=fp_sim_meas
          )
        }
        # Get the fingerprinting ranks.
        rank1 <- rank2 <- vector("numeric", nSubjects)
        for (ss in seq(nSubjects)) {
          rank1[ss] <- rank(-sim[ss,], ties.method="min")[ss]
          rank2[ss] <- rank(-sim[,ss], ties.method="min")[ss]
        }

        fp[cmask_name, vv, ] <- c(c(rank1 == 1), c(rank2 == 1))
      }
    }
  }
  fp
}

qc_fc <- function(x){
  # match scan index
  x <- x[c(1,3,5,7,2,4,6,8),,,,drop=FALSE]
  n1 <- dimnames(x[[1]])
  n2 <- apply(meanFD[seq(8),], 1, function(x){paste0(ifelse(x["test"], "test", "retest"), "_", x["acquisition"], x["visit"])})

  # match subjects
  x <- x[,as.character(meanFD[seq(42)*8,]$subject),,,drop=FALSE]

  # calculate correlation
  r <- apply(x, c(3,4), function(a,b){cor(c(a), b)}, b=meanFD$FDmean)

  # calculate p-value
  t_from_r <- function(r, n){r * sqrt((42*8 - 2) / (1 - (r^2)))}
  tvals <- t_from_r(r, 42*8-2)
  pvals <- pt(-abs(tvals), 42*8-2) * 2
  # tp[] <- p.adjust(tp, "BH")

  return(list(r=r, pvals=pvals))
}

pbeta_make <- function(FC, what, parallel=FALSE) {
  
  candidate_alphas <- c(0, .25, .5, .75, 1)
  
  # Folds
  folds <- vector("numeric", nrow(dg))
  n_folds <- 8
  for (sex in c("F", "M")) {
    for (age_group in levels(dg$Age)) {
      idx <- which(dg$Gender == sex & dg$Age == age_group)
      n_group <- length(idx)
      folds[idx] <- c(rep(seq(n_folds), ceiling(n_group/n_folds)))[seq(length(idx))]
    }
  }
  cat("Folds:\n")
  print(table(folds))
  
  split <- "subjects"
  FC <- matrix(FC, ncol=dim(FC)[3]) # as matrix
  
  # Sex ------------------------------------------------------------------------
  if (what == "Sex") {
    cat("Sex", "\n")
    alpha_error_best <- Inf; q_best <- NULL; alpha_best <- NULL
    for (jj in seq(length(candidate_alphas))) {
      alpha <- candidate_alphas[jj]
      cat("Trying alpha=", alpha, "\n")
      q <- cv.glmnet(
        FC, dg2$Gender, 
        family="binomial", foldid=rep(folds, each=nVisits), 
        alpha=alpha, intercept=TRUE,
        parallel=parallel
      )
      alpha_error <- q$cvm[which(q$lambda == q$lambda.1se)]
      if (alpha_error < alpha_error_best) { 
        q_best <- q; alpha_error_best <- alpha_error; alpha_best <- alpha
      }
    }
    cat("Best alpha:", alpha_best, "\n")
    lambda <- q_best$lambda.1se
    return(list(
      betas=q_best$glmnet.fit$beta[,which(q_best$lambda == lambda)], 
      lambda=lambda,
      alpha=alpha_best
    ))
  }
  
  # Cog_* ----------------------------------------------------------------------
  else {
    candidate_alphas=c(.25, .5, .75, 1)
    cat(what, "\n")
    sex <- switch(what, Cog_M="M", Cog_F="F")
    yvar <- paste0("Cog_", sex)
    cat(yvar, "\n")
    sex_mask <- dg$Gender == sex & (!is.na(dg$CogTotalComp_Unadj))
    sex_mask2 <- rep(sex_mask, each=nVisits)
    alpha_error_best <- Inf; q_best <- NULL
    ti <- Sys.time()
    for (jj in seq(length(candidate_alphas))) {
      alpha <- candidate_alphas[jj]
      cat("Trying alpha=", alpha, "\n")
      q <- cv.glmnet(
        FC[sex_mask2,], dg2[[what]][sex_mask2], parallel=parallel,
        family="gaussian", foldid=rep(folds[sex_mask], each=nVisits), alpha=alpha
      )
      alpha_error <- q$cvm[which(q$lambda == q$lambda.1se)]
      if (alpha_error < alpha_error_best) { 
        q_best <- q; alpha_error_best <- alpha_error; alpha_best <- alpha
      }
    }
    cat("Best alpha:", alpha_best, "\n")
    lambda <- q_best$lambda.1se
    return(list(
      betas=q_best$glmnet.fit$beta[,which(q_best$lambda == lambda)], 
      lambda=lambda,
      alpha=alpha_best
    ))
  }
}

pred_make <- function(FC, ii=1, what, parallel=FALSE) {
  
  # Folds
  folds <- foldss <- vector("numeric", nrow(dg))
  n_folds <- 8
  for (sex in c("F", "M")) {
    for (age_group in levels(dg$Age)) {
      idx <- which(dg$Gender == sex & dg$Age == age_group)
      n_group <- length(idx)
      folds[idx] <- c(rep(seq(n_folds), ceiling(n_group/n_folds)))[seq(length(idx))]
      foldss[idx] <- with(set.seed(ii), sample(folds[idx]))
    }
  }
  
  FC <- matrix(FC, ncol=dim(FC)[3]) # as matrix
  
  preds <- vector("numeric", nSubjects*nVisits) * NA
  alpha_best <- vector("numeric", n_folds) * NA
  #colnames(preds) <- names(alpha_best) <- c("Sex", "Cog_F", "Cog_M")
  
  # Sub-folds
  for (fold_id in seq(n_folds)) {
    cat("--------------------\n")
    cat("Subfold\t", fold_id, "\n")
    cat("--------------------\n")
    idx_train <- which(foldss != fold_id)
    idx_test <- which(foldss == fold_id)
    # use same folds for lambda CV
    foldss2 <- rep(as.numeric(factor(foldss[idx_train])), each=nVisits)
    idx_train2 <- which(rep(foldss, each=4) != fold_id)
    idx_test2 <- which(rep(foldss, each=4) == fold_id)
    
    candidate_alphas <- c(0, .25, .5, .75, 1)
    
    # Sex ----------------------------------------------------------------------
    if (what == "Sex") {
      cat("Sex\n")
      alpha_error_best <- Inf; q_best <- NULL
      ti <- Sys.time()
      for (jj in seq(length(candidate_alphas))) {
        alpha <- candidate_alphas[jj]
        cat("\tTrying alpha=", alpha, "\n")
        q <- cv.glmnet(
          FC[idx_train2,,drop=FALSE], dg2$Gender[idx_train2],
          family="binomial", foldid=foldss2, alpha=alpha, intercept=TRUE,
          parallel=parallel
        )
        alpha_error <- q$cvm[which(q$lambda == q$lambda.1se)]
        if (alpha_error < alpha_error_best) { 
          q_best <- q; alpha_error_best <- alpha_error; alpha_best[fold_id] <- alpha
        }
        print(Sys.time() - ti); ti <- Sys.time()
      }
      cat("Best alpha:", alpha_best[fold_id], "\n")
      preds[idx_test2] <- predict(
        q_best, FC[idx_test2,,drop=FALSE], type="response"
      )
    }
    
    # Cog_* --------------------------------------------------------------------
    else if (what %in% c("Cog_F", "Cog_M")) {
      candidate_alphas=c(.25, .5, .75, 1)
      cat(what, "\n")
      sex <- switch(what, Cog_M="M", Cog_F="F")
      yvar <- paste0("Cog_", sex)
      cat(yvar, "\n")
      sex_mask <- dg$Gender == sex & (!is.na(dg$CogTotalComp_Unadj))
      sex_mask2 <- rep(sex_mask, each=nVisits)
      alpha_error_best <- Inf; q_best <- NULL
      ti <- Sys.time()
      for (jj in seq(length(candidate_alphas))) {
        mask_jj <- intersect(idx_train, which(sex_mask))
        mask2_jj <- intersect(idx_train2, which(sex_mask2))
        foldss2jj <- rep(as.numeric(factor(foldss[mask_jj])), each=nVisits)
        alpha <- candidate_alphas[jj]
        cat("\tTrying alpha=", alpha, "\n")
        q <- cv.glmnet(
          FC[mask2_jj,], dg2$CogTotalComp_Unadj[mask2_jj],
          family="gaussian", foldid=foldss2jj,
          alpha=alpha, parallel=parallel
        )
        alpha_error <- q$cvm[which(q$lambda == q$lambda.1se)]
        if (alpha_error < alpha_error_best) {
          q_best <- q; alpha_error_best <- alpha_error; alpha_best[fold_id] <- alpha
        }
        print(Sys.time() - ti); ti <- Sys.time()
      }
      cat("Best alpha:", alpha_best[fold_id], "\n")
      mask_pred <- intersect(idx_test2, which(sex_mask2))
      preds[mask_pred] <- predict(
        q_best, FC[mask_pred,,drop=FALSE], type="response"
      )
    }
    
    else if (what %in% c("Vocab")) {
      candidate_alphas <- .5
      cat("Vocab\n")
      alpha_error_best <- Inf; q_best <- NULL
      ti <- Sys.time()
      for (jj in seq(length(candidate_alphas))) {
        alpha <- candidate_alphas[jj]
        cat("\tTrying alpha=", alpha, "\n")
        q <- cv.glmnet(
          FC[idx_train2,,drop=FALSE], dg2$PicVocab_Unadj[idx_train2],
          family="gaussian", foldid=foldss2, alpha=alpha, intercept=TRUE,
          parallel=parallel
        )
        alpha_error <- q$cvm[which(q$lambda == q$lambda.1se)]
        if (alpha_error < alpha_error_best) { 
          q_best <- q; alpha_error_best <- alpha_error; alpha_best[fold_id] <- alpha
        }
        print(Sys.time() - ti); ti <- Sys.time()
      }
      cat("Best alpha:", alpha_best[fold_id], "\n")
      preds[idx_test2] <- predict(
        q_best, FC[idx_test2,,drop=FALSE], type="response"
      )
    }
    
    else { stop("What is ", what, "?") }
  }
  list(preds=preds, alpha_best=alpha_best)
}
