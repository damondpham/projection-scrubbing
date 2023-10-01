## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("0_SharedCode.R")
stopifnot(SharedCode_version == c(11,0))

baseNames <- c("CC2MP6")

iters2 <- expand.grid(
  visit=seq(2),
  test=c(TRUE, FALSE)
)

wmeanFC <- function(x){
  # Compute the number of volumes after scrubbing for each session.
  x <- lapply(x, function(y){y$nVol <- y$nScrub[2]-y$nScrub[1]; y})
  # For each session, the number of volumes will be the weight. Sum the weights.
  wsum_ii <- sum(do.call(c, lapply(x, '[[', "nVol")))
  # Compute the weighted sum of the FC values.
  x <- lapply(x, function(x){x$FC * x$nVol})
  x <- rowSums(do.call(cbind, x)) / wsum_ii
}

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
time <- Sys.time()

for (baseName in baseNames) {
  cat(baseName, "~~~~~~~~~~~~~~~~~~~~~~\n")
  
  all_FCval_fname <- file.path(dir_FCval, "../3_FCval2_union", paste0("FC_", baseName, ".rds"))
  if (file.exists(all_FCval_fname)) { next }

  agg_nVol <- array(NA, dim=c(length(subjects), nrow(iters2), 9))
  agg <- array(NA, dim=c(length(subjects), nrow(iters2), 9, 87571))

  for (nn in seq(length(subjects))) {
    subject <- subjects[nn]
    cat(paste0(baseName, ", Subject ", subject, "\n"))

    # Get session results.
    iters_nn <- iters[iters$subject==subject,]
    sess_nn <- paste0(
      iters_nn$subject, "_v", iters_nn$visit + 
      (!iters_nn$test)*2, "_", iters_nn$acquisition
    )
    gt_nn <- file.path(dir_FCval, "../3_FCval2", baseName, paste0(sess_nn, ".rds"))
    if (!all(file.exists(gt_nn))) { stop("Missing data.") }
    sess_nn <- file.path(dir_FCval, "../3_FCval2_union", baseName, paste0(sess_nn, ".rds"))
    if (!all(file.exists(sess_nn))) { stop("Missing data.") }
    sess_nn <- lapply(sess_nn, readRDS)

    # Get FC estimates from full sessions using the ground truth method.
    gt_nn <- lapply(gt_nn, readRDS)
    gt_name <- "FD___og_nfc_l4___0.2"
    gt_nn <- lapply(sess_nn, '[[', gt_name)

    for (ii in seq(nrow(iters2))) {
      visit <- iters2[ii, "visit"]
      test <- iters2[ii, "test"]
      pair_ii <- iters_nn$visit==visit & iters_nn$test == test

      # Compute ground truth.
      gt_ii <- wmeanFC(c(
        lapply(gt_nn[!pair_ii], '[[', "full"),
        lapply(gt_nn[pair_ii], '[[', "out")
      ))

      # Compute LR/RL pair estimates for each scrubbing method. 
      est_ii <- setNames(
        vector("list", length(sess_nn[[1]])), 
        names(sess_nn[[1]])
      )
      agg_nVol[nn,ii,] <- lapply(est_ii, function(q){
        q$mid$nScrub[2] - q$mid$nScrub[1]
      })

      for (s_name in names(sess_nn[[1]])) {
        sess_nn_s <- lapply(sess_nn, '[[', s_name)
        est_ii[[s_name]] <- wmeanFC(lapply(sess_nn_s[pair_ii], '[[', "mid"))
      }

      # Compute error.
      est_ii <- lapply(est_ii, function(x){ psych::fisherz(x) - psych::fisherz(gt_ii) })
      est_ii <- do.call(rbind, est_ii)
      
      # Add to results.
      agg[nn,ii,,] <- est_ii
    }
  }
  
  iters2$test <- ifelse(iters2$test, "test", "retest")
  dimnames(agg) <- list(
    subject=paste0("s", subjects),
    pair=apply(iters2, 1, paste, collapse=", "),
    scrubbing=rownames(est_ii),
    FC=NULL
  )
  dimnames(agg_nVol) <- list(
    subject=paste0("s", subjects),
    pair=apply(iters2, 1, paste, collapse=", ")
  ) 

  saveRDS(agg, all_FCval_fname)
  saveRDS(agg_nVol, gsub("FC", "nVol", all_FCval_fname))
  agg2 <- apply(agg^2, c(3,4), mean)
  saveRDS(agg2, file.path(dir_FCval, "../3_FCval2_union", paste0("FC_", baseName, "_mean.rds")))
}
