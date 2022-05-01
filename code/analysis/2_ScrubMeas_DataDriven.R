## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("0_SharedCode.R")
stopifnot(SharedCode_version == c(11,0))

# Run this script twice, first without all subjects, then with.
# The latter is only used for the prediction analysis.
Use_S1200 <- TRUE
if (Use_S1200){
  subjects <- readRDS(file.path(dir_data_misc, "HCP_S1200_LargeBalancedSample.rds"))
  iters <- expand.grid(
    visit=seq(2), 
    acquisition=c("LR", "RL"), 
    subject=subjects,
    test=TRUE
  )
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
baseNames <- c("CC2MP6", "CC2", "CC5", "MPP", "DCT4", "9P", "36P")
if (Use_S1200) { baseNames <- c("CC2MP6", "36P") }
no_ms <- Use_S1200 && COMPUTER=="MPro"

for (ii in seq(nrow(iters))) {
  # Get iteration info ---------------------------------------------------------
  subject <- iters[ii, "subject"]
  acquisition <- iters[ii, "acquisition"]
  test <- iters[ii, "test"]
  visit <- iters[ii, "visit"]
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat(paste0(
    "Subject ", subject, ", ", as.character(acquisition), " ", 
    ifelse(test, "test", "retest"), " ", visit, 
    " (", ii, " of ", nrow(iters), ")", "\n"
  ))
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

  # Get nuisance regressors and CIFTI data. ------------------------------------
  suffix <- paste0(subject, "_v", visit + (!test)*2, "_", acquisition)

  # Mean signals
  if (!no_ms) {
    ms <- data.frame(readRDS(file.path(dir_meanSignals, paste0(suffix, ".rds"))))
    ms <- as.matrix(ms[,c("wm", "csf", "wholebrain")]) 
  }

  fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)
  data_dir <- switch(COMPUTER,
     RED = file.path(subject, "MNINonLinear", "Results", fname_prefix),
     MPro = file.path(subject, fname_prefix)
  )
  
  # CompCor
  cc <- readRDS(switch(COMPUTER,
    RED=file.path(dir_CompCor, paste0(suffix, ".rds")),
    MPro=file.path(dir_HCP_test, data_dir, paste0(fname_prefix, "_CompCor.rds"))
  ))
  cc <- cbind(cc$PCs$wm_cort[,seq(5)], cc$PCs$csf[,seq(5)], cc$PCs$wm_cblm[,seq(5)])
  cc <- scale(cc)
  
  # CIFTI
  fnames <- list(
    CIFTI = file.path(data_dir, paste0(fname_prefix, "_Atlas.dtseries.nii")),
    RP = file.path(data_dir, "Movement_Regressors.txt")
  )
  
  # If retest, the data needs to be loaded.
  if (!test) {
    data_zip_MPP <- file.path(
      dir_HCP_retest_archive, paste0(subject, "_3T_rfMRI_REST", visit, "_preproc.zip")
    )
    
    for (fname in fnames) {
      if (!file.exists(file.path(dir_HCP_retest, fname))) {
        cmd <- paste("unzip", data_zip_MPP, fname, "-d", dir_HCP_retest)
        system(cmd)
        stopifnot(file.exists(file.path(dir_HCP_retest, fname)))
      }
    }
  }

  fnames <- lapply(fnames, function(x){file.path(ifelse(test, dir_HCP_test, dir_HCP_retest), x)})

  if (Use_S1200) {
    if (!file.exists(fnames$RP)) { next }
    if (!file.exists(fnames$CIFTI)) { next }
  }
    
  rp <- as.matrix(read.table(fnames$RP))
  p6 <- scale(rp[seq(nDrop+1, hcp_T),seq(6)])
  if (!no_ms) {
    p36 <- scale(cbind(ms, rbind(0, diff(ms)), rp))
    p36 <- scale(cbind(p36, p36^2))[seq(nDrop+1, hcp_T),] 
  }
    
  # Read CIFTI. Drop first 15 frames. ------------------------------------------
  cat("\tReading data.\n")
  cii0 <- t(as.matrix(read_cifti(fnames$CIFTI, brainstructures="all")))[seq(nDrop+1, hcp_T),]
  cii_mean <- mean(cii0)
  #cii_mode <- fMRIscrub:::Mode(as.vector(cii))
  #cii_mode100 <- fMRIscrub:::Mode(round(as.vector(cii)/100)*100)

  time <- Sys.time()

  for (baseName in baseNames) {

    # Get output files; skip if done -------------------------------------------
    lev_fname <- file.path( 
      file.path(dir_scrubMeas, "Lev", baseName),
      paste0(
        "LEV_", subject, 
        "_v", visit + (!test)*2, "_", acquisition, ".rds"
      )
    )
    if (!dir.exists(dirname(lev_fname))) { dir.create(dirname(lev_fname)) }
    dvars_fname <- file.path(
      file.path(dir_scrubMeas, "DVARS", baseName), 
      paste0(
        "DVARS_", subject, 
        "_v", visit + (!test)*2, "_", acquisition, ".rds"
      )
    )
    if (file.exists(lev_fname) && file.exists(dvars_fname)) { next }

    # Nuisance regression ------------------------------------------------------
    cat(baseName, "\n")
    nreg <- switch(baseName,
      CC2 = cbind(1, dct4, cc[,c(1,2,6,7,11,12)]),
      CC5 = cbind(1, dct4, cc),
      MPP = as.matrix(rep(1, hcp_T - nDrop)),
      DCT4 = cbind(1, dct4),
      `9P` = cbind(1, dct4, p36[,c(seq(3), seq(7,12))]),
      `36P` = cbind(1, dct4, p36),
      CC2MP6 = cbind(1, dct4, cc[,c(1,2,6,7,11,12)], p6)
    )
    cii <- nuisance_regression(cii0, nreg)
    
    # Projection scrubbing -----------------------------------------------------
    if (!file.exists(lev_fname) && baseName != "MPP") {
      if (Use_S1200) {
        projections <- c("ICA", "ICA_kurt")
      } else {
        projections <- c(
          "PCA", "PCA_kurt",
          "fusedPCA", "fusedPCA_kurt",
          "ICA", "ICA_kurt"
        )
      }

      cat("\tProjection scrubbing.\n")
      pscrub <- fMRIscrub:::pscrub_multi(
        cii, projection=projections,
        nuisance = NULL, get_dirs=FALSE, get_outliers=FALSE, verbose=TRUE
      )
      # Do not keep unnecessary metadata.
      pscrub$measure_info <- NULL
      # Save projection scrubbing results.
      saveRDS(pscrub, lev_fname)
    }

    # DVARS --------------------------------------------------------------------
    if (!file.exists(dvars_fname)) {
      cat("\tDVARS.\n")
      cii <- t((cii + cii_mean) / cii_mean * 1000)
      cii <- t(cii - apply(cii, 1, mean))
      dv <- DVARS(cii, normalize=FALSE)$measure
      saveRDS(list(dv=dv, grand=c(mean=cii_mean)), dvars_fname)
    }

    print(Sys.time() - time)
    time <- Sys.time()
    cat("\n")
  }
  
  # Unload retest data.
  if (!test) {
    unlink(file.path(dir_HCP_retest, fnames$CIFTI))
    unlink(file.path(dir_HCP_retest, fnames$RP))
  }
  cat("\n")
}

