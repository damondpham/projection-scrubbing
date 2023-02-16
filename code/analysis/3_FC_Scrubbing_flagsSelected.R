## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("0_SharedCode.R")
stopifnot(SharedCode_version == c(11,0))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Parcellation -----------------------------------------------------------------
ParcMat <- readRDS(parc_fname)
cor_mask <- upper.tri(diag(parc_res2))

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
baseNames <- c("CC2MP6")

time <- Sys.time()

for (baseName in baseNames) {
  cat(baseName, "~~~~~~~~~~~~~~~~~~~~~~\n")
  
  all_flag_fname <- file.path(dir_FC, paste0("flag_", baseName, ".rds"))
  if (file.exists(all_flag_fname)) { next }
  
  # Directories --------------------------------
  fd_dir <- file.path(dir_scrubMeas, "FD")
  lev_dir <- file.path(dir_scrubMeas, "Lev", baseName)
  dvars_dir <- file.path(dir_scrubMeas, "DVARS", baseName)

  flag_fname <- file.path(dir_FC, baseName, "flag.rds")
  if (file.exists(flag_fname)) { next }
  flag <- vector("list")

  for (ii in seq(nrow(iters))) {
    # Get iteration info -------------------------------------------
    subject <- iters[ii, "subject"]
    acquisition <- as.character(iters[ii, "acquisition"])
    test <- iters[ii, "test"]
    visit <- iters[ii, "visit"]
  
    cat(paste0(
      baseName, ", Subject ", subject, ", ", acquisition, " ", 
      ifelse(test, "test", "retest"), " ", visit, 
      " (", ii, " of ", nrow(iters), ")", "\n"
    ))
        
    # Get files; skip if done --------------------------------
    suffix <- paste0(subject, "_v", visit + (!test)*2, "_", acquisition)

    # CompCor
    cc <- readRDS(file.path(dir_CompCor, paste0(suffix, ".rds")))
    cc <- scale(cbind(cc$PCs$wm_cort[,seq(5)], cc$PCs$csf[,seq(5)], cc$PCs$wm_cblm[,seq(5)]))
    
    # CIFTI
    fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)
    data_dir <- file.path(subject, "MNINonLinear", "Results", fname_prefix)
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

    read_dir <- ifelse(test, dir_HCP_test, dir_HCP_retest)
    
    rp <- read.table(file.path(ifelse(test, dir_HCP_test, dir_HCP_retest), fnames$RP))
    p6 <- scale(rp[seq(nDrop+1, hcp_T),seq(6)])

    # Read CIFTI. Drop first 15 frames. Parcellate.
    # (Is equivalent to parcellating after nuisance regression.)
    cii <- as.matrix(read_xifti(file.path(read_dir, fnames$CIFTI), brainstructures="all"))[,seq(nDrop+1, hcp_T)]
    cii <- t(cii - rowMeans(cii)) %*% ParcMat

    nreg <- switch(baseName,
      CC2MP6 = cbind(1, dct4, cc[,c(1,2,6,7,11,12)], p6)
    )
    
    # Scrubbing measures
    fd <- readRDS(file.path(fd_dir, paste0("FD_", suffix, ".rds")))
    fd <- fd[c("og", "og_nfc_l4")]
    lev <- readRDS(file.path(lev_dir, paste0("LEV_", suffix, ".rds")))
    dvars <- readRDS(file.path(dvars_dir, paste0("DVARS_", suffix, ".rds")))
  
    flag_ii <- vector("list")
    
    # Get FC values ----------------------------------------------------------  
    # FD
    cat("\n\tFD")
    FD_cuts <- seq(.1, .8, .1)
    for (FD_cut in FD_cuts) {
      for (FD_type in names(fd)) {
        this_name <- paste("FD", FD_type, FD_cut, sep="___")
        flag_ii[[this_name]] <- fd[[FD_type]] > FD_cut
      }
    }
    
    # Median leverage
    cat("\n\tLev")
    medlev_cuts <- seq(1, 8)
    for (proj in "ICA_kurt") {
      for (medlev_cut in medlev_cuts) {
        this_name <- paste("Lev", proj, medlev_cut, sep="___")
        flag_ii[[this_name]] <- lev$measure[[proj]] > medlev_cut * median(lev$measure[[proj]])
      } 
    }
    
    # DVARS dual (DPD & ZD)
    cat("\n\tDVARS dual")
    flag_ii[["DVARSdual"]] <- (dvars$dv$DPD > 5) & (dvars$dv$ZD > qnorm(1-.05/(hcp_T-nDrop))) 
    cat("\n")
    
    # Unload retest data.
    if (!test) {
      unlink(file.path(dir_HCP_retest, fnames$CIFTI))
      unlink(file.path(dir_HCP_retest, fnames$RP))
    }

    flag[[suffix]] <- flag_ii
    
    print(Sys.time() - time)
    time <- Sys.time()
  }

  saveRDS(flag, flag_fname)
}
