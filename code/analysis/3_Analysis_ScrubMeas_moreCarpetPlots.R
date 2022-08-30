## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("0_SharedCode.R")
source("0_Analysis_helper.R")
stopifnot(SharedCode_version == c(11,0))

cat("Newest version :)")

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
splits <- c(
  round(seq(0, 10)/10*nrow(meanFD)),
  round((seq(0, 9)*2+1)/20*nrow(meanFD))
)
splits <- pmax(1, splits)
subs <- meanFD[round(splits),]
baseNames <- c("MPP", "DCT4", "36P")

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------

for (ii in seq(nrow(subs))) {
  # Get iteration info ---------------------------------------------------------
  subject <- subs[ii, "subject"]
  acquisition <- subs[ii, "acquisition"]
  test <- subs[ii, "test"]
  visit <- subs[ii, "visit"]

  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat(paste0(
    "Subject ", subject, ", ", as.character(acquisition), " ", 
    ifelse(test, "test", "retest"), " ", visit, 
    " (", ii, " of ", nrow(subs), ")", "\n"
  ))
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  
  gp_fname <- paste0(
    subject, "_v", visit + (!test)*2, "_", acquisition, 
    "_carpetplot_", baseNames, ".pdf"
  )
  if (all(file.exists(file.path(dir_carpetPlots, gp_fname)))) { next }

  # Get nuisance regressors and CIFTI data. ------------------------------------
  cat("Input files.\n")
  # Mean signals
  ms_fname <- file.path(
    dir_meanSignals, 
    paste0(
      subject, "_v", visit + (!test)*2, "_", acquisition, ".rds"
    )
  )
  ms <- data.frame(readRDS(ms_fname))
  ms <- as.matrix(ms[,c("wm", "csf", "wholebrain")])

  # CompCor
  cc_fname <- file.path(
    dir_CompCor, 
    paste0(
      subject, "_v", visit + (!test)*2, "_", acquisition, ".rds"
    )
  )
  cc <- readRDS(cc_fname)
  cc <- cbind(cc$PCs$wm_cort[,seq(5)], cc$PCs$csf[,seq(5)], cc$PCs$wm_cblm[,seq(5)])
  cc <- scale(cc)
  
  # CIFTI
  fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)
  
  if (COMPUTER == "RED") {
    data_dir <- file.path(subject, "MNINonLinear", "Results", fname_prefix)
    fnames <- list(
      CIFTI = file.path(data_dir, paste0(fname_prefix, "_Atlas.dtseries.nii")),
      RP = file.path(data_dir, "Movement_Regressors.txt")
    )
    cii_fname <- file.path(ifelse(test, dir_HCP_test, dir_HCP_retest), fnames$CIFTI)
  
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
  
    rp <- read.table(file.path(ifelse(test, dir_HCP_test, dir_HCP_retest), fnames$RP))
  
  # on Personal Computer -------------------------------------------------------
  } else {
    cii_fname <- file.path(dir_data_misc, "SelectedScans", paste0(subject, "_", fname_prefix, "_Atlas.dtseries.nii"))
    rp <- read.table(file.path(dir_data_misc, "SelectedScans", paste0(subject, "_", fname_prefix, "_Movement_Regressors.txt")))
  }
  
  # ----------------------------------------------------------------------------
  
  p6 <- scale(rp[seq(nDrop+1, hcp_T),seq(6)])
  p36 <- scale(cbind(ms, rbind(0, diff(ms)), rp))
  p36 <- scale(cbind(p36, p36^2))[seq(nDrop+1, hcp_T),]
    
  cii_list <- setNames(vector("list", length(baseNames)), baseNames)
  cii0 <- t(as.matrix(read_cifti(cii_fname, brainstructures="all")))[seq(nDrop+1, hcp_T),]

  time <- Sys.time()
  
  # Get cleaned CIFTI data and leverage from each denoising method -------------
  cat("CIFTI data and leverage.\n")
  for (baseName in baseNames) {

    # Nuisance regression 
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
    cii_b <- nuisance_regression(cii0, nreg)

    gp_fname <- paste0(
      subject, "_v", visit + (!test)*2, "_", acquisition, 
      "_carpetplot_", baseName, ".pdf"
    )
    if (!file.exists(file.path(dir_carpetPlots, gp_fname))) {
      carpetplot(
        cii_b, 
        fname=file.path(dir_carpetPlots, gp_fname),
        width=16, height=3
      ) 
    }

    print(Sys.time() - time)
    time <- Sys.time()
    cat("\n") 
  }
  
  # Unload retest data.
  if (!test && COMPUTER=="RED") {
    unlink(file.path(dir_HCP_retest, fnames$CIFTI))
    unlink(file.path(dir_HCP_retest, fnames$RP))
  }
}