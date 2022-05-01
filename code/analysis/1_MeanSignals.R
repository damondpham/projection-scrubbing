## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
time <- Sys.time()

for (ii in seq(nrow(iters))) {
  subject <- iters[ii, "subject"]
  acquisition <- as.character(iters[ii, "acquisition"])
  test <- iters[ii, "test"]
  visit <- iters[ii, "visit"]

  # Get output file.
  ms_fname <- file.path(
    dir_meanSignals, 
    paste0(subject, "_v", visit + (!test)*2, "_", acquisition, ".rds")
  )
  if (file.exists(ms_fname)) { next }
  
  cat(paste0(
    "Subject ", subject, ", ", as.character(acquisition), " ", 
    ifelse(test, "test", "retest"), " ", visit, 
    " (", ii, " of ", nrow(iters), ")", "\n"
  ))
  
  # Get input files.
  fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)
  data_dir <- file.path(subject, "MNINonLinear", "Results", fname_prefix)
  fnames <- list(
    CIFTI = file.path(data_dir, paste0(fname_prefix, "_Atlas.dtseries.nii")),
    NIFTI = file.path(data_dir, paste0(fname_prefix, ".nii.gz")),
    NIFTI_labs = file.path(data_dir, "../../ROIs/Atlas_wmparc.2.nii.gz")
  )
  
  # If retest, the data needs to be loaded.
  if (!test) {
    data_zip_MPP <- file.path(
      dir_HCP_retest_archive, paste0(subject, "_3T_rfMRI_REST", visit, "_preproc.zip")
    )
    
    for (fname in fnames[c("CIFTI", "NIFTI")]) {
      if (!file.exists(file.path(dir_HCP_retest, fname))) {
        cmd <- paste("unzip", data_zip_MPP, fname, "-d", dir_HCP_retest)
        system(cmd)
        stopifnot(file.exists(file.path(dir_HCP_retest, fname)))
      }
    }
  }
  
  nii <- file.path(ifelse(test, dir_HCP_test, dir_HCP_retest), fnames$NIFTI)
  if (Use_S1200 && !file.exists(nii)) { cat("missing NIFTI\n"); stop }
  nii <- RNifti::readNifti(nii)
  niiLabs <- file.path(dir_HCP_test, fnames$NIFTI_labs)
  if (Use_S1200 && !file.exists(niiLabs)) { cat("missing NIFTI labels\n"); stop }
  niiLabs <- RNifti::readNifti(niiLabs)
  noiseROIs <- fMRIscrub:::get_NIFTI_ROI_masks(file.path(dir_HCP_test, fnames$NIFTI_labs), c("wm_cort", "csf", "wm_cblm"))
  noiseROIs <- list(wm=noiseROIs$wm_cort | noiseROIs$wm_cblm, csf=noiseROIs$csf, wholebrain=niiLabs>0)
  noise_erosion <- list(wm=2, csf=1, wholebrain=0)
  for (rr in seq(length(noiseROIs))) {
    noiseROIs[[rr]][,,] <- as.logical(noiseROIs[[rr]]) * 1
    noiseROIs[[rr]][,,] <- erode_vol(noiseROIs[[rr]], noise_erosion[[rr]], c(-1, 0, NA, NaN))
    noiseROIs[[rr]] <- apply(matrix(nii[noiseROIs[[rr]] > 0], ncol=1200), 2, mean)
  }
  rm(nii)
  
  cii <- file.path(ifelse(test, dir_HCP_test, dir_HCP_retest), fnames$CIFTI)
  if (Use_S1200 && !file.exists(cii)) { cat("missing CIFTI\n"); stop }
  cii <- tryCatch({read_xifti(cii, brainstructures="all", flat=TRUE)}, error=function(cond){NA})
  if (identical(cii, NA)) { cat("\tBad CIFTI"); stop }
  noiseROIs$cii <- colMeans(as.matrix(cii))
  
  saveRDS(noiseROIs, ms_fname)
  rm(cii)
  
  # Unload retest data.
  if (!test) {
    unlink(file.path(dir_HCP_retest, fnames$NIFTI))
    unlink(file.path(dir_HCP_retest, fnames$CIFTI))
  }
  
  print(Sys.time() - time)
  time <- Sys.time()
}

