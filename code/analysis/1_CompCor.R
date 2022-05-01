## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("0_SharedCode.R")
stopifnot(SharedCode_version == c(11,0))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
time <- Sys.time()

for (ii in seq(nrow(iters))) {
  subject <- iters[ii, "subject"]
  acquisition <- as.character(iters[ii, "acquisition"])
  test <- iters[ii, "test"]
  visit <- iters[ii, "visit"]

  # Get output file.
  cc_fname <- file.path(
    dir_CompCor, 
    paste0(subject, "_v", visit + (!test)*2, "_", acquisition, ".rds")
  )
  if (file.exists(cc_fname)) { next }
  
  cat(paste0(
    "Subject ", subject, ", ", as.character(acquisition), " ", 
    ifelse(test, "test", "retest"), " ", visit, 
    " (", ii, " of ", nrow(iters), ")", "\n"
  ))
  
  # Get input files.
  fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)
  data_dir <- file.path(subject, "MNINonLinear", "Results", fname_prefix)
  fnames <- list(
    NIFTI = file.path(data_dir, paste0(fname_prefix, ".nii.gz")),
    NIFTI_labs = file.path(data_dir, "../../ROIs/Atlas_wmparc.2.nii.gz")
  )
  
  # If retest, the data needs to be loaded.
  if (!test) {
    data_zip_MPP <- file.path(
      dir_HCP_retest_archive, paste0(subject, "_3T_rfMRI_REST", visit, "_preproc.zip")
    )
    
    if (!file.exists(file.path(dir_HCP_retest, fnames$NIFTI))) {
      cmd <- paste("unzip", data_zip_MPP, fnames$NIFTI, "-d", dir_HCP_retest)
      system(cmd)
      stopifnot(file.exists(file.path(dir_HCP_retest, fnames$NIFTI)))
    }
  }
  
  # Compute and save CompCor results.
  saveRDS(
    fMRIscrub::CompCor_HCP(
      nii=file.path(ifelse(test, dir_HCP_test, dir_HCP_retest), fnames$NIFTI), 
      nii_labels=file.path(dir_HCP_test, fnames$NIFTI_labs),
      ROI_noise=c("wm_cort", "csf", "wm_cblm"),
      idx=seq(nDrop+1, hcp_T), noise_nPC=10, noise_erosion=c(2,1,2),
      verbose=TRUE
    ), 
  cc_fname)
  
  # If retest, the data needs to be unloaded.
  if (!test) {
    unlink(file.path(dir_HCP_retest, fnames$NIFTI))
  }
  
  print(Sys.time() - time)
  time <- Sys.time()
}

