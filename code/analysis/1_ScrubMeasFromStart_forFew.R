## --------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This script shows how to compute scrubbing measures directly from the HCP data
#   (without using pre-computed CompCor, etc.). It was just used for a few scans
#   and is not required to replicate the results. 

source("0_SharedCode.R")
stopifnot(SharedCode_version == c(11,0))

subjects <- c(114318, 118528, 122620)
iters <- expand.grid(
  visit=seq(1), 
  acquisition=c("LR"), 
  subject=subjects,
  test=TRUE
)

nfb <- c(.2, .5) * .72 * 2
nfc <- c(.31, .43) * .72 * 2
nfiltb <- gsignal::butter(n=10, w=nfb, type="stop")
nfiltc <- gsignal::cheby2(2, Rs=20, w=nfc, type="stop")

time <- Sys.time()

for (ii in seq(nrow(iters))) {
  subject <- iters[ii, "subject"]
  acquisition <- as.character(iters[ii, "acquisition"])
  test <- iters[ii, "test"]
  visit <- iters[ii, "visit"]

  out_fname <- file.path(dir_data_misc, "scrubMeas", paste0(subject, ".rds"))
  if (file.exists(out_fname)) { next }

  cat(paste0(
    "Subject ", subject, ", ", as.character(acquisition), " ", 
    ifelse(test, "test", "retest"), " ", visit, 
    " (", ii, " of ", nrow(iters), ")", "\n"
  ))
  
  # Get input files.
  fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)
  data_dir <- file.path(subject, "MNINonLinear", "Results", fname_prefix)
  prefix_dir <- ifelse(test, dir_HCP_test, dir_HCP_retest)
  fnames <- list(
    CIFTI = file.path(prefix_dir, data_dir, paste0(fname_prefix, "_Atlas.dtseries.nii")),
    NIFTI = file.path(prefix_dir, data_dir, paste0(fname_prefix, ".nii.gz")),
    NIFTI_labs = file.path(dir_HCP_test, data_dir, "../../ROIs/Atlas_wmparc.2.nii.gz"),
    RP = file.path(prefix_dir, data_dir, "Movement_Regressors.txt"),
    RP_dt = file.path(prefix_dir, data_dir, "Movement_Regressors_dt.txt")
  )

  # # If retest, the data needs to be loaded.
  # if (!test) {
  #   data_zip_MPP <- file.path(
  #     dir_HCP_retest_archive, paste0(subject, "_3T_rfMRI_REST", visit, "_preproc.zip")
  #   )
    
  #   for (fname in fnames[c("CIFTI", "NIFTI", "RP", "RP_dt")]) {
  #     if (!file.exists(fname)) {
  #       cmd <- paste("unzip", data_zip_MPP, fname, "-d", dir_HCP_retest)
  #       system(cmd)
  #       stopifnot(file.exists(fname))
  #     }
  #   }
  # }

  # CompCor
  cat("CompCor.\n")
  cc <- fMRIscrub::CompCor_HCP(
    nii=fnames$NIFTI, 
    nii_labels=fnames$NIFTI_labs,
    ROI_noise=c("wm_cort", "csf", "wm_cblm"),
    idx=seq(nDrop+1, hcp_T), 
    noise_nPC=10, noise_erosion=c(2,1,2),
    verbose=TRUE
  )
  cc <- cbind(cc$PCs$wm_cort[,seq(5)], cc$PCs$csf[,seq(5)], cc$PCs$wm_cblm[,seq(5)])
  cc <- scale(cc)

  # Mean Signals and read in CIFTI
  cat("Reading NIFTI & processing mean signals.\n")
  nii <- RNifti::readNifti(fnames$NIFTI)
  niiLabs <- RNifti::readNifti(fnames$NIFTI_labs)
  ms <- fMRIscrub:::get_NIFTI_ROI_masks(fnames$NIFTI_labs, c("wm_cort", "csf", "wm_cblm"))
  ms <- list(wm=ms$wm_cort | ms$wm_cblm, csf=ms$csf, wholebrain=niiLabs>0)
  noise_erosion <- list(wm=2, csf=1, wholebrain=0)
  for (rr in seq(length(ms))) {
    ms[[rr]][,,] <- as.logical(ms[[rr]]) * 1
    ms[[rr]][,,] <- erode_vol(ms[[rr]], noise_erosion[[rr]], c(-1, 0, NA, NaN))
    ms[[rr]] <- apply(matrix(nii[ms[[rr]] > 0], ncol=1200), 2, mean)
  }
  rm(nii); gc()
  ms <- as.matrix(data.frame(ms)[,c("wm", "csf", "wholebrain")]) 

  cat("Reading CIFTI.\n")
  cii0 <- t(as.matrix(read_cifti(fnames$CIFTI, brainstructures="all")))
  cii0 <- cii0[seq(nDrop+1, hcp_T),]
  cii_mean <- mean(cii0)

  # FD
  cat("Motion Scrubbing.\n")
  rp <- as.matrix(read.table(fnames$RP)[,seq(6)])
  p6 <- scale(rp[seq(nDrop+1, hcp_T),])
  p36 <- scale(cbind(ms, rbind(0, diff(ms)), rp))
  p36 <- scale(cbind(p36, p36^2))[seq(nDrop+1, hcp_T),] 
  rp_nfc <- gsignal::filtfilt(nfiltc, rp)
  fd <- list(
    og = fMRIscrub::FD(rp),
    og_nfc_l4 = fMRIscrub::FD(rp_nfc, lag=4)
  )

  # Data-driven scrubbing
  cat("Data-driven scrubbing.\n")
  dd <- list(CC2MP6=list(pscrub=NULL, dv=NULL), `36P`=list(pscrub=NULL, dv=NULL))
  for (baseName in c("CC2MP6", "36P")) {
    cat(paste0("\t", baseName, "."))
    nreg <- switch(baseName,
      CC2 = cbind(1, dct4, cc[,c(1,2,6,7,11,12)]),
      CC5 = cbind(1, dct4, cc),
      MPP = as.matrix(rep(1, hcp_T - nDrop)),
      DCT4 = cbind(1, dct4),
      `9P` = cbind(1, dct4, p36[,c(seq(3), seq(7,12))]),
      `36P` = cbind(1, dct4, p36),
      CC2MP6 = cbind(1, dct4, cc[,c(1,2,6,7,11,12)], p6)
    )
    cat("\t\tNuisance regression.")
    cii <- nuisance_regression(cii0, nreg)

    # Projection scrubbing -----------------------------------------------------
    cat("\t\tProjection scrubbing.")
    dd[[baseName]]$pscrub <- fMRIscrub:::pscrub_multi(
      cii, projection=c("PCA", "PCA_kurt", "fusedPCA", "fusedPCA_kurt", "ICA", "ICA_kurt"),
      nuisance = NULL, get_dirs=FALSE, get_outliers=FALSE, verbose=TRUE
    )

    # DVARS --------------------------------------------------------------------
    cat("\t\tDVARS.")
    cii <- t((cii + cii_mean) / cii_mean * 1000)
    cii <- t(cii - apply(cii, 1, mean))
    dd[[baseName]]$dv <- DVARS(cii, normalize=FALSE)$measure
    rm(cii); gc()
  }

  saveRDS(list(fd=fd, dd=dd), out_fname)
  
  print(Sys.time() - time)
  time <- Sys.time()

  # if (!test) {
  #   unlink(fnames$CIFTI)
  #   unlink(fnames$NIFTI)
  #   unlink(fnames$RP)
  #   unlink(fnames$RP_dt)
  # }
}
