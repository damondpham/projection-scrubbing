## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("0_SharedCode.R")
stopifnot(SharedCode_version == c(11,0))

# Run this script twice, first without all subjects, then with.
# The latter is only used for the prediction analysis.
Use_S1200 <- TRUE
if (Use_S1200){
  subjects <- readRDS(file.path(dir_data_misc, "HCP_subjects.rds"))
  iters <- expand.grid(
    visit=seq(2), 
    acquisition=c("LR", "RL"), 
    subject=subjects,
    test=TRUE
  )
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

nfb <- c(.2, .5) * .72 * 2
nfc <- c(.31, .43) * .72 * 2
nfiltb <- gsignal::butter(n=10, w=nfb, type="stop")
nfiltc <- gsignal::cheby2(2, Rs=20, w=nfc, type="stop")

time <- Sys.time()

for (ii in seq(nrow(iters))) {
  # Get iteration info -------------------------------------------
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

  
  # Get output file; skip if done --------------------------------
  fd_fname <- file.path(
    dir_scrubMeas, "FD", 
    paste0(
      "FD_", subject, "_v", visit + (!test)*2, 
      "_", acquisition, ".rds"
    )
  )
  if (file.exists(fd_fname)) { next }
  
  # Get input files. ---------------------------------------------
  # RPs
  fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)
  data_dir <- switch(COMPUTER,
    RED = file.path(subject, "MNINonLinear", "Results", fname_prefix),
    MPro = file.path(subject, fname_prefix)
  )
  fnames <- list(
    RP = file.path(data_dir, "Movement_Regressors.txt"),
    RP_dt = file.path(data_dir, "Movement_Regressors_dt.txt")
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
    if (!file.exists(fnames$RP_dt)) { next }
  }
  
  rp <- as.matrix(read.table(fnames$RP)[,seq(6)])
  rp_nfc <- gsignal::filtfilt(nfiltc, rp)
  if (!Use_S1200) {
    rp_nfb <- gsignal::filtfilt(nfiltb, rp)
    rp_dt <- as.matrix(read.table(fnames$RP_dt)[,seq(6)])
    rp_dt_nfb <- gsignal::filtfilt(nfiltb, rp_dt)
    rp_dt_nfc <- gsignal::filtfilt(nfiltc, rp_dt) 
  }

  fd <- if (Use_S1200) {
    list(
      og = fMRIscrub::FD(rp),
      og_nfc_l4 = fMRIscrub::FD(rp_nfc, lag=4)
    )
  } else {
    list(
      og = fMRIscrub::FD(rp),
      og_nfb = fMRIscrub::FD(rp_nfb),
      og_nfc = fMRIscrub::FD(rp_nfc),
      og_l4 = fMRIscrub::FD(rp, lag=4),
      og_nfb_l4 = fMRIscrub::FD(rp_nfb, lag=4),
      og_nfc_l4 = fMRIscrub::FD(rp_nfc, lag=4),
      dt = fMRIscrub::FD(rp_dt),
      dt_nfb = fMRIscrub::FD(rp_dt_nfb),
      dt_nfc = fMRIscrub::FD(rp_dt_nfc),
      dt_l4 = fMRIscrub::FD(rp_dt, lag=4),
      dt_nfb_l4 = fMRIscrub::FD(rp_dt_nfb, lag=4),
      dt_nfc_l4 = fMRIscrub::FD(rp_dt_nfc, lag=4)
    ) 
  }
  
  fd <- lapply(fd, function(x){x$measure[seq(nDrop+1, hcp_T)]})

  # Save FD.
  saveRDS(fd, fd_fname)
  
  # Unload retest data.
  if (!test) {
    unlink(fnames$RP)
    unlink(fnames$RP_dt)
  }
  
  print(Sys.time() - time)
  time <- Sys.time()
}

