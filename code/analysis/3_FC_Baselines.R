## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("0_SharedCode.R")
stopifnot(SharedCode_version == c(11,0))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Parcellation -----------------------------------------------------------------
ParcMat <- readRDS(parc_fname)
cor_mask <- upper.tri(diag(parc_res2))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Edit this function to control skipping iterations.
skip_ii <- function(FC_fname) {
  # Will skip this iteration if the FC file already exists.
  file.exists(FC_fname)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nreg_base <- c(
  "None__DCT0",
  "None__DCT4",
  "CompCor__PC01_withCblm_DCT0",
  "CompCor__PC01_withCblm_DCT4",
  "CompCor__PC02_withCblm_DCT0",
  "CompCor__PC02_withCblm_DCT4",
  "CompCor__PC05_withCblm_DCT0",
  "CompCor__PC05_withCblm_DCT4",
  "CompCor__PC10_withCblm_DCT0",
  "CompCor__PC10_withCblm_DCT4",
  "CompCor__PC01_woutCblm_DCT0",
  "CompCor__PC01_woutCblm_DCT4",
  "CompCor__PC02_woutCblm_DCT0",
  "CompCor__PC02_woutCblm_DCT4",
  "CompCor__PC05_woutCblm_DCT0",
  "CompCor__PC05_woutCblm_DCT4",
  "CompCor__PC10_woutCblm_DCT0",
  "CompCor__PC10_woutCblm_DCT4",
  "HCP__ICA-FIX",
  "ICA-FIX__seq_24P_soft_DCT0", # replicate HCP
  "ICA-FIX__seq_24P_soft_DCT4",
  "ICA-FIX__sim_24P_soft_DCT0",
  "ICA-FIX__sim_24P_soft_DCT4",
  "ICA-FIX__seq_00P_soft_DCT4",
  "ICA-FIX__xxx_00P_soft_DCT0",
  "ICA-FIX__sim_00P_soft_DCT4",
  "ICA-FIX__xxx_24P_hard_DCT0",
  "ICA-FIX__xxx_24P_hard_DCT4",
  "ICA-FIX__xxx_00P_hard_DCT0",
  "ICA-FIX__xxx_00P_hard_DCT4",
  "9P__grayo_DCT0",
  "9P__whole_DCT0",
  "9P__grayo_DCT4",
  "9P__whole_DCT4",
  "36P__grayo_DCT0",
  "36P__whole_DCT0",
  "36P__grayo_DCT4",
  "36P__whole_DCT4",
  "CSF-WM__DCT0",
  "CSF-WM__DCT4",
  "GSR__grayo_DCT0",
  "GSR__whole_DCT0",
  "GSR__grayo_DCT4",
  "GSR__whole_DCT4",
  "CC2MP__06P_DCT0",
  "CC2MP__06P_DCT4",
  "CC2MP__24P_DCT0",
  "CC2MP__24P_DCT4"
)

baseNames <- baseName <- "Baselines"

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
get_FC <- function(this_flag, this_cii) {
  this_T <- length(this_flag)
  nScrub <- sum(this_flag)

  if (nScrub < this_T - 2) {
    return(cor(this_cii[!this_flag,])[cor_mask])
  } else {
    # Too many frames were scrubbed to calculate FC
    return(matrix(NA, nrow=parc_res2, ncol=parc_res2)[cor_mask])
  }
}

get_randFC <- function(that_flag, that_cii, nrand=10){
  out <- vector("list", length=nrand)
  for (rr in seq(length(nrand))) {
    out[[rr]] <- get_FC(sample(that_flag), that_cii)
  }
  out <- abind::abind(out, along=2)
  out <- rowMeans(out)
}

flag_wrap1 <- function(my_flag, my_cii){
  list(FC=get_FC(my_flag, my_cii), randFC=get_randFC(my_flag, my_cii), nScrub=sum(my_flag))
}

vols_from_nT <- function(nT){
  v0 <- ceiling((hcp_T-15)/2) - floor(nT/2)
  seq(v0, v0+nT-1)
}

flag_wrap2 <- function(flag=NULL, nT=c(400, 800, 1185), seed=NULL) {
  if (is.null(flag)) { 
    flag <- rep(FALSE, nT[length(nT)]) 
  } else { 
    flag <- as.logical(flag)
    stopifnot(length(flag) == nT[length(nT)])
  }
  out <- vector("list", length=length(nT))
  names(out) <- paste0("T_", nT)
  
  for (tt in seq(nT)) {
    vols_tt <- vols_from_nT(nT[tt])
    flag_tt <- flag[vols_tt]
    
    T_ <- length(flag_tt)
    nScrub <- sum(flag_tt)
    
    if (nScrub < T_ - 2) {
      if (nScrub > 0) {
        # One-hot encode outlier flags
        spikes <- matrix(0, nrow=T_, ncol=nScrub)
        spikes[seq(0, nScrub-1)*T_ + which(flag_tt)] <- 1
      } else {
        # No scrubbing
        spikes <- NULL
      }
      cii2 <- nuisance_regression(
        cii_bb[vols_tt,], 
        cbind(nreg[vols_tt,], spikes)
      )
    } else {
      spikes <- NULL
      cii2 <- cii_bb * NA
    }
    
    if (!is.null(seed)) {
      out[[tt]] <- with(set.seed(seed), flag_wrap1(flag_tt, cii2))
    } else {
      out[[tt]] <- flag_wrap1(flag_tt, cii2)
    }
  }
  out
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
time <- Sys.time()

for (ii in seq(nrow(iters))) {
  
  subject <- iters[ii, "subject"]
  acquisition <- as.character(iters[ii, "acquisition"])
  test <- iters[ii, "test"]
  visit <- iters[ii, "visit"]

  cat(paste0(
    "Subject ", subject, ", ", acquisition, " ", 
    ifelse(test, "test", "retest"), " ", visit, 
    " (", ii, " of ", nrow(iters), ")", "\n"
  ))

  # Get files; skip if done --------------------------------
  suffix <- paste0(subject , "_v", visit + (!test)*2, "_", acquisition)
  FC_fname <- file.path(dir_FC, baseName, paste0(suffix, ".rds"))
  if (skip_ii(FC_fname)) { next }  
  
  # -------------------------------------------------------------------------------------------------------------------
  cat("Reading files.\n")
  # Get input files.
  fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)
  data_dir <- file.path(subject, "MNINonLinear", "Results", fname_prefix)
  fnames <- list(
    RP = file.path(data_dir, "Movement_Regressors.txt"),
    CIFTI = file.path(data_dir, paste0(fname_prefix, "_Atlas.dtseries.nii")),
    CIFTI_FIX = file.path(data_dir, paste0(fname_prefix, "_Atlas_hp2000_clean.dtseries.nii")),
    FIX_IC = file.path(data_dir, paste0(fname_prefix, "_hp2000.ica"), "filtered_func_data.ica/melodic_mix"),
    FIX_labs = file.path(data_dir, paste0(fname_prefix, "_hp2000.ica"), "Noise.txt")
  )
  
  precomputed_fname <- paste0(subject, "_v", visit + (!test)*2, "_", acquisition, ".rds")
  cc <- readRDS(file.path(dir_CompCor, precomputed_fname))
  cc$PCs <- lapply(cc$PCs, scale)
  ms <- scale(as.data.frame(readRDS(file.path(dir_meanSignals, precomputed_fname)))[seq(nDrop+1, hcp_T),])
  
  # If retest, the data needs to be loaded.
  if (!test) {
    data_zip_MPP <- file.path(
      dir_HCP_retest_archive, paste0(subject, "_3T_rfMRI_REST", visit, "_preproc.zip")
    )
    data_zip_FIX <- file.path(
      dir_HCP_retest_archive, paste0(subject, "_3T_rfMRI_REST", visit, "_fixextended.zip")
    )
    for (fname in fnames[c("RP", "CIFTI")]) {
      if (!file.exists(file.path(dir_HCP_retest, fname))) {
        cmd <- paste("unzip", data_zip_MPP, fname, "-d", dir_HCP_retest)
        system(cmd)
        stopifnot(file.exists(file.path(dir_HCP_retest, fname)))
      }
    }
    for (fname in fnames[c("CIFTI_FIX", "FIX_IC", "FIX_labs")]) {
      if (!file.exists(file.path(dir_HCP_retest, fname))) {
        cmd <- paste("unzip", data_zip_FIX, fname, "-d", dir_HCP_retest)
        system(cmd)
        stopifnot(file.exists(file.path(dir_HCP_retest, fname)))
      }
    }
  }
  
  read_dir <- ifelse(test, dir_HCP_test, dir_HCP_retest)
  
  nBase <- length(nreg_base)
  FC_ii <- setNames(vector("list", length=nBase), nreg_base)
  
  # -------------------------------------------------------------------------------------------------------------------
  cat("HCP ICA-FIX\n")
  cii_bb <- as.matrix(read_xifti(file.path(read_dir, fnames$CIFTI_FIX), brainstructures="all"))[,seq(nDrop+1, hcp_T)]
  cii_bb <- t(cii_bb - rowMeans(cii_bb)) %*% ParcMat
  nreg <- as.matrix(rep(sqrt(1/(hcp_T-nDrop)), hcp_T-nDrop))
  FC_ii[["HCP__ICA-FIX"]] <- lapply(flag_wrap2(seed=ii), `[[`, "FC")
  
  # -------------------------------------------------------------------------------------------------------------------
  # CIFTI
  cii <- as.matrix(read_xifti(file.path(read_dir, fnames$CIFTI), brainstructures="all"))
  cii <- cii - rowMeans(cii)
  cat("Highpassing the CIFTI data (to use only for some FCs).\n")
  ciiHP <- t(fMRIscrub::fsl_bptf(t(cii), 0.5*2000/.72))
  cii <- cii[,seq(nDrop+1, hcp_T)]
  ciiHP <- ciiHP[,seq(nDrop+1, hcp_T)]
  cii <- t(cii - rowMeans(cii))
  ciiHP <- t(ciiHP - rowMeans(ciiHP))

  # FIX ICA timecourses
  ic <- scale(as.matrix(read.table(file.path(read_dir, fnames$FIX_IC)))[seq(nDrop+1, hcp_T),])
  icid <- as.numeric(read.table(file.path(read_dir, fnames$FIX_labs)))
  
  # Load RPs. Compute their highpass and first-order differences
  rp <- as.matrix(read.table(file.path(read_dir, fnames$RP)))
  rp <- scale(cbind(rp, rp^2))
  rpHP <- scale(fMRIscrub:::fsl_bptf(rp, 0.5*2000/.72))
  rp <- scale(rp[seq(nDrop+1, hcp_T),])
  rpHP <- scale(rpHP[seq(nDrop+1, hcp_T),])
  
  for (bb in seq(length(nreg_base))) {
    base <- nreg_base[bb]
    if (base == "HCP__ICA-FIX") { next }
    cat("\t", base, "\n")
    base_type <- gsub("__.*", "", base)
    
    if (base_type == "ICA-FIX") {
      cii_bb <- ciiHP
      rp_bb <- rpHP
    } else {
      cii_bb <- cii
      rp_bb <- rp
    }
    
    use <- function(x){grepl(x, base)}
    
    # Intercept
    nreg <- as.matrix(rep(sqrt(1/(hcp_T-nDrop)), hcp_T-nDrop))
    if (use("DCT4")) { nreg <- cbind(nreg, dct4) }
    if (use("24P|36P")) { nreg <- cbind(nreg, rp_bb) }

    ciibbp <- cii_bb %*% ParcMat
    
    # Different types
    if (base_type == "ICA-FIX") {
      # ICA-FIX
      if (use("soft")) {
        if (use("seq")) {
          ic_bb <- nuisance_regression(ic, nreg)
          cii_bb <- nuisance_regression(cii_bb, nreg) 
          beta <- pracma::pinv(ic_bb) %*% cii_bb
          cii_bb <- cii_bb - ic_bb[,icid,drop=FALSE] %*% beta[icid,,drop=FALSE] # TxP
        } else {
          rmid <- c(seq(ncol(nreg)), icid+ncol(nreg))
          nreg <- cbind(nreg, ic)
          beta <- pracma::pinv(nreg) %*% cii_bb
          cii_bb <- cii_bb - nreg[,rmid,drop=FALSE] %*% beta[rmid,,drop=FALSE] # TxP 
        }
      } else {
        nreg <- cbind(nreg, ic[,icid])
        cii_bb <- nuisance_regression(cii_bb, nreg)
      }
      cii_bb <- cii_bb %*% ParcMat
      nreg <- as.matrix(rep(sqrt(1/(hcp_T-nDrop)), hcp_T-nDrop))
      FC_ii[[base]] <- lapply(flag_wrap2(seed=ii), `[[`, "FC")
    } else {
      cii_bb <- ciibbp
      if (base_type == "None") {
        NULL
      } else if (base_type == "CompCor") {
        nPC_bb <- as.numeric(substr(base, 12, 13))
        stopifnot(nPC_bb %in% c(1,2,5,10))
        nreg <- cbind(nreg, cc$PCs$wm_cort[,seq(nPC_bb)], cc$PCs$csf[,seq(nPC_bb)])
        if (use("withCblm")) { nreg <- cbind(nreg, cc$PCs$wm_cblm[,seq(nPC_bb)]) }
      } else if (base_type %in% c("9P", "36P")) {
        ms_bb <- ms[,c("wm", "csf", ifelse(use("whole"), "wholebrain", "cii"))]
        if (base_type == "36P") {
          # Derivative (first difference)
          ms_bb <- cbind(ms_bb, rbind(0, diff(ms_bb)))
          # Square
          ms_bb <- cbind(ms_bb, ms_bb^2)
        }
        nreg <- cbind(nreg, ms_bb)
        if (base_type == "36P") {
          nreg <- cbind(nreg, rp)
        } else {
          nreg <- cbind(nreg, rp[,seq(6)])
        }
      } else if (base_type == "CSF-WM") {
        nreg <- cbind(nreg, ms[,c("wm", "csf")])
      } else if (base_type == "GSR") {
        nreg <- cbind(nreg, ms[,ifelse(use("whole"), "wholebrain", "cii")])
      } else if (base_type == "CC2MP") {
        nPC_bb <- 2
        nreg <- cbind(nreg, cc$PCs$wm_cort[,seq(nPC_bb)], cc$PCs$csf[,seq(nPC_bb)], cc$PCs$wm_cblm[,seq(nPC_bb)])
        if (use("24P")) {
          nreg <- cbind(nreg, rp)
        } else {
          nreg <- cbind(nreg, rp[,seq(6)])
        }
      } else {
        stop()
      }
      FC_ii[[base]] <- lapply(flag_wrap2(seed=ii), `[[`, "FC")
    }
  }
  
  saveRDS(list(FC=FC_ii, randFC=NULL, nScrub=NULL), FC_fname)
    
  # Unload retest data.
  if (!test) {
    cat("Removing data.\n")
    unlink(file.path(dir_HCP_retest, fnames$RP))
    unlink(file.path(dir_HCP_retest, fnames$CIFTI))
    unlink(file.path(dir_HCP_retest, fnames$CIFTI_FIX))
    unlink(file.path(dir_HCP_retest, fnames$FIX_IC))
    unlink(file.path(dir_HCP_retest, fnames$FIX_labs))
  }
  
  print(Sys.time() - time)
  time <- Sys.time()
}
