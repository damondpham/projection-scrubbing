## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("0_SharedCode.R")
stopifnot(SharedCode_version == c(11,0))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Parcellation -----------------------------------------------------------------
ParcMat <- readRDS(parc_fname)
cor_mask <- upper.tri(diag(parc_res2))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Edit this function to control skipping iterations.
skip_ii <- function(FCval_fname) {
  # Will skip this iteration if the FC file already exists.
  file.exists(FCval_fname)
}

baseNames <- c("CC2MP6")

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

flag_wrap1_val <- function(my_flag, my_cii){
  list(FC=get_FC(my_flag, my_cii), nScrub=c(sum(my_flag), length(my_flag)))
}

vols_from_nT <- function(nT){
  v0 <- ceiling((hcp_T-15)/2) - floor(nT/2)
  seq(v0, v0+nT-1)
}

flag_wrap2_val <- function(flag=NULL, tmasks="default") {
  nT <- hcp_T-15
  if (is.null(flag)) { 
    flag <- rep(FALSE, nT[length(nT)]) 
  } else { 
    flag <- as.logical(flag)
    stopifnot(length(flag) == nT[length(nT)])
  }

  if (tmasks=="default") {
    tmasks <- list(
      mid = vols_from_nT(417), # 5*60/.72
      out = setdiff(seq(nT), vols_from_nT(417))
    )
  } else if (tmasks=="modFD_strict") {
    tmasks <- list(
      full = seq(nT),
      mid = vols_from_nT(417), # 3*60/.72
      out = setdiff(seq(nT), vols_from_nT(417))
    )
  } else { stop() }

  flag <- lapply(tmasks, function(x){flag[x]})

  out <- vector("list", length=length(tmasks))
  names(out) <- names(tmasks)
  
  for (tt in seq(length(tmasks))) {
    vols_tt <- tmasks[[tt]]
    flag_tt <- flag[[tt]]
    
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
        cii[vols_tt,], 
        cbind(nreg[vols_tt,], spikes)
      )
    } else {
      spikes <- NULL
      cii2 <- cii * NA
    }
    
    out[[tt]] <- flag_wrap1_val(flag_tt, cii2)
  }
  out
}

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
time <- Sys.time()

for (baseName in baseNames) {
  cat(baseName, "~~~~~~~~~~~~~~~~~~~~~~\n")
  
  all_FCval_fname <- file.path(dir_FCval, "../3_FCval2_union", paste0("FC_", baseName, ".rds"))
  if (file.exists(all_FCval_fname)) { next }
  
  # Directories --------------------------------
  fd_dir <- file.path(dir_scrubMeas, "FD")
  lev_dir <- file.path(dir_scrubMeas, "Lev", baseName)
  dvars_dir <- file.path(dir_scrubMeas, "DVARS", baseName)

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

    FCval_fname <- file.path(dir_FCval, "../3_FCval2_union", baseName, paste0(suffix, ".rds"))
    #if (skip_ii(FCval_fname)) { next }

    # # Mean signals
    # ms <- data.frame(readRDS(file.path(dir_meanSignals, paste0(suffix, ".rds"))))
    # ms <- as.matrix(ms[,c("wm", "csf", "wholebrain")])

    # CompCor
    cc <- readRDS(file.path(dir_CompCor, paste0(suffix, ".rds")))
    cc <- scale(cbind(cc$PCs$wm_cort[,seq(5)], cc$PCs$csf[,seq(5)], cc$PCs$wm_cblm[,seq(5)]))
    
    # CIFTI
    fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)
    data_dir <- if (COMPUTER=="MPro") {
      file.path(subject, fname_prefix)
    } else {
      file.path(subject, "MNINonLinear", "Results", fname_prefix)
    }
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
    # p36 <- scale(cbind(ms, rbind(0, diff(ms)), rp))
    # p36 <- scale(cbind(p36, p36^2))[seq(nDrop+1, hcp_T),]

    # Read CIFTI. Drop first 15 frames. Parcellate.
    # (Is equivalent to parcellating after nuisance regression.)
    cat("\tReading CIFTI and parcellating.\n")
    cii <- as.matrix(read_xifti(file.path(read_dir, fnames$CIFTI), brainstructures="all"))[,seq(nDrop+1, hcp_T)]
    cii <- t(cii - rowMeans(cii)) %*% ParcMat

    nreg <- switch(baseName,
      # CC2 = cbind(1, dct4, cc[,c(1,2,6,7,11,12)]),
      # CC5 = cbind(1, dct4, cc),
      # MPP = as.matrix(rep(1, hcp_T - nDrop)),
      # DCT4 = cbind(1, dct4),
      # `9P` = cbind(1, dct4, p36[,c(seq(3), seq(7,12))]),
      # `36P` = cbind(1, dct4, p36),
      CC2MP6 = cbind(1, dct4, cc[,c(1,2,6,7,11,12)], p6)
    )
    
    # Scrubbing measures
    cat("\tReading scrubbing measures.\n")
    fd <- readRDS(file.path(fd_dir, paste0("FD_", suffix, ".rds")))
    if (baseName != "CC2MP6") { fd <- fd[c("og", "og_nfc_l4")] } # post-hoc, saves memory.
    lev <- readRDS(file.path(lev_dir, paste0("LEV_", suffix, ".rds")))
    dvars <- readRDS(file.path(dvars_dir, paste0("DVARS_", suffix, ".rds")))
  
    FCval_ii <- vector("list")
    
    # Get FC values ----------------------------------------------------------  
    # Nothing
    cat("\tBase")
    FCval_ii[["Base"]] <- flag_wrap2_val()

    flag <- list(
      modFD = fd[[FD_type]] > .5,
      proj = lev$measure[["ICA_kurt"]] > 3 * median(lev$measure[["ICA_kurt"]]),
      DVARS = (dvars$dv$DPD > 5) & (dvars$dv$ZD > qnorm(1-.05/(hcp_T-nDrop))) 
    )

    # for verification
    FCval_ii[["modFD"]] <- flag_wrap2_val(flag$modFD)
    FCval_ii[["DVARS"]] <- flag_wrap2_val(flag$proj)
    FCval_ii[["proj"]] <- flag_wrap2_val(flag$DVARS)

    # unions
    FCval_ii[["modFD_DVARS"]] <- flag_wrap2_val(flag$modFD | flag$DVARS)
    FCval_ii[["modFD_proj"]] <- flag_wrap2_val(flag$modFD | flag$proj)
    FCval_ii[["DVARS_proj"]] <- flag_wrap2_val(flag$DVARS | flag$proj)
    FCval_ii[["modFD_DVARS_proj"]] <- flag_wrap2_val(flag$modFD | flag$DVARS | flag$proj)
    
    # Unload retest data.
    if (!test) {
      unlink(file.path(dir_HCP_retest, fnames$CIFTI))
      unlink(file.path(dir_HCP_retest, fnames$RP))
    }
    
    print(Sys.time() - time)
    time <- Sys.time()
  }
}
