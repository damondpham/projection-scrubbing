## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
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

baseNames <- c("CC2MP6", "CC2", "CC5", "DCT4", "9P", "36P")

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
  out <- matrix(NA, nrow=sum(cor_mask), ncol=nrand)
  for (rr in seq(nrand)) {
    out[,rr] <- get_FC(sample(that_flag), that_cii)
  }
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
        cii[vols_tt,], 
        cbind(nreg[vols_tt,], spikes)
      )
    } else {
      spikes <- NULL
      cii2 <- cii * NA
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


for (baseName in baseNames) {
  cat(baseName, "~~~~~~~~~~~~~~~~~~~~~~\n")
  
  all_FC_fname <- file.path(dir_FC, paste0("FC_", baseName, ".rds"))
  all_randFC_fname <- file.path(dir_FC, paste0("FC_", baseName, "_random.rds"))
  if (file.exists(all_FC_fname)) { next }
  
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

    FC_fname <- file.path(dir_FC, baseName, paste0(suffix, ".rds"))
    if (skip_ii(FC_fname)) { next }

    # Mean signals
    ms <- data.frame(readRDS(file.path(dir_meanSignals, paste0(suffix, ".rds"))))
    ms <- as.matrix(ms[,c("wm", "csf", "wholebrain")])

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
    p36 <- scale(cbind(ms, rbind(0, diff(ms)), rp))
    p36 <- scale(cbind(p36, p36^2))[seq(nDrop+1, hcp_T),]

    # Read CIFTI. Drop first 15 frames. Parcellate.
    # (Is equivalent to parcellating after nuisance regression.)
    cii <- as.matrix(read_xifti(file.path(read_dir, fnames$CIFTI), brainstructures="all"))[,seq(nDrop+1, hcp_T)]
    cii <- t(cii - rowMeans(cii)) %*% ParcMat

    nreg <- switch(baseName,
      CC2 = cbind(1, dct4, cc[,c(1,2,6,7,11,12)]),
      CC5 = cbind(1, dct4, cc),
      MPP = as.matrix(rep(1, hcp_T - nDrop)),
      DCT4 = cbind(1, dct4),
      `9P` = cbind(1, dct4, p36[,c(seq(3), seq(7,12))]),
      `36P` = cbind(1, dct4, p36),
      CC2MP6 = cbind(1, dct4, cc[,c(1,2,6,7,11,12)], p6)
    )
    
    # Scrubbing measures
    fd <- readRDS(file.path(fd_dir, paste0("FD_", suffix, ".rds")))
    if (baseName != "CC2MP6") { fd <- fd[c("og", "og_nfc_l4")] } # post-hoc, saves memory.
    lev <- readRDS(file.path(lev_dir, paste0("LEV_", suffix, ".rds")))
    dvars <- readRDS(file.path(dvars_dir, paste0("DVARS_", suffix, ".rds")))
  
    FC_ii <- randFC_ii <- nScrub_ii <- vector("list")
    
    # Get FC values ----------------------------------------------------------  
    # Nothing
    cat("\tBase")
    x <- flag_wrap2(seed=ii)
    FC_ii[["Base"]] <- lapply(x, `[[`, "FC")
    randFC_ii[["Base"]] <- lapply(x, `[[`, "randFC")
    nScrub_ii[["Base"]] <- lapply(x, `[[`, "nScrub")

    # FD
    cat("\n\tFD")
    FD_cuts <- seq(.1, .8, .1)
    for (FD_cut in FD_cuts) {
      for (FD_type in names(fd)) {
        this_name <- paste("FD", FD_type, FD_cut, sep="___")
        flag <- fd[[FD_type]] > FD_cut
        x <- flag_wrap2(flag, seed=ii)
        FC_ii[[this_name]] <- lapply(x, `[[`, "FC")
        randFC_ii[[this_name]] <- lapply(x, `[[`, "randFC")
        nScrub_ii[[this_name]] <- lapply(x, `[[`, "nScrub")
      }
    }
    
    # Median leverage
    cat("\n\tLev")
    medlev_cuts <- seq(1, 8)
    for (proj in colnames(lev$measure)) {
      for (medlev_cut in medlev_cuts) {
        this_name <- paste("Lev", proj, medlev_cut, sep="___")
        flag <- lev$measure[[proj]] > medlev_cut * median(lev$measure[[proj]])
        x <- flag_wrap2(flag, seed=ii)
        FC_ii[[this_name]] <- lapply(x, `[[`, "FC")
        randFC_ii[[this_name]] <- lapply(x, `[[`, "randFC")
        nScrub_ii[[this_name]] <- lapply(x, `[[`, "nScrub")
      } 
    }
    
    # DVARS dual (DPD & ZD)
    cat("\n\tDVARS dual")
    flag <- (dvars$dv$DPD > 5) & (dvars$dv$ZD > qnorm(1-.05/(hcp_T-nDrop))) 
    x <- flag_wrap2(flag, seed=ii)
    FC_ii[["DVARSdual"]] <- lapply(x, `[[`, "FC")
    randFC_ii[["DVARSdual"]] <- lapply(x, `[[`, "randFC")
    nScrub_ii[["DVARSdual"]] <- lapply(x, `[[`, "nScrub")      
    cat("\n")
    
    saveRDS(list(FC=FC_ii, randFC=randFC_ii, nScrub=nScrub_ii), FC_fname)
    
    # Unload retest data.
    if (!test) {
      unlink(file.path(dir_HCP_retest, fnames$CIFTI))
      unlink(file.path(dir_HCP_retest, fnames$RP))
    }
    
    print(Sys.time() - time)
    time <- Sys.time()
  }
}
