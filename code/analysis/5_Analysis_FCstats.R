## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("0_SharedCode.R")
source("0_Analysis_helper.R")
stopifnot(SharedCode_version == c(11,0))

library(abind)
library(tidyr)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
FC_fnames <- list.files("../../analysis-results/3_FC", "rds")

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
subjects <- sort(subjects)
nSubjects <- length(subjects)
visits <- sapply(
  seq(8), function(x){paste0(
    ifelse((x+1)%%4 >1, "test", "retest"), "_",
    ifelse(x>4, "RL", "LR"),
    ifelse(x%%2==1, 1, 2)
  )}
)
nVisits <- length(visits)

n_FC <- (parc_res2^2 - parc_res2)/2

dg <- read.csv(hcp_dg_fname)
dg <- subset(dg, Subject %in% subjects)
dg$Gender <- factor(dg$Gender)

check_FC <- function(FC) {
  stopifnot(!any(is.na(FC)))
  #n_parcels <- 119#419
  #n_parcpairs <- (n_parcels * n_parcels - n_parcels) / 2 #7021
  stopifnot(dim(FC)[seq(2)] == c(nVisits, nSubjects)) #n_parcpairs
  stopifnot(all(dimnames(FC)[[1]] == visits))
  stopifnot(all(dimnames(FC)[[2]] == as.character(dg$Subject)))
  stopifnot(all(dg$Subject == as.numeric(gsub("subject_", "", dimnames(FC)[[3]]))))
}

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load variables used for fingerprinting.
q <- readRDS(file.path(dir_Parc, "ParcLabels.rds"))
PLabs <- q$PLabs
SLabs <- q$SLabs
Labs <- q$Labs
rm(q)

# Organize parcellation pairs by networks/subcortex, and off/on-diagonal -------
labs <- unique(as.character(Labs$network))
nlabs <- length(labs)
labs <- setNames(vector("list", nlabs+1), c(labs, "All"))
for (ll in seq(nlabs + 1)) {
  if (ll == nlabs + 1) {
    labs[[ll]] <- list( on = rep(TRUE, n_FC), off = rep(TRUE, n_FC) ); next
  }
  q <- outer(Labs$network==names(labs)[ll], Labs$network==names(labs)[ll], `+`)
  q <- q[upper.tri(q)]
  labs[[ll]] <- list( on = q > 1, off = q > 0 )
}

# Load mean FD, for QC-FC
meanFD <- readRDS(file.path(dir_data_misc, "HCP_meanFD.rds"))

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Filter FC files to process here.
FC_fnames <- FC_fnames[!grepl("random|Baselines|S1200|_dt_", FC_fnames)]
#FC_fnames <- FC_fnames[grepl("CC2MP6", FC_fnames)]
FC_fnames <- FC_fnames[grepl("400", FC_fnames)]

for (ff in seq(length(FC_fnames))) {
  FC_fname <- file.path(dir_FC, FC_fnames[ff])
  cat(basename(FC_fname), "\n")
  # Load FC values.
  FC <- readRDS(FC_fname)
  FC <- FC[,order(dimnames(FC)[[2]]),,]
  
  # Format and check.
  if (length(dim(FC)) == 3) {
    q <- dimnames(FC)
    dim(FC) <- c(dim(FC), 1)
    dimnames(FC)[seq(3)] <- q
    dimnames(FC)[[4]] <- gsub("^FC_|\\..*", "", basename(FC_fname))
  }
  # Skip these baselines: CompCor w/o cerebellum, and ICA-FIX w/ hard regression.
  #   for the purpose of reducing memory requirements.
  my_methods <- dimnames(FC)[[4]][!grepl("woutCblm|hard", dimnames(FC)[[4]])]
  FC <- FC[,,,my_methods,drop=FALSE]
  FC <- FC[,,,!apply(is.na(FC), 4, any),drop=FALSE]
  my_methods <- dimnames(FC)[[4]]
  check_FC(FC)
  
  # Mean FC.
  meanFC_fname <- file.path(dir_analysis, "meanFC", gsub("^FC_", "", basename(FC_fname)))
  if (!file.exists(meanFC_fname)) {
    cat("\tMean FC")
    saveRDS(colMeans(FC, dims=2), meanFC_fname)
  }
  
  # z-transform.
  FCz <- psych::fisherz(FC)
  rm(FC) # Not used.
  gc()
  
  # Mean z-transformed FC.
  meanFCz_fname <- file.path(dir_analysis, "meanFCz", gsub("^FC_", "", basename(FC_fname)))
  if (!file.exists(meanFCz_fname)) {
    cat("\tMean z-transformed FC")
    saveRDS(colMeans(FCz, dims=2), meanFCz_fname)
  }
  
  # Reliability
  icc_fname <- file.path(dir_analysis, "reliability", gsub("^FC_", "", basename(FC_fname)))
  if (!file.exists(icc_fname)) {
    cat("\tReliability")
    
    # Vectorize across FC pairs & methods.
    FCz2 <- FCz
    d <- dim(FCz2)
    dim(FCz2) <- c(d[1], d[2], d[3]*d[4])
    # Compute ICC (and MSB, and MSR).
    icc <- get_ICC31_decomp(FCz2)
    # Re-format and unvectorize.
    icc2 <- vector("list", length=ncol(icc))
    names(icc2) <- colnames(icc)
    for (jj in seq(ncol(icc))) {
      z <- matrix(icc[,jj], nrow=d[3], ncol=d[4])
      colnames(z) <- my_methods
      icc2[[jj]] <- z
    }
    # Save.
    saveRDS(icc2, icc_fname)
    rm(FCz2); rm(icc)
  }
  
  # Do the other analyses only for some files.
  more_analysis <- TRUE
  more_analysis <- more_analysis && (!grepl("random", FC_fname))
  more_analysis <- more_analysis && (grepl("1185", FC_fname))
  if (more_analysis) {
    # Fingerprinting
    fp_fname <- file.path(dir_analysis, "fingerprinting", gsub("^FC_", "", basename(FC_fname)))
    if (!file.exists(fp_fname)) {
      cat("\tFingerprinting")
      fp <- setNames(vector("list", dim(FCz)[4]), dimnames(FCz)[[4]])
      for (mm in seq(dim(FCz)[4])) {
        fp[[mm]] <- fingerprint(FCz[,,,mm])
      }
      saveRDS(fp, fp_fname)
    }
    
    # QC-FC
    qcfc_fname <- file.path(dir_analysis, "qcfc", gsub("^FC_", "", basename(FC_fname)))
    if (!file.exists(qcfc_fname)) {
      cat("\tQC-FC")
      q <- qc_fc(FCz)
      saveRDS(q, qcfc_fname)
    }
    
    # MAC-rsFC
    # "All RSFC correlations were Fisher's r-to-Z transformed immediately after
    #   their calculation, prior to their use in any further computations or
    #   analyses, and are reported in this form throughout this manuscript. In
    #   all cases, (Z-transformed) RSFC correlations for each subject were
    #   calculated separately over runs and averaged together."
    mac_fname <- file.path(dir_analysis, "mac", gsub("^FC_", "", basename(FC_fname)))
    if (!file.exists(mac_fname)) {
      cat("\tMAC")
      FCzR <- readRDS(file.path(
        dir_FC, gsub(".rds", "_random.rds", basename(FC_fname), fixed=TRUE)
      ))
      FCzR <- psych::fisherz(FCzR[,order(dimnames(FCzR)[[2]]),,,drop=FALSE])
      saveRDS(colMeans(abs(colMeans(FCz-FCzR)), dims=2), mac_fname)
      rm(FCzR)
    }
    
    # DMN-Motion
    cat("\n")
  }
  
  next # TEMP
  
  
  
  cat("\n")
}
