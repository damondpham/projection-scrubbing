## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("0_SharedCode.R")
stopifnot(SharedCode_version == c(11,0))

subjects <- readRDS(file.path(dir_data_misc, "HCP_S1200_LargeBalancedSample.rds"))
iters <- expand.grid(
  visit=seq(2), 
  acquisition=c("LR", "RL"), 
  subject=subjects,
  test=TRUE
)

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Parcellation -----------------------------------------------------------------
ParcMat <- readRDS(parc_fname)
cor_mask <- upper.tri(diag(parc_res2))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
baseNames <- list(
  `36P`=c("Base", "Lev___ICA", "DVARSdual", "FD___og___", "FD___og_nfc_l4___"),
  `CC2MP6`=c("Base", "Lev___ICA", "DVARSdual", "FD___og___", "FD___og_nfc_l4___")
)
ParcMat <- readRDS(parc_fname)
nT <- 1185

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

visits <- sapply(
  seq(4), function(x){paste0(
    "test_",
    ifelse(x>2, "RL", "LR"),
    ifelse(x%%2==1, 1, 2)
  )}
)

stopifnot(all(
  paste0(ifelse(iters$test, "test", "retest"), "_", iters$acquisition, iters$visit) == rep(visits, length(subjects))
))

n_FC <- (parc_res2^2 - parc_res2)/2

time <- Sys.time()

tname <- paste0("T_", nT)
for (bb in seq(length(baseNames))) {
  baseName <- names(baseNames)[bb]
  cat(baseName, "\n\n")
  
  for (ss in seq(length(baseNames[[bb]]))) {
    
    subsetName <- baseNames[[bb]][ss]
    cat(subsetName, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    FCagg_fname <- file.path(dir_FC, paste0("FC_", gsub("___$", "", baseName), "_", subsetName, "_", nT, "_S1200.rds"))
    if (file.exists(FCagg_fname)) { next }
    print(basename(FCagg_fname))
    
    for (ii in seq(nrow(iters))) {
      # Get iters info
      subject <- iters[ii, "subject"]
      acquisition <- as.character(iters[ii, "acquisition"])
      test <- iters[ii, "test"]
      visit <- iters[ii, "visit"]
      
      # Get FC values
      FC_fname <- file.path(
        dir_FC, baseName,
        paste0(subject, "_v", visit + (!test)*2, "_", acquisition, ".rds")
      )
      cat("\t", basename(FC_fname), "\n")
      if (!file.exists(FC_fname)) { stop("Missing file: ", FC_fname) }
      FC <- readRDS(FC_fname)
      
      # Initialize big FC array, if first iteration
      if (ii == 1) {
        methods <- names(FC$FC)[grepl(subsetName, names(FC$FC))]
        n_methods <- length(methods)
        FCagg <- array(NA, c(nrow(iters), n_FC, n_methods))
      }
      FCagg[ii,,] <- do.call(cbind, lapply(FC$FC[methods], function(x){x[[tname]]}))
    }
    
    dim(FCagg) <- c(length(visits), length(subjects), n_FC, n_methods)
    dimnames(FCagg) <- list(visits, subjects, NULL, methods)
    
    cat("Saving...\n")
    saveRDS(FCagg, FCagg_fname)
  }
}
