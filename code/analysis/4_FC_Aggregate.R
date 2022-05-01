## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("0_SharedCode.R")
stopifnot(SharedCode_version == c(11,0))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
baseNames <- list(
  Baselines=c("None", "CompCor", "HCP__ICA-FIX", "ICA-FIX", "9P", "36P", "CSF-WM", "GSR", "CC2MP"),
  CC2MP6=c("Base", "Lev___PCA", "Lev___fusedPCA", "Lev___ICA", "DVARSdual"),
  CC2MP6=c("FD___og___", "FD___og_nfb___", "FD___og_nfc___", "FD___og_l4___", "FD___og_nfb_l4___", "FD___og_nfc_l4___"),
  CC2MP6=c("FD___dt___", "FD___dt_nfb___", "FD___dt_nfc___", "FD___dt_l4___", "FD___dt_nfb_l4___", "FD___dt_nfc_l4___"),
  DCT4=c("Base", "Lev___PCA", "Lev___fusedPCA", "Lev___ICA", "DVARSdual", "FD___og_nfc_l4", "FD___og___"),
  CC5=c("Base", "Lev___PCA", "Lev___fusedPCA", "Lev___ICA", "DVARSdual", "FD___og_nfc_l4", "FD___og___"),
  `9P`=c("Base", "Lev___PCA", "Lev___fusedPCA", "Lev___ICA", "DVARSdual", "FD___og_nfc_l4", "FD___og___"),
  `36P`=c("Base", "Lev___PCA", "Lev___fusedPCA", "Lev___ICA", "DVARSdual", "FD___og_nfc_l4", "FD___og___")
)
ParcMat <- readRDS(parc_fname)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

visits <- sapply(
  seq(8), function(x){paste0(
    ifelse((x+1)%%4 >1, "test", "retest"), "_",
    ifelse(x>4, "RL", "LR"),
    ifelse(x%%2==1, 1, 2)
  )}
)

stopifnot(all(
  paste0(ifelse(iters$test, "test", "retest"), "_", iters$acquisition, iters$visit) == rep(visits, length(subjects))
))

n_FC <- (parc_res2^2 - parc_res2)/2

time <- Sys.time()

for (nT in c(400, 800, 1185)) {
  tname <- paste0("T_", nT)
  for (bb in seq(length(baseNames))) {
    baseName <- names(baseNames)[bb]
    cat(baseName, "\n\n")
    
    #if ((baseName %in% c("DCT4", "CC5", "9P", "36P")) && (nT != 1185)) { next }
    useNScrub <- baseName != "Baselines"

    for (ss in seq(length(baseNames[[bb]]))) {
      
      subsetName <- baseNames[[bb]][ss]
      cat(subsetName, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
      FCagg_fname <- file.path(dir_FC, paste0("FC_", gsub("___$", "", baseName), "_", subsetName, "_", nT, ".rds"))
      if (file.exists(FCagg_fname)) { next }
      print(basename(FCagg_fname))
      nScrub_fname <- file.path(dir_analysis, "nScrub", paste0(baseName, "_", subsetName, "_", nT, ".rds"))
      
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
        if (useNScrub) { ns <- readRDS(FC_fname)$nScrub }
        
        # Initialize big FC array, if first iteration
        if (ii == 1) {
          methods <- names(FC$FC)[grepl(subsetName, names(FC$FC))]
          n_methods <- length(methods)
          FCagg <- array(NA, c(nrow(iters), n_FC, n_methods))
          if (useNScrub) { nScrub <- array(NA, c(nrow(iters), n_methods)) }
          useRand <- ("randFC" %in% names(FC)) && !is.null(FC$randFC)
          if (useRand) { randFCagg <- FCagg }
        }
        FCagg[ii,,] <- do.call(cbind, lapply(FC$FC[methods], function(x){x[[tname]]}))
        if (useNScrub) {
          nScrub[ii,] <- do.call(c, lapply(ns[methods], function(y){y[[tname]]}))
        }
        if (useRand) {
          randFCagg[ii,,] <- do.call(cbind, lapply(FC$randFC[methods], function(x){x[[tname]]}))
        }
      }
      
      dim(FCagg) <- c(length(visits), length(subjects), n_FC, n_methods)
      dimnames(FCagg) <- list(visits, subjects, NULL, methods)

      dim(nScrub) <- c(length(visits), length(subjects), n_methods)
      dimnames(nScrub) <- list(visits, subjects, methods)

      if (useRand) {
        dim(randFCagg) <- c(length(visits), length(subjects), n_FC, n_methods)
        dimnames(randFCagg) <- list(visits, subjects, NULL, methods)
      }
      
      cat("Saving...\n")
      saveRDS(FCagg, FCagg_fname)
      if (useNScrub) {
        saveRDS(nScrub, nScrub_fname)
      }
      if (useRand) {
        saveRDS(randFCagg, gsub(".rds", "_random.rds", FCagg_fname, fixed=TRUE))
      }
    }
  }
}
