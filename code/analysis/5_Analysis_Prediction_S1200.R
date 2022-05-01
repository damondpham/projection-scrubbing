## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("0_SharedCode.R")
source("0_Analysis_helper.R")
stopifnot(SharedCode_version == c(11,0))

subjects <- readRDS(file.path(dir_data_misc, "HCP_S1200_LargeBalancedSample.rds"))
iters <- expand.grid(
  visit=seq(2), 
  acquisition=c("LR", "RL"), 
  subject=subjects,
  test=TRUE
)

library(abind)
library(tidyr)
library(glmnet)

library(parallel)
library(doParallel)

n_cores <- 7
stopifnot(n_cores < parallel::detectCores() - 1)
print(n_cores)
my_cluster <- parallel::makeCluster(
  n_cores, type="PSOCK", outfile=""
)
doParallel::registerDoParallel(cl=my_cluster)

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
subjects <- sort(subjects)
nSubjects <- length(subjects)
visits <- sapply(
  seq(4), function(x){paste0(
    "test_",
    ifelse(x>2, "RL", "LR"),
    ifelse(x%%2==1, 1, 2)
  )}
)
nVisits <- length(visits)

dg <- read.csv(hcp_dg_fname)
dg <- subset(dg, Subject %in% subjects)
dg$Gender <- factor(dg$Gender)
dg$Age <- factor(dg$Age)
dg2 <- dg[rep(seq(nrow(dg)), each=nVisits),]

check_FC <- function(FC) {
  stopifnot(!any(is.na(FC)))
  stopifnot(dim(FC)[seq(2)] == c(nVisits, nSubjects))
  stopifnot(all(dimnames(FC)[[1]] == visits))
  stopifnot(all(dimnames(FC)[[2]] == as.character(dg$Subject)))
  stopifnot(all(dg$Subject == as.numeric(gsub("subject_", "", dimnames(FC)[[3]]))))
}

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define FC files to process here.
FC_fnames <- paste0("CC2MP6_", c(
  "Base_1185_S1200",
  "DVARSdual_1185_S1200",
  "FD___og____1185_S1200",
  "FD___og_nfc_l4____1185_S1200",
  "Lev___ICA_1185_S1200"
))

for (ff in seq(length(FC_fnames))) {
  FC_fname <- file.path(dir_FC, paste0("FC_", FC_fnames[ff], ".rds"))
  cat(basename(FC_fname), "\n")
  # Load FC values.
  FCz <- readRDS(FC_fname)
  FCz <- FCz[,order(dimnames(FCz)[[2]]),,]
  
  # Format.
  if (length(dim(FCz)) == 3) {
    q <- dimnames(FCz)
    dim(FCz) <- c(dim(FCz), 1)
    dimnames(FCz)[seq(3)] <- q
    dimnames(FCz)[[4]] <- gsub("^FC_|\\..*", "", basename(FC_fname))
  }
  FCz <- FCz[,,,!apply(is.na(FCz), 4, any),drop=FALSE]
  # Skip non-kurtosis.
  if (grepl("Lev___ICA", FC_fname)) {
    FCz <- FCz[,,,c("Lev___ICA_kurt___4", "Lev___ICA_kurt___5")]
  }
  
  print(dimnames(FCz)[[4]])
  
  if (dim(FCz)[[4]] < 1) { next }
  # Check.
  check_FC(FCz)
  # z-transform.
  FCz <- psych::fisherz(FCz)
  
  # Prediction betas
  for (rvar in c("Sex") {
    time <- Sys.time()
    pbeta_fname <- file.path(dir_analysis, "pbetas", rvar, gsub("^FC_", "", basename(FC_fname)))
    if (!file.exists(pbeta_fname)) {
      cat("\tPrediction\n")
      pred <- vector("list", dim(FCz)[4])
      names(pred) <- dimnames(FCz)[[4]]
      for (pp in seq(dim(FCz)[4])) {
        cat("\tMethod ", pp, "\n")
        cat("===========================================\n")
        pred[[pp]] <- pbeta_make(FCz[,,,pp,drop=FALSE], what=rvar, parallel=TRUE)
      }
      saveRDS(pred, pbeta_fname)
    }
    print(time - Sys.time())
  }
  
  # Prediction 
  for (ii in seq(20)) {
    cat("Replication: ", ii, "\n")
    for (rvar in c("Sex")) {
      time <- Sys.time()
      pred_fname <- file.path(dir_analysis, "pred", rvar, gsub("^FC_", paste0(ii, "_"), basename(FC_fname)))
      if (!file.exists(pred_fname)) {
        cat("\tPrediction\n")
        pred <- vector("list", dim(FCz)[4])
        names(pred) <- dimnames(FCz)[[4]]
        for (pp in seq(dim(FCz)[4])) {
          cat("\tMethod ", pp, "\n")
          cat("===========================================\n")
          pred[[pp]] <- pred_make(FCz[,,,pp,drop=FALSE], ii, what=rvar, parallel=TRUE)
        }
        saveRDS(pred, pred_fname)
      }
      print(Sys.time() - time)
    }
  }
}

parallel::stopCluster(cl = my_cluster)
