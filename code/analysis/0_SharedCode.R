SharedCode_version = c(11, 0)

# Change these depending on where files are in your computer -------------------
COMPUTER <- ifelse(
  grepl("Carbonate|/N", getwd()), 
  "RED", 
  ifelse(grepl("Volume", getwd()), "MPro", "Personal")
)

secret_fnames <- "0_Secret-FilePaths.txt"
if (!file.exists(secret_fnames)) { 
  if (file.exists(file.path("code/analysis", secret_fnames))) {
    secret_fnames <- file.path("code/analysis", secret_fnames)
  } else if (file.exists(file.path("../analysis", secret_fnames))) {
    secret_fnames <- file.path("../analysis", secret_fnames)
  }
}
secret_fnames <- readLines(secret_fnames)

# HCP data, for RED
dir_HCP_test <- switch(COMPUTER,
  RED = secret_fnames[1],
  MPro = secret_fnames[10]
)
dir_HCP_retest_archive <- secret_fnames[2]
# where to unzip HCP retest files for temporary use
dir_HCP_retest <- secret_fnames[3]

# Project directory
dir_project <- switch(COMPUTER,
  RED = secret_fnames[4],
  MPro = secret_fnames[5],
  Personal = secret_fnames[6]
)

# Connectome Workbench
wb_path <- switch(COMPUTER,
  RED = secret_fnames[7],
  MPro = secret_fnames[8],
  Personal = secret_fnames[9]
)

# Below should be the same -----------------------------------------------------
dir_data_misc <- file.path(dir_project, "data")
dir_Parc <- file.path(dir_data_misc, "Parc")
dir_CompCor <- file.path(dir_project, "analysis-results/1_CompCor")
dir_meanSignals <- file.path(dir_project, "analysis-results/1_MeanSignals")
dir_scrubMeas <- file.path(dir_project, "analysis-results/2_ScrubMeas")
dir_FC <- file.path(dir_project, "analysis-results/3_FC")
dir_FCval <- file.path(dir_project, "analysis-results/3_FCval")
dir_FC_shorter <- file.path(dir_project, "analysis-results/3_FC_shorter")
dir_analysis <- file.path(dir_project, "analysis-results/4_Analysis")
dir_carpetPlots <- file.path(dir_analysis, "carpetPlots")
dir_ciiDN <- file.path(dir_analysis, "ciiDN")
dir_scrubMeasPlots <- file.path(dir_analysis, "scrubMeasPlots")
dir_plots <- file.path(dir_project, "plots")

parc_cii_fname <- file.path(dir_Parc, "Schaefer2018_400Parcels_Kong2022_17Networks_order.dlabel.nii")
hcp_dg_fname <- file.path(dir_project, "data/unrestricted_HCP_demographics.csv")

# subjects <- list.files(dir_HCP_test)
# Only use retest subjects
subjects <- c(
  103818, 187547,
  105923, 192439,
  111312, 194140,
  114823, 195041,
  115320, 200109,
  122317, 200614,
  125525, 204521,
  130518, 250427,
  135528, 287248,
  137128, 341834,
  139839, 433839,
  143325, 562345,
  144226, 599671,
  146129, 601127,
  149337, 627549,
  149741, 660951,
  151526, 662551,
  158035, 783462,
  169343, 859671,
  172332, 861456,
  175439, 877168,
  177746, 917255,
  185442
)
omit_subjects <- c(
  627549, # does not have retest zip files
  341834, # does not have retest FIX files in the zip (REST1 RL)
  143325  # retest MPP CIFTI for REST2_RL is truncated to 939 timepoints
)
subjects <- subjects[!(subjects %in% omit_subjects)]

# Packages
#install.packages("fMRIscrub")
library(fMRIscrub)
stopifnot(utils::packageVersion("fMRIscrub") >= "0.8.6")
#install.packages("ciftiTools")
library(ciftiTools)
stopifnot(utils::packageVersion("ciftiTools") >= "0.8.0")
ciftiTools.setOption("wb_path", wb_path)

hcp_T <- 1200
nDrop <- 15
dct4 <- scale(dct_bases(hcp_T-nDrop, 4))
iters <- expand.grid(
  visit=seq(2),
  test=c(TRUE, FALSE),
  acquisition=c("LR", "RL"),
  subject=subjects
)

parc_fname <- file.path(dir_Parc, "ParcMat.rds")
parc_res <- 400
parc_res2 <- parc_res + 19

color_m <- c(
  Baseline = "#999999",
  modFD = "#7f2811",# "#7c532e",
  FD = "#cb3b02",
  ICA = "#1555bb",
  PCA = "#cd9bea",
  FusedPCA = "#8fbffd",
  DVARS = "#e6ab02",
  DVARS2 = "#bcb09c",
  ICAFIX = "#4c7a27",
  `36P` = "#914aa0",
  lightGray= "#cccccc",
  darkGray = "#333333",
  altGreen = "#95af5c",
  altTan = "#d8b999",
  altPink = "#E78AC3"
)
