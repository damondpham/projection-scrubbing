## -------------------------------------------------------------------------------------------------------------------------------------
source("0_SharedCode.R")
stopifnot(SharedCode_version == c(11,0))

stopifnot(COMPUTER == "MPro")

# Load list of S1200 subjects.
subjects <- readRDS(file.path(dir_data_misc, "HCP_subjects.rds"))

# Load demographic data.
hcp_dg <- read.csv(hcp_dg_fname, sep=",")
hcp_dg$Gender <- factor(hcp_dg$Gender)
hcp_dg$Age <- factor(hcp_dg$Age)

# Load iters + mean FDs for each scan.
# Mean FD is `NA` for scans that were missing RP data.
iters <- readRDS(file.path(dir_data_misc, "HCP_S1200_meanFD.rds"))
iters$subject <- as.numeric(as.character(iters$subject))

# Option: decrease the number of subjects in the study with this variable
N <- NULL
N <- 48 # must be < 74, maybe even stricter limit.
smallerSample <- !is.null(N)

## -------------------------------------------------------------------------------------------------------------------------------------

# Identify scans that are missing fMRI data.
iters$hasFiles <- FALSE
for (ii in seq(nrow(iters))) {
  subject <- as.numeric(iters[ii, "subject"])
  acquisition <- as.character(iters[ii, "acquisition"])
  test <- iters[ii, "test"]
  visit <- iters[ii, "visit"]
  
  prefix_ii <- paste0("rfMRI_REST", visit + (!test)*2, "_", acquisition)
  cat(subject, prefix_ii, "\n")
  
  # CIFTI data.
  cii_fname <- file.path(dir_HCP_test, subject, prefix_ii, paste0(prefix_ii, "_Atlas.dtseries.nii"))
  if (!file.exists(cii_fname)) { cat("\tMissing CIFTI data.\n"); next }
  
  # CompCor: requires NIFTI fMRI data & labels.
  # Has been computed for an older project, using the same script.
  cc_fname <- file.path(dir_HCP_test, subject, prefix_ii, paste0(prefix_ii, "_CompCor.rds"))
  if (!file.exists(cc_fname)) { cat("\tMissing CompCor data.\n"); next }
  
  # Realignment parameters.
  cc_fname <- file.path(dir_HCP_test, subject, prefix_ii, paste0("Movement_Regressors.txt"))
  if (!file.exists(cc_fname)) { cat("\tMissing RP data.\n"); next }
  
  iters[ii, "hasFiles"] <- TRUE
}

# Flag bad scans: any missing data, or mean FD > .3
iters$bad <- ifelse(is.na(iters$FDmean), TRUE, iters$FDmean > .3) | (!iters$hasFiles)

# Get list of subjects with all good scans.
iters2 <- subset(iters, !bad)
subjects2 <- table(iters2$subject)
subjects2 <- as.numeric(names(subjects2[subjects2==4]))
hcp_dg2 <- subset(hcp_dg, hcp_dg$Subject %in% subjects2)

# Get mean FD across all four scans for each subject, and add to `hcp_dg2`.
q <- aggregate(iters2$FDmean, list(iters2$subject), mean)
hcp_dg2$FDmean <- q[match(hcp_dg2$Subject, q[,1]),2]
  
## -------------------------------------------------------------------------------------------------------------------------------------

# Balance age and sex.
age_groups <- levels(hcp_dg2$Age)
to_rm <- setNames(vector("list", length(age_groups)), age_groups)
to_rmN <- vector("numeric")
# For each age group...
for (aa in seq(length(age_groups))) {
  age_group <- age_groups[aa]
  q <- sort(table(subset(hcp_dg2, Age==age_group)$Gender))
  # Determine which sex is overrepresented, and by how many people.
  overrep_sex <- names(q)[2]
  overrep_dif <- diff(q)
  # If making the sample even smaller, sample the underrepresented sex.
  if (smallerSample && aa < 4) {
    q <- hcp_dg2[which(hcp_dg2$Age==age_group & hcp_dg2$Gender!=overrep_sex),"Subject"]
    rmN <- length(q)-N
    with(set.seed(aa), to_rmN <- c(to_rmN, sample(q, rmN)))
    overrep_dif <- overrep_dif + rmN
  }
  # Get the median mean FD of the underrepresented sex.
  # Count how many subjects of overrepresented sex are above vs below this value.
  split <- median(subset(hcp_dg2, Age==age_group & Gender!=overrep_sex)$FDmean)
  split_idx <- list(
    under = which(hcp_dg2$Age==age_group & hcp_dg2$Gender==overrep_sex & hcp_dg2$FDmean < split),
    over = which(hcp_dg2$Age==age_group & hcp_dg2$Gender==overrep_sex & hcp_dg2$FDmean > split)
  )
  split_idx <- lapply(split_idx, function(x){hcp_dg2[x,"Subject"]})
  # Determine how many subjects of overrepresented sex to sample from the above
  #   group and the below group, to have the resulting median mean FD be similar
  #   to that of the underrepresented sex.
  split_dif <- diff(vapply(split_idx, length, 0))
  split_dif_half <- floor(abs(split_dif)/2)
  nsamp <- c(floor(abs(overrep_dif)/2), ceiling(abs(overrep_dif)/2))
  nsamp <- nsamp + sign(split_dif) * c(-split_dif_half, split_dif_half)
  # Adjust if the best choice of `nsamp` is beyond the available data.
  if (nsamp[1] > length(split_idx$under)) {
    nsamp <- c(length(split_idx$under), overrep_dif - length(split_idx$under))
  } else if (nsamp[1] <0) {
    nsamp <- c(0, overrep_dif)
  }
  # Sample the overrepresented sex.
  set.seed(aa+100)
  to_rm[[aa]] <- c(
    sample(split_idx$under, nsamp[1]),
    sample(split_idx$over, nsamp[2])
  )
}

# Get revised sample of subjects.
subjects3 <- subjects2[!(subjects2 %in% as.numeric(do.call(c, to_rm)))]
if (smallerSample) { subjects3 <- subjects3[!(subjects3 %in% to_rmN)] }
hcp_dg3 <- subset(hcp_dg2, hcp_dg2$Subject %in% subjects3)

# Remove subjects over 36 years old: too few in this category.
hcp_dg4 <- subset(hcp_dg3, Age != "36+")
subjects4 <- subjects3[subjects3 %in% hcp_dg4$Subject]

## -------------------------------------------------------------------------------------------------------------------------------------
table(hcp_dg4$Gender, hcp_dg4$Age)

## ----fig.width=6, fig.height=5--------------------------------------------------------------------------------------------------------
library(ggplot2); library(cowplot)
ggplot(hcp_dg4, aes(x=Age, fill=Gender, y=FDmean)) +
  geom_boxplot() +
  theme_cowplot() + ylab("Mean FD out/scan") +
  scale_fill_brewer(palette="Set1") + ggtitle("Motion amounts for large balanced sample")
ggsave(file.path(dir_data_misc, "HCP_S1200_LargeBalancedSample_FD.png"))

## -------------------------------------------------------------------------------------------------------------------------------------
saveRDS(subjects4, file.path(dir_data_misc, "HCP_S1200_LargeBalancedSample.rds"))

