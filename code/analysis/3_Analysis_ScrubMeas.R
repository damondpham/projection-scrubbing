## Prelim ----------------------------------------------------------------------
source("0_SharedCode.R")
source("0_Analysis_helper.R")
stopifnot(SharedCode_version == c(11,0))

## Load mean FD-----------------------------------------------------------------
meanFD_fname <- file.path(dir_data_misc, "HCP_meanFD.rds")
if (!file.exists(meanFD_fname)) {

  meanFD <- iters
  meanFD$FDmean <- NA

  for (ii in seq(nrow(meanFD))) {
    subject <- iters[ii, "subject"]
    acquisition <- iters[ii, "acquisition"]
    test <- iters[ii, "test"]
    visit <- iters[ii, "visit"]

    meanFD[ii, "FDmean"] <- mean(readRDS(file.path(
      dir_scrubMeas, "FD",
      paste0(
        "FD_", subject, "_v", visit + (!test)*2, "_", acquisition, ".rds"
      )
    ))$og)
  }
  meanFD <- meanFD[order(meanFD$FDmean),]

  saveRDS(meanFD, meanFD_fname)

}
meanFD <- readRDS(meanFD_fname)
#meanFD <- meanFD[order(as.numeric(rownames(meanFD))),]

meanFD2_fname <- file.path(dir_data_misc, "HCP_S1200_meanFD.rds")
if (!file.exists(meanFD2_fname)) {

  meanFD2 <- expand.grid(
    visit=seq(2),
    acquisition=c("LR", "RL"),
    subject=readRDS(file.path(dir_data_misc, "HCP_subjects.rds")),
    test=TRUE
  )
  meanFD2$FDmean <- NA

  for (ii in seq(nrow(meanFD2))) {
    subject <- meanFD2[ii, "subject"]
    acquisition <- meanFD2[ii, "acquisition"]
    visit <- meanFD2[ii, "visit"]

    FD_fname <- file.path(
      dir_scrubMeas, "FD",
      paste0(
        "FD_", subject, "_v", visit, "_", acquisition, ".rds"
      )
    )
    if (!file.exists(FD_fname)) { next }
    meanFD2[ii, "FDmean"] <- mean(readRDS(FD_fname)$og)
  }
  meanFD2 <- meanFD2[order(meanFD2$FDmean),]

  saveRDS(meanFD2, meanFD2_fname)
}
meanFD2 <- readRDS(meanFD2_fname)
meanFD2 <- meanFD2[order(as.numeric(rownames(meanFD2))),]

## Do overlap ------------------------------------------------------------------
base <- "CC2MP6"
overlap_fname <- file.path(dir_analysis, "overlap", paste0(base, ".rds"))
if (!file.exists(overlap_fname)) {
  overlap <- scrub_flag_overlap(base)
  saveRDS(overlap, overlap_fname)
}

overlap2_fname <- file.path(dir_analysis, "overlap", paste0(base, "_2.rds"))
if (!file.exists(overlap2_fname)) {
  overlap <- scrub_flag_overlap2(base)
  saveRDS(overlap, overlap2_fname)
}

## Load labels -----------------------------------------------------------------
q <- readRDS(file.path(dir_Parc, "ParcLabels.rds"))
PLabs <- q$PLabs
SLabs <- q$SLabs
Labs <- q$Labs
rm(q)

## Load splits -----------------------------------------------------------------
splits <- c(round(seq(0, 10)/10*nrow(meanFD)))
splits <- pmax(1, splits)
subs <- meanFD[round(splits),]
baseNames <- "CC2MP6" #c("MPP", "DCT4", "CC2MP6")

## 9 Sessions ------------------------------------------------------------------
library(ggplot2)
library(cowplot)

subs_9 <- subs[seq(2,10),]
dfgg <- vector("list", nrow(subs_9))
for (ii in seq(nrow(subs_9))) {
  subject <- subs_9[ii, "subject"]
  acquisition <- subs_9[ii, "acquisition"]
  test <- subs_9[ii, "test"]
  visit <- subs_9[ii, "visit"]
  id <- paste0(subject, "_v", visit + (!test)*2, "_", acquisition)
  scrubMeas <- list(
    lev = readRDS(file.path(dir_scrubMeas, "Lev/CC2MP6", paste0("LEV_", id, ".rds"))),
    #dv = readRDS(file.path(dir_scrubMeas, "DVARS/CC2MP6", paste0("DVARS_", id, ".rds"))),
    fd = readRDS(file.path(dir_scrubMeas, "FD", paste0("FD_", id, ".rds")))
  )
  dfgg[[ii]] <- data.frame(
    idx = ii,
    v = seq(1185),
    FD = pmin(1, scrubMeas$fd$og),
    modFD = pmin(1, scrubMeas$fd$og_nfc_l4),
    proj_ICA = pmin(10, scrubMeas$lev$measure$ICA_kurt / median(scrubMeas$lev$measure$ICA_kurt))
  )
}

dfgg <- do.call(rbind, dfgg)
dfgg <- tidyr::pivot_longer(dfgg, seq(3,5))
dfgg$idx <- factor(dfgg$idx, levels=seq(1, 9), labels=paste0((seq(1, 9))/10*100, "%"))
dfgg$name <- factor(dfgg$name, levels=c("FD", "modFD", "proj_ICA"))

p1 <- ggplot(subset(dfgg, name=="FD"), aes(x=v, y=value)) +
  scale_y_continuous(limits=c(0, 1), expand=expansion(mult=c(0, .02)), breaks=c(0, .2, .5)) +
  geom_hline(yintercept=c(.2, .3, .4, .5), linetype="dashed") +
  geom_point(size=.5, col=color_m["FD"]) + facet_grid(idx~.) +
  geom_hline(yintercept=c(.2, .3, .4, .5), linetype="dashed", alpha=.5) +
  geom_hline(yintercept=c(0)) +
  xlab("Index (Time Point)") + ylab("Value") +
  coord_cartesian(clip="off") +
  theme_cowplot()

p2 <- ggplot(subset(dfgg, name=="modFD"), aes(x=v, y=value)) +
  scale_y_continuous(limits=c(0, 1), expand=expansion(mult=c(0, .02)), breaks=c(0, .2, .5)) +
  geom_hline(yintercept=c(.2, .3, .4, .5), linetype="dashed") +
  geom_point(size=.5, col=color_m["modFD"]) + facet_grid(idx~.) +
  geom_hline(yintercept=c(.2, .3, .4, .5), linetype="dashed", alpha=.5) +
  geom_hline(yintercept=c(0)) +
  xlab("Index (Time Point)") + ylab("Value") +
  coord_cartesian(clip="off") +
  theme_cowplot()

p3 <- ggplot(subset(dfgg, name=="proj_ICA"), aes(x=v, y=value)) +
  scale_y_continuous(limits=c(0, 10), expand=expansion(mult=c(0, .02)), breaks=c(0, 3, 5)) +
  geom_hline(yintercept=c(3, 4, 5), linetype="dashed") +
  geom_point(size=.5, col=color_m["ICA"]) + facet_grid(idx~.) +
  geom_hline(yintercept=c(3, 4, 5), linetype="dashed", alpha=.5) +
  geom_hline(yintercept=c(0)) +
  xlab("Index (Time Point)") + ylab("Value") +
  coord_cartesian(clip="off") +
  theme_cowplot()

pdf(file.path(dir_analysis, "scrubMeasPlots/9sub.pdf"), width=14, height=17)
cowplot::plot_grid(p1, p2, p3, nrow=1)
dev.off()

# Denoised data and carpetplot -------------------------------------------------
splits <- c(
  round(seq(0, 10)/10*nrow(meanFD)),
  round((seq(0, 9)*2+1)/20*nrow(meanFD))
)
splits <- pmax(1, splits)
subs <- meanFD[round(splits),]
baseNames <- c("MPP", "DCT4", "36P", "CC2MP6", "CC2MP24", "32P")

get_intercept <- function(Y, design, int_idx=1){
  stopifnot(var(design[,int_idx])==0)
  Z <- design
	if(nrow(Y) != nrow(Z)) stop('Y and Z must have same number of rows')
 	invZtZ <- solve(crossprod(Z))
	betahat <- invZtZ %*% t(Z) %*% Y
	betahat[int_idx,]
}

useFIX <- identical(baseNames, "ICAFIX")
do <- c("ciiDN", "carpet", "pscrub")

# subs <- subs[rownames(subs) %in% c(34, 122, 126, 326, 97, 318),]

for (ii in seq(nrow(subs))) {
  # Get iteration info ---------------------------------------------------------
  subject <- subs[ii, "subject"]
  acquisition <- subs[ii, "acquisition"]
  test <- subs[ii, "test"]
  visit <- subs[ii, "visit"]

  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat(paste0(
    "Subject ", subject, ", ", as.character(acquisition), " ",
    ifelse(test, "test", "retest"), " ", visit,
    " (", ii, " of ", nrow(subs), ")", "\n"
  ))
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

  prefix <- paste0(subject, "_v", visit + (!test)*2, "_", acquisition)

  skip <- TRUE
  if ("saveDN" %in% do) {
    all_ciiDN_fname <- paste0(prefix, "_ciiDN_", baseNames, ".rds")
    if (!all(file.exists(file.path(dir_ciiDN, all_ciiDN_fname)))) {
      skip <- FALSE
    }
  }
  if ("carpet" %in% do) {
    all_carpet_fname <- paste0(prefix, "_carpetplot_", baseNames, ".pdf")
    if (!all(file.exists(file.path(dir_carpetPlots, all_carpet_fname)))) {
      skip <- FALSE
    }
  }
  if ("pscrub" %in% do) {
    all_lev_fname <- file.path(
      dir_scrubMeas, "Lev_withDirs", baseNames,
      paste0("LEV_", prefix, ".rds")
    )
    if (!all(file.exists(all_lev_fname))) {
      skip <- FALSE
    }
  }
  if (skip) { cat("Skipping.\n"); next }

  # Get nuisance regressors and CIFTI data. ------------------------------------
  cat("Input files.\n")
  # Mean signals
  ms_fname <- file.path(
    dir_meanSignals,
    paste0(
      subject, "_v", visit + (!test)*2, "_", acquisition, ".rds"
    )
  )
  ms <- data.frame(readRDS(ms_fname))
  ms <- as.matrix(ms[,c("wm", "csf", "wholebrain")])

  # CompCor
  cc_fname <- file.path(
    dir_CompCor,
    paste0(
      subject, "_v", visit + (!test)*2, "_", acquisition, ".rds"
    )
  )
  cc <- readRDS(cc_fname)
  cc <- cbind(cc$PCs$wm_cort[,seq(5)], cc$PCs$csf[,seq(5)], cc$PCs$wm_cblm[,seq(5)])
  cc <- scale(cc)

  # CIFTI
  fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)

  if (COMPUTER == "RED") {
    data_dir <- file.path(subject, "MNINonLinear", "Results", fname_prefix)
    fnames <- list(
      CIFTI = ifelse(useFIX,
        file.path(data_dir, paste0(fname_prefix, "_Atlas_hp2000_clean.dtseries.nii")),
        file.path(data_dir, paste0(fname_prefix, "_Atlas.dtseries.nii"))
      ),
      RP = file.path(data_dir, "Movement_Regressors.txt")
    )
    cii_fname <- file.path(ifelse(test, dir_HCP_test, dir_HCP_retest), fnames$CIFTI)

    # If retest, the data needs to be loaded.
    if (!test) {
      if (useFIX) {
        data_zip_MPP <- file.path(
          dir_HCP_retest_archive, paste0(subject, "_3T_rfMRI_REST", visit, "_fixextended.zip")
        )
        for (fname in fnames$CIFTI) {
          if (!file.exists(file.path(dir_HCP_retest, fname))) {
            cmd <- paste("unzip", data_zip_MPP, fname, "-d", dir_HCP_retest)
            system(cmd)
            stopifnot(file.exists(file.path(dir_HCP_retest, fname)))
          }
        }
        data_zip_MPP <- file.path(
          dir_HCP_retest_archive, paste0(subject, "_3T_rfMRI_REST", visit, "_preproc.zip")
        )
        for (fname in fnames$RP) {
          if (!file.exists(file.path(dir_HCP_retest, fname))) {
            cmd <- paste("unzip", data_zip_MPP, fname, "-d", dir_HCP_retest)
            system(cmd)
            stopifnot(file.exists(file.path(dir_HCP_retest, fname)))
          }
        }
      } else {
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
    }

    rp <- read.table(file.path(ifelse(test, dir_HCP_test, dir_HCP_retest), fnames$RP))

  # on Personal Computer -------------------------------------------------------
  } else {
    cii_fname <- ifelse(useFIX,
      file.path(dir_data_misc, "SelectedScans", paste0(subject, "_", fname_prefix, "_Atlas_hp2000_clean.dtseries.nii")),
      file.path(dir_data_misc, "SelectedScans", paste0(subject, "_", fname_prefix, "_Atlas.dtseries.nii"))
    )
    rp <- read.table(file.path(dir_data_misc, "SelectedScans", paste0(subject, "_", fname_prefix, "_Movement_Regressors.txt")))
  }

  # ----------------------------------------------------------------------------

  p6 <- scale(rp[seq(nDrop+1, hcp_T),seq(6)])
  p36 <- scale(cbind(ms, rbind(0, diff(ms)), rp))
  p36 <- scale(cbind(p36, p36^2)[seq(nDrop+1, hcp_T),])

  cii_list <- setNames(vector("list", length(baseNames)), baseNames)
  cii0 <- as.matrix(read_cifti(cii_fname, brainstructures="all"))[,seq(nDrop+1, hcp_T)]

  time <- Sys.time()

  # Get cleaned CIFTI data and leverage from each denoising method -------------
  cat("CIFTI data and leverage.\n")
  for (baseName in baseNames) {

    # Nuisance regression
    cat(baseName, "\n")
    nreg <- switch(baseName,
      CC2 = cbind(1, dct4, cc[,c(1,2,6,7,11,12)]),
      CC5 = cbind(1, dct4, cc),
      ICAFIX = as.matrix(rep(1, hcp_T - nDrop)),
      MPP = as.matrix(rep(1, hcp_T - nDrop)),
      DCT4 = cbind(1, dct4),
      `9P` = cbind(1, dct4, p36[,c(seq(3), seq(7,12))]),
      `36P` = cbind(1, dct4, p36),
      `32P` = cbind(1, dct4, p36[,setdiff(seq(36), c(3, 6, 21, 24))]),
      CC2MP6 = cbind(1, dct4, cc[,c(1,2,6,7,11,12)], p6),
      CC2MP24 = cbind(1, dct4, cc[,c(1,2,6,7,11,12)], p36[,c(seq(7,18), seq(25, 36))])
    )
    cii_b <- nuisance_regression(cii0, nreg)
    cii_i <- get_intercept(cii0, nreg)

    if ("ciiDN" %in% do) {
      ciiDN_fname <- paste0(prefix, "_ciiDN_", baseName, ".rds")
      if (!file.exists(ciiDN_fname)) {
        saveRDS(cii_b + cii_i, file.path(dir_ciiDN, ciiDN_fname))
      }
    }

    if ("carpet" %in% do) {
      carpet_fname <- paste0(prefix, "_carpetplot_", baseName, ".pdf")
      if (!file.exists(carpet_fname)) {
        carpetplot(t(cii_b),
          fname=file.path(dir_carpetPlots, carpet_fname),
          width=16, height=3
        )
      }
    }

    if ("pscrub" %in% do) {
      lev_fname <- file.path(
        dir_scrubMeas, "Lev_withDirs", baseName,
        paste0("LEV_", prefix, ".rds")
      )
      if (!file.exists(lev_fname)) {
        projections <- c(
          #"PCA", "PCA_kurt",
          #"fusedPCA", "fusedPCA_kurt",
          #"ICA",
          "ICA_kurt"
        )
        lev <- fMRIscrub:::pscrub_multi(
          t(cii_b), projection=projections,
          nuisance = NULL, get_dirs=TRUE, get_outliers=FALSE, verbose=TRUE
        )
        # Do not keep unnecessary metadata.
        lev$measure_info <- NULL
        # Save clever.
        saveRDS(lev, lev_fname)
      }
    }

    print(Sys.time() - time)
    time <- Sys.time()
    cat("\n")
  }

  # Unload retest data.
  if (!test && COMPUTER=="RED") {
    unlink(file.path(dir_HCP_retest, fnames$CIFTI))
    unlink(file.path(dir_HCP_retest, fnames$RP))
  }
}

# Scrub meas plots -------------------------------------------------------------

xii <- read_xifti(file.path(dir_data_misc, "rfMRI_empty.dtseries.nii"), brainstructures="all")

library(ggplot2)

for (ii in seq(nrow(subs))) {
  # Get iteration info ---------------------------------------------------------
  subject <- subs[ii, "subject"]
  acquisition <- subs[ii, "acquisition"]
  test <- subs[ii, "test"]
  visit <- subs[ii, "visit"]

  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat(paste0(
    "Subject ", subject, ", ", as.character(acquisition), " ",
    ifelse(test, "test", "retest"), " ", visit,
    " (", ii, " of ", nrow(subs), ")", "\n"
  ))
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

  # Load scrubbing measures ----------------------------------------------------
  cat("Scrubbing measures.\n")
  # FD
  fd_fname <- file.path(
    file.path(dir_scrubMeas, "FD"),
    paste0(
      "FD_", subject, "_v", visit + (!test)*2, "_", acquisition, ".rds"
    )
  )
  fd <- readRDS(fd_fname)
  # DVARS
  dv_fnames <- file.path(
    file.path(dir_scrubMeas, "DVARS", baseNames),
    paste0(
      "DVARS_", subject, "_v", visit + (!test)*2, "_", acquisition, ".rds"
    )
  )
  names(dv_fnames) <- baseNames
  dv <- lapply(dv_fnames, readRDS)
  # Leverage
  lev_fnames <- file.path(
    dir_scrubMeas, "Lev_withDirs", baseNames,
    paste0(
      "LEV_", subject,
      "_v", visit + (!test)*2, "_", acquisition, ".rds"
    )
  )
  names(lev_fnames) <- baseNames
  lev <- lapply(lev_fnames, readRDS)

  # Visualize the scrubbing measures. ------------------------------------------
  cat("Scrubbing measures visualizations.\n")
  plt <- list()
  plt$fd <- fMRIscrub:::scrub_plot(fd$og, cut = .3, colors = color_m["FD"], ylab="FD\n", ylim_max=1
  ) + guides(color="none") + scale_y_continuous(breaks=c(0, .3, .5, 1.0))
  plt$modfd <- fMRIscrub:::scrub_plot(fd$og_nfc_l4, cut = .2, colors = color_m["modFD"], ylab="modFD\n", ylim_max=1
  ) + guides(color="none") + scale_y_continuous(breaks=c(0, .2, .5, 1.0))
  plt$dv <- fMRIscrub:::scrub_plot(
    dv$CC2MP6$dv[c("DPD", "ZD")], cut = c(5, qnorm(1 - .05 / 1185)),
    colors = color_m[c("DVARS", "DVARS2")], ylab="DVARS", ylim_max=120,
    flag_intersect=TRUE
  ) + guides(color="none") + scale_y_continuous(breaks=c(0, 40, 80))
  plt$lev <- fMRIscrub:::scrub_plot(
    lev$CC2MP6$measure$ICA_kurt / median(lev$CC2MP6$measure$ICA_kurt), cut = 3,
    colors = color_m[c("ICA")], ylab="Leverage/Median", ylim_max=10
  ) + guides(color="none") + scale_y_continuous(breaks=c(0, 3, 5, 10))

  print(cowplot::plot_grid(plotlist=plt[c("fd", "modfd", "dv", "lev")], ncol=1, align="v"))
  ggplot2::ggsave(file.path(
    dir_scrubMeasPlots,
    paste0(subject, "_v", visit + (!test)*2, "_", acquisition, "_scrubMeas_xticks.pdf")
  ), width=5, height=4.5)

  for (pp in seq(length(plt))) {
    plt[[pp]] <- plt[[pp]] + theme(
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    )
  }

  print(cowplot::plot_grid(plotlist=plt[c("fd", "modfd", "dv", "lev")], ncol=1, align="v"))
  ggplot2::ggsave(file.path(
    dir_scrubMeasPlots,
    paste0(subject, "_v", visit + (!test)*2, "_", acquisition, "_scrubMeas.pdf")
  ), width=5, height=4.5)

  # ICs ------------------------------------------------------------------------
  cat("ICs.\n")
  ic <- lev$CC2MP6$ICA
  ic_sign <- colMeans(ic$S) > apply(ic$S, 2, median)
  ic$S[,!ic_sign] <- -ic$S[,!ic_sign]
  ic$M[,!ic_sign] <- -ic$M[,!ic_sign]

  # ickplot --------------------------------------------------------------------
  pcols <- as.character(
    ciftiTools::expand_color_pal(ciftiTools::make_color_pal(c("white", color_m["ICA"])), 255)$color
  )
  pcols <- c("#FFFFFF", pcols)
  add_sep <- function(x, q=5){
    n <- ncol(x); x <- cbind(x, 0);
    x <- x[,ifelse(seq(n*q)/q ==rep(seq(n), each=q), n+1, rep(seq(n), each=q))]
  }
  fMRIscrub::carpetplot(
    add_sep(abs(ic$M[,which(ic$highkurt)])), center=FALSE, qcut=.001, colors=pcols,
    width=15, height=6,
    fname=file.path(dir_scrubMeasPlots, paste0(subject, "_v", visit + (!test)*2, "_", acquisition, "_ickplot.pdf"))
  )

  # Component signal time series -----------------------------------------------
  cat("Component signal time series.\n")
  mask <- lev$CC2MP6$measure$ICA_kurt > median(lev$CC2MP6$measure$ICA_kurt)*5

  squish <- -.95
  for (squareVals in c(FALSE, TRUE)) {
    plt <- list()
    plt_kurt <- vector("numeric")
    ylim <- ifelse(squareVals, .4, .6)
    for (jj in seq(sum(ic$highkurt))) {
      jj2 <- which(ic$highkurt)[jj]
      jj2s <- gsub(" ", "0", format(jj2, width=3))
      v <- ic$M[,jj2]
      if (squareVals) { v <- v * abs(v) }
      plt[[paste0("ic", jj2)]] <- fMRIscrub:::scrub_plot(v, geom = "line", color=color_m["ICA"],
        ylab=paste0("IC\n#", jj2),
        ylim_min=-ylim, ylim_max=ylim
      ) + guides(color="none") +
        coord_cartesian(ylim=c(-ylim, ylim)) + scale_y_continuous(breaks=c(-ylim, 0, ylim), labels=c("120", "0.0", "120")) +
        geom_hline(yintercept=0, linetype="dashed", color="black")
      if (squareVals) {
        plt[[paste0("ic", jj2)]] <- plt[[paste0("ic", jj2)]] + theme_void()
      } else {
        plt_kurt[[paste0("ic", jj2)]] <- e1071::kurtosis(ic$M[,jj2], type=1)
        plt[[paste0("ic", jj2)]] <- plt[[paste0("ic", jj2)]] + theme(
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + annotate(
            "text", label=paste0("kurt: ", round(plt_kurt[[paste0("ic", jj2)]], 1)), x=5, hjust=0, y=.5)
      }

    }
    for (jj in seq(min(3, sum(ic$highkurt)))) {
      if (jj == 1) { hiVarLoKurt <- which(!ic$highkurt)[jj] }
      jj2 <- which(!ic$highkurt)[jj]
      jj2s <- gsub(" ", "0", format(jj2, width=3))
      v <- ic$M[,jj2]
      if (squareVals) { v <- v * abs(v) }
      plt[[paste0("ic", jj2)]] <- fMRIscrub:::scrub_plot(v, geom = "line", color=color_m["altGreen"],
        ylab=paste0("IC\n#", jj2),
        ylim_min=-ylim, ylim_max=ylim
      ) + guides(color="none") +
        coord_cartesian(ylim=c(-ylim, ylim)) + scale_y_continuous(breaks=c(-ylim, 0, ylim), labels=c("120", "0.0", "120")) +
        geom_hline(yintercept=0, linetype="dashed", color="black")
      if (squareVals) {
        plt[[paste0("ic", jj2)]] <- plt[[paste0("ic", jj2)]] + theme_void()
      } else {
        plt[[paste0("ic", jj2)]] <- plt[[paste0("ic", jj2)]] + theme(
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + annotate(
            "text", label=paste0("kurt: ", round(e1071::kurtosis(ic$M[,jj2], type=1), 1)), x=5, hjust=0, y=.5)
      }
    }

    for (kk in seq(2)) {
      if (!squareVals) {
        if (kk == 2) { plt <- plt[c(names(rev(sort(plt_kurt))[seq(min(length(plt_kurt), 4))]), paste0("ic", hiVarLoKurt))]  }
        print(cowplot::plot_grid(plotlist=plt, ncol=1, align="v"))
      } else {
        if (kk == 2) { next }
        plt2 <- vector("list", length=length(plt)*2-1)
        plt2[seq(1, length(plt2), 2)] <- plt
        print(cowplot::plot_grid(plotlist=plt2, ncol=1, align="v", rel_heights=c(rep(c(1, squish), length(plt)), 1)))
      }

      ggplot2::ggsave(file.path(
        dir_scrubMeasPlots,
        paste0(
          subject, "_v", visit + (!test)*2, "_", acquisition,
          "_ICmixing", ifelse(squareVals, "_sqVals", ""),
          ifelse(kk==2, "_subset", ""), ".pdf"
        )
      ), width=5, height=length(plt)*ifelse(squareVals, .5, 1))
    }
  }

  # ICs on surface -------------------------------------------------------------
  prfx <- paste0(subject, "_v", visit + (!test)*2, "_", acquisition, "_IC")
  cat("ICs on surface.\n")
  for (jj in seq(sum(ic$highkurt))) {
    jj2 <- which(ic$highkurt)[jj]
    jj2s <- gsub(" ", "0", format(jj2, width=3))
    y <- plot(newdata_xifti(xii, ic$S[,jj2]), zlim=c(-5, 5), what="surf", legend_embed=FALSE, height=1000, width=1400, fname=tempfile())
    z <- plot(newdata_xifti(xii, ic$S[,jj2]), zlim=c(-5, 5), what="vol", legend_embed=FALSE, height=1000, width=750, fname=tempfile())
    png(file.path(dir_scrubMeasPlots, paste0(prfx, jj2s, "_high.png")), width=1900, height=1000)
    ciftiTools::view_comp(c(y[1], z[1]), nrow=1, widths=c(1.7, 1))
    dev.off()
  }

  for (jj in seq(min(sum(!ic$highkurt), 3))) {
    jj2 <- which(!ic$highkurt)[jj]
    jj2s <- gsub(" ", "0", format(jj2, width=3))
    y <- plot(newdata_xifti(xii, ic$S[,jj2]), zlim=c(-5, 5), what="surf", legend_embed=FALSE, height=1000, width=1400, fname=tempfile())
    z <- plot(newdata_xifti(xii, ic$S[,jj2]), zlim=c(-5, 5), what="vol", legend_embed=FALSE, height=1000, width=750, fname=tempfile())
    png(file.path(dir_scrubMeasPlots, paste0(prfx, jj2s, "_low.png")), width=1900, height=1000)
    ciftiTools::view_comp(c(y[1], z[1]), nrow=1, widths=c(1.7, 1))
    dev.off()
  }
}
