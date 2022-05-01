# Organize (re-order, rename, group) the cortical and subcortical labels.

source("0_SharedCode.R")
stopifnot(SharedCode_version == c(11,0))

### `ParcMat` ------------------------------------------------------------------
# Read in Kong parcellation
parc_cii <- read_cifti(parc_cii_fname)
nVox <- nrow(parc_cii)
stopifnot(nrow(parc_cii$meta$cifti$labels[[1]])-1 == parc_res)
parc_cols <- parc_cii$meta$cifti$labels$parcels
parc_labs <- rownames(parc_cols)

# Begin to make `ParcMat`: vertices along rows, parcels along columns
ParcMat <- matrix(0, nrow=parc_res+1, ncol=nVox)
ParcMat[c(1+as.matrix(parc_cii)) + seq(0, nVox-1)*(parc_res+1)] <- 1
# table(colSums(ParcMat)) all ones

# Get CIFTI medial wall mask and subcortex labels.
# These masks are the same for all CIFTI files in the HCP.
hcp_cii <- read_cifti(file.path(
  dir_data_misc, "rfMRI_empty.dtseries.nii"
), brainstructures="all", idx=1)
hcp_mwall <- do.call(c, hcp_cii$meta$cortex$medial_wall_mask)
hcp_subLabs <- hcp_cii$meta$subcort$labels

# Remove medial wall vertices; add subcortical "parcels" (each brainstructure)
ParcMat <- ParcMat[,hcp_mwall,drop=FALSE]
ParcMatSub <- t(model.matrix(~0+hcp_subLabs))
ParcMatSub <- ParcMatSub[seq(3,21),] # Remove empty left & right cortex parcels

# Combine cortical parcels and subcortical parcels into one multiplication matrix
ParcMat <- cbind(ParcMat, matrix(0, nrow=nrow(ParcMat), ncol=ncol(ParcMatSub)))
ParcMatSub <- cbind(matrix(0, nrow=nrow(ParcMatSub), ncol=ncol(ParcMat)-ncol(ParcMatSub)), ParcMatSub)
ParcMat <- rbind(ParcMat, ParcMatSub)

# By dividing each row by its sum, we can use matrix multiplication to compute the mean fMRI signal
ParcMat <- t(ParcMat / apply(ParcMat, 1, sum))
dimnames(ParcMat) <- NULL

# Number of vertices/voxels per parcel
#   summary(apply(ParcMat > 0, 2, sum))
# 806 vertices not included in any parcel
#   (location: below medial wall)
#   sum(ParcMat[,1])
ParcMat <- ParcMat[,seq(2,ncol(ParcMat))]
saveRDS(ParcMat, parc_fname)

### `PLabs` --------------------------------------------------------------------
PLabs <- gsub("17networks_", "", parc_labs)
PLabs <- PLabs[PLabs != "???"] # medial wall
PLabs <- strsplit(PLabs, "_")
PLabs <- data.frame(
  label = vapply(PLabs, paste, collapse="_", ""),
  bs = ifelse(vapply(PLabs, '[', i=1, "")=="LH", "left", "right"),
  network2 = vapply(PLabs, '[', i=2, ""),
  group = vapply(PLabs, '[', i=3, ""),
  group2 = vapply(PLabs, '[', i=2, ""),
  gidx = as.numeric(vapply(PLabs, '[', i=4, "")) # formerly gidx
)
PLabs$network <- gsub("A$|B$|C$", "", PLabs$network2)
# PLabs$network <- gsub("Language|Aud", "LanAud", PLabs$network)

PLabs$idx <- seq(parc_res)
# `idx2`: # network, bs/hemisphere, network2, group, gidx
PLabs$idx2 <- order(order(PLabs$network))
PLabs$label_color <- apply(
  parc_cols[
    match(paste0("17networks_", PLabs$label), parc_labs),
    c("Red", "Green", "Blue")
  ], 
  1, 
  function(x){rgb(x[1], x[2], x[3])}
)

PLabs$group_color <- NA
PLabs$network_color <- NA

# Organize the colors for each subcortical structure
# Source: https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT
# Modification: cerebellum was lightened and thalamus was darkened, to enhance contrast between the two

### `SLabs` --------------------------------------------------------------------
SLabs <- c(
  ACCUMBENS_LEFT = "#007b7d",
  ACCUMBENS_RIGHT = "#007b7d",
  AMYGDALA_LEFT = "#dcd814",
  AMYGDALA_RIGHT = "#dcd814",
  BRAIN_STEM = "#ee9349",
  CAUDATE_LEFT = "#30998f",
  CAUDATE_RIGHT = "#30998f",
  CEREBELLUM_LEFT = "#d27eb9",
  CEREBELLUM_RIGHT = "#d27eb9",
  DIENCEPHALON_VENTRAL_LEFT = "#f4be91",
  DIENCEPHALON_VENTRAL_RIGHT = "#f4be91",
  HIPPOCAMPUS_LEFT = "#fefb7d",
  HIPPOCAMPUS_RIGHT = "#fefb7d",
  PALLIDUM_LEFT = "#7fd4d6",
  PALLIDUM_RIGHT = "#7fd4d6",
  PUTAMEN_LEFT = "#adf4fe",
  PUTAMEN_RIGHT = "#adf4fe",
  THALAMUS_LEFT = "#6e74c9",
  THALAMUS_RIGHT = "#6e74c9"
)
SLabs2 <- ciftiTools:::substructure_table()
names(SLabs) <- SLabs2[match(names(SLabs), SLabs2$Original_Name), "ciftiTools_Name"]
rm(SLabs2)

# Subcortical labels
SLabs <- data.frame(
  label = names(SLabs),
  bs = "subcort",
  network2 = "Subcort",
  group = c("Basal Ganglia", "Basal Ganglia", "Hippocampus-Amygdala", 
            "Hippocampus-Amygdala", "Brainstem-Diencephalon", "Basal Ganglia", 
            "Basal Ganglia", "Cerebellum", "Cerebellum", "Brainstem-Diencephalon", 
            "Brainstem-Diencephalon", "Hippocampus-Amygdala", 
            "Hippocampus-Amygdala", "Basal Ganglia", "Basal Ganglia", 
            "Basal Ganglia", "Basal Ganglia", "Thalamus", "Thalamus"),
  group2 = NA,
  gidx = NA,
  network = "Subcort", #ifelse(grepl("Cerebellum", names(SLabs)), "Cerebellum", "Subcort"),
  idx = parc_res + seq(1, 19),
  idx2 = NA,
  label_color = SLabs,
  group_color = NA,
  network_color = "#999999"
)
SLabs$group <- factor(
  SLabs$group, 
  levels=c("Cerebellum", "Basal Ganglia", "Hippocampus-Amygdala", "Brainstem-Diencephalon", "Thalamus")
)
SLabs$idx2 <- parc_res + order(order(SLabs$group))

### Merge `PLabs` and `SLabs` --------------------------------------------------
NetworkAlt <- data.frame(
  network = c("Cont", "Default", "DorsAttn", "Language", "Aud", 
              "SalVenAttn", "SomMot", "Visual", "Subcort"),
  full = c("Frontoparietal Control", "Default", "Dorsal Attention", "Language", "Auditory", 
           "Salience/Ventral Attention", "Somatomotor", "Visual", "Subcortical"),
  abrev = c("FP", "DMN", "dATTN", "L", "A", 
            "vATTN", "SMOT", "V", "SC")
)
GroupAlt <- data.frame(
  group = c("Cerebellum", "Basal Ganglia", "Hippocampus-Amygdala", "Brainstem-Diencephalon", "Thalamus"),
  abrev= c("CBLM", "BG", "H&A", "B&D", "TH")
)

# Merge the labels
Labs <- rbind(PLabs, SLabs)
rownames(Labs) <- NULL
Labs$bs <- factor(Labs$bs, levels=c("left", "right", "subcort"))
Labs$network <- factor(Labs$network, levels=NetworkAlt$network)

Labs$network_first[order(Labs$idx2)] <- c(
  TRUE, 
  Labs$network[order(Labs$idx2)][seq(parc_res2-1)] != Labs$network[order(Labs$idx2)][seq(2,parc_res2)]
)
Labs$group_first[order(Labs$idx2)] <- as.logical(c(
  rep(FALSE, parc_res), TRUE,
  Labs$group[order(Labs$idx2)][seq(parc_res+1, parc_res2-1)] != Labs$group[order(Labs$idx2)][seq(parc_res+2, parc_res2)]
))

Labs2 <- subset(Labs, network_first | group_first)

Labs$network2 <- NetworkAlt$full[match(Labs$network, NetworkAlt$network)]
Labs$network3 <- NetworkAlt$abrev[match(Labs$network, NetworkAlt$network)]
Labs$group2[Labs$bs=="subcort"] <- GroupAlt$abrev[match(Labs$group[Labs$bs=="subcort"], GroupAlt$group)]

Labs$network_color[seq(parc_res)] <- Labs2$label_color[match(Labs$network[seq(parc_res)], Labs2$network)]
Labs$network_color[grepl("Default", Labs$network2)] <- Labs$label_color[which(Labs$group2=="DefaultB")[1]] # Manually change this one.
Labs$network_color[grepl("Auditory", Labs$network2)] <- Labs$label_color[which(Labs$group2=="DefaultA")[1]] # Manually change this one.
Labs$group_color[seq(parc_res+1, parc_res2)] <- Labs2$label_color[match(Labs$group[seq(parc_res+1, parc_res2)], Labs2$group)]
# Lighter subcortical shades
Labs$group_color[Labs$group_color=="#007b7d"] <- "#30998f"
Labs$group_color[Labs$group_color=="#dcd814"] <- "#fefb7d"
Labs$group_color[Labs$group_color=="#ee9349"] <- "#f4be91"

saveRDS(list(Labs=Labs, PLabs=PLabs, SLabs=SLabs), file.path(dir_Parc, "ParcLabels.rds"))

parc_cii2 <- parc_cii
parc_cii2$meta$cifti$labels$parcels[seq(parc_res)+1,c("Red", "Green", "Blue")] <- c(t(col2rgb(Labs$network_color[seq(parc_res)], alpha=FALSE))/255)
plot(parc_cii2, borders=TRUE, fname=file.path(dir_Parc, "Kong2022_recolor"), title="Kong 2022 colored by network")
