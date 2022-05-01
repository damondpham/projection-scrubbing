#' Create a qualitative color legend.
qual_cleg <- function(labels, colors, ncol, title=NULL, bgcol="#ffffff") {
  scale_color <- setNames(as.list(colors), labels)
  df_cleg <- unique(data.frame(value=1, label=labels, stringsAsFactors=FALSE))
  plt <- ggplot(data = df_cleg, aes(x=value, y=value, color=labels)) +
    geom_point(size=5, shape=15) + theme_bw() +
    scale_color_manual(breaks=df_cleg$label, values=scale_color, name=title) +
    guides(color=guide_legend(label.theme=element_text(color="black"), ncol=ncol, override.aes=list(fill=bgcol))) +
    theme(
      legend.title=element_text(size=rel(1.5)), 
      legend.text=element_text(color="black", size=rel(1.2)),
      legend.background=element_rect(fill=bgcol),
      legend.box.background = element_rect(fill=bgcol)
    )
  ggpubr::as_ggplot(ggpubr::get_legend(plt))
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------
plot(
NetMat, 
 fname=plt_fmt("Ruby17Net_onSurf"), idx=2, borders=TRUE, legend_fname=FALSE
)

plot(
NetMat, 
 fname=plt_fmt("Ruby17Net_onSurf_reColor"), idx=1, borders=TRUE, legend_fname=FALSE
)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
PLabs2 <- subset(Labs[seq(parc_res),], Labs$network_first[seq(parc_res)])
PLabs2$network4 <- paste0(PLabs2$network2, " (", PLabs2$network3, ")")
PLabs2 <- PLabs2[order(PLabs2$idx2),]

pdf(plt_fmt("Ruby17Net_ColorLegend.pdf"))

print(
  qual_cleg(PLabs2$network4, PLabs2$network_color, ncol=1)
)

dev.off()

rm(PLabs2)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# adapted from ciftiTools:::get_data_meta_from_cifti_xml
get_subcort_coordlist <- function(cifti_fname, intent=3000) {
  intent <- 3006
  xml <- ciftiTools:::xml_cifti(cifti_fname)
  ii <- which(names(xml$CIFTI$Matrix) == "MatrixIndicesMap")[2]
  xml <- xml$CIFTI$Matrix[[ii]]

  bs_names <- as.character(lapply(xml, function(x){attr(x, "BrainStructure")}))
  bs_names <- gsub("CIFTI_STRUCTURE_", "", bs_names)

  meta <- list(
    brainstructures=NULL,
    cortex_left_mwall=NULL,
    cortex_right_mwall=NULL,
    subcort_trans_mat=NULL,
    subcort_labs=NULL,
    subcort_dims=NULL
  )

  # Subcortical Transformation Matrix
  if ("Volume" %in% names(xml)) {
    meta$brainstructures <- c(meta$brainstructures, "subcortical")
    vol_tmat <- strsplit(xml$Volume$TransformationMatrixVoxelIndicesIJKtoXYZ[[1]], "\n")[[1]]
    vol_tmat <- vol_tmat[vol_tmat != ""]
    vol_tmat <- do.call(rbind, lapply(vol_tmat, function(x){as.numeric(strsplit(x, " ")[[1]])}))
    meta$subcort_trans_mat <- vol_tmat

    meta$subcort_dims <- as.numeric(strsplit(attr(xml$Volume, "VolumeDimensions"), ",")[[1]])
  }

  # Subcortical
  sub_names <- substructure_table()$Original_Name[3:nrow(substructure_table())]
  bs_subcort <- bs_names %in% sub_names
  if (any(bs_subcort)) {
    subcort_dat <- vector("list", length(sub_names))
    for (ii in seq_len(length(sub_names))) {
      sub_name <- sub_names[ii]
      xml_idx <- which(bs_names == sub_name)
      # stopifnot(length(xml_idx) > 1) # By CIFTI definition, only one match should be present.
      if (length(xml_idx) > 0) {
        subcort_dat[[ii]] <- cbind(do.call(
          rbind,
          lapply(strsplit(
            strsplit(xml[[xml_idx]]$VoxelIndicesIJK[[1]], "\n")[[1]], " "
          ), as.numeric)
          # Add two to skip cortex L/R. Line after already adds the other 1.
        ), ii+1)
      }
    }
    # Add one to start indexing at 1 instead of zero
    subcort_dat <- do.call(rbind, subcort_dat) + 1
  }
  
  colnames(subcort_dat) <- c("x", "y", "z", "s")
  subcort_dat <- as.data.frame(subcort_dat)
  subcort_dat
}

q <- get_subcort_coordlist(ciftiTools:::ciftiTools.files()$cifti["dscalar_ones"])
q$s <- factor(SLabs$label[q$s-2], levels=SLabs$label)
#qlims <- as.data.frame(apply(q[,seq(3)], 2, function(x){c(min(x), max(x))}))

# MNI T1 template
mni <- RNifti::readNifti(
  system.file("extdata", "MNI152_T1_2mm_crop.nii.gz", package="ciftiTools")
)
mni2 <- do.call(expand.grid, lapply(as.list(dim(mni)), seq))
colnames(mni2) <- c("x", "y", "z")
mni2$v <- as.vector(mni)
# Crop
mni2$v[mni2$x < 17] <- 0
mni2$v[mni2$x > 75] <- 0
mni2$v[mni2$z > 60] <- 0
mni <- mni2; rm(mni2)

plot_slice <- function(q, slicedim="z", sliceval=20, mni=NULL){
  ctable <- as.list(SLabs$label_color)
  names(ctable) <- SLabs$label
  
  dim_x <- switch(slicedim, x="y", y="x", z="x")
  dim_y <- switch(slicedim, x="z", y="z", z="y")
  
  if (!is.null(mni)) {
    mni_x <- mni[mni$v != 0,dim_x]
    mni_y <- mni[mni$v != 0,dim_y]
    xlim <- c(max(min(mni_x), min(q[,dim_x]) - 2), min(max(mni_x), max(q[,dim_x]) + 2))
    ylim <- c(max(min(mni_y), min(q[,dim_y]) - 2), min(max(mni_y), max(q[,dim_y]) + 2))
  } else {
    xlim <- c(min(q[,dim_x]) - 2, max(q[,dim_x]) + 2)
    ylim <- c(min(q[,dim_y]) - 2, max(q[,dim_y]) + 2)
  }
  
  if (!is.null(mni)) {
    mni <- subset(mni, mni[[slicedim]]==sliceval)
    q1 <- quantile(mni$v, .02, na.rm=TRUE)
    mni$c <- (mni$v - q1) / ( quantile(mni$v, .98, na.rm=TRUE) - q1 )
    # # Increase contrast
    # mni$c <- 1/(1 + exp(-2*(2*mni$c - 1)))
    # Darken
    mni$c <- (.8 * mni$c) / (.8 - mni$c + 1)
    # Raise blackpoint
    mni$c <- mni$c * .9 + .1
    mni$c <- scales::colour_ramp(c("black", "white"))(mni[mni[,slicedim]==sliceval,"c"])
    q <- subset(q, q[[slicedim]]==sliceval)
    # for (ii in seq(nrow(q))) {
    #   mni[mni[[dim_x]]==q[ii,dim_x] & mni[[dim_y]]==q[ii,dim_y],"c"] <- ctable[[q[ii,"s"]]]
    # }
    
    plt <- ggplot() + 
      annotate("raster", x=mni[,dim_x], y=mni[,dim_y], fill=mni$c, interpolate=TRUE) +
      annotate("raster", x=q[,dim_x], y=q[,dim_y], fill=ctable[q$s])
  } else {
    plt <- ggplot() +
      geom_raster(aes_string(x=dim_x, y=dim_y, fill="s"), data = subset(q, q[[slicedim]]==sliceval)) 
  }
    
  plt <- plt + theme_void() + 
    coord_equal(xlim=xlim, ylim=ylim, expand=FALSE) + 
    scale_fill_manual(values=ctable) +
    theme(legend.position="none") + ggtitle(paste0(" ", paste(slicedim, "=", sliceval)))
  
  return(plt)
}

pdf(plt_fmt("SubcortSlices.pdf"))
for (slicedim in c("x", "y", "z")) {
  if (slicedim!="y") { next }
  cat(slicedim ,"\n")
  sliceseq <- seq(min(q[[slicedim]]), max(q[[slicedim]]))
  for (ss in sliceseq) {
    print(plot_slice(q, slicedim, ss))
  }
} 
dev.off()

pdf(plt_fmt("SubcortSlices_withMNI.pdf"))
for (slicedim in c("x", "y", "z")) {
  if (slicedim!="y") { next }
  cat(slicedim ,"\n")
  sliceseq <- seq(min(q[[slicedim]]), max(q[[slicedim]]))
  for (ss in sliceseq) {
    print(plot_slice(q, slicedim, ss, mni=mni))
  }
} 
dev.off()

pdf(plt_fmt("SubcortSlices_withMNI2.pdf"))
ykeep <- c(42, 59, 66)
q <- subset(q, y %in% ykeep)
sliceseq <- ykeep
plts <- vector("list", length(sliceseq))
for (ss in seq(length(sliceseq))) {
  plts[[ss]] <- plot_slice(q, "y", sliceseq[ss], mni=mni) + 
    theme(
      plot.margin = unit(c(0, .2, 0, .2), "lines"),
      plot.title = element_text(size=11, margin = margin(b = -1.4, unit = "lines"))
    )
}
print(
  cowplot::plot_grid(plotlist=plts, nrow=1, aligh="v", axis="tb")
)
dev.off()


## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------
## view_subcort_struc <- function(){
##   ## make rgl window as big as you want it
##   rgl::open3d()
##   rgl::par3d(windowRect = c(20, 20, 2500, 2500))
## 
##   ##################################################################
##   # Draw each brainstructure separately
##   cii_ones <- read_cifti(ciftiTools:::demo_files()$cifti["dscalar_ones"], brainstructures="all")
##   sub_vol <- unmask_vol(as.numeric(cii_ones$meta$subcort$labels)-2, cii_ones$meta$subcort$mask, 0)
##   for(ll in SLabs$idx-parc_res){
##     misc3d::contour3d(sub_vol==ll, color=SLabs$color2[ll], level=0, smooth=TRUE, material="dull", specular="black", add=TRUE)
##   }
## 
##   #################################################################
##   # Add "cap" at the bottom of the brainstem to close the mesh.
##   # (Keep RGL window open.)
## 
## 
##   # First, convert to verts and faces ------------------------
## 
##   mesh <- misc3d::contour3d(cii_ones$meta$subcort$mask, level=0, smooth=TRUE, material="dull", specular="black", draw=FALSE)
## 
##   tri <- matrix(nrow=nrow(mesh$v1)*3, ncol=ncol(mesh$v1))
##   tri[seq(1, nrow(tri), 3),] <- mesh$v1
##   tri[seq(1, nrow(tri), 3)+1,] <- mesh$v2
##   tri[seq(1, nrow(tri), 3)+2,] <- mesh$v3
## 
##   #tri <- apply(tri, 2, function(x){x - min(x) + 1})
## 
##   verts = unique(tri)
##   get_ID <- function(x){x[,1] + x[,2]*(max(x[,1])+1) + x[,3]*(max(x[,1]+1) * max(x[,2]+1))}
##   vertID <- get_ID(verts)
##   faces <- matrix(match(get_ID(tri), vertID), nrow=3)
## 
##   # Then, manually add that brainstem bottom -------------------
## 
##   # It's located where the third coordiante column == 1
##   # verts[verts[,3]==1,]
##   # I wrote down the coordinates and manually came up with a mesh to "patch" the hole
##   v1 = which(verts[,3]==1)
## 
##   new_verts <- rbind(
##     c(26+18,26+16),
##     cbind(27:31+17,26+16),
##     cbind(27:31+17,27+16)
##   )
##   v2 <- 1:nrow(new_verts) + nrow(verts)
##   verts <- rbind(verts, cbind(new_verts[,1], new_verts[,2], 1))
## 
##   new_faces <- rbind(
##     c(v1[1],v1[2],v2[1]),
##     c(v1[2],v1[9],v2[1]),
##     c(v1[9],v1[11],v2[7]),
##     c(v1[7],v1[8],v2[6]),
##     c(v1[10],v1[15],v2[11]),
##     c(v1[8],v1[10],v2[6]),
##     c(v2[6],v2[11],v1[10])
##   )
## 
##   gaps <- list(
##     cbind(v1[3:7], v2[2:6], v2[7:11], v1[11:15]),
##     cbind(v1[c(1,3)], v2[c(1,2)], c(v1[9],v2[7]))
##   )
##   for(gap in gaps){
##     for(i in 2:nrow(gap)){
##       for(j in 2:ncol(gap)){
##         new_faces <- rbind(new_faces, c(gap[i,j], gap[i-1,j], gap[i-1,j-1]))
##         new_faces <- rbind(new_faces, c(gap[i,j], gap[i,j-1], gap[i-1,j-1]))
##       }
##     }
##   }
## 
## 
##   # Plot -------------------------------------------------------
##   # (Keep RGL window open.)
## 
##   #faces <- cbind(faces, t(new_faces))
##   #colors <- c(mesh$color, rep("#779FB000", nrow(new_faces)))
##   #mesh2 <- tmesh3d(t(cbind(verts, 1)), faces, meshColor = "faces")
##   #mesh2 <- addNormals(mesh2)
## 
##   mesh3 <- tmesh3d(t(cbind(verts, 1)), t(new_faces), meshColor = "faces")
##   # mesh3 <- addNormals(mesh3) # Doesn't work
## 
##   # Just the patch
##   shade3d(
##     mesh3,
##     color="#779FB000",
##     specular="black",
##     alpha=1
##   )
## 
##   ###############################################################
##   # Photoshoot
##   #
## 
##   photoshoot <- function(){
##     ZOOM = .8
##     FOV = 0
## 
##     rgl::rgl.viewpoint(userMatrix=diag(4), zoom=ZOOM, fov=FOV)
##     rgl::rgl.snapshot("subcort_top.png", top=TRUE)
## 
##     rot <- list(
##       bottom = rbind(
##         c(-1, 0, 0, 0),
##         c( 0, 1, 0, 0),
##         c( 0, 0,-1, 0),
##         c( 0, 0, 0, 1)
##       ),
##       back = rbind(
##         c(-1, 0, 0, 0),
##         c( 0, 0,-1, 0),
##         c( 0,-1, 0, 0),
##         c( 0, 0, 0, 1)
##       )
##     )
## 
##     rgl::rgl.viewpoint(userMatrix=rot$bottom, zoom=ZOOM, fov=FOV)
##     rgl::rgl.snapshot(plt_fmt("Subcort_view_bottom.png"), top=TRUE)
## 
##     rgl::rgl.viewpoint(userMatrix=rot$back, zoom=ZOOM, fov=FOV)
##     rgl::rgl.snapshot(plt_fmt("Subcort_view_back.png"), top=TRUE)
## 
##     rgl::rgl.viewpoint(theta = 0, phi = 90, zoom=ZOOM, fov=FOV)
##     rgl::rgl.snapshot(plt_fmt("Subcort_view_front.png"), top=TRUE)
## 
##     rgl::rgl.viewpoint(theta = 90, phi = 0, zoom=ZOOM, fov=FOV)
##     rgl::rgl.snapshot(plt_fmt("Subcort_view_right.png"), top=TRUE)
## 
##     rgl::rgl.viewpoint(theta = -90, phi = 0, zoom=ZOOM, fov=FOV)
##     rgl::rgl.snapshot(plt_fmt("Subcort_view_left.png"), top=TRUE)
##   }
## 
##   photoshoot()
## }
## view_subcort_struc()


## -----------------------------------------------------------------------------------------------------------------------------------------------------
pdf(plt_fmt("Subcort_ColorLegend.pdf"))

SLabs2 <- data.frame(
  label_color=Labs[order(Labs$idx2)[seq(19)+parc_res],"label_color"], 
  label=Labs[order(Labs$idx2)[seq(19)+parc_res],"label"], 
  stringsAsFactors = FALSE
)
SLabs2$label <- gsub("-.*", "", SLabs2$label)
SLabs2 <- unique(SLabs2)
print(
  qual_cleg(SLabs2$label, SLabs2$label_color, ncol=1) #title="Subcortical structures"
)

SLabs2 <- data.frame(
  group_color=Labs[order(Labs$idx2)[seq(19)+parc_res],"group_color"], 
  group=Labs[order(Labs$idx2)[seq(19)+parc_res],"group"], 
  group2=Labs[order(Labs$idx2)[seq(19)+parc_res],"group2"], 
  stringsAsFactors = FALSE
)
SLabs2$group <- paste0(SLabs2$group, " (", SLabs2$group2, ")")
SLabs2 <- unique(SLabs2)
print(
  qual_cleg(SLabs2$group, SLabs2$group_color, ncol=1) #title="Subcortical structures"
)

dev.off()

rm(SLabs2)


## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------
## pdf(plt_fmt("ColorStrip_All.pdf"))
## ggplot() +
##   geom_tile(aes(x=x, y=0, width=1, height=1), data=data.frame(x=seq(nrow(Labs))), fill=Labs$color[order(Labs$idx2)]) +
##   geom_tile(aes(x=x, y=0, width=1, height=1), data=data.frame(x=SLabs$idx), fill=SLabs$color2[order(SLabs$idx2)]) +
##   theme(legend.position="none") +
##   coord_cartesian(expand=FALSE) + theme_void()
## dev.off()
## 
## pdf(plt_fmt("ColorStrip_All2.pdf"))
## ggplot() +
##   geom_tile(aes(x=x, y=0, width=1, height=1), data=data.frame(x=seq(nrow(Labs))), fill=Labs$color[order(Labs$idx2)]) +
##   geom_tile(aes(x=x, y=.2, width=1, height=.6), data=data.frame(x=SLabs$idx), fill=SLabs$color2[order(SLabs$idx2)]) +
##   theme(legend.position="none") +
##   coord_cartesian(expand=FALSE) + theme_void()
## dev.off()

