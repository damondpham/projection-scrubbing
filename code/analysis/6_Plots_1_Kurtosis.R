## --------------------------------------------------------------------------------------------------------------------------------------------------------------
library(splines) # bs in est_trend
library(MASS) # mvrnorm
library(e1071) # kurtosis

color_phi <- c("#000000", "#fd68ad", "#95db3a", "#257cd1", "#491ebd") #1e88bd
color_th <- c("#ebe5d4", "#a39d8d")

## --------------------------------------------------------------------------------------------------------------------------------------------------------------
kurt_sample_var <- function(n){
  return( (24*n*(n-1)^2) / ((n-3)*(n-2)*(n+3)*(n+5)) )
}

# Returns the trend estimated in `ts`.
est_trend <- function(ts, n_knots=5){
	t = length(ts)
	i = 1:t

	loc_knots <- function(t, n_knots){
		if(n_knots >= t){stop('n_knots >= length(ts)')}
		return( t/(n_knots+1)*c(1:(n_knots)) )
	}

	# Do not use outliers (MAD > 3).
	out = abs(ts/stats::mad(ts)) > 3
	i[out] = NA

	# Use cubic spline with five knots.
	model = lm(ts~bs(i, knots=loc_knots(t, n_knots),
					 degree=3, Boundary.knots=c(1, t))
	)

	est = predict(model, newdata=data.frame(i=1:t))

	return(est)
}

set.seed(0)
ts1 <- rnorm(30)
ts1[c(1:4)*4] <- ts1[c(1:4)*4] + 1
ts1[24:25] <- ts1[24:25] + 3.5
ts1[3:18] <- ts1[c(11:18,7:10,3:6)]

ts2 <- ts1[order(ts1)]
set.seed(0)
for(k in 1:20){
  ts2[k:(k+4)] <- sample(ts2[k:(k+4)])
}
ts2[c(11:21)] <- ts2[c(16:21, 11:15)]

ts3 <- ts1[order(ts1)]
ts3[20:30] <- ts3[c(25,27,26,30,29,28,24,23, 20:22)]

df <- data.frame(
  idx = c(1:30),
  ts1 = ts1,
  ts2 = ts2,
  ts3 = ts3,
  ts1.t = est_trend(ts1),
  ts2.t = est_trend(ts2),
  ts3.t = est_trend(ts3),
  ts1.dt = ts1 - est_trend(ts1),
  ts2.dt = ts2 - est_trend(ts2),
  ts3.dt = ts3 - est_trend(ts3)
)

df2 <- reshape::melt(df, id.vars='idx')
df2$Timecourse <- as.numeric(
  gsub('\\..*', '', gsub('ts', '', df2$variable)))
df2$what <- 'Data'
df2$what[grep('\\.t', df2$variable)] <- 'Trend'
df2$what[grep('\\.dt', df2$variable)] <- 'Detrended'
df2 <- df2[c('idx', 'value', 'Timecourse', 'what')]

df2$Timecourse <- paste0('Timecourse #', df2$Timecourse)
p1 <- ggplot(mapping=aes(x=idx,y=value)) + 
  geom_line(data=df2[df2$what=='Trend',], aes(color='Trendline')) + 
  geom_point(data=df2[df2$what=='Data',]) + 
  facet_wrap(Timecourse ~ ., scales='free', 
             dir='h', strip.position = 'top') +
  scale_color_manual(values=c('blue')) +
  labs(colour='') + 
  theme_cowplot() +
  theme(legend.position = 'bottom') + 
  xlab('Timepoint') + ylab('Value') + ylim(c(min(df$ts1), max(df$ts1)))

p1 <- p1 + geom_text(
  data=data.frame(
    label=paste0('Kurtosis: ', round(by(df2[df2$what=='Data','value'], 
                                    df2[df2$what=='Data','Timecourse'], kurtosis, type=1), 3)),
    Timecourse=unique(df2$Timecourse)),
  mapping = aes(x = 10, y = 3.5, label = label)
)

pleg <- cowplot::get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

df2$Timecourse <- gsub('Timecourse', 'Detrended', df2$Timecourse)
p2 <- ggplot(mapping=aes(x=idx,y=value)) + 
  geom_hline(yintercept=0, color='blue') + 
  geom_point(data=df2[df2$what=='Detrended',]) + 
  facet_wrap(Timecourse ~ ., scales='free', 
             dir='h', strip.position = 'top') +
  scale_color_manual(values=c('blue')) +
  theme_cowplot() +
  labs(colour='Trendline') + 
  xlab('Timepoint') + ylab('Value') + ylim(c(min(df$ts1), max(df$ts1)))

p2 <- p2 + geom_text(
  data=data.frame(
    label=paste0('Kurtosis: ', round(by(df2[df2$what=='Detrended','value'], 
                                    df2[df2$what=='Detrended','Timecourse'], kurtosis, type=1), 3)),
    Timecourse=unique(df2$Timecourse)),
  mapping = aes(x = 10, y = 3.5, label = label)
)

pdf(plt_fmt("kurt_detrending.pdf"), width=8, height=2.5)
print(p1)
print(p2)
plot_grid(NULL, pleg, ncol=1)
dev.off()

rm(df, df2, p1, p2, pleg, k, ts1, ts2, ts3, est_trend)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
n_rep <- 100000
t_levels <- c(100, 200, 500, 1000)
t_max <- max(t_levels)
n_t <- length(t_levels)
phi_levels <- c(0, .2, .4, .6, .8)
names(color_phi) <- as.character(phi_levels)
n_phi <- length(phi_levels)

x_ar1_fname <- file.path(dir_plots, "helper", "x_ar1.rds")
if (!file.exists(x_ar1_fname)) {
  x_ar1 <- vector("list", length=n_phi)
  names(x_ar1) <- paste("phi", phi_levels)
  for (ii in seq(n_phi)) {
    cat(phi_levels[ii], "\n")
    x_ar1[[ii]] <- t(mvrnorm(
      n_rep, mu=rep(0, t_max), 
      Sigma=toeplitz(ARMAacf(ar=phi_levels[ii], lag.max=t_max-1))
    ))
  } 
  saveRDS(x_ar1, x_ar1_fname)
}

x_kurt_fname <- file.path(dir_plots, "helper", "x_kurt.rds")
if (!file.exists(x_kurt_fname)) {
  x_ar1 <- readRDS(x_ar1_fname)
  x_kurt <- vector("list", length=n_t)
  names(x_kurt) <- paste("T =", t_levels)
  for (ii in seq(n_t)) {
    cat(t_levels[ii], "\n")
    x_kurt[[ii]] <- as.data.frame(lapply(
      x_ar1, function(m){apply(m[seq(t_levels[ii]),], 2, kurtosis, type=1)}
    ))
  }
  saveRDS(x_kurt, x_kurt_fname)
  rm(x_ar1)
}
x_kurt <- readRDS(x_kurt_fname)

get_hist <- function(d, breaks=100){
  x = hist(d, breaks=breaks, plot=FALSE)
  n = length(x$breaks)
  b = (x$breaks[c(2:n)] + x$breaks[c(1:(n-1))]) / 2
  out = data.frame(breaks=b, density=x$density)
  return(out)
}

x_kurt2 <- x_kurt
for (ii in seq(n_t)) {
  cat(t_levels[ii], "\n")
  x_kurt2[[ii]] <- lapply(x_kurt[[ii]], get_hist)
  for (jj in seq(n_phi)) {
    x_kurt2[[ii]][[jj]]$AR <- phi_levels[jj]
    x_kurt2[[ii]][[jj]]$type <- "empirical"
  }
  x_kurt2[[ii]][[n_phi+1]] <- data.frame(
    breaks = x_kurt2[[ii]][["phi.0"]]$breaks,
    density = dnorm(x_kurt2[[ii]][["phi.0"]]$breaks, sd=sqrt(kurt_sample_var(t_levels[ii]))),
    AR = 0,
    type="theoretical"
  )
  x_kurt2[[ii]] <- do.call(rbind, x_kurt2[[ii]])
  x_kurt2[[ii]]$t <- names(x_kurt2)[ii]
}
x_kurt2 <- do.call(rbind, x_kurt2)
x_kurt2$AR <- factor(x_kurt2$AR, levels=as.character(phi_levels))
x_kurt2$type <- factor(x_kurt2$type)
x_kurt2$t <- factor(x_kurt2$t, levels=paste("T =", t_levels))
colnames(x_kurt2)[colnames(x_kurt2)=="breaks"] <- "value"

x_99q <- data.frame(sapply(x_kurt, apply, 2, quantile, .99))
colnames(x_99q) <- paste("T =", t_levels)
x_99q <- rbind(x_99q, qnorm(0.99) * sqrt(kurt_sample_var(t_levels)))
x_99q <- stack(x_99q)
colnames(x_99q) <- c("value", "t")
x_99q$AR <- factor(rep(c(phi_levels, 0), n_t), levels=as.character(phi_levels))
x_99q$type <- factor(rep(c(rep("empirical", n_phi), "theoretical"), n_t))

pdf(plt_fmt("kurt_ar_ntime.pdf"), width=6.7, height=8)

p1 <- ggplot(x_kurt2, aes(x=value)) +
  # Densities
  geom_area(aes(y=density), fill=color_th[1], data=subset(x_kurt2, type=="theoretical")) +
  geom_line(aes(y=density, color=AR), data=subset(x_kurt2, type=="empirical" & AR=="0"), size=1.5*.8, position="identity") +
  geom_line(aes(y=density, color=AR), data=subset(x_kurt2, type=="empirical" & AR!="0"), size=.8*.8, position="identity") +
  # Quantiles
  geom_vline(aes(xintercept=value, color=AR), data=subset(x_99q, type=="empirical" & AR=="0"), size=1.5, show.legend=FALSE) +
  geom_vline(aes(xintercept=value), color=color_th[2], data=subset(x_99q, type=="theoretical"), size=1.5, show.legend=FALSE) +
  geom_vline(aes(xintercept=value, color=AR), data=subset(x_99q, type=="empirical" & AR!="0"), size=.8, linetype="dashed", show.legend=FALSE) +
  # Colors
  scale_color_manual(values=color_phi, name="AR(1) coef.") +
  # Formatting
  facet_wrap(t ~ ., scales='free', dir='h', strip.position = 'top', ncol=1) +
  theme_cowplot() +
  theme(
    legend.position = 'bottom', 
    strip.background = element_rect(fill="#ffffff"), strip.text=element_text(hjust=0, face="bold"),
    panel.spacing = unit(1.15, "lines")
  ) +
  xlab('Sample kurtosis') + ylab('Density') +
  labs(sub='Densities') +
  #coord_fixed() +
  coord_cartesian(xlim=c(-1, 1.5)) +
  scale_x_continuous(breaks=c(-1, 0, 1.5)) +
  scale_y_continuous(expand=expansion(add=c(0,.2))) +
  #ggtitle('Sampling distributions of kurtosis for AR(1) data') +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
print(p1)
dev.off()

rm(
  n_rep, t_levels, t_max, n_t, phi_levels, n_phi, 
  x_99q, x_kurt, x_kurt2, x_ar1_fname, x_kurt_fname
)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
lev_examples_fnames <- file.path(dir_scrubMeas, "Lev/CC2MP6", c(
  "LEV_103818_v1_LR.rds", "LEV_139839_v1_LR.rds", "LEV_662551_v2_RL.rds"
))

pdf(plt_fmt("kurtosis_examples.pdf"), width=8.5, height=3)

for (ll in seq(length(lev_examples_fnames))) {
  lev <- readRDS(lev_examples_fnames[ll])

  hk <- "highkurt"
  comps <- "M"
  
  stopifnot(sum(lev$ICA[[hk]]) >= 3)
  ICs_lowkurt <- which(!lev$ICA[[hk]])[seq(3)]
  ICs_highkurt <- rev(order(apply(lev$ICA[[comps]], 2, e1071::kurtosis, type=1)))[seq(3)]
  
  df_k <- data.frame(
    Value = c(
      as.vector(apply(lev$ICA[[comps]][,ICs_lowkurt], 2, scale)), 
      as.vector(apply(lev$ICA[[comps]][,ICs_highkurt], 2, scale))
    ), 
    col = factor(paste(c("First", "Second", "Third")[rep(rep(seq(3), each=hcp_T-nDrop), 2)], "IC")),
    Timepoint = rep(rep(seq(hcp_T-nDrop), 3), 2),
    kurt = factor(rep(c("Low kurtosis", "High kurtosis"), each=(hcp_T-nDrop)*3))
  )
  
  df_kv <- rbind(
    data.frame(
      KurtVal=apply(lev$ICA[[comps]][,ICs_lowkurt], 2, e1071::kurtosis, type=1), 
      col=factor(paste(c("First", "Second", "Third"), "IC")), kurt="Low kurtosis"
    ),
    data.frame(
      KurtVal=apply(lev$ICA[[comps]][,ICs_highkurt], 2, e1071::kurtosis, type=1), 
      col=factor(paste(c("First", "Second", "Third"), "IC")), kurt="High kurtosis"
    )
  )
  df_kv$KurtVal <- paste("Kurtosis:", round(df_kv$KurtVal, 2))
  
  plt_k <- ggplot(subset(df_k, kurt=="High kurtosis")) + 
    #geom_vline(data=df_kl, aes(xintercept=Timepoint), linetype="dotted", color="#3366FF") +
    geom_text(data=subset(df_kv, kurt=="High kurtosis"), aes(x=565, y=12, label=KurtVal), color="blue", hjust=0) +
    geom_line(aes(x=Timepoint, y=Value)) + facet_grid(.~col) + 
    cowplot::theme_minimal_grid() + 
    cowplot::panel_border(size=3) + 
    scale_x_continuous(breaks=c(0, 300, 600, 900, 1185)) +
    coord_cartesian(xlim=c(0,1200), ylim=c(-15, 15), expand=FALSE) + 
    theme(panel.spacing.x = unit(1.7, "lines"), panel.spacing.y = unit(0, "lines"), strip.background.x = element_blank(), strip.text.x = element_blank(), plot.margin = grid::unit(c(2,8,2,2), "mm")) +
    xlab("Timepoint") + ylab("Value") + ggtitle(paste0("High-kurtosis ICs"))
  
  print(plt_k)
  
  plt_k <- ggplot(subset(df_k, kurt=="Low kurtosis")) + 
    #geom_vline(data=df_kl, aes(xintercept=Timepoint), linetype="dotted", color="#3366FF") +
    geom_text(data=subset(df_kv, kurt=="Low kurtosis"), aes(x=565, y=12, label=KurtVal), color="blue", hjust=0) +
    geom_line(aes(x=Timepoint, y=Value)) + facet_grid(.~col) + 
    cowplot::theme_minimal_grid() + 
    cowplot::panel_border(size=3) + 
    scale_x_continuous(breaks=c(0, 300, 600, 900, 1185)) +
    coord_cartesian(xlim=c(0,1200), ylim=c(-15, 15), expand=FALSE) + 
    theme(panel.spacing.x = unit(1.7, "lines"), panel.spacing.y = unit(0, "lines"), strip.background.x = element_blank(), strip.text.x = element_blank(), plot.margin = grid::unit(c(2,8,2,2), "mm")) +
    xlab("Timepoint") + ylab("Value") + ggtitle(paste0("Low-kurtosis ICs"))
  
  print(plt_k)
}

dev.off()

rm(
  df_k, df_kv, hk, comps, plt_k, ICs_highkurt, ICs_lowkurt,
  lev_examples_fnames
)

