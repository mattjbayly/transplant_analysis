### LOAD LIBRARIES ---------------------------
library(tidyverse)
library(cowplot)

### LOAD AND PREP DATA FRAMES -------------------------
climvars_sites <- read_csv("Data/climate_enm_variables.csv")
streamvars_sites <- read_csv("Data/stream_enm_variables.csv")
#climvars_presabs <- read_csv("Data/climate_enm_presabs_records.csv") #this file has bad values for bio12 and bio14 (mostly NAs, bio14 wrong scale)
climvars_presabs <- read_csv("Data/pres.records.aug.31.csv")
streamvars_presabs <- read_csv("Data/stream_enm_presabs_records.csv")
lams <- read_csv("Robjects/site.lambdas.bootstrap.csv")

dat <- left_join(lams,climvars_sites, by=c("Site"="site", "ID1"="ID1", "lat"="lat", "long"="long", "el"="elev_m", "region"="region"))
dat <- left_join(dat,streamvars_sites,by=c("Site"="site", "ID1"="ID1", "lat"="lat", "long"="long", "el"="elev_m", "region"="region", "site_label"="site_label"))

# variables used in climate ENM: precip seasonality (bio15), log precip driest month (log bio14), log annual precip (log bio12), mean temp coldest quarter (bio11), log mean temp warmest quarter (log bio10), temperature seasonality (bio4), log isothermality (log bio3), mean diurnal range (bio2); these are already log-transformed as needed in input file

# variables used in stream ENM: log annual discharge (log bio12), discharge seasonality (bio15_stream), slope of stream reach, log drainage area, topographic roughness

dat <- dat %>% 
  select(X1, Site, region, lat, long, lambda, upper, lower, bio15_clim, bio14_clim, bio12_clim, bio11_clim, bio10_clim, bio4_clim, bio3_clim, bio2_clim, logbio12_stream, bio15_stream, SLOPE_stream, logDrainAre_stream, terrough20C_stream)

dat2 <- climvars_presabs %>% 
  filter(DATASET=="herb") %>% 
  select(Latitude, Longitude, Elevation, bio15, bio14, bio12, bio11, bio10, bio4, bio3, bio2)

dat3 <- streamvars_presabs %>% 
  filter(presabs==1) %>% 
  select(-modeldataset, -presabs, -id, -bio12_stream)


### UNIVARIATE MODELS OF LAM ~ ENM VARS -------------------------
moda <- lm(lambda ~ bio15_clim, data=dat)
summary(moda)

modb <- lm(lambda ~ bio14_clim, data=dat)
summary(modb)

modc <- lm(lambda ~ bio12_clim, data=dat)
summary(modc)

modd <- lm(lambda ~ bio11_clim, data=dat)
summary(modd)

mode <- lm(lambda ~ bio10_clim, data=dat)
summary(mode)

modf <- lm(lambda ~ bio4_clim, data=dat)
summary(modf)

modg <- lm(lambda ~ bio3_clim, data=dat)
summary(modg)

modh <- lm(lambda ~ bio2_clim, data=dat)
summary(modh)

modi <- lm(lambda ~ logbio12_stream, data=dat)
summary(modi)

modj <- lm(lambda ~ bio15_stream, data=dat)
summary(modj)

modk <- lm(lambda ~ SLOPE_stream, data=dat)
summary(modk)

modl <- lm(lambda ~ logDrainAre_stream, data=dat)
summary(modl)

modm <- lm(lambda ~ terrough20C_stream, data=dat)
summary(modm)

### UNIVARIATE PLOTS OF LAM ~ ENM VARS -------------------------
a <- ggplot(dat, aes(bio15_clim, lambda)) +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Precip. seasonality") + #clim bio15
  xlim(53,80) +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,60,70), "pt"))

b <- ggplot(dat, aes(bio14_clim, lambda)) +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Precip. driest month") + #clim bio14
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,60,70), "pt"))

c <- ggplot(dat, aes(bio12_clim, lambda)) +
  geom_smooth(method=lm, se=FALSE, color="black", linetype="dashed") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Annual precip.") + #clim bio12
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,60,70), "pt"))

d <- ggplot(dat, aes(bio11_clim, lambda)) +
  geom_smooth(method=lm, se=FALSE, color="black", linetype="dashed") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Temp. coldest quarter") + #clim bio11
  xlim(4.3,6) +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,60,70), "pt"))

e <- ggplot(dat, aes(bio10_clim, lambda)) +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Temp. warmest quarter") + #clim bio10
  xlim(17.4, 19.2) +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,60,70), "pt"))

f <- ggplot(dat, aes(bio4_clim, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Temp. seasonality") + #clim bio4
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,60,70), "pt"))

g <- ggplot(dat, aes(bio3_clim, lambda)) +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Isothermality") + #clim bio3
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,60,70), "pt"))

h <- ggplot(dat, aes(bio2_clim, lambda)) +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Diurnal Range") + #clim bio2
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,60,70), "pt"))
  
i <- ggplot(dat, aes(logbio12_stream, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Annual discharge") + #stream bio12
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,60,70), "pt"))

j <- ggplot(dat, aes(bio15_stream, lambda)) +
  geom_smooth(method=lm, se=FALSE, color="black", linetype="dashed") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Discharge seasonality") + #stream bio15
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,60,70), "pt"))

k <- ggplot(dat, aes(SLOPE_stream, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Stream slope") + 
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,60,70), "pt"))

l <- ggplot(dat, aes(logDrainAre_stream, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Drainage area") + 
  xlim(0.44,0.68) +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,60,70), "pt"))

m <- ggplot(dat, aes(terrough20C_stream, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Drainage roughness") + 
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,60,70), "pt"))

lamclim <- plot_grid(a + theme(legend.position = "none"), 
                     b + theme(legend.position = "none"), 
                     c + theme(legend.position = "none"), 
                     d + theme(legend.position = "none"), 
                     e + theme(legend.position = "none"), 
                     f + theme(legend.position = "none"), 
                     g + theme(legend.position = "none"),
                     h + theme(legend.position = "none"),
                     nrow=2, labels="AUTO", label_x=0.9)

lamstream <- plot_grid(i + theme(legend.position = "none"),
                     j + theme(legend.position = "none"),
                     k + theme(legend.position = "none"),
                     l + theme(legend.position = "none"),
                     m + theme(legend.position = "none"),
                     nrow=2, ncol=3, labels="AUTO", label_x=0.9)

legend <- get_legend(a)
left_label <- expression(paste("Population growth rate (", lambda, ")"))
bottom_label1 <- "Climate ENM variables"
bottom_label2 <- "Stream ENM variables"

lamclim2 <- plot_grid(lamclim, legend, rel_widths = c(5, 0.5))
lamclim3 <- ggdraw(lamclim2) + draw_label(left_label, angle=90, x=0.05, size=24)
lamclim4 <- ggdraw(lamclim3) + draw_label(bottom_label1, angle=0, y=0.05, size=24)
save_plot("Figures/Lambda_vs_ClimENMVars.png", lamclim4, base_width=22, base_height=11)

lamstream2 <- plot_grid(lamstream, legend, rel_widths = c(5, 0.5))
lamstream3 <- ggdraw(lamstream2) + draw_label(left_label, angle=90, x=0.05, size=24)
lamstream4 <- ggdraw(lamstream3) + draw_label(bottom_label2, angle=0, y=0.05, size=24)
save_plot("Figures/Lambda_vs_StreamENMVars.png", lamstream4, base_width=22, base_height=11)


### SITES IN MULTIVARIATE ENVIRO SPACE -------------------------

##PCA with transplant sites only

clim.dat <- dat %>% 
  select(Pseas=bio15_clim, 
         Pdri=bio14_clim, 
         Pann=bio12_clim, 
         Tcold=bio11_clim, 
         Twarm=bio10_clim, 
         Tseas=bio4_clim, 
         Tiso=bio3_clim,
         Tdiurn=bio2_clim) 

stream.dat <- dat %>% 
  select(Dann=logbio12_stream, 
         Dseas=bio15_stream, 
         Wslope=SLOPE_stream, 
         Warea=logDrainAre_stream, 
         Wrough=terrough20C_stream) 

all.dat <- cbind(stream.dat, clim.dat)

PCbiplot <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  region <- c(rep("in",4),rep("out",4))
  data <- cbind(data.frame(obsnames=row.names(as.data.frame(PC$x)), PC$x),region)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- ggplot(data, aes_string(x=x, y=y)) + 
    geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, hjust="outward", vjust="outward", color="grey") + 
    geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="grey") + 
    geom_point(aes(fill=region), size=6, shape=21) + 
    scale_fill_manual(values=c("grey", "white")) +
    geom_text(size=6, aes(label=obsnames))
  theme_classic() + 
    theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,0,0), "pt"))
  plot
}

pcclim <- prcomp(clim.dat, center=TRUE, scale=TRUE)
#summary(pcclim)
#biplot(pcclim)
plotclim <- PCbiplot(pcclim)

pcstream <- prcomp(stream.dat, center=TRUE, scale=TRUE)
#summary(pcstream)
#biplot(pcstream)
plotstream <- PCbiplot(pcstream)

pcall <- prcomp(all.dat, center=TRUE, scale=TRUE)
#summary(pcall)
#biplot(pcall)
plotall <- PCbiplot(pcall)

pc_plots <- plot_grid(plotclim + theme(legend.position = "none"),
                      plotstream + theme(legend.position = "none"),
                      plotall + theme(legend.position = "none"),
                      nrow=3, labels="AUTO", label_x=0.9)
legend <- get_legend(plotclim)
pc_plots2 <- plot_grid(legend, pc_plots, nrow=2, rel_heights = c(1,5,5))
save_plot("Figures/PCA_ENMVars.png", pc_plots2, base_width=8, base_height=11)


## PCA with presence records

clim.pres <- dat2 %>% 
  select(-Latitude, -Longitude, -Elevation, 
         ppt_seasonal=bio15, 
         ppt_drimonth=bio14, 
         ann_ppt=bio12, 
         temp_coldquart=bio11, 
         temp_warmquart=bio10, 
         temp_seasonal=bio4, 
         isotherm=bio3,
         diurn_range=bio2) 

clim.all <- bind_rows(clim.dat, clim.pres)

PC <- prcomp(clim.all, center=TRUE, scale=TRUE)
#summary(pcclim2)
#biplot(pcclim2)
plotclim2 <- PCbiplot(pcclim2)

PCbiplot.climall <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  region <- c(rep("in",4),rep("out",4),rep(NA,431))
  data <- cbind(data.frame(obsnames=row.names(as.data.frame(PC$x)), PC$x),region)
  data.sites <- data %>% filter(region=="in" | region=="out")
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- ggplot(data, aes_string(x=x, y=y)) + 
    geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, hjust="outward", vjust="outward", color="black") +  
    geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="grey") + 
    geom_point(color="grey") + 
    geom_point(data=data.sites, aes(fill=region), size=6, shape=21) + 
    scale_fill_manual(values=c("grey", "white")) +
    geom_text(data=data.sites, size=6, aes(label=obsnames)) + 
    xlim(-4.5,4.5) + 
    ylim(-4.5,4.5) +
    theme_classic() + 
    theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,0,0), "pt"))
  plot
}
plotclim2 <- PCbiplot.climall(pcclim2)


stream.pres <- dat3 %>% 
  select(-lat, -long, -elev_m,
         ann_discharge=logbio12_stream, 
         discharge_seasonal=bio15_stream, 
         slope=SLOPE_stream, 
         drain_area=logDrainAre_stream, 
         roughness=terrough20C_stream) %>% 
  na.omit()

stream.all <- bind_rows(stream.dat, stream.pres)

pcstream2 <- prcomp(stream.all, center=TRUE, scale=TRUE)
#summary(pcstream2)
#biplot(pcstream2)

PCbiplot.streamall <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  region <- c(rep("in",4),rep("out",4),rep(NA,243))
  data <- cbind(data.frame(obsnames=row.names(as.data.frame(PC$x)), PC$x),region)
  data.sites <- data %>% filter(region=="in" | region=="out")
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- ggplot(data, aes_string(x=x, y=y)) + 
    geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust="outward", hjust="outward", color="grey") + 
    geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="grey") + 
    geom_point(color="grey") + 
    geom_point(data=data.sites, aes(fill=region), size=6, shape=21) + 
    scale_fill_manual(values=c("grey", "white")) +
    geom_text(data=data.sites, size=6, aes(label=obsnames)) + 
    xlim(-4,6) + 
    ylim(-3,3) +
  theme_classic() + 
    theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,0,0), "pt"))
  plot
}
plotstream2 <- PCbiplot.streamall(pcstream2)

pc_plots_pres <- plot_grid(plotclim2 + theme(legend.position = "none"),
                      plotstream2 + theme(legend.position = "none"),
                      nrow=1, labels="AUTO", label_x=0.9)
legend <- get_legend(plotclim2)
pc_plots_pres2 <- plot_grid( legend, pc_plots_pres, nrow=2, rel_heights = c(1,5))
save_plot("Figures/PCA_ENMVars_Presences.png", pc_plots_pres2, base_width=8, base_height=5)
