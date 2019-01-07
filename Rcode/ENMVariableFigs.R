### LOAD LIBRARIES ---------------------------
library(tidyverse)
library(cowplot)

### LOAD AND PREP DATA FRAMES -------------------------
climvars <- read_csv("Data/climate_enm_variables.csv")
streamvars <- read_csv("Data/stream_enm_variables.csv")
lams <- read_csv("Robjects/site.lambdas.bootstrap.csv")

dat <- left_join(lams,climvars, by=c("Site"="site", "ID1"="ID1", "lat"="lat", "long"="long", "el"="elev_m", "region"="region"))
dat <- left_join(dat,streamvars,by=c("Site"="site", "ID1"="ID1", "lat"="lat", "long"="long", "el"="elev_m", "region"="region", "site_label"="site_label"))

# variables used in climate ENM: precip seasonality (bio15), log precip driest month (log bio14), log annual precip (log bio12), mean temp coldest quarter (bio11), log mean temp warmest quarter (log bio10), temperature seasonality (bio4), log isothermality (log bio3), mean diurnal range (bio2); these are already log-transformed as needed in input file

# variables used in stream ENM: log annual discharge (log bio12), discharge seasonality (bio15_stream), slope of stream reach, log drainage area, topographic roughness

dat <- dat %>% 
  select(X1, Site, region, lat, long, lambda, upper, lower, bio15_clim, bio14_clim, bio12_clim, bio11_clim, bio10_clim, bio4_clim, bio3_clim, bio2_clim, logbio12_stream, bio15_stream, SLOPE_stream, logDrainAre_stream, terrough20C_stream)

### UNIVARIATE PLOTS OF LAM ~ ENM VARS -------------------------
a <- ggplot(dat, aes(bio15_clim, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Precip. seasonality") + #clim bio15
  xlim(53,80) +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,0,80), "pt"))

b <- ggplot(dat, aes(bio14_clim, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Precip. driest month") + #clim bio14
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,0,0), "pt"))

c <- ggplot(dat, aes(bio12_clim, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Annual precip.") + #clim bio12
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,0,0), "pt"))

d <- ggplot(dat, aes(bio11_clim, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Temp. coldest quarter") + #clim bio11
  xlim(4.3,6) +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,0,0), "pt"))

e <- ggplot(dat, aes(bio10_clim, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Temp. warmest quarter") + #clim bio10
  xlim(17.4, 19.2) +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,0,80), "pt"))

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
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,0,0), "pt"))

g <- ggplot(dat, aes(bio3_clim, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Isothermality") + #clim bio3
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,0,0), "pt"))

h <- ggplot(dat, aes(logbio12_stream, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Annual discharge") + #stream bio12
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,0,0), "pt"))

i <- ggplot(dat, aes(bio15_stream, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Discharge seasonality") + #stream bio15
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,80,80), "pt"))

j <- ggplot(dat, aes(SLOPE_stream, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Stream slope") + 
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,80,0), "pt"))

k <- ggplot(dat, aes(logDrainAre_stream, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Drainage area") + 
  xlim(0.44,0.68) +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,80,0), "pt"))

l <- ggplot(dat, aes(terrough20C_stream, lambda)) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_point(aes(fill=region), size=6, shape=21) +
  scale_fill_manual(values=c("grey", "white")) +
  geom_text(aes(label=X1), size=4) + 
  xlab("Drainage roughness") + 
  #xlim() +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,80,0), "pt"))

lamclim <- plot_grid(a + theme(legend.position = "none"), 
                     b + theme(legend.position = "none"), 
                     c + theme(legend.position = "none"), 
                     d + theme(legend.position = "none"), 
                     e + theme(legend.position = "none"), 
                     f + theme(legend.position = "none"), 
                     g + theme(legend.position = "none"),
                     h + theme(legend.position = "none"),
                     i + theme(legend.position = "none"),
                     j + theme(legend.position = "none"),
                     k + theme(legend.position = "none"),
                     l + theme(legend.position = "none"),
                     nrow=3, labels="AUTO", label_x=0.9)
legend <- get_legend(a)
lamclim2 <- plot_grid(lamclim, legend, rel_widths = c(5, 0.5))
left_label <- expression(paste("Population growth rate (", lambda, ")"))
lamclim3 <- ggdraw(lamclim2) + draw_label(left_label, angle=90, x=0.05, size=24)
bottom_label <- "ENM variables"
lamclim4 <- ggdraw(lamclim3) + draw_label(bottom_label, angle=0, y=0.05, size=24)
save_plot("Figures/Lambda_vs_ENMVars2.png", lamclim4, base_width=22, base_height=11)

### SITES IN MULTIVARIATE ENVIRO SPACE -------------------------

PC <- prcomp(dat[,c(9:21)], center=TRUE, scale=TRUE)
summary(pcall)
biplot(pcall)

pcclim <- prcomp(dat[,c(9:16)], center=TRUE, scale=TRUE)
summary(pcclim)
biplot(pcclim)

pcstream <- prcomp(dat[,c(17:21)], center=TRUE, scale=TRUE)
summary(pcstream)
biplot(pcstream)


PCbiplot <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=3, aes(label=obsnames))
  #plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  plot
}

fit <- prcomp(USArrests, scale=T)
PCbiplot(fit)

PCbiplot(pcclim)

