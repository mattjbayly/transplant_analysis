library(tidyverse)
library(MuMIn)
library(cowplot)

# lambdas
lams <- read_csv("Robjects/site.lambdas.bootstrap.csv")
head(lams)

# suitability (redone by Amy March 2018, using final Am Nat objects)
# ENM models applied to long-term average climate of the sites
# climate values from sites from single spot
suit <- read_csv("Data/site_preds_average.csv") %>% 
  filter(ID1=="trans") %>% droplevels() %>% 
  dplyr::select(ID2, LRavg, GAMavg, RFavg, BRTavg, MAXavg) %>% 
  mutate(Site = if_else(ID2=="Rock Creek", "ROCK", ifelse(ID2=="Coast Fork Willamette", "COAST", ifelse(ID2=="Mosby Creek", "MOSBY", ifelse(ID2=="Calapooia Creek", "CALAPOOIA", ifelse(ID2=="Wiley Creek", "WILEY", ifelse(ID2=="Thomas", "THOMAS", ifelse(ID2=="Hunter", "HUNTER", ifelse(ID2=="Looking Glass", "LOOK", NA)))))))))

dat <- left_join(lams, suit)
head(dat)

for (i in 1:dim(dat)[1]) {
  dat$Ens[i] = mean(c(dat$LRavg[i], dat$GAMavg[i], dat$RFavg[i], dat$BRTavg[i], dat$MAXavg[i]))
}

write_csv(dat, "Robjects/site.lambdas.suitability.csv")

# matthew's ENM predictions (redone March 2018, using final Am Nat objects)
# ENM models applied to experimental conditions during transplant
# climate values from zonal average of all plots at a site
matt.preds <- read_csv("data/Site Level ENM Preds.csv")

dat <- left_join(dat, matt.preds)
check.merge <- dat %>% dplyr::select(Site, site) # good

write_csv(dat, "Robjects/site.lambdas.suitability.csv")

# add SD among models (for each ensemble)
dat <- dat %>% 
  group_by(Site) %>%
  mutate(clim_amy_sd = sd(c(LRavg,GAMavg,RFavg,BRTavg,MAXavg)),
         clim_matt_sd = sd(c(glm_climate,gam_climate,rf_climate,brt_climate,max_climate)),
         stream_sd = sd(c(glm_stream,gam_stream,rf_stream,brt_stream,max_stream)),
         clim_amy_se = sd(c(LRavg,GAMavg,RFavg,BRTavg,MAXavg))/sqrt(5),
         clim_matt_se = sd(c(glm_climate,gam_climate,rf_climate,brt_climate,max_climate))/sqrt(5),
         stream_se = sd(c(glm_stream,gam_stream,rf_stream,brt_stream,max_stream))/sqrt(5)) %>% 
  ungroup()

plot(dat$Ens, dat$mean_pred_climate) #WTF
plot(dat$LRavg, dat$glm_climate) #pretty good
plot(dat$GAMavg, dat$gam_climate) #one bad outlier
plot(dat$RFavg, dat$rf_climate) #not good
plot(dat$BRTavg, dat$brt_climate) #terrible
plot(dat$MAXavg, dat$max_climate) #good

plot_ens <- ggplot(dat, aes(Ens, mean_pred_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Amy climate ENS") + 
  ylab("Matt climate ENS") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right") 

plot_glm <- ggplot(dat, aes(LRavg, glm_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Amy climate GLM") + 
  ylab("Matt climate GLM") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right") 

plot_gam <- ggplot(dat, aes(GAMavg, gam_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Amy climate GAM") + 
  ylab("Matt climate GAM") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right") 

plot_rf <- ggplot(dat, aes(RFavg, rf_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Amy climate RF") + 
  ylab("Matt climate RF") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right") 

plot_brt <- ggplot(dat, aes(BRTavg, brt_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Amy climate BRT") + 
  ylab("Matt climate BRT") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right") 

plot_max <- ggplot(dat, aes(MAXavg, max_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Amy climate MAX") + 
  ylab("Matt climate MAX") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none") 

model_comp <- plot_grid(plot_ens + theme(legend.position = "none"), 
                        plot_glm + theme(legend.position = "none"), 
                        plot_gam + theme(legend.position = "none"), 
                        plot_rf + theme(legend.position = "none"), 
                        plot_brt + theme(legend.position = "none"), 
                        plot_max + theme(legend.position = "none"))
legend <- get_legend(plot_ens)
model_comp2 <- plot_grid(model_comp, legend)
save_plot("Figures/ClimateENMComparions1.png", model_comp, base_width=11, base_height=8.5)

plot(dat$mean_pred_climate, dat$mean_pred_stream)

### linear models of suitability vs latitude
## ave suitability
# Amy's values, climate long-term
mod0.ens <- lm(Ens ~ 1, data=dat)
mod1.ens <- lm(Ens ~ lat, data=dat)
model.sel(mod0.ens, mod1.ens)
summary(mod1.ens)

# Matt's values, climate short-term
mod0.ens.matt <- lm(mean_pred_climate ~ 1, data=dat)
mod1.ens.matt <- lm(mean_pred_climate ~ lat, data=dat)
model.sel(mod0.ens.matt, mod1.ens.matt)
summary(mod1.ens.matt)

# Matt's values, stream
mod0.ens.stream <- lm(mean_pred_stream ~ 1, data=dat)
mod1.ens.stream <- lm(mean_pred_stream ~ lat, data=dat)
model.sel(mod0.ens.stream, mod1.ens.stream)
summary(mod1.ens.stream)

## individual ENM algorithms
# Amy's values, climate long-term
mod0.lr <- lm(LRavg ~ 1, data=dat)
mod1.lr <- lm(LRavg ~ lat, data=dat)
model.sel(mod0.lr, mod1.lr)
summary(mod1.lr)

mod0.gam <- lm(GAMavg ~ 1, data=dat)
mod1.gam <- lm(GAMavg ~ lat, data=dat)
model.sel(mod0.gam, mod1.gam)
summary(mod1.gam)

mod0.rf <- lm(RFavg ~ 1, data=dat)
mod1.rf <- lm(RFavg ~ lat, data=dat)
model.sel(mod0.rf, mod1.rf)
summary(mod1.rf)

mod0.brt <- lm(BRTavg ~ 1, data=dat)
mod1.brt <- lm(BRTavg ~ lat, data=dat)
model.sel(mod0.brt, mod1.brt)
summary(mod1.brt)

mod0.max <- lm(MAXavg ~ 1, data=dat)
mod1.max <- lm(MAXavg ~ lat, data=dat)
model.sel(mod0.max, mod1.max)
summary(mod1.max)

# Matt's values, climate short-term
mod0.lr.matt <- lm(glm_climate ~ 1, data=dat)
mod1.lr.matt <- lm(glm_climate ~ lat, data=dat)
model.sel(mod0.lr.matt, mod1.lr.matt)
summary(mod1.lr.matt)

mod0.gam.matt <- lm(gam_climate ~ 1, data=dat)
mod1.gam.matt <- lm(gam_climate ~ lat, data=dat)
model.sel(mod0.gam.matt, mod1.gam.matt)
summary(mod1.gam.matt)

mod0.rf.matt <- lm(rf_climate ~ 1, data=dat)
mod1.rf.matt <- lm(rf_climate ~ lat, data=dat)
model.sel(mod0.rf.matt, mod1.rf.matt)
summary(mod1.rf.matt)

mod0.brt.matt <- lm(brt_climate ~ 1, data=dat)
mod1.brt.matt <- lm(brt_climate ~ lat, data=dat)
model.sel(mod0.brt, mod1.brt.matt)
summary(mod1.brt.matt)

mod0.max.matt <- lm(max_climate ~ 1, data=dat)
mod1.max.matt <- lm(max_climate ~ lat, data=dat)
model.sel(mod0.max.matt, mod1.max.matt)
summary(mod1.max.matt)

# Matt's values, stream
mod0.lr.stream <- lm(glm_stream ~ 1, data=dat)
mod1.lr.stream <- lm(glm_stream ~ lat, data=dat)
model.sel(mod0.lr.stream, mod1.lr.stream)
summary(mod1.lr.stream)

mod0.gam.stream <- lm(gam_stream ~ 1, data=dat)
mod1.gam.stream <- lm(gam_stream ~ lat, data=dat)
model.sel(mod0.gam.stream, mod1.gam.stream)
summary(mod1.gam.stream)

mod0.rf.stream <- lm(rf_stream ~ 1, data=dat)
mod1.rf.stream <- lm(rf_stream ~ lat, data=dat)
model.sel(mod0.rf.stream, mod1.rf.stream)
summary(mod1.rf.stream)

mod0.brt.stream <- lm(brt_stream ~ 1, data=dat)
mod1.brt.stream <- lm(brt_stream ~ lat, data=dat)
model.sel(mod0.brt, mod1.brt.stream)
summary(mod1.brt.stream)

mod0.max.stream <- lm(max_stream ~ 1, data=dat)
mod1.max.stream <- lm(max_stream ~ lat, data=dat)
model.sel(mod0.max.stream, mod1.max.stream)
summary(mod1.max.stream)

### anova of suitability by range position
# long-term average climate, Amy
mod2.ens <- lm(Ens ~ region, data=dat)
summary(mod2.ens)

# short-term average climate, Matt
mod2.ens.matt <- lm(mean_pred_climate ~ region, data=dat)
summary(mod2.ens.matt)

# average stream
mod2.ens.stream <- lm(mean_pred_stream ~ region, data=dat)
summary(mod2.ens.stream)

# individual long-term climate models, amy
mod2.lr <- lm(LRavg ~ region, data=dat)
summary(mod2.lr)

mod2.gam <- lm(GAMavg ~ region, data=dat)
summary(mod2.gam)

mod2.rf <- lm(RFavg ~ region, data=dat)
summary(mod2.rf)

mod2.brt <- lm(BRTavg ~ region, data=dat)
summary(mod2.brt)

mod2.max <- lm(MAXavg ~ region, data=dat)
summary(mod2.max)

plot_lat_suit1 <- ggplot(dat, aes(lat, Ens)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  #xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlab("") + 
  xlim(43,45.5) +
  ylab("Climate (1980-2010)") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
#ggplot2::ggsave("Figures/LongTermClimateENM_vs_Latitude.png", width=8, height=8)

plot_lat_suit2 <- ggplot(dat, aes(lat, mean_pred_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  #xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlab("") + 
  xlim(43,45.5) +
  ylab("Climate (2014-15)") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt"))
#ggplot2::ggsave("Figures/ShortTermClimateENM_vs_Latitude.png", width=8, height=8)

plot_lat_suit3 <- ggplot(dat, aes(lat, mean_pred_stream)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlim(43,45.5) +
  ylab("Physical habitat") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt"))
#ggplot2::ggsave("Figures/StreamENM_vs_Latitude.png", width=8, height=8)

latsuit <- plot_grid(plot_lat_suit1 + theme(legend.position = "none"), 
                     plot_lat_suit2 + theme(legend.position = "none"), 
                     plot_lat_suit3 + theme(legend.position = "none"),
                     nrow=3)
legend <- get_legend(plot_lat_suit1)
latsuit2 <- plot_grid(legend, latsuit, ncol=1, rel_heights = c(0.3, 3))
left_label <- "ENM suitability"
latsuit3 <- ggdraw(latsuit2) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/Suitability_vs_Latitude_3Panel.png", latsuit3, base_width=5, base_height=11)

# individual models
LRlat <- ggplot(dat, aes(lat, LRavg)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("LR") +
  ylim(0,0.8) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
GAMlat <- ggplot(dat, aes(lat, GAMavg)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("GAM") +
  ylim(0,0.8) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
RFlat <- ggplot(dat, aes(lat, RFavg)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("RF") +
  ylim(0,0.8) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
BRTlat <- ggplot(dat, aes(lat, BRTavg)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("BRT") +
  ylim(0,0.8) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
MAXlat <- ggplot(dat, aes(lat, MAXavg)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlim(43,45.5) +
  ylab("MAX") +
  ylim(0,0.8) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 

latsuit_ind_long <- plot_grid(LRlat + theme(legend.position = "none"), 
                   GAMlat + theme(legend.position = "none"), 
                   RFlat + theme(legend.position = "none"), 
                   BRTlat + theme(legend.position = "none"), 
                   MAXlat + theme(legend.position = "none"), 
                   ncol=1, labels="AUTO", label_x=0.9)
legend <- get_legend(LRlat)
latsuit_ind_long2 <- plot_grid(legend, latsuit_ind_long, ncol=1, rel_heights = c(0.3, 5))
left_label <- "ENM suitability: long-term climate"
latsuit_ind_long3 <- ggdraw(latsuit_ind_long2) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/LongTermClimateENM_vs_Latitude_IndModels.png", latsuit_ind_long3, base_width=5, base_height=11)

LRlat2 <- ggplot(dat, aes(lat, glm_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("LR") +
  ylim(0.7,0.95) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
GAMlat2 <- ggplot(dat, aes(lat, gam_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("GAM") +
  ylim(0.7,0.95) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
RFlat2 <- ggplot(dat, aes(lat, rf_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("RF") +
  ylim(0.2,0.8) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
BRTlat2 <- ggplot(dat, aes(lat, brt_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("BRT") +
  ylim(0.1,0.9) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
MAXlat2 <- ggplot(dat, aes(lat, max_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlim(43,45.5) +
  ylab("MAX") +
  ylim(0.4,0.6) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 

latsuit_ind_short <- plot_grid(LRlat2 + theme(legend.position = "none"), 
                         GAMlat2 + theme(legend.position = "none"), 
                         RFlat2 + theme(legend.position = "none"), 
                         BRTlat2 + theme(legend.position = "none"), 
                         MAXlat2 + theme(legend.position = "none"), 
                         ncol=1, labels="AUTO", label_x=0.9)
legend <- get_legend(LRlat2)
latsuit_ind_short2 <- plot_grid(legend, latsuit_ind_short, ncol=1, rel_heights = c(0.3, 5))
left_label <- "ENM suitability: short-term climate"
latsuit_ind_short3 <- ggdraw(latsuit_ind_short2) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/ShortTermClimateENM_vs_Latitude_IndModels.png", latsuit_ind_short3, base_width=5, base_height=11)

LRlat3 <- ggplot(dat, aes(lat, glm_stream)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("LR") +
  ylim(0,0.1) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
GAMlat3 <- ggplot(dat, aes(lat, gam_stream)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("GAM") +
  ylim(0,0.1) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
RFlat3 <- ggplot(dat, aes(lat, rf_stream)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("RF") +
  ylim(0,0.7) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
BRTlat3 <- ggplot(dat, aes(lat, brt_stream)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("BRT") +
  ylim(0,0.2) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
MAXlat3 <- ggplot(dat, aes(lat, max_stream)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlim(43,45.5) +
  ylab("MAX") +
  ylim(0.3,0.8) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 

latsuit_ind_stream <- plot_grid(LRlat3 + theme(legend.position = "none"), 
                               GAMlat3 + theme(legend.position = "none"), 
                               RFlat3 + theme(legend.position = "none"), 
                               BRTlat3 + theme(legend.position = "none"), 
                               MAXlat3 + theme(legend.position = "none"), 
                               ncol=1, labels="AUTO", label_x=0.9)
legend <- get_legend(LRlat3)
latsuit_ind_stream2 <- plot_grid(legend, latsuit_ind_stream, ncol=1, rel_heights = c(0.3, 5))
left_label <- "ENM suitability: physical habitat"
latsuit_ind_stream3 <- ggdraw(latsuit_ind_stream2) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/StreamENM_vs_Latitude_IndModels.png", latsuit_ind_stream3, base_width=5, base_height=11)

### linear models of lambda vs suitability 
# climate suitability (Amy March 2018)
mod0b.ens <- lm(lambda ~ 1, data=dat)
mod1b.ens <- lm(lambda ~ Ens, data=dat)
model.sel(mod0b.ens, mod1b.ens)
summary(mod1b.ens)

# climate suitability (Matt March 2018)
mod0b.ens.matt <- lm(lambda ~ 1, data=dat)
mod1b.ens.matt <- lm(lambda ~ mean_pred_climate, data=dat)
model.sel(mod0b.ens.matt, mod1b.ens.matt)
summary(mod1b.ens.matt)

# Stream suitability 
mod0b.ens.stream <- lm(lambda ~ 1, data=dat)
mod1b.ens.stream <- lm(lambda ~ mean_pred_stream, data=dat)
model.sel(mod0b.ens.stream, mod1b.ens.stream)
summary(mod1b.ens.stream)

# individual models, Amy climate
mod0b.lr <- lm(lambda ~ 1, data=dat)
mod1b.lr <- lm(lambda ~ LRavg, data=dat)
model.sel(mod0b.lr, mod1b.lr)
summary(mod1b.lr)

mod0b.gam <- lm(lambda ~ 1, data=dat)
mod1b.gam <- lm(lambda ~ GAMavg, data=dat)
model.sel(mod0b.gam, mod1b.gam)
summary(mod1b.gam)

mod0b.rf <- lm(lambda ~ 1, data=dat)
mod1b.rf <- lm(lambda ~ RFavg, data=dat)
model.sel(mod0b.rf, mod1b.rf)
summary(mod1b.rf)

mod0b.brt <- lm(lambda ~ 1, data=dat)
mod1b.brt <- lm(lambda ~ BRTavg, data=dat)
model.sel(mod0b.brt, mod1b.brt)
summary(mod1b.brt)

mod0b.max <- lm(lambda ~ 1, data=dat)
mod1b.max <- lm(lambda ~ MAXavg, data=dat)
model.sel(mod0b.max, mod1b.max)
summary(mod1b.max)

# individual models, Matt climate
mod0b.lr <- lm(lambda ~ 1, data=dat)
mod1b.lr <- lm(lambda ~ glm_climate, data=dat)
model.sel(mod0b.lr, mod1b.lr)
summary(mod1b.lr)

mod0b.gam <- lm(lambda ~ 1, data=dat)
mod1b.gam <- lm(lambda ~ gam_climate, data=dat)
model.sel(mod0b.gam, mod1b.gam)
summary(mod1b.gam)

mod0b.rf <- lm(lambda ~ 1, data=dat)
mod1b.rf <- lm(lambda ~ rf_climate, data=dat)
model.sel(mod0b.rf, mod1b.rf)
summary(mod1b.rf)

mod0b.brt <- lm(lambda ~ 1, data=dat)
mod1b.brt <- lm(lambda ~ brt_climate, data=dat)
model.sel(mod0b.brt, mod1b.brt)
summary(mod1b.brt)

mod0b.max <- lm(lambda ~ 1, data=dat)
mod1b.max <- lm(lambda ~ max_climate, data=dat)
model.sel(mod0b.max, mod1b.max)
summary(mod1b.max)

# individual models, stream
mod0b.lr.stream <- lm(lambda ~ 1, data=dat)
mod1b.lr.stream <- lm(lambda ~ glm_stream, data=dat)
model.sel(mod0b.lr.stream, mod1b.lr.stream)
summary(mod1b.lr.stream)

mod0b.gam.stream <- lm(lambda ~ 1, data=dat)
mod1b.gam.stream <- lm(lambda ~ gam_stream, data=dat)
model.sel(mod0b.gam.stream, mod1b.gam.stream)
summary(mod1b.gam.stream)

mod0b.rf.stream <- lm(lambda ~ 1, data=dat)
mod1b.rf.stream<- lm(lambda ~ rf_stream, data=dat)
model.sel(mod0b.rf.stream, mod1b.rf.stream)
summary(mod1b.rf)

mod0b.brt.stream <- lm(lambda ~ 1, data=dat)
mod1b.brt.stream <- lm(lambda ~ brt_stream, data=dat)
model.sel(mod0b.brt.stream, mod1b.brt.stream)
summary(mod1b.brt.stream)

mod0b.max.stream <- lm(lambda ~ 1, data=dat)
mod1b.max.stream <- lm(lambda ~ max_stream, data=dat)
model.sel(mod0b.max.stream, mod1b.max.stream)
summary(mod1b.max.stream)

# lambda vs ensemble suitability
plot_lam_suit1 <- ggplot(dat, aes(Ens, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0.01) + 
  geom_errorbarh(aes(xmin=Ens+clim_amy_se, xmax=Ens-clim_amy_se),height=0.1)+
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("Climate (1980-2010)") + 
  #xlim(0.2,0.57) + 
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(1.85)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(6,0,40,0), "pt")) 
#ggplot2::ggsave("Figures/Lambda_vs_ClimateENM_8010.png", width=8, height=8)

plot_lam_suit2 <- ggplot(dat, aes(mean_pred_climate, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0.01) + 
  geom_errorbarh(aes(xmin=mean_pred_climate+clim_matt_se, xmax=mean_pred_climate-clim_matt_se),height=0.1)+
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("Climate (2014-15)") + 
  #xlim(0.44,0.75) + 
  ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  #ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(1.85)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(6,0,40,0), "pt")) 
#ggplot2::ggsave("Figures/Lambda_vs_ClimateENM_1415.png", width=8, height=8)

plot_lam_suit3 <- ggplot(dat, aes(mean_pred_stream, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0.01) + 
  geom_errorbarh(aes(xmin=mean_pred_stream+stream_se, xmax=mean_pred_stream-stream_se),height=0.1)+
  geom_smooth(method=lm, se=FALSE, color="black", linetype="dashed") + 
  xlab("Physical habitat") + 
  #xlim(0.04,0.36) + 
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(1.85)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(6,0,40,0), "pt")) 
#ggplot2::ggsave("Figures/Lambda_vs_StreamENM.png", width=8, height=8)

lamsuit <- plot_grid(plot_lam_suit2 + theme(legend.position = "none"), 
                   plot_lam_suit1 + theme(legend.position = "none"), 
                   plot_lam_suit3+ theme(legend.position = "none"), 
                   nrow=1, labels="AUTO", label_x=0.9)
legend <- get_legend(plot_lam_suit1)
lamsuit2 <- plot_grid(lamsuit, legend, rel_widths = c(3, 0.3))
bottom_label <- "ENM Suitability"
lamsuit3 <- ggdraw(lamsuit2) + draw_label(bottom_label, angle=0, y=0.05, size=24)
save_plot("Figures/Lambda_vs_Suitability_3Panel_error.png", lamsuit3, base_width=11, base_height=5)

# individual models, Amy climate
LRlam <- ggplot(dat, aes(LRavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  #xlab("LR") + 
  xlab("") + 
  xlim(0.1,0.7) +
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("Climate 1980") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,0,100), "pt"))
GAMlam <- ggplot(dat, aes(GAMavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  #xlab("GAM") + 
  xlab("") + 
  xlim(0.1,0.7) +
  #ylab("Climate 1980-2010") +
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", plot.margin = unit(c(0,0,0,0), "pt")) 
RFlam <- ggplot(dat, aes(RFavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + 
  #xlab("RF") + 
  xlab("") + 
  xlim(0,0.7) +
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", plot.margin = unit(c(0,0,0,0), "pt")) 
BRTlam <- ggplot(dat, aes(BRTavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + 
  #xlab("BRT") + 
  xlab("") + 
  #xlim(0,0.85) +
  scale_x_continuous(limits=c(0, 0.85), breaks=c(0.1, 0.4, 0.7)) + 
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", plot.margin = unit(c(0,0,0,0), "pt")) 
MAXlam <- ggplot(dat, aes(MAXavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  #xlab("MAX") + 
  xlab("") + 
  #xlim(0.27,0.5) +
  scale_x_continuous(limits=c(0.27, 0.5), breaks=c(0.3, 0.4, 0.5)) + 
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", plot.margin = unit(c(0,0,0,0), "pt")) 

suitlam_long_ind <- plot_grid(LRlam + theme(legend.position="none"), 
                               GAMlam + theme(legend.position="none"), 
                               RFlam + theme(legend.position="none"), 
                               BRTlam + theme(legend.position="none"), 
                               MAXlam + theme(legend.position="none"), 
                               nrow=1, labels="AUTO", label_x=0.9)
legend <- get_legend(LRlam)
suitlam_long_ind2 <- plot_grid(suitlam_long_ind, legend, nrow=1, rel_widths = c(5, 0.5))
bottom_label = "ENM suitability: long-term climate"
suitlam_long_ind3 <- ggdraw(suitlam_long_ind2) + draw_label(bottom_label, angle=0, y=0.05, size=24)
save_plot("Figures/Lambda_vs_ClimateENM_8010_IndModels.png", suitlam_long_ind3, base_width=11, base_height=5)

# individual models, Matt climate
LRlam2 <- ggplot(dat, aes(glm_climate, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  #xlab("LR") + 
  xlab("") + 
  #xlim(0.7, 0.95) +
  scale_x_continuous(limits=c(0.7, 0.95), breaks=c(0.7, 0.8, 0.9)) + 
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("Climate 2014") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,0,100), "pt"))
GAMlam2 <- ggplot(dat, aes(gam_climate, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  #xlab("GAM") + 
  xlab("") + 
  #xlim(0.7, 0.95) +
  scale_x_continuous(limits=c(0.7, 0.95), breaks=c(0.7, 0.8, 0.9)) + 
  ylab("") +
  #ylab("Climate 2014-15") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", plot.margin = unit(c(0,0,0,0), "pt")) 
RFlam2 <- ggplot(dat, aes(rf_climate, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  #xlab("RF") + 
  xlab("") + 
  xlim(0.2,0.85) +
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", plot.margin = unit(c(0,0,0,0), "pt")) 
BRTlam2 <- ggplot(dat, aes(brt_climate, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + 
  #xlab("BRT") + 
  xlab("") + 
  #xlim(0.1,0.85) +
  #scale_x_continuous(limits=c(0, 0.85), breaks=c(0.1, 0.4, 0.7)) + 
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", plot.margin = unit(c(0,0,0,0), "pt")) 
MAXlam2 <- ggplot(dat, aes(max_climate, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  #xlab("MAX") + 
  xlab("") + 
  xlim(0.37,0.55) +
  #scale_x_continuous(limits=c(0.37, 0.55), breaks=c(0.4, 0.5)) + 
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", plot.margin = unit(c(0,0,0,0), "pt")) 

suitlam_short_ind <- plot_grid(LRlam2 + theme(legend.position="none"), 
                              GAMlam2 + theme(legend.position="none"), 
                              RFlam2 + theme(legend.position="none"), 
                              BRTlam2 + theme(legend.position="none"), 
                              MAXlam2 + theme(legend.position="none"), 
                              nrow=1, labels="AUTO", label_x=0.9)
legend <- get_legend(LRlam2)
suitlam_short_ind2 <- plot_grid(suitlam_short_ind, legend, nrow=1, rel_widths = c(5, 0.5))
bottom_label = "ENM suitability: short-term climate"
suitlam_short_ind3 <- ggdraw(suitlam_short_ind2) + draw_label(bottom_label, angle=0, y=0.05, size=24)
save_plot("Figures/Lambda_vs_ClimateENM_1415_IndModels.png", suitlam_short_ind3, base_width=11, base_height=5)

# individual models, stream
LRlam3 <- ggplot(dat, aes(glm_stream, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + 
  xlab("LR") + 
  xlim(0, 0.07) +
  #scale_x_continuous(limits=c(0.7, 0.95), breaks=c(0.7, 0.8, 0.9)) + 
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("Habitat") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,0,80,100), "pt"))
GAMlam3 <- ggplot(dat, aes(gam_stream, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + 
  xlab("GAM") + 
  xlim(0, 0.07) +
  #scale_x_continuous(limits=c(0.7, 0.95), breaks=c(0.7, 0.8, 0.9)) + 
  ylab("") +
  #ylab("Physical habitat") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", plot.margin = unit(c(0,0,80,0), "pt")) 
RFlam3 <- ggplot(dat, aes(rf_stream, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("RF") + 
  xlim(0,0.7) +
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", plot.margin = unit(c(0,0,80,0), "pt")) 
BRTlam3 <- ggplot(dat, aes(brt_stream, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + 
  xlab("BRT") + 
  #xlim(0.1,0.85) +
  #scale_x_continuous(limits=c(0, 0.85), breaks=c(0.1, 0.4, 0.7)) + 
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", plot.margin = unit(c(0,0,80,0), "pt")) 
MAXlam3 <- ggplot(dat, aes(max_stream, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + 
  xlab("MAX") + 
  #xlim(0.37,0.55) +
  #scale_x_continuous(limits=c(0.37, 0.55), breaks=c(0.4, 0.5)) + 
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", plot.margin = unit(c(0,0,80,0), "pt")) 

suitlam_stream_ind <- plot_grid(LRlam3 + theme(legend.position="none"), 
                               GAMlam3 + theme(legend.position="none"), 
                               RFlam3 + theme(legend.position="none"), 
                               BRTlam3 + theme(legend.position="none"), 
                               MAXlam3 + theme(legend.position="none"), 
                               nrow=1, labels="AUTO", label_x=0.9)
legend <- get_legend(LRlam3)
suitlam_stream_ind2 <- plot_grid(suitlam_stream_ind, legend, nrow=1, rel_widths = c(5, 0.5))
bottom_label = "ENM suitability: physical habitat"
suitlam_stream_ind3 <- ggdraw(suitlam_stream_ind2) + draw_label(bottom_label, angle=0, y=0.05, size=24)
save_plot("Figures/Lambda_vs_StreamENM_IndModels.png", suitlam_stream_ind3, base_width=11, base_height=5)

suitlam_all_ind <- plot_grid(LRlam + theme(legend.position="none"), 
                             GAMlam + theme(legend.position="none"), 
                             RFlam + theme(legend.position="none"), 
                             BRTlam + theme(legend.position="none"), 
                             MAXlam + theme(legend.position="none"),
                             LRlam2 + theme(legend.position="none"), 
                             GAMlam2 + theme(legend.position="none"), 
                             RFlam2 + theme(legend.position="none"), 
                             BRTlam2 + theme(legend.position="none"), 
                             MAXlam2 + theme(legend.position="none"),
                             LRlam3 + theme(legend.position="none"), 
                             GAMlam3 + theme(legend.position="none"),
                             RFlam3 + theme(legend.position="none"),
                             BRTlam3 + theme(legend.position="none"),
                             MAXlam3 + theme(legend.position="none"),
                             nrow=3, labels="AUTO", label_x=0.9, align="hv")
legend <- get_legend(LRlam3)
suitlam_all_ind2 <- plot_grid(suitlam_all_ind, legend, nrow=1, rel_widths = c(5, 0.5))
left_label <- expression(paste("Population growth rate (", lambda, ")"))
bottom_label = "ENM algorithm"
suitlam_all_ind3 <- ggdraw(suitlam_all_ind2) + draw_label(bottom_label, angle=0, y=0.05, size=24)
suitlam_all_ind4 <- ggdraw(suitlam_all_ind3) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/Lambda_vs_AllENM_IndModels.png", suitlam_all_ind4, base_width=22, base_height=11)

suitlam_best_ind <- plot_grid(GAMlam + theme(legend.position="none"), 
                             RFlam + theme(legend.position="none"), 
                             GAMlam2 + theme(legend.position="none"), 
                             RFlam2 + theme(legend.position="none"), 
                             GAMlam3 + theme(legend.position="none"),
                             RFlam3 + theme(legend.position="none"),
                             nrow=3, labels="AUTO", label_x=0.9, align="hv")
legend <- get_legend(GAMlam3)
suitlam_best_ind2 <- plot_grid(suitlam_best_ind, legend, nrow=1, rel_widths = c(2, 0.5))
left_label <- expression(paste("Population growth rate (", lambda, ")"))
bottom_label = "ENM algorithm"
suitlam_best_ind3 <- ggdraw(suitlam_best_ind2) + draw_label(bottom_label, angle=0, y=0.05, size=24)
suitlam_best_ind4 <- ggdraw(suitlam_best_ind3) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/Lambda_vs_BestENM_IndModels.png", suitlam_best_ind4, base_width=8.5, base_height=11)

# exploring individual vital rates
vitals <- read_csv("Robjects/vital_rate_coefficients.csv")

dat2 <- left_join(dat, vitals, by=c("Site"="site"))

# linear models of vitals vs latitude
surv.lat <- lm(surv.siteint ~ lat, data=dat2)
summary(surv.lat)

growth.lat <- lm(growth.siteint ~ lat, data=dat2)
summary(growth.lat)

flower.lat <- lm(flowering.siteint ~ lat, data=dat2)
summary(flower.lat)

fruit.lat <- lm(fruits.siteint ~ lat, data=dat2)
summary(fruit.lat)

# linear models of vitals vs region
surv.reg <- lm(surv.siteint ~ region, data=dat2)
summary(surv.reg)
visreg(surv.reg)

growth.reg <- lm(growth.siteint ~ region, data=dat2)
summary(growth.reg)

flower.reg <- lm(flowering.siteint ~ region, data=dat2)
summary(flower.reg)

fruit.reg <- lm(fruits.siteint ~ region, data=dat2)
summary(fruit.reg)

# individual models
survlat <- ggplot(dat2, aes(lat, surv.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("Survival") +
  ylim(-2.5, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.8)), plot.margin = unit(c(0,6,0,40), "pt")) 
growthlat <- ggplot(dat2, aes(lat, growth.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("Growth") +
  ylim(-2, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", plot.margin = unit(c(0,6,0,40), "pt"))
flowerlat <- ggplot(dat2, aes(lat, flowering.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("") + 
  xlim(43,45.5) +
  ylab("Flowering") +
  ylim(-4, 1) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", plot.margin = unit(c(0,6,0,40), "pt"))
fruitlat <- ggplot(dat2, aes(lat, fruits.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlim(43,45.5) +
  ylab("Fruits") +
  ylim(-1, 0.3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", plot.margin = unit(c(0,6,0,40), "pt"))

vitallat <- plot_grid(survlat + theme(legend.position = "none"), 
                      growthlat + theme(legend.position = "none"), 
                      flowerlat + theme(legend.position = "none"), 
                      fruitlat + theme(legend.position = "none"), 
                      nrow=4, ncol=1, labels="AUTO", label_x=0.9)
legend <- get_legend(survlat)
vitallat2 <- plot_grid(legend, vitallat, nrow=2, rel_heights=c(0.2, 4))
left_label <- "Vital rate intercepts"
vitallat3 <- ggdraw(vitallat2) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/Vitals_vs_Latitude.png", vitallat3, base_width=5, base_height=11)

## linear models of vitals vs suitability
# Amy climate
surv.suit <- lm(surv.siteint ~ Ens, data=dat2)
summary(surv.suit)

growth.suit <- lm(growth.siteint ~ Ens, data=dat2)
summary(growth.suit)

flower.suit <- lm(flowering.siteint ~ Ens, data=dat2)
summary(flower.suit)

fruit.suit <- lm(fruits.siteint ~ Ens, data=dat2)
summary(fruit.suit)

# Matt climate
surv.suit.matt <- lm(surv.siteint ~ mean_pred_climate, data=dat2)
summary(surv.suit.matt)

growth.suit.matt <- lm(growth.siteint ~ mean_pred_climate, data=dat2)
summary(growth.suit.matt)

flower.suit.matt <- lm(flowering.siteint ~ mean_pred_climate, data=dat2)
summary(flower.suit.matt)

fruit.suit.matt <- lm(fruits.siteint ~ mean_pred_climate, data=dat2)
summary(fruit.suit.matt)

# stream
surv.suit.stream <- lm(surv.siteint ~ mean_pred_stream, data=dat2)
summary(surv.suit.stream)

growth.suit.stream <- lm(growth.siteint ~ mean_pred_stream, data=dat2)
summary(growth.suit.stream)

flower.suit.stream <- lm(flowering.siteint ~ mean_pred_stream, data=dat2)
summary(flower.suit.stream)

fruit.suit.stream <- lm(fruits.siteint ~ mean_pred_stream, data=dat2)
summary(fruit.suit.stream)

# individual vitals
survsuit.amy <- ggplot(dat2, aes(Ens, surv.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black", linetype="dashed") +
  xlab("") + 
  #ylab("Survival") +
  ylab("") +
  ylim(-3, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.8)), plot.margin = unit(c(0,0,0,0), "pt"))
growthsuit.amy <- ggplot(dat2, aes(Ens, growth.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black") +
  xlab("") + 
  #ylab("Growth") +
  ylab("") +
  ylim(-2, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,0,0,0), "pt"))
flowersuit.amy <- ggplot(dat2, aes(Ens, flowering.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black", linetype="dashed") +
  xlab("") + 
  #ylab("Flowering") +
  ylab("") +
  ylim(-4, 1) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,0,0,0), "pt"))
fruitsuit.amy <- ggplot(dat2, aes(Ens, fruits.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("Climatic suitability \n(1980-2010)") + 
  ylab("") +
  #ylab("Fruits") +
  #ylim(-1, 0.3) + 
  scale_y_continuous(limit=c(-1.5,0.5), breaks=c(-1,0)) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,0,0,0), "pt"))

vitalsuit_long <- plot_grid(survsuit.amy + theme(legend.position = "none"), 
                            growthsuit.amy + theme(legend.position = "none"), 
                            flowersuit.amy + theme(legend.position = "none"), 
                            fruitsuit.amy + theme(legend.position = "none"), 
                            ncol=1, labels="AUTO", label_x=0.9)
legend <- get_legend(survsuit.amy)
vitalsuit_long2 <- plot_grid(legend, vitalsuit_long, ncol=1, rel_heights = c(0.2, 4))
left_label <- "Vital rate intercepts"
vitalsuit_long3 <- ggdraw(vitalsuit_long2) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/Vitals_vs_Suit_Clim8010.png", vitalsuit_long3, base_width=5, base_height=11)

survsuit.matt <- ggplot(dat2, aes(mean_pred_climate, surv.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black", linetype="dashed") +
  xlab("") + 
  scale_x_continuous(limits=c(0.4, 0.8), breaks=seq(0.4, 0.8, by=0.1)) + 
  ylab("Survival") +
  ylim(-3, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.8)), plot.margin = unit(c(0,0,0,40), "pt"))
growthsuit.matt <- ggplot(dat2, aes(mean_pred_climate, growth.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, linetype="dashed", colour="black") +
  xlab("") + 
  scale_x_continuous(limits=c(0.4, 0.8), breaks=seq(0.4, 0.8, by=0.1)) +
  ylab("Growth") +
  ylim(-2, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,0,0,40), "pt"))
flowersuit.matt <- ggplot(dat2, aes(mean_pred_climate, flowering.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black", linetype="dashed") +
  xlab("") + 
  scale_x_continuous(limits=c(0.4, 0.8), breaks=seq(0.4, 0.8, by=0.1)) +
  ylab("Flowering") +
  ylim(-4, 1) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,0,0,40), "pt"))
fruitsuit.matt <- ggplot(dat2, aes(mean_pred_climate, fruits.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + 
  xlab("Climatic suitability \n(2014-2015)") + 
  scale_x_continuous(limits=c(0.4, 0.8), breaks=seq(0.4, 0.8, by=0.1)) +
  ylab("Fruits") +
  #ylim(-1, 0.3) + 
  scale_y_continuous(limit=c(-1.5,0.5), breaks=c(-1,0)) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,0,0,40), "pt"))

vitalsuit_short <- plot_grid(survsuit.matt + theme(legend.position = "none"), 
                            growthsuit.matt + theme(legend.position = "none"), 
                            flowersuit.matt + theme(legend.position = "none"), 
                            fruitsuit.matt + theme(legend.position = "none"), 
                            ncol=1, labels="AUTO", label_x=0.9)
legend <- get_legend(survsuit.matt)
vitalsuit_short2 <- plot_grid(legend, vitalsuit_short, ncol=1, rel_heights = c(0.2, 4))
left_label <- "Vital rate intercepts"
vitalsuit_short3 <- ggdraw(vitalsuit_short2) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/Vitals_vs_Suit_Clim1415.png", vitalsuit_short3, base_width=5, base_height=11)

vitalsuit_clim <- plot_grid(survsuit.matt + theme(legend.position="none"),
                            survsuit.amy + theme(legend.position="none"), 
                            growthsuit.matt + theme(legend.position="none"),
                            growthsuit.amy + theme(legend.position="none"),
                            flowersuit.matt + theme(legend.position="none"),
                            flowersuit.amy + theme(legend.position="none"),
                            fruitsuit.matt + theme(legend.position="none"),
                            fruitsuit.amy + theme(legend.position="none"),
                            ncol=2, nrow=4, labels="AUTO", label_x=0.9)
left_label <- "Vital rate intercept"
vitalsuit_clim2 <- ggdraw(vitalsuit_clim) + draw_label(left_label, angle=90, size=24, x=0.02)
save_plot("Figures/Vitals_vs_Suit_Clim.png", vitalsuit_clim2, base_width=8.5, base_height=11)

survsuit.stream <- ggplot(dat2, aes(mean_pred_stream, surv.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black", linetype="dashed") +
  xlab("") + 
  ylab("Survival") +
  ylim(-2.5, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.8)), plot.margin = unit(c(0,10,0,40), "pt"))
growthsuit.stream <- ggplot(dat2, aes(mean_pred_stream, growth.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, linetype="dashed", colour="black") +
  xlab("") + 
  ylab("Growth") +
  ylim(-2, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,10,0,40), "pt"))
flowersuit.stream <- ggplot(dat2, aes(mean_pred_stream, flowering.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black", linetype="dashed") +
  xlab("") + 
  ylab("Flowering") +
  ylim(-4, 1) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,10,0,40), "pt"))
fruitsuit.stream <- ggplot(dat2, aes(mean_pred_stream, fruits.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + 
  xlab("Physical habitat suitability") + 
  ylab("Fruits") +
  ylim(-1, 0.3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,10,40,40), "pt"))

vitalsuit_stream <- plot_grid(survsuit.stream + theme(legend.position = "none"), 
                             growthsuit.stream + theme(legend.position = "none"), 
                             flowersuit.stream + theme(legend.position = "none"), 
                             fruitsuit.stream + theme(legend.position = "none"), 
                             ncol=1, labels="AUTO", label_x=0.9)
legend <- get_legend(survsuit.stream)
vitalsuit_stream2 <- plot_grid(legend, vitalsuit_stream, ncol=1, rel_heights = c(0.2, 4))
left_label <- "Vital rate intercepts"
vitalsuit_stream3 <- ggdraw(vitalsuit_stream2) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/Vitals_vs_Suit_Stream.png", vitalsuit_stream3, base_width=5, base_height=11)

# correlations among vital rates
vitals.slim <- dat2 %>% dplyr::select(surv.siteint, growth.siteint, flowering.siteint, fruits.siteint, lambda)
vitals.slim <- as.data.frame(vitals.slim)

panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}
pairs(vitals.slim[1:5], lower.panel=panel.smooth, upper.panel=panel.cor)


########## Weighted ensembles

#### Calculate weighted ensemble scores
AUC <- read_csv("Data/AUCscores.csv") %>% select(Model, ValidationMethod,AUC_LR,AUC_GAM,AUC_RF,AUC_BRT,AUC_MAX) %>% filter(Model=="PhysicalHabitat"|Model=="Climate")

stream_AUC_ext <- AUC %>% filter(Model=="PhysicalHabitat", ValidationMethod=="ExternalValidation") %>% select(AUC_LR.stream = AUC_LR, AUC_GAM.stream=AUC_GAM, AUC_RF.stream=AUC_RF, AUC_BRT.stream=AUC_BRT, AUC_MAX.stream=AUC_MAX) 

clim_AUC_ext <- AUC %>% filter(Model=="Climate", ValidationMethod=="ExternalValidation") %>% select(AUC_LR.clim = AUC_LR, AUC_GAM.clim=AUC_GAM, AUC_RF.clim=AUC_RF, AUC_BRT.clim=AUC_BRT, AUC_MAX.clim=AUC_MAX) 

joint_AUC_ext <- cbind(stream_AUC_ext, clim_AUC_ext)

joint_scores <- dat %>% select(site, LRavg, GAMavg, RFavg, BRTavg, MAXavg, glm_climate, gam_climate, rf_climate, brt_climate, max_climate, glm_stream, gam_stream, rf_stream, brt_stream, max_stream)

joint <- cbind(joint_AUC_ext, joint_scores)

joint_wtd_ens <- joint %>% 
  group_by(site) %>%
  mutate(
    stream_wtd_ens = AUC_LR.stream*glm_stream + AUC_GAM.stream*gam_stream + AUC_RF.stream*rf_stream + AUC_BRT.stream*brt_stream + AUC_MAX.stream*max_stream,
    stream_wtd_ens_sd = sd(c(AUC_LR.stream*glm_stream, AUC_GAM.stream*gam_stream, AUC_RF.stream*rf_stream, AUC_BRT.stream*brt_stream, AUC_MAX.stream*max_stream)),
    clim_wtd_ens_shortterm = AUC_LR.clim*glm_climate + AUC_GAM.clim*gam_climate + AUC_RF.clim*rf_climate + AUC_BRT.clim*brt_climate + AUC_MAX.clim*max_climate,
    clim_wtd_ens_shortterm_sd = sd(c(AUC_LR.clim*glm_climate,AUC_GAM.clim*gam_climate,AUC_RF.clim*rf_climate,AUC_BRT.clim*brt_climate,AUC_MAX.clim*max_climate)),
    clim_wtd_ens_longterm = AUC_LR.clim*LRavg + AUC_GAM.clim*GAMavg + AUC_RF.clim*RFavg + AUC_BRT.clim*BRTavg + AUC_MAX.clim*MAXavg,
    clim_wtd_ens_longterm_sd = sd(c(AUC_LR.clim*LRavg,AUC_GAM.clim*GAMavg,AUC_RF.clim*RFavg,AUC_BRT.clim*BRTavg,AUC_MAX.clim*MAXavg)),
    joint_wtd_ens_longterm = (stream_wtd_ens + clim_wtd_ens_longterm)/2,
    joint_wtd_ens_longterm_sd = sd(c(AUC_LR.clim*LRavg, AUC_GAM.clim*GAMavg, AUC_RF.clim*RFavg, AUC_BRT.clim*BRTavg, AUC_MAX.clim*MAXavg, AUC_LR.stream*glm_stream, AUC_GAM.stream*gam_stream, AUC_RF.stream*rf_stream, AUC_BRT.stream*brt_stream, AUC_MAX.stream*max_stream)),
    joint_wtd_ens_shortterm = (stream_wtd_ens + clim_wtd_ens_shortterm)/2,
    joint_wtd_ens_shortterm_sd = sd(c(AUC_LR.clim*glm_climate, AUC_GAM.clim*gam_climate, AUC_RF.clim*rf_climate, AUC_BRT.clim*brt_climate, AUC_MAX.clim*max_climate, AUC_LR.stream*glm_stream, AUC_GAM.stream*gam_stream, AUC_RF.stream*rf_stream, AUC_BRT.stream*brt_stream, AUC_MAX.stream*max_stream))
  ) %>% 
  ungroup()

dat <- left_join(dat, joint_wtd_ens)

#### Redo key stats and plots with new weighted ensembles

### linear models of suitability vs latitude
## ave suitability
# Amy's values, climate long-term
mod0.ens <- lm(clim_wtd_ens_longterm ~ 1, data=dat)
mod1.ens <- lm(clim_wtd_ens_longterm ~ lat, data=dat)
model.sel(mod0.ens, mod1.ens)
summary(mod1.ens)

# Matt's values, climate short-term
mod0.ens.matt <- lm(clim_wtd_ens_shortterm ~ 1, data=dat)
mod1.ens.matt <- lm(clim_wtd_ens_shortterm ~ lat, data=dat)
model.sel(mod0.ens.matt, mod1.ens.matt)
summary(mod1.ens.matt)

# Matt's values, stream
mod0.ens.stream <- lm(stream_wtd_ens ~ 1, data=dat)
mod1.ens.stream <- lm(stream_wtd_ens ~ lat, data=dat)
model.sel(mod0.ens.stream, mod1.ens.stream)
summary(mod1.ens.stream)

# joint values, climate long-term
mod0.ens.joint.long <- lm(joint_wtd_ens_longterm ~ 1, data=dat)
mod1.ens.joint.long <- lm(joint_wtd_ens_longterm ~ lat, data=dat)
model.sel(mod0.ens.joint.long, mod1.ens.joint.long)
summary(mod1.ens.joint.long)

# joint values, climate short-term
mod0.ens.joint.short <- lm(joint_wtd_ens_shortterm ~ 1, data=dat)
mod1.ens.joint.short <- lm(joint_wtd_ens_shortterm ~ lat, data=dat)
model.sel(mod0.ens.joint.short, mod1.ens.joint.short)
summary(mod1.ens.joint.short)

### anova of suitability by range position
# long-term average climate, Amy
mod2.ens <- lm(clim_wtd_ens_longterm ~ region, data=dat)
summary(mod2.ens)

# short-term average climate, Matt
mod2.ens.matt <- lm(clim_wtd_ens_shortterm ~ region, data=dat)
summary(mod2.ens.matt)

# average stream
mod2.ens.stream <- lm(stream_wtd_ens ~ region, data=dat)
summary(mod2.ens.stream)

# joint, climate long-term
mod2.ens.joint.long <- lm(joint_wtd_ens_longterm ~ region, data=dat)
summary(mod2.ens.joint.long)

# joint, climate short-term
mod2.ens.joint.short <- lm(joint_wtd_ens_shortterm ~ region, data=dat)
summary(mod2.ens.joint.short)

plot_lat_suit1 <- ggplot(dat, aes(lat, clim_wtd_ens_longterm)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  #xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlab("") + 
  xlim(43,45.5) +
  ylab("Climate (1980-2010)") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt")) 
plot_lat_suit2 <- ggplot(dat, aes(lat, clim_wtd_ens_shortterm)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  #xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlab("") + 
  xlim(43,45.5) +
  ylab("Climate (2014-15)") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt"))

plot_lat_suit3 <- ggplot(dat, aes(lat, stream_wtd_ens)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlim(43,45.5) +
  ylab("Physical habitat") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt"))

latsuit <- plot_grid(plot_lat_suit1 + theme(legend.position = "none"), 
                     plot_lat_suit2 + theme(legend.position = "none"), 
                     plot_lat_suit3 + theme(legend.position = "none"),
                     nrow=3)
legend <- get_legend(plot_lat_suit1)
latsuit2 <- plot_grid(legend, latsuit, ncol=1, rel_heights = c(0.3, 3))
left_label <- "ENM suitability, weighted ensemble"
latsuit3 <- ggdraw(latsuit2) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/WtdSuitability_vs_Latitude_3Panel.png", latsuit3, base_width=5, base_height=11)

plot_lat_suit4 <- ggplot(dat, aes(lat, joint_wtd_ens_longterm)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlim(43,45.5) +
  ylab("Climate + habitat ensemble") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt"))

plot_lat_suit5 <- ggplot(dat, aes(lat, joint_wtd_ens_shortterm)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlim(43,45.5) +
  ylab("Climate + habitat ensemble") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(0,6,0,40), "pt"))

# climate suitability (Amy March 2018)
mod0b.ens <- lm(lambda ~ 1, data=dat)
mod1b.ens <- lm(lambda ~ clim_wtd_ens_longterm, data=dat)
model.sel(mod0b.ens, mod1b.ens)
summary(mod1b.ens)

# climate suitability (Matt March 2018)
mod0b.ens.matt <- lm(lambda ~ 1, data=dat)
mod1b.ens.matt <- lm(lambda ~ clim_wtd_ens_shortterm, data=dat)
model.sel(mod0b.ens.matt, mod1b.ens.matt)
summary(mod1b.ens.matt)

# Stream suitability 
mod0b.ens.stream <- lm(lambda ~ 1, data=dat)
mod1b.ens.stream <- lm(lambda ~ stream_wtd_ens, data=dat)
model.sel(mod0b.ens.stream, mod1b.ens.stream)
summary(mod1b.ens.stream)

# Joint suitability, long-term climate
mod0b.joint.long <- lm(lambda ~ 1, data=dat)
mod1b.joint.long <- lm(lambda ~ joint_wtd_ens_longterm, data=dat)
model.sel(mod0b.joint.long, mod1b.joint.long)
summary(mod1b.joint.long)

# Joint suitability, short-term climate
mod0b.joint.short <- lm(lambda ~ 1, data=dat)
mod1b.joint.short <- lm(lambda ~ joint_wtd_ens_shortterm, data=dat)
model.sel(mod0b.joint.short, mod1b.joint.short)
summary(mod1b.joint.short)

# lambda vs ensemble suitability
plot_lam_suit1 <- ggplot(dat, aes(clim_wtd_ens_longterm, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_errorbarh(aes(xmin=clim_wtd_ens_longterm+clim_wtd_ens_longterm_sd, xmax=clim_wtd_ens_longterm-clim_wtd_ens_longterm_sd),height=0)+
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("Climate (1980-2010)") + 
  #xlim(0.7,2.05) + 
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(1.85)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(6,6,40,30), "pt")) 
#ggplot2::ggsave("Figures/Lambda_vs_ClimateENM_8010.png", width=8, height=8)

plot_lam_suit2 <- ggplot(dat, aes(clim_wtd_ens_shortterm, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_errorbarh(aes(xmin=clim_wtd_ens_shortterm+clim_wtd_ens_shortterm_sd, xmax=clim_wtd_ens_shortterm-clim_wtd_ens_shortterm_sd),height=0)+
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("Climate (2014-15)") + 
  #xlim(1.5,3.1) + 
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(1.85)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(6,6,40,30), "pt")) 
#ggplot2::ggsave("Figures/Lambda_vs_ClimateENM_1415.png", width=8, height=8)

plot_lam_suit3 <- ggplot(dat, aes(stream_wtd_ens, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_errorbarh(aes(xmin=stream_wtd_ens+stream_wtd_ens_sd, xmax=stream_wtd_ens-stream_wtd_ens_sd),height=0)+
  geom_smooth(method=lm, se=FALSE, color="black", linetype="dashed") + 
  xlab("Physical habitat") + 
  #xlim(0.1,1.0) + 
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(1.85)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(6,6,40,30), "pt")) 
#ggplot2::ggsave("Figures/Lambda_vs_StreamENM.png", width=8, height=8)

lamsuit <- plot_grid(plot_lam_suit2 + theme(legend.position = "none"), 
                     plot_lam_suit1 + theme(legend.position = "none"), 
                     plot_lam_suit3+ theme(legend.position = "none"), 
                     nrow=1, labels="AUTO", label_x=0.9)
legend <- get_legend(plot_lam_suit1)
lamsuit2 <- plot_grid(lamsuit, legend, rel_widths = c(3, 0.3))
bottom_label <- "ENM suitability, weighted ensemble"
lamsuit3 <- ggdraw(lamsuit2) + draw_label(bottom_label, angle=0, y=0.05, size=24)
save_plot("Figures/Lambda_vs_WtdSuitability_3Panel.png", lamsuit3, base_width=11, base_height=5)

plot_lam_suit4 <- ggplot(dat, aes(joint_wtd_ens_longterm, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_errorbarh(aes(xmin=joint_wtd_ens_longterm+joint_wtd_ens_longterm_sd, xmax=joint_wtd_ens_longterm-joint_wtd_ens_longterm_sd),height=0)+
  geom_smooth(method=lm, se=FALSE, color="black", linetype="dashed") + 
  xlab("Climate + habitat ensemble") + 
  #xlim(0.1,1.0) + 
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(1.85)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(6,6,40,30), "pt"))

plot_lam_suit5 <- ggplot(dat, aes(joint_wtd_ens_shortterm, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0) + 
  geom_errorbarh(aes(xmin=joint_wtd_ens_shortterm+joint_wtd_ens_shortterm_sd, xmax=joint_wtd_ens_shortterm-joint_wtd_ens_shortterm_sd),height=0)+
  geom_smooth(method=lm, se=FALSE, color="black", linetype="dashed") + 
  xlab("Climate + habitat ensemble") + 
  #xlim(0.1,1.0) + 
  #ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylab("") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(1.85)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), plot.margin = unit(c(6,6,40,6), "pt"))

lamsuit2 <- plot_grid(plot_lam_suit2 + theme(legend.position = "none"), 
                     plot_lam_suit1 + theme(legend.position = "none"), 
                     plot_lam_suit3 + theme(legend.position = "none"),
                     plot_lam_suit4 + theme(legend.position = "none"),
                     nrow=2, labels="AUTO", label_x=0.9)
legend <- get_legend(plot_lam_suit1)
lamsuit3 <- plot_grid(lamsuit2, legend, rel_widths = c(3, 0.5))
bottom_label <- "ENM suitability"
lamsuit4 <- ggdraw(lamsuit3) + draw_label(bottom_label, angle=0, y=0.05, size=24)
left_label <- expression(paste("Population growth rate (", lambda, ")"))
lamsuit5 <- ggdraw(lamsuit4) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/Lambda_vs_WtdSuitability_4Panel.png", lamsuit5, base_width=10, base_height=8)


## linear models of vitals vs suitability
# Amy climate
surv.suit <- lm(surv.siteint ~ clim_wtd_ens_longterm, data=dat2)
summary(surv.suit)

growth.suit <- lm(growth.siteint ~ clim_wtd_ens_longterm, data=dat2)
summary(growth.suit)

flower.suit <- lm(flowering.siteint ~ clim_wtd_ens_longterm, data=dat2)
summary(flower.suit)

fruit.suit <- lm(fruits.siteint ~ clim_wtd_ens_longterm, data=dat2)
summary(fruit.suit)

# Matt climate
surv.suit.matt <- lm(surv.siteint ~ clim_wtd_ens_shortterm, data=dat2)
summary(surv.suit.matt)

growth.suit.matt <- lm(growth.siteint ~ clim_wtd_ens_shortterm, data=dat2)
summary(growth.suit.matt)

flower.suit.matt <- lm(flowering.siteint ~ clim_wtd_ens_shortterm, data=dat2)
summary(flower.suit.matt)

fruit.suit.matt <- lm(fruits.siteint ~ clim_wtd_ens_shortterm, data=dat2)
summary(fruit.suit.matt)

# stream
surv.suit.stream <- lm(surv.siteint ~ stream_wtd_ens, data=dat2)
summary(surv.suit.stream)

growth.suit.stream <- lm(growth.siteint ~ stream_wtd_ens, data=dat2)
summary(growth.suit.stream)

flower.suit.stream <- lm(flowering.siteint ~ stream_wtd_ens, data=dat2)
summary(flower.suit.stream)

fruit.suit.stream <- lm(fruits.siteint ~ stream_wtd_ens, data=dat2)
summary(fruit.suit.stream)

# individual vitals
survsuit.amy <- ggplot(dat2, aes(clim_wtd_ens_longterm, surv.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black", linetype="dashed") +
  xlab("") + 
  xlim(0.7, 2.1) +
  #ylab("Survival") +
  ylab("") +
  ylim(-3, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.8)), plot.margin = unit(c(0,0,0,0), "pt"))
growthsuit.amy <- ggplot(dat2, aes(clim_wtd_ens_longterm, growth.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black") +
  xlab("") + 
  xlim(0.7, 2.1) +
  #ylab("Growth") +
  ylab("") +
  ylim(-2, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,0,0,0), "pt"))
flowersuit.amy <- ggplot(dat2, aes(clim_wtd_ens_longterm, flowering.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black", linetype="dashed") +
  xlab("") + 
  xlim(0.7, 2.1) +
  #ylab("Flowering") +
  ylab("") +
  ylim(-4, 1) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,0,0,0), "pt"))
fruitsuit.amy <- ggplot(dat2, aes(clim_wtd_ens_longterm, fruits.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("Climatic suitability \n(1980-2010)") + 
  xlim(0.7, 2.1) +
  ylab("") +
  #ylab("Fruits") +
  #ylim(-1, 0.3) + 
  scale_y_continuous(limit=c(-1.5,0.5), breaks=c(-1,0)) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,0,0,0), "pt"))

vitalsuit_long <- plot_grid(survsuit.amy + theme(legend.position = "none"), 
                            growthsuit.amy + theme(legend.position = "none"), 
                            flowersuit.amy + theme(legend.position = "none"), 
                            fruitsuit.amy + theme(legend.position = "none"), 
                            ncol=1, labels="AUTO", label_x=0.9)
legend <- get_legend(survsuit.amy)
vitalsuit_long2 <- plot_grid(legend, vitalsuit_long, ncol=1, rel_heights = c(0.2, 4))
left_label <- "Vital rate intercepts"
vitalsuit_long3 <- ggdraw(vitalsuit_long2) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/Vitals_vs_WtdSuit_Clim8010.png", vitalsuit_long3, base_width=5, base_height=11)

survsuit.matt <- ggplot(dat2, aes(clim_wtd_ens_shortterm, surv.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black", linetype="dashed") +
  xlab("") + 
  scale_x_continuous(limits=c(1.7, 2.8), breaks=seq(1.7, 2.8, by=0.5)) + 
  ylab("Survival") +
  ylim(-3, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.8)), plot.margin = unit(c(0,0,0,40), "pt"))
growthsuit.matt <- ggplot(dat2, aes(clim_wtd_ens_shortterm, growth.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, linetype="dashed", colour="black") +
  xlab("") + 
  scale_x_continuous(limits=c(1.7, 2.8), breaks=seq(1.7, 2.8, by=0.5)) + 
  ylab("Growth") +
  ylim(-2, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,0,0,40), "pt"))
flowersuit.matt <- ggplot(dat2, aes(clim_wtd_ens_shortterm, flowering.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black", linetype="dashed") +
  xlab("") + 
  scale_x_continuous(limits=c(1.7, 2.8), breaks=seq(1.7, 2.8, by=0.5)) + 
  ylab("Flowering") +
  ylim(-4, 1) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,0,0,40), "pt"))
fruitsuit.matt <- ggplot(dat2, aes(clim_wtd_ens_shortterm, fruits.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + 
  xlab("Climatic suitability \n(2014-2015)") + 
  scale_x_continuous(limits=c(1.7, 2.8), breaks=seq(1.7, 2.8, by=0.5)) + 
  ylab("Fruits") +
  #ylim(-1, 0.3) + 
  scale_y_continuous(limit=c(-1.5,0.5), breaks=c(-1,0)) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,0,0,40), "pt"))

vitalsuit_short <- plot_grid(survsuit.matt + theme(legend.position = "none"), 
                             growthsuit.matt + theme(legend.position = "none"), 
                             flowersuit.matt + theme(legend.position = "none"), 
                             fruitsuit.matt + theme(legend.position = "none"), 
                             ncol=1, labels="AUTO", label_x=0.9)
legend <- get_legend(survsuit.matt)
vitalsuit_short2 <- plot_grid(legend, vitalsuit_short, ncol=1, rel_heights = c(0.2, 4))
left_label <- "Vital rate intercepts"
vitalsuit_short3 <- ggdraw(vitalsuit_short2) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/Vitals_vs_WtdSuit_Clim1415.png", vitalsuit_short3, base_width=5, base_height=11)

vitalsuit_clim <- plot_grid(survsuit.matt + theme(legend.position="none"),
                            survsuit.amy + theme(legend.position="none"), 
                            growthsuit.matt + theme(legend.position="none"),
                            growthsuit.amy + theme(legend.position="none"),
                            flowersuit.matt + theme(legend.position="none"),
                            flowersuit.amy + theme(legend.position="none"),
                            fruitsuit.matt + theme(legend.position="none"),
                            fruitsuit.amy + theme(legend.position="none"),
                            ncol=2, nrow=4, labels="AUTO", label_x=0.9)
left_label <- "Vital rate intercept"
vitalsuit_clim2 <- ggdraw(vitalsuit_clim) + draw_label(left_label, angle=90, size=24, x=0.02)
save_plot("Figures/Vitals_vs_WtdSuit_Clim.png", vitalsuit_clim2, base_width=8.5, base_height=11)

survsuit.stream <- ggplot(dat2, aes(stream_wtd_ens, surv.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black", linetype="dashed") +
  xlab("") + 
  ylab("Survival") +
  #ylim(-2.5, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="top", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.8)), plot.margin = unit(c(0,10,0,40), "pt"))
growthsuit.stream <- ggplot(dat2, aes(stream_wtd_ens, growth.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, linetype="dashed", colour="black") +
  xlab("") + 
  ylab("Growth") +
  ylim(-2, 1.5) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,10,0,40), "pt"))
flowersuit.stream <- ggplot(dat2, aes(stream_wtd_ens, flowering.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black", linetype="dashed") +
  xlab("") + 
  ylab("Flowering") +
  ylim(-4, 1) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,10,0,40), "pt"))
fruitsuit.stream <- ggplot(dat2, aes(stream_wtd_ens, fruits.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + 
  xlab("Physical habitat suitability") + 
  ylab("Fruits") +
  ylim(-1, 0.3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none", plot.margin = unit(c(0,10,40,40), "pt"))

vitalsuit_stream <- plot_grid(survsuit.stream + theme(legend.position = "none"), 
                              growthsuit.stream + theme(legend.position = "none"), 
                              flowersuit.stream + theme(legend.position = "none"), 
                              fruitsuit.stream + theme(legend.position = "none"), 
                              ncol=1, labels="AUTO", label_x=0.9)
legend <- get_legend(survsuit.stream)
vitalsuit_stream2 <- plot_grid(legend, vitalsuit_stream, ncol=1, rel_heights = c(0.2, 4))
left_label <- "Vital rate intercepts"
vitalsuit_stream3 <- ggdraw(vitalsuit_stream2) + draw_label(left_label, angle=90, x=0.05, size=24)
save_plot("Figures/Vitals_vs_WtdSuit_Stream.png", vitalsuit_stream3, base_width=5, base_height=11)

