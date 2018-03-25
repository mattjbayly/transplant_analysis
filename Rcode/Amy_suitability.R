library(tidyverse)
library(MuMIn)
library(cowplot)

# lambdas
lams <- read_csv("Robjects/site.lambdas.bootstrap.csv")
head(lams)

# suitability (redone by Amy March 2018)
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

# matthew's ENM predictions (redone March 2018)
matt.preds <- read_csv("data/Site Level ENM Preds.csv")

dat <- left_join(dat, matt.preds)
check.merge <- dat %>% dplyr::select(Site, site) # good

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
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none") 

plot_glm <- ggplot(dat, aes(LRavg, glm_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Amy climate GLM") + 
  ylab("Matt climate GLM") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="nonw") 

plot_gam <- ggplot(dat, aes(GAMavg, gam_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Amy climate GAM") + 
  ylab("Matt climate GAM") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none") 

plot_rf <- ggplot(dat, aes(RFavg, rf_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Amy climate RF") + 
  ylab("Matt climate RF") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none") 

plot_brt <- ggplot(dat, aes(BRTavg, brt_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Amy climate BRT") + 
  ylab("Matt climate BRT") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none") 

plot_max <- ggplot(dat, aes(MAXavg, max_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Amy climate MAX") + 
  ylab("Matt climate MAX") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none") 

model_comp <- plot_grid(plot_ens, plot_glm, plot_gam, plot_rf, plot_brt, plot_max)
save_plot("Figures/ClimateENMComparions1.png", model_comp, base_width=11, base_height=8.5)

plot(dat$mean_pred_climate, dat$mean_pred_stream)

### linear models of suitability vs latitude
# ave climate suitability, Amy's values
mod0.ens <- lm(Ens ~ 1, data=dat)
mod1.ens <- lm(Ens ~ lat, data=dat)
model.sel(mod0.ens, mod1.ens)
summary(mod1.ens)

# ave climate suitability, Matt's values
mod0.ens.matt <- lm(mean_pred_climate ~ 1, data=dat)
mod1.ens.matt <- lm(mean_pred_climate ~ lat, data=dat)
model.sel(mod0.ens.matt, mod1.ens.matt)
summary(mod1.ens.matt)

# ave stream suitability, Matt's values
mod0.ens.stream <- lm(mean_pred_stream ~ 1, data=dat)
mod1.ens.stream <- lm(mean_pred_stream ~ lat, data=dat)
model.sel(mod0.ens.stream, mod1.ens.stream)
summary(mod1.ens.stream)

# individual ENM algorithms
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

## anova of suitability by range position
mod2.ens <- lm(Ens ~ region, data=dat)
summary(mod2.ens)

mod2.ens.matt <- lm(mean_pred_climate ~ region, data=dat)
summary(mod2.ens.matt)

mod2.ens.stream <- lm(mean_pred_stream ~ region, data=dat)
summary(mod2.ens.stream)

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
  xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlim(43,45.5) +
  ylab("Climate ENM suitability") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2))) 
ggplot2::ggsave("Figures/ClimateENM_vs_Latitude.png", width=8, height=8)

plot_lat_suit2 <- ggplot(dat, aes(lat, mean_pred_climate)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlim(43,45.5) +
  ylab("Climate ENM suitability (Matt)") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)))

plot_lat_suit3 <- ggplot(dat, aes(lat, mean_pred_stream)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlim(43,45.5) +
  ylab("Stream ENM suitability") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)))

plot_grid(plot_lat_suit1, plot_lat_suit2)

# individual models
LRlat <- ggplot(dat, aes(lat, LRavg)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("LF suitability") +
  ylim(0,0.8) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="none") 
GAMlat <- ggplot(dat, aes(lat, GAMavg)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("GAM suitability") +
  ylim(0,0.8) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="none") 
RFlat <- ggplot(dat, aes(lat, RFavg)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("RF suitability") +
  ylim(0,0.8) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="none") 
BRTlat <- ggplot(dat, aes(lat, BRTavg)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("BRT suitability") +
  ylim(0,0.8) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="none") 
MAXlat <- ggplot(dat, aes(lat, MAXavg)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlim(43,45.5) +
  ylab("MAX suitability") +
  ylim(0,0.8) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="none") 

multi <- plot_grid(LRlat, GAMlat, RFlat, BRTlat, MAXlat, ncol=1, labels="AUTO", label_x=0.9)
save_plot("Figures/ClimateENM_vs_Latitude_IndModels.png", multi, base_width=5, base_height=11)


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

# Stream suitability (Amy March 2018)
mod0b.ens.stream <- lm(lambda ~ 1, data=dat)
mod1b.ens.stream <- lm(lambda ~ mean_pred_stream, data=dat)
model.sel(mod0b.ens.stream, mod1b.ens.stream)
summary(mod1b.ens.stream)

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

# lambda vs ensemble suitability
plot_lam_suit1 <- ggplot(dat, aes(Ens, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("Climatic suitability") + 
  ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2))) 
ggplot2::ggsave("Figures/Lambda_vs_ClimateENM.png", width=8, height=8)

plot_lam_suit2 <- ggplot(dat, aes(mean_pred_climate, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("Climatic suitability (Matt)") + 
  ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2))) 

plot_lam_suit3 <- ggplot(dat, aes(mean_pred_stream, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black", linetype="dashed") + 
  xlab("Stream suitability") + 
  ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2))) 

multi2 <- plot_grid(plot_lam_suit1, plot_lam_suit2)
multi4 <- plot_grid(plot_lat_suit1, plot_lat_suit2, plot_lam_suit1, plot_lam_suit2)
save_plot("Figures/ClimateENMComparions.png", multi4, base_width=8, base_height=8)

# individual models
LRlam <- ggplot(dat, aes(LRavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("LR suitability") + 
  xlim(0,0.7) +
  ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none") 
GAMlam <- ggplot(dat, aes(GAMavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("GAM suitability") + 
  xlim(0,0.7) +
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none") 
RFlam <- ggplot(dat, aes(RFavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + 
  xlab("RF suitability") + 
  xlim(0,0.7) +
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="top", legend.text=element_text(size=rel(1.2)), legend.title=element_text(size=rel(1.8))) 
BRTlam <- ggplot(dat, aes(BRTavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + 
  xlab("BRT suitability") + 
  xlim(0,0.8) +
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none") 
MAXlam <- ggplot(dat, aes(MAXavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("MAX suitability") + 
  xlim(0.2,0.5) +
  ylab("") +
  ylim(0,3) + 
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none") 

multi2 <- plot_grid(LRlam, GAMlam, RFlam, BRTlam, MAXlam, nrow=1, labels="AUTO", label_x=0.9)
save_plot("Figures/Lambda_vs_ClimateENM_IndModels.png", multi2, base_width=11, base_height=5)


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
  ylab("Survival intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="top", legend.text=element_text(size=rel(1.2)), legend.title=element_text(size=rel(1.8))) 
growthlat <- ggplot(dat2, aes(lat, growth.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("Growth intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none")
flowerlat <- ggplot(dat2, aes(lat, flowering.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("") + 
  xlim(43,45.5) +
  ylab("Flowering intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none")
fruitlat <- ggplot(dat2, aes(lat, fruits.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab(expression(paste("Latitude (", degree, "N)"))) + 
  xlim(43,45.5) +
  ylab("Fruits intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none")

multi3 <- plot_grid(survlat, growthlat, flowerlat, fruitlat, nrow=4, ncol=1, labels="AUTO", label_x=0.9)
save_plot("Figures/Vitals_vs_Latitude.png", multi3, base_width=5, base_height=11)

# linear models of vitals vs suitability
surv.suit <- lm(surv.siteint ~ Ens, data=dat2)
summary(surv.suit)

growth.suit <- lm(growth.siteint ~ Ens, data=dat2)
summary(growth.suit)

flower.suit <- lm(flowering.siteint ~ Ens, data=dat2)
summary(flower.suit)

fruit.suit <- lm(fruits.siteint ~ Ens, data=dat2)
summary(fruit.suit)

# individual models
survsuit <- ggplot(dat2, aes(Ens, surv.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  ylab("Survival intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="top", legend.text=element_text(size=rel(1.2)), legend.title=element_text(size=rel(1.8))) 
growthsuit <- ggplot(dat2, aes(Ens, growth.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, colour="black") +
  xlab("") + 
  ylab("Growth intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none")
flowersuit <- ggplot(dat2, aes(Ens, flowering.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  ylab("Flowering intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none")
fruitsuit <- ggplot(dat2, aes(Ens, fruits.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("Climatic suitability") + 
  ylab("Fruits intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none")

multi4 <- plot_grid(survsuit, growthsuit, flowersuit, fruitsuit, ncol=1, labels="AUTO", label_x=0.9)
save_plot("Figures/Vitals_vs_Suitability.png", multi4, base_width=5, base_height=11)

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
