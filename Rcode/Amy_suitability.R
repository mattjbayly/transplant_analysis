library(tidyverse)
library(MuMIn)
library(cowplot)

lams <- read_csv("Robjects/site.lambdas.bootstrap.csv")
head(lams)

suit <- read_csv("Data/site_preds_average.csv") %>% 
  filter(ID1=="trans") %>% droplevels() %>% 
  dplyr::select(ID2, LRavg, GAMavg, RFavg, BRTavg, MAXavg) %>% 
  mutate(Site = if_else(ID2=="Rock Creek", "ROCK", ifelse(ID2=="Coast Fork Willamette", "COAST", ifelse(ID2=="Mosby Creek", "MOSBY", ifelse(ID2=="Calapooia Creek", "CALAPOOIA", ifelse(ID2=="Wiley Creek", "WILEY", ifelse(ID2=="Thomas", "THOMAS", ifelse(ID2=="Hunter", "HUNTER", ifelse(ID2=="Looking Glass", "LOOK", NA)))))))))

dat <- left_join(lams, suit)
head(dat)

for (i in 1:dim(dat)[1]) {
  dat$Ens[i] = mean(c(dat$LRavg[i], dat$GAMavg[i], dat$RFavg[i], dat$BRTavg[i], dat$MAXavg[i]))
}

write.csv(dat, "Robjects/site.lambdas.suitability.csv")

# linear models of suitability vs latitude
mod0.ens <- lm(Ens ~ 1, data=dat)
mod1.ens <- lm(Ens ~ lat, data=dat)
model.sel(mod0.ens, mod1.ens)
summary(mod1.ens)

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

# anova of suitability by range position
mod2.ens <- lm(Ens ~ region, data=dat)
summary(mod2.ens)

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

ggplot(dat, aes(lat, Ens)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Latitude") + 
  xlim(43,45.5) +
  ylab("Climate ENM suitability") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2))) 
ggplot2::ggsave("Figures/ClimateENM_vs_Latitude.png", width=8, height=8)

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
  xlab("Latitude") + 
  xlim(43,45.5) +
  ylab("MAX suitability") +
  ylim(0,0.8) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.8)), legend.position="none") 

multi <- plot_grid(LRlat, GAMlat, RFlat, BRTlat, MAXlat, ncol=1, labels="AUTO", label_x=0.9)
save_plot("Figures/ClimateENM_vs_Latitude_IndModels.png", multi, base_width=5, base_height=11)


# linear models of lambda vs suitability 
mod0b.ens <- lm(lambda ~ 1, data=dat)
mod1b.ens <- lm(lambda ~ Ens, data=dat)
model.sel(mod0b.ens, mod1b.ens)
summary(mod1b.ens)

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
ggplot(dat, aes(Ens, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, color="black") + 
  xlab("Climate ENM suitability") + 
  ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="right", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2))) 
gplot2::ggsave("Figures/Lambda_vs_ClimateENM.png", width=8, height=8)

# individual models
LRlam <- ggplot(dat, aes(LRavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, aes(color="black")) + 
  xlab("Climate ENM suitability: LR") + 
  ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none") 
GAMlam <- ggplot(dat, aes(GAMavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, aes(color="black")) + 
  xlab("Climate ENM suitability: GAM") + 
  ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none") 
RFlam <- ggplot(dat, aes(RFavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", aes(color="black")) + 
  xlab("Climate ENM suitability: RF") + 
  ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none") 
BRTlam <- ggplot(dat, aes(BRTavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", aes(color="black")) + 
  xlab("Climate ENM suitability: BRT") + 
  ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none") 
MAXlam <- ggplot(dat, aes(MAXavg, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, aes(color="black")) + 
  xlab("Climate ENM suitability: MAX") + 
  ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none") 

multi2 <- plot_grid(LRlam, GAMlam, RFlam, BRTlam, MAXlam, labels="AUTO")
save_plot("Figures/Lambda_vs_ClimateENM_IndModels.png", multi2, base_width=8.5, base_height=11)


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

# individual models
survlat <- ggplot(dat2, aes(lat, surv.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("") + 
  xlim(43,45.5) +
  ylab("Survival intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none") 
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
  scale_color_manual(values=c("black", "black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, aes(color="black")) + 
  xlab("") + 
  xlim(43,45.5) +
  ylab("Flowering intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none")
fruitlat <- ggplot(dat2, aes(lat, fruits.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_manual(values=c("black", "grey")) +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Latitude") + 
  xlim(43,45.5) +
  ylab("Fruits intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none")

multi3 <- plot_grid(survlat, growthlat, flowerlat, fruitlat, nrow=4, ncol=1, labels="AUTO")
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
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Climatic suitability") + 
  ylab("Survival intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none") 
growthsuit <- ggplot(dat2, aes(Ens, growth.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method="lm", se=FALSE, aes(colour="black")) +
  xlab("Climatic suitability") + 
  ylab("Growth intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none")
flowersuit <- ggplot(dat2, aes(Ens, flowering.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Climatic suitability") + 
  ylab("Flowering intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none")
fruitsuit <- ggplot(dat2, aes(Ens, fruits.siteint)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  geom_smooth(method=lm, se=FALSE, aes(color="black")) + 
  xlab("Climatic suitability") + 
  ylab("Fruits intercept") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1)), axis.title=element_text(size=rel(1.2)), legend.position="none")

multi4 <- plot_grid(survsuit, growthsuit, flowersuit, fruitsuit, labels="AUTO")
save_plot("Figures/Vitals_vs_Suitability.png", multi4, base_width=8.5, base_height=11)

