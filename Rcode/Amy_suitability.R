library(tidyverse)
library(MuMIn)

lams <- read_csv("Robjects/site.lambdas.bootstrap.csv")
head(lams)

suit <- read_csv("Data/site_preds_average.csv") %>% 
  filter(ID1=="trans") %>% droplevels() %>% 
  select(ID2, LRavg, GAMavg, RFavg, BRTavg, MAXavg) %>% 
  mutate(Site = if_else(ID2=="Rock Creek", "ROCK", ifelse(ID2=="Coast Fork Willamette", "COAST", ifelse(ID2=="Mosby Creek", "MOSBY", ifelse(ID2=="Calapooia Creek", "CALAPOOIA", ifelse(ID2=="Wiley Creek", "WILEY", ifelse(ID2=="Thomas", "THOMAS", ifelse(ID2=="Hunter", "HUNTER", ifelse(ID2=="Looking Glass", "LOOK", NA)))))))))

dat <- left_join(lams, suit)
head(dat)

for (i in 1:dim(dat)[1]) {
  dat$Ens[i] = mean(c(dat$LRavg[i], dat$GAMavg[i], dat$RFavg[i], dat$BRTavg[i], dat$MAXavg[i]))
}

write.csv(dat, "Robjects/site.lambdas.suitability.csv")

# linear models of lambda vs latitude
mod0 <- lm(Ens ~ 1, data=dat)
mod1 <- lm(Ens ~ lat, data=dat)
mod2 <- lm(Ens ~ poly(lat, 2), data=dat)
AIC(mod0, mod1, mod2)
model.sel(mod0, mod1, mod2)
summary(mod0)
summary(mod1)
summary(mod2)

# without northernmost site
mod0b <- lm(Ens ~ 1, data=dat[dat$lat<44.8,])
mod1b <- lm(Ens ~ lat, data=dat[dat$lat<44.8,])
mod2b <- lm(Ens ~ poly(lat, 2), data=dat[dat$lat<44.8,])
AIC(mod0b, mod1b, mod2b)
model.sel(mod0b, mod1b, mod2b)
summary(mod0b)
summary(mod1b)
summary(mod2b)

# anova of lambda by range position
mod1c <- lm(Ens ~ region, data=dat)
summary(mod1c)

# anova of lambda by range position without northernmost site
mod1d <- lm(Ens ~ region, data=dat[dat$lat<44.8,])
summary(mod1d)

ggplot(dat, aes(lat, Ens)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Latitude") + 
  ylab("Climate ENM suitability") +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none") 
ggsave("Figures/ClimateENM_vs_Latitude.png", width=8, height=8)

ggplot(dat, aes(Ens, lambda)) +
  geom_point(aes(colour=region), size=5) +
  scale_color_grey() +
  geom_point(shape=1, size=5, colour="black") +
  xlab("Climate ENM suitability") + 
  ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  theme_classic() + 
  theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(2)), legend.position="none") 
ggsave("Figures/Lambda_vs_ClimateENM.png", width=8, height=8)