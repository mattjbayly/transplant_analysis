library(tidyverse)
library(ggplot2)

Thomas.vitals <- read_csv("Robjects/vital_rate_coefficients.csv") %>% 
  filter(site=="THOMAS") %>% 
  droplevels()

Thomas.boot.vitals <- read_csv("Robjects/Mcard_transplant_BOOTSTRAP_params.csv") %>% 
  filter(Site=="THOMAS") %>% 
  droplevels()

ggplot(Thomas.boot.vitals, aes(surv.int)) +
  geom_histogram() + 
  geom_vline(data=Thomas.vitals, aes(xintercept=surv.int))
# real value is a bit right-shifted but not dramatically

ggplot(Thomas.boot.vitals, aes(growth.int)) +
  geom_histogram() + 
  geom_vline(data=Thomas.vitals, aes(xintercept=growth.int))
# real value is near center of bootstrapped distribution

ggplot(Thomas.boot.vitals, aes(flowering.int)) +
  geom_histogram() + 
  geom_vline(data=Thomas.vitals, aes(xintercept=flowering.int))
# real value is near right tail of bootstrapped distribution, which leads to observed lambda being higher than most bootstrap values.

ggplot(Thomas.boot.vitals, aes(fruits.int)) +
  geom_histogram() + 
  geom_vline(data=Thomas.vitals, aes(xintercept=fruits.int))
# real value is near center of bootstrapped distribution

ggplot(Thomas.boot.vitals, aes(seeds.per.fruit)) +
  geom_histogram() + 
  geom_vline(data=Thomas.vitals, aes(xintercept=seeds.per.fruit))
# real value is near center of bootstrapped distribution

ggplot(Thomas.boot.vitals, aes(establishment.prob)) +
  geom_histogram() + 
  geom_vline(data=Thomas.vitals, aes(xintercept=establishment.prob))
# real value is higher than tail of bootstrapped distribution, meaning observed is larger than most bootstrap replicates, which leads to observed lambda being greater than most bootstraps. this is a vital rate with high sensitivity.

ggplot(Thomas.boot.vitals, aes(recruit.logSize.mean)) +
  geom_histogram() + 
  geom_vline(data=Thomas.vitals, aes(xintercept=recruit.logSize.mean))
# real value is near center of bootstrapped distribution

ggplot(Thomas.boot.vitals, aes(recruit.logSize.sd)) +
  geom_histogram() + 
  geom_vline(data=Thomas.vitals, aes(xintercept=recruit.logSize.sd))
# real value is near center of bootstrapped distribution


# based on these, seems like establishment probability is the most likely culprit of the extreme mismatch between observed and bootstrapped lambdas. 

# flowering probability could make a secondary contribution.

# but establishment prob is fixed for all sites, is it having the same effect everywhere?

Coast.vitals <- read_csv("Robjects/vital_rate_coefficients.csv") %>% 
  filter(site=="COAST") %>% 
  droplevels()

Coast.boot.vitals <- read_csv("Robjects/Mcard_transplant_BOOTSTRAP_params.csv") %>% 
  filter(Site=="COAST") %>% 
  droplevels()

ggplot(Coast.boot.vitals, aes(surv.int)) +
  geom_histogram() + 
  geom_vline(data=Coast.vitals, aes(xintercept=surv.int))
# real value is a bit right-shifted but not dramatically

ggplot(Coast.boot.vitals, aes(growth.int)) +
  geom_histogram() + 
  geom_vline(data=Coast.vitals, aes(xintercept=growth.int))
# real value is near center of bootstrapped distribution

ggplot(Coast.boot.vitals, aes(flowering.int)) +
  geom_histogram() + 
  geom_vline(data=Coast.vitals, aes(xintercept=flowering.int))
# real value is near right tail of bootstrapped distribution, which leads to observed lambda being higher than most bootstrap values.

ggplot(Coast.boot.vitals, aes(fruits.int)) +
  geom_histogram() + 
  geom_vline(data=Coast.vitals, aes(xintercept=fruits.int))
# real value is near center of bootstrapped distribution

ggplot(Coast.boot.vitals, aes(seeds.per.fruit)) +
  geom_histogram() + 
  geom_vline(data=Coast.vitals, aes(xintercept=seeds.per.fruit))
# real value is near center of bootstrapped distribution

ggplot(Coast.boot.vitals, aes(establishment.prob)) +
  geom_histogram() + 
  geom_vline(data=Coast.vitals, aes(xintercept=establishment.prob))
# real value is higher than tail of bootstrapped distribution, meaning observed is larger than most bootstrap replicates, which leads to observed lambda being greater than most bootstraps. this is a vital rate with high sensitivity.

ggplot(Coast.boot.vitals, aes(recruit.logSize.mean)) +
  geom_histogram() + 
  geom_vline(data=Coast.vitals, aes(xintercept=recruit.logSize.mean))
# real value is near center of bootstrapped distribution

ggplot(Coast.boot.vitals, aes(recruit.logSize.sd)) +
  geom_histogram() + 
  geom_vline(data=Coast.vitals, aes(xintercept=recruit.logSize.sd))
# real value is near center of bootstrapped distribution

# ok, so this is a problem for all sites, but for some reason it is having an exaggerated effect at Thomas (mabye sensitivity is higher for its life history structure, with highest observed lambda)