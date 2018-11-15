############################################################################################
##  ORIGIONAL CODE FROM: 
##		Rees, M., Childs, D. Z., Ellner, S. P. (2014)
## 		Building integral projection models: a user's guide. Journal of Animal Ecology, 83: 528–545
##		doi: 10.1111/1365-2656.12178
##
##      M A I N     I P M    S C R I P T 
##		
##		Modified for Mimulus cardinalis dataset (April 2015)
###########################################################################################

## INDEX: 
	## Section 1 - Fit statistical models to simulated data and plot these where appropriate
		## 1 - SURVIVAL (binary indicator 'Surv' = 1 if survived)
		## 1.2 - growth (continuous variable 'z1')
		## 1.3A - reproduction (binary indicator 'Repr' = 1 if reproduced)
		## 1.3B - recruitment (binary indicator 'Recr' = 1 if recruitment)
		## 1.4 - recruit size (continuous variable 'Rcsz')




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Ungulate IBM to illustrate the construction of an IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(doBy)
require(car)
require(mgcv)

rm(list=ls(all=TRUE))

set.seed(270875)

## working directory must be set here, so the source()'s below run
	setwd("C:/Users/DW/Desktop/transplant_analysis/Planning_Docs/2.IPM_tutorials/Rees_2014_how to IPM/Reese example")

## run the utility functions
	source("./Standard Graphical Pars.R")

## run the ungulate IBM
	source("./Ungulate Demog Funs.R")
	cheat_sheet

## set the simulation parameters
	init.pop.size <- 500
	n.yrs <- 400
	m.par <- m.par.true # fixed values from "real population"

## run the ungulate IBM
	source("./Ungulate Simulate IBM.R")
	cat(pop.size,"\n")
	head(sim.data); dim(sim.data) # example dataset

## take a random sample of 3000 observations
## we need to take a random sample because new recruits are placed at the start of the data frame 
rows.to.keep <- sample.int(nrow(sim.data), 3000)
sim.data <- sim.data[rows.to.keep,]; dim(sim.data)
# summary(sim.data)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Fit statistical models to simulated data and plot these where appropriate
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## sort by size and print a sample to the screen
	sim.data <- sim.data[order(sim.data$z),]					# sort smallest -> largest (z) 
	round(sim.data[sample(1:nrow(sim.data),8),],2) 				# trim down extra decimals 

## set the plotting region
	## dev.new(6, 6)
	set_graph_pars(ptype = "panel4")							# custom plot par() arguments
	plot.range <- c(2.0, 3.5)

## 1.1 - SURVIVAL (binary indicator 'Surv' = 1 if survived)

surv.plot.data <- within(sim.data, {
    z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
    z.quantiles[1] <- z.quantiles[1]-0.1 				# make sure we include the smallest
    z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
    rm(z.quantiles)
})

# survival model
	mod.surv.glm <- glm(Surv ~ z  , family = binomial, data = surv.plot.data)
	summary(mod.surv.glm)

# calculate groupwise summary statistics (estimates across size range - 16 classes) 
		# (actual values for plotting vs. predicted) 
	surv.ps <- summaryBy(z + Surv ~ z.classes, data = surv.plot.data)
	surv.ps

	plot(Surv.mean ~ z.mean,
		 data = surv.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1),
		 xlab  =expression("Size , "*italic(z)), ylab = "Probability of surviving")

	# Could be fitted line from MEM
	lines(fitted(mod.surv.glm) ~ z, data = surv.plot.data, col = "red")
	add_panel_label(ltype="a") # A.) B.) C.) ... 

	
#----------------------------------------
## 1.2 - growth (continuous variable 'z1')

	# Growth (given survival), conditional on survival 
	grow.plot.data <- subset(sim.data, !is.na(z1))

	# Growth regression, simple linear 
	mod.grow <- lm(z1 ~ z, data = grow.plot.data)

	plot(z1 ~ z,
		 data = grow.plot.data,
		 xlim = plot.range, ylim = plot.range, pch = 16, cex = 0.25,
		 xlab = expression("Initial size, "*italic(z)),
		 ylab = expression("Final size, "*italic(z)*"'"))

	# fitted line from growth regression 
	abline(mod.grow, col="red")
	abline(0, 1, lty=2) # line of zero growth (not positive or negative), useful for negative growth at large size classes. 
	add_panel_label(ltype="b")

	
#----------------------------------------
## 1.3A - reproduction (binary indicator 'Repr' = 1 if reproduced)

# given survival & final size 
	repr.plot.data <- within(subset(sim.data, Surv==1), {
		z.quantiles <- quantile(z, probs = seq(0, 1, length = 16)) # make quantile classes of size range 
		z.quantiles[1] <- z.quantiles[1]-0.1 # make sure we include the smallest
		z.classes <- cut(z, z.quantiles)
		rm(z.quantiles)
	})

# pFlower linear regression - did the plant become reproductive (0/1) given survival & final size 
# IN REESE TUTORIAL STATE VARIABLE IS INITIAL SSIZE Z RATHER THAN FINAL SIZE Z' 
	mod.repr <- glm(Repr ~ z, family = binomial, data = repr.plot.data)

	repr.ps <- summaryBy(z + Repr ~ z.classes, data = repr.plot.data)
	repr.ps

	plot(Repr.mean ~ z.mean, data = repr.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1),
		 xlab=expression("Initial size, "*italic(z)), ylab="Probability of reproducing")

	lines(fitted(mod.repr) ~ z, data = repr.plot.data, col = "red")
	add_panel_label(ltype="c")


#----------------------------------------
## 1.3B - recruitment (binary indicator 'Recr' = 1 if recruitment)

	## < no plot > - spring-> summer recruitment (binary indicator 'Recr' = 1 if recruitment)
	summary(mod.recr <- glm(Recr ~ 1, family=binomial, data=sim.data))


#----------------------------------------
## 1.4 - recruit size (continuous variable 'Rcsz')

	# conditional on recruitment = 1, alternatively could use either Recr or Rcsz here. 
	cheat_sheet
	rcsz.plot.data <- subset(sim.data, !is.na(Rcsz))

	summary(mod.rcsz <- lm(Rcsz ~ z, data=rcsz.plot.data)) # simple linear relationship to initial size

	# would show strong evidence of maternal effect - but for cards won't be able to track identity of  
	# offspring. 
	plot(Rcsz ~ z,
		 data = rcsz.plot.data,
		 xlim = plot.range, ylim = plot.range, pch = 16, cex = 0.25,
		 xlab = expression("Maternal mass t, "*italic(z)),
		 ylab = expression("Offspring mass t+1, "*italic(z)*"'"))

abline(mod.rcsz, col="red")
add_panel_label(ltype="d")

## close the graphics device here 


#----------------------------------------
## Finally, store the estimated parameters
##

m.par.est <- c(## survival
               surv      = coef(mod.surv.glm),
               ## growth 
               grow      =  coef(mod.grow),
               grow.sd   =  summary(mod.grow)$sigma,
               ## reproduce or not
               repr      =  coef(mod.repr),
               ## recruit or not
               recr      =  coef(mod.recr),
               ## recruit size
               rcsz      =  coef(mod.rcsz),
               rcsz.sd   =  summary(mod.rcsz)$sigma)

names(m.par.est) <- names(m.par.true)
cheat_sheet
m.par.true # parameters stored for IPM 



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Construct kernels and compute population growth rates  
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	nBigMatrix <- 100

## size range - check the min/max sizes ever seen
	with(sim.data, min(z))
	with(sim.data, max(z))
## so let's set the lower and upper limits for the size range at 1.6 and 3.6 so slightly smaller/bigger
## than actually observed
# !!!! WILL NEED TO CHANGE THIS FOR CARDINALIS DATA ????
# !!!! IS THIS SUFFICENT TO ACCOUNT FOR EVICTION????
	min.size <- 1.6
	max.size <- 3.7

## make our true and estimated projection kernels
	IPM.true <- mk_K(nBigMatrix, m.par.true, min.size, max.size)
	IPM.est <- mk_K(nBigMatrix, m.par.est, min.size, max.size)
		str(IPM.true); class(IPM.true) # List with K vals, mesh pts & P & F kernels
		# example 
			head(IPM.est$meshpts) # mesh points
			
## calculate the population growth rate
	lam.est <- Re(eigen(IPM.est$K)$values[1])
	lam.est
	lam.true <- Re(eigen(IPM.true$K)$values[1])
	lam.true

## estimate the growth rate of the simulated population
	fit.pop.growth <- lm(log(pop.size.t) ~ seq.int(along = pop.size.t))
	fit.pop.growth
	exp(coef(fit.pop.growth)[2])

	
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 3 - Construct kernels, compute useful summary quantities and plot size by age distribution
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## extract the mesh points
	meshpts <- IPM.true$meshpts

## normalised stable size distribution
	# for estimate/ simulated data 
		w.est <- Re(eigen(IPM.est$K)$vectors[, 1])
		stable.z.dist.est <- w.est/sum(w.est)
		stable.z.repr.dist.est <- s_z(meshpts, m.par.est) * pb_z(meshpts, m.par.est) * w.est/sum(s_z(meshpts, 
				m.par.est) * pb_z(meshpts, m.par.est) * w.est)
	# for true fixed values 
		w.true <- Re(eigen(IPM.true$K)$vectors[, 1])
		stable.z.dist.true <- w.true/sum(w.true)
		stable.z.repr.dist.true <- s_z(meshpts, m.par.true) * pb_z(meshpts, m.par.true) * w.true/sum(s_z(meshpts, 
				m.par.true) * pb_z(meshpts, m.par.true) * w.true)

## mean size - log size scale (i.e. geometric mean taken on original scale size, then logged)
	mean.z.est <- sum(stable.z.dist.est * meshpts)
	mean.z.est
	mean.z.true <- sum(stable.z.dist.true * meshpts)
	mean.z.true

## Average size @ flowering
## mean size at reproduction - log size scale (i.e. geometric mean taken on original scale size, then logged)
	mean.z.repr.est <- sum(stable.z.repr.dist.est * meshpts)
	mean.z.repr.est
	mean.z.repr.true <- sum(stable.z.repr.dist.true * meshpts)
	mean.z.repr.true

## variance in size
	var.z.est <- sum(stable.z.dist.est * meshpts^2) - mean.z.est^2
	var.z.est
	var.z.true <- sum(stable.z.dist.true * meshpts^2) - mean.z.true^2
	var.z.true

## mean and variance in size - original untransformed scale
	mean.z.ari.est <- sum(stable.z.dist.est * exp(meshpts))
	mean.z.ari.est
	var.z.ari.est <- sum(stable.z.dist.true * exp(2 * meshpts)) - mean.z.ari.est^2
	var.z.ari.est

## compute the size distribution for each age class with repetitious code...
	a0.z.dist.est <- IPM.est$F %*% stable.z.dist.est/lam.est
	a1.z.dist.est <- IPM.est$P %*% a0.z.dist.est/lam.est
	a2.z.dist.est <- IPM.est$P %*% a1.z.dist.est/lam.est
	a3.z.dist.est <- IPM.est$P %*% a2.z.dist.est/lam.est

## ...or we can do this with a loop
	z.dist.by.age <- list()
	z.dist.by.age[[1]] <- IPM.est$F %*% stable.z.dist.est/lam.est
	for (i in 2:150) z.dist.by.age[[i]] <- IPM.est$P %*% z.dist.by.age[[i - 1]]/lam.est

## build a little helper function to compute the means & variances
		mk_moments <- function(z.dist, meshpts) {
			z.dist <- z.dist/sum(z.dist)
			mean.z <- sum(z.dist * meshpts)
			var.z <- sum(z.dist * meshpts^2) - mean.z^2
			return(c(mean = mean.z, sd = sqrt(var.z)))
		}

## Plot of expected size distribution in relation to age 
	set_graph_pars(ptype = "panel1")
	## age = 0 (new recruits)
	z.dist <- z.dist.by.age[[1]]
	plot(meshpts, z.dist, type = "n", xlab = "Mass, z", ylab = "Density", xlim = c(1.5, 4.2))
	## age = 1, 2, 3 and 4
		for (A in 0:4) {
			z.dist <- z.dist.by.age[[A + 1]]
			lines(meshpts, z.dist)
			moments <- round(mk_moments(z.dist, meshpts), 2)
			text(x = moments["mean"], y = max(z.dist) + 5e-05, pos = 4, cex = 0.75, labels = paste("A = ", 
				A, " (mean = ", moments["mean"], ", s.d. = ", moments["sd"], ")", sep = ""))
		}

		
		
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 4 - make plots comparing the behaviour of the IPM and the IBM 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dev.off(); set_graph_pars(ptype = "panel4")

## 1 - plot population density versus time...

	plot(seq.int(along = pop.size.t), log(pop.size.t), type = "l", xlab = "Time", ylab = "Population size")
	mod.pop.growth <- lm(log(pop.size.t) ~ seq.int(along = pop.size.t))
	abline(mod.pop.growth, col = "blue")

	#mtext(side = 3, line=0.5, adj = 0, text = "A)")
	add_panel_label(ltype = "a")
	## ...roughly linear for log(Nt) vs time so exponential growth


## 2 - plot mean size versus time...
	plot(seq.int(along = pop.size.t), mean.z.t[seq.int(along = pop.size.t)], type = "l", xlab = "Time", 
		ylab = "Mean sheep mass")
	abline(h = mean.z.true, col = "red")
	add_panel_label(ltype = "b")
	## ...mean size seems to settle down after a while

## 3 - plot mean flowering size versus time...
	plot(seq.int(along = pop.size.t), mean.z.repr.t[seq.int(along = pop.size.t)], type = "l", xlab = "Time", 
		ylab = "Mean reproductive mass")
	abline(h = mean.z.repr.true, col = "red")
	add_panel_label(ltype = "c")
	## ...mean flowering size seems to settle down after a while

## 4 - plot of density estimates at time 50 and the end 
	plot(density(sim.data$z), main = "", ylim = c(0, 3), xlab="Mass")
	#lines(density(sim.data$z))
	lines(IPM.est$meshpts, stable.z.dist.est/diff(IPM.est$meshpts)[1], col = "red")
	add_panel_label(ltype = "d")


	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 4 - Model diagnostics for the regressions underpinning the basic growth kernel and the survival
## function
## April 2014: Will need to revise this from Zurr for MEMs 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pick.data <- sample(seq(1, nrow(sim.data)))

#----------------------------------------
## GROWTH
	## Construct a data set of plausible size (subsampling is fudged so that we can see patterns discussed)
	strt <- 50
	end <- strt + 240 - 1
	test.data <- sim.data[pick.data[strt:end], ]
	test.data <- na.omit(subset(test.data, select = c(z, z1, Surv)))
	test.data <- test.data[order(test.data$z), ]
	summary(test.data); dim(test.data)
	
## refit models to the reduced data set 
	# basic growth model
		mod.grow <- lm(z1 ~ z, data = test.data)
		cat(length(mod.grow$fitted))

## setup the plot region for dignostics

dev.off(); set_graph_pars(ptype = "panel4")

## plot residuals versus fitted for growth model (mod.grow) 
	zhat <- fitted(mod.grow)
	resid <- residuals(mod.grow)
	plot(zhat, resid, xlab = "Fitted values", ylab = "Residuals")
	gam.resid <- gam(resid ~ s(zhat), method="REML")
	rhat <- predict(gam.resid, type = "response")
	points(zhat, rhat, type = "l")
	add_panel_label(ltype = "a")
	
## Normal qq-plot for growth model
	sresid <- rstandard(mod.grow)
	qqPlot(sresid, main = "", xlab = "Normal quantiles", ylab = "Residual quantiles",
		   col.lines = "black", lwd = 1)
	add_panel_label(ltype = "b")

## Absolute residuals versus fitted
	plot(zhat, sqrt(abs(sresid)), xlab = "Fitted values", ylab = "sqrt(|Residuals|)")
	gam.sresid <- gam(sqrt(abs(sresid)) ~ s(zhat), method = "REML")
	rhat <- predict(gam.sresid, type = "response")
	points(zhat, rhat, type = "l")
	add_panel_label(ltype = "c")

## compare to a gam fit
	gam.grow <- gam(z1 ~ s(z), data = test.data, method = "REML")
	gam.grow.fitted <- predict(gam.grow, type = "response")
	matplot(test.data$z, cbind(fitted(mod.grow), gam.grow.fitted), type = "l", lty = c(1, 2), lwd = 2, 
		xlab = "Size t", ylab = "Fitted size t+1", col=rep(1,2))
	add_panel_label(ltype = "d")

## close plot device

## miscellaneous checks of the growth function
	shapiro.test(sresid)
	cor.test(zhat, sqrt(abs(sresid)), method="k")
	AIC(gam.grow, mod.grow)
	anova(gam.grow, gam(z1 ~ z, data = test.data, method = "REML"), test="Chi")

#----------------------------------------	
## SURVIVAL functions

## Construct a data set of plausible size
	strt <- 50
	end <- strt + 340 - 1
	test.data <- sim.data[pick.data[strt:end], ]
	test.data <- na.omit(subset(test.data, select = c(z, Surv)))
	test.data <- test.data[order(test.data$z), ]
		summary(test.data); str(test.data); dim(test.data)
	
## setup the plot region for dignostics
dev.off(); set_graph_pars(ptype = "panel2")

# fit models to the reduced data set 
	mod.surv <- glm(Surv ~ z, family = binomial, data = test.data)
	gam.surv <- gam(Surv ~ s(z), family = binomial, data = test.data, method = "REML")

# compare fitted glm, grouped data, and fitted gam 
		sizerank <- seq_along(test.data$z)/nrow(test.data)
		test.data$sizeclass <- round(9 * sizerank) + 1 # make size classes based on quantiles
		surv.ps <- summaryBy(z + Surv ~ sizeclass, data = test.data, na.rm = TRUE)
		with(test.data, {
			plot(surv.ps$z.mean, surv.ps$Surv.mean, xlab = "Size z", ylab = "Survival probability",
				 pch = 1, xlim = range(z), ylim = c(0, 1))
			points(test.data$z, fitted(mod.surv), type = "l", lty = 1, lwd = 2)
			svals <- seq(min(z), max(z), length = 12)
			ghat <- predict(gam.surv, newdata = data.frame(z = svals), type = "response")
			points(svals, ghat, type = "p", col = "red", pch = 16)
		})
		add_panel_label(ltype = "a")

summary(gam.surv)

# Show the fitted GAM's linear predictor
	plot(gam.surv, seWithMean = TRUE, xlab = "Size z", ylab = "Spline(z, edf=1)")
	add_panel_label(ltype = "b")

## close plot device
