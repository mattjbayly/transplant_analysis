############################################################################################
##  ORIGIONAL CODE FROM: 
##		Rees, M., Childs, D. Z., Ellner, S. P. (2014)
## 		Building integral projection models: a user's guide. Journal of Animal Ecology, 83: 528–545
##		doi: 10.1111/1365-2656.12178
##
##      D E M O G R A P H Y    F U N C T I O N S
##		
##		Modified for Mimulus cardinalis dataset (April 2015)
############################################################################################

## IBM to generate the data for the simple ungulate example
## Size is on a log scale so don't worry about the -ive ones!

## The parameters and life cycle correspond to the Soay sheep of St Kilda, estimated from 26 years of data
## (1985-2010)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Define the demographic functions and parameters
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the true parameter vector, each parameter is given a name so the formulae, below, are easier to
## read. We'll use 'model.term' scheme to name elements of the vector
## the intercepts & slopes here are listed directly here as true values for the population 
## they are then inserted directly into the ipm, rather than estimates

m.par.true <- c(## survival
                surv.int  = -9.65e+0,
                surv.z    =  3.77e+0,
                ## growth 
                grow.int  =  1.41e+0,
                grow.z    =  5.57e-1,
                grow.sd   =  7.99e-2,
                ## reproduce or not
                repr.int  = -7.23e+0,
                repr.z    =  2.60e+0,                
                ## recruit or not
                recr.int  =  1.93e+0,
                ## recruit size
                rcsz.int  =  3.62e-1,
                rcsz.z    =  7.09e-1,
                rcsz.sd   =  1.59e-1)

				
##
## Define the various demographic functions, we pass a vector of parameters "m.par" to each function
## this makes it easier to compare the results of different paramter sets, say the true values and
## estimated ones
##

# GROWTH - PROBABILITY DENSITY FUNCTION (g_z1z)
## Growth function, given you are size z now returns the pdf (probability density function) of size z1 next time

g_z1z <- function(z1, z, m.par)
{
    mu <- m.par["grow.int"] + m.par["grow.z"] * z           # mean size next year (mu = intercept + B(initial size))
    sig <- m.par["grow.sd"]                                 # sd about mean (sig = stdev of growth)
    p.den.grow <- dnorm(z1, mean = mu, sd = sig)            # pdf that you are size z1 given you were size z
    return(p.den.grow)										# returns probability density function of growth 
}

# SURVIVAL - FUNCTION logit scale (s_z)
## Survival function, logistic regression

s_z <- function(z, m.par)
{
    linear.p <- m.par["surv.int"] + m.par["surv.z"] * z       # linear predictor (linear.p = intercept + B(initial size))
    p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability (p - logit scale) 
    return(p)												  # returns p on logit scale 
}

# FECUNDITY - FUNCTION logit scale (pb_z)
## Reproduction function, logistic regression

pb_z <- function(z, m.par)
{
    linear.p <- m.par["repr.int"] + m.par["repr.z"] * z       # linear predictor
    p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability
    return(p)												  # response 
}

# RECRUITMENT - function (just an intercept?)
## Recruitment function (N.B - from birth in spring to first summer), logistic regression

pr_z <- function(m.par)
{
    linear.p <- m.par["recr.int"]                             # linear predictor
    p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability
    return(p)
}

## OFFSPRING SIZE - PROBABILITY DENSITY FUNCTION (p.den.rcsz)
## Recruit size function

c_z1z <- function(z1, z, m.par)
{
    mu <- m.par["rcsz.int"] + m.par["rcsz.z"] * z           # mean size of recuits next year
    sig <- m.par["rcsz.sd"]                                 # sd about mean
    p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)            # pdf that offspring are size z1 given you were size z
    return(p.den.rcsz)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build IPM kernels P, F, and K
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## P kernel - survivorship & growth 
## F kernel - fecundity kernel
## K kernel - 

## Define the survival kernel (fxn * pdf)
	# P KERNEL - [size dependant survival, s(z)]	X	[G(z',z) probability density of size z′ at the next census for survivors having initial size z.]
	P_z1z <- function (z1, z, m.par) {
		# P(z', z)	=	s(z)	x	G(z',z)
		return( s_z(z, m.par) * g_z1z(z1, z, m.par) )
	}

## Define the reproduction kernel
F_z1z <- function (z1, z, m.par) {
	# F KERNEL - survive & grow to size z*, flower, have b(z*) offspring, some of which are size z*
		# need to sum over all possible sizes of z (upper to lower limit)
		# F(z',z) = s(z) integral[G(z*,z)pb(z*)b(z*)C0(z',z*)]dz*
    return( s_z(z, m.par) * pb_z(z, m.par) * (1/2) * pr_z(m.par) * c_z1z(z1, z, m.par) )
	# Survival probability * Flowering probability * 1/2 * Recruitment fxn * Offspring p.d.f
}


## Build the discretized kernel
# K - KERNEL  K(z',z) = P(z',z) + F(z',z)
mk_K <- function(m, m.par, L, U) {

	# mesh points 
	h <- (U - L)/m		# range 
	meshpts <- L + ((1:m) - 1/2) * h	
	P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
	F <- h * (outer(meshpts, meshpts, F_z1z, m.par = m.par))
	K <- P + F
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 3 - Make a cheat sheet for quick future reference
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my_notes <- c(
"Survival Int = surv.int",
"Survival Slope = surv.z",
"Growth Int = growt.int",
"Growth Slope = growt.z",
"Growth Stdev = growt.sd",
"Flower Int = repr.int",
"Flower Slope = repr.int",
"Seed recruit Int = recr.int",
"Recruit size Int = rcsz.int",
"Recruit size Slope = rcsz.z",
"Recruit size Stdev = rcsz.sd")

my_notes2 <- c(
"G(z',z) Growth p.d.f = g_z1z",
"s(z) Survival probability = s_z", 
"Flowering probability = pb_z",
"Recruitment function = pr_z",
"Offspring p.d.f = c_z1z"
)

my_notes3 <- c(
"P-kernel P(z',z) ~ s(z)G(z',z) == P_z1z",
"F-kernel F(z',z)  == F_z1z",
"K-kernel K(z',z) ~ P(z',z) + F(z',z) == mk_K"
)

my_notes4 <- c(
"IPM.est = estimate matrix"
)


cheat_sheet <- list(my_notes, my_notes2, my_notes3, my_notes4)
names(cheat_sheet) <- c("Parameters", "Vital Rate Functions", "Kernels", "R objects")
