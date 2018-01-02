#######################################
setwd("C:/Users/DW/Desktop/transplant_analysis/Data/raw_data/WunderGround")
dir()


dat <- read.csv(file="StationData.csv")
head(dat)
str(dat)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}


levels(dat$Var)


temp <- dat[which(dat$Var=='Temperature'), ]
pre <- dat[which(dat$Var=='Precipitation'), ]

head(temp)
head(pre)

test <- strsplit(as.character(temp$high), " °C")
test <- unlist(test)
temp$high2 <- test
temp$high2 <- as.factor(temp$high2)
levels(temp$high2)[levels(temp$high2)=="°C"] <- "NA"
temp$high2 <- as.numeric.factor(temp$high2)





temp



as.numeric.factor(temp$high2) 
 
 
"°C"
temp
 
16.1 °C
70.4 mm


x <- c(as = "asfef", qu = "qwerty", "yuiop[", "b", "stuff.blah.yech")
# split x on the letter e
strsplit(x, "e")
Var    high      low     avg

