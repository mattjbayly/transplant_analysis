# Xpairs function
panel.line<- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
cex = 1, ...)
{
points(x, y, pch = pch, col = col, bg = bg, cex = cex)
ok <- is.finite(x) & is.finite(y)
if (any(ok))
md<-step(lm(y~poly(x,2)))
xx<- seq(min(x),max(x),length=100)
yy<-predict(md,data.frame(x=xx),se=T,type="response")
lines(xx,yy$fit,col=1)
lines(xx,yy$fit+2*yy$se.fit,col=3,lty=2)
lines(xx,yy$fit-2*yy$se.fit,col=2,lty=2)
}
panel.hist <- function(x, ...)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(usr[1:2], 0, 1.5) )
h <- hist(x, plot = FALSE)
breaks <- h$breaks; nB <- length(breaks)
y <- h$counts; y <- y/max(y)
rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
r <- abs(cor(x, y))
txt <- format(c(r, 0.123456789), digits=digits)[1]
txt <- paste(prefix, txt, sep="")
if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
text(0.5, 0.5, txt, cex = cex * r*2)
}
Xpairs <- function(...)pairs(...,lower.panel=panel.line, upper.panel=panel.cor,diag.panel=panel.hist)
