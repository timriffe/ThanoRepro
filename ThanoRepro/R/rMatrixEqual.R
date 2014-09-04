source("/home/tim/git/ThanoRepro/ThanoRepro/R/Functions.R")

# I finally got the discrete (pseudo continuous) equation to match
# perfectly, so now it seems relevant to try again with the matrix.
# this will require making some compromises on what seems to be best 
# practice for discretizing actual projection matrices so that things line
# up exactly. Let's see how far we can take this.

Data                        <- local(get(load("/home/tim/git/ThanoRepro/ThanoRepro/Data/DataAll.Rdata")))
Data                        <- Data[Data$Sex == "f", ]
Data$Fy[is.na(Data$Fy)]     <- 0
Data$Fyf[is.nan(Data$Fyf)]  <- 0
fa <- with(Data, Fxf[Year == 2000 & Code == "SWE"])
# 1) use l(a) instead of L(a)
# further make age strictly integer...
la2 <- with(Data, lx[Year == 2000 & Code == "SWE"])
# 2) force d(a) to entail l(a) and vice versa.
da2 <- -diff(c(la2,0))
a2      <- 0:110
(rL2    <- rLotkaCoale(fa,la2,x=a2))
# 5) recalc stable age and births using these 'clean' inputs
ca2     <- getca(la2,rL2,a2)
ba2     <- fa*ca2 
# 6) and convert these to thanatological equivalents *without*
# staggering the d(a) vector.
.cy2    <- rowSums(Thano(ca2,da2,stagger=FALSE))
.by2    <- rowSums(Thano(ba2,da2,stagger=FALSE))
# these is the 'clean' stable thanatological fertility schedule
.fy2    <- .by2 / .cy2
# ok, now make simplified Leslie:

L <- rbind(fa, cbind(diag(la2[2:111] / la2[1:110]), 0))
Y       <- outer(da2, .fy2, "*") + rbind(cbind(0, diag(110)), 0)

# 
log(Re(eigen(L)$values[1]))
log(Re(eigen(Y)$values[1]))

# 5th decimal. Can it be the close-out?

# same first eigenvalue(roughly), but different first eigenvectors...
# these are the stable age structures
plot(Re(eigen(L)$vectors[, 1]),type='l')
lines(abs(Re(eigen(Y)$vectors[, 1])),col="green