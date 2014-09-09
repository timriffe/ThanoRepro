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

L <- rbind(fa, cbind(diag(la2[2:111] / la2[1:110]), 0))

# 2) force d(a) to entail l(a) and vice versa.
da2 <- -diff(c(la2,0))
a2      <- 0:110
#(rL2    <- rLotkaCoale(fa,la2,x=a2)
(rL2 <- log(Re(eigen(L)$values[1])))
# 5) recalc stable age and births using these 'clean' inputs
ca2     <- getca(la2,rL2,a2)

ca2 <- Re(eigen(L)$vectors[, 1])
ca2 <- ca2/sum(ca2)

ba2     <- fa*ca2 
# 6) and convert these to thanatological equivalents *without*
# staggering the d(a) vector.
.cy2    <- rowSums(Thano(ca2,da2,stagger=FALSE))
.by2    <- rowSums(Thano(ba2,da2,stagger=FALSE))
# these is the 'clean' stable thanatological fertility schedule
.fy2    <- .by2 / .cy2
# ok, now make simplified Leslie:

Y       <- outer(da2, .fy2, "*") + rbind(cbind(0, diag(110)), 0)

# wowsa
log(Re(eigen(L)$values[1]))
log(Re(eigen(Y)$values[1]))

## 5th decimal. Can it be the close-out?
#
## same first eigenvalue(roughly), but different first eigenvectors...
## these are the stable age structures
plot(Re(eigen(L)$vectors[, 1]),type='l')
lines(abs(Re(eigen(Y)$vectors[, 1])),col="green")

ca2L <- Re(eigen(L)$vectors[, 1])
ca2L <- ca2/sum(ca2)
#ca2L-ca2
#
## can I get it exactly on a 4x4 matrix?
#
#la3 <- c(1,.8,.4,.1)
#da3 <- c(.2,.4,.3,.1)
## fa3 <- c(0,1,.5,0)
#fa3 <- c(0,1,.6,0)
#a3 <- 0:3
#L3 <- rbind(fa3, cbind(diag(la3[2:4] / la3[1:3]), 0))
#log(Re(eigen(L3)$values[1]))
#
##(rL3    <- rLotkaCoale(fa3,la3,x=a3,T.guess=1.6,maxit = 1e3))
##sum(exp(-a3*rL3)*fa3*la3
#rL3      <- log(Re(eigen(L3)$values[1]))
#ca3      <- getca(la3,rL3,a3)
#ba3      <- fa3*ca3 
#cy3      <- rowSums(Thano(ca3,da3,stagger=FALSE))
#by3      <- rowSums(Thano(ba3,da3,stagger=FALSE))
## these is the 'clean' stable thanatological fertility schedule
#fy3      <- by3 / cy3
#Y3       <- outer(da3, fy3, "*") + rbind(cbind(0, diag(3)), 0)
#
## interesting ... both agree about 0 when fa3 = c(0,1,.5,0)
#log(Re(eigen(L3)$values[1]))
#log(Re(eigen(Y3)$values[1]))

da2
M <- diag(da2)
col(M)
M <- matrix(da2[col(M)],ncol=ncol(M))
M[lower.tri(M)] <- 0


M <- matrix(0,ncol=length(da2),nrow=length(da2))
M[upper.tri(M,TRUE)] <- 1
day <- diag(1/la2) %*% M %*% diag(da2)
plot(rowSums(day*ca2)-ca2)

matplot(day %*% diag(ca2),type='l')

# try igraph stuff
install.packages("igraph")



