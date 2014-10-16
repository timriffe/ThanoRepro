# for Tim, this will choke
if (system("hostname",intern=TRUE)=="triffe-N80Vm"){
  # if I'm on the laptop
  setwd("/home/tim/git/ThanoRepro/ThanoRepro")
} else {
  # in that case I'm on Berkeley system, and other people in the dept can run this too
  setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ThanoRepro/ThanoRepro"))
}

source("R/Functions.R")

# I finally got the discrete (pseudo continuous) equation to match
# perfectly, so now it seems relevant to try again with the matrix.
# this will require making some compromises on what seems to be best 
# practice for discretizing actual projection matrices so that things line
# up exactly. Let's see how far we can take this.

Data                        <- local(get(load("Data/DataAll.Rdata")))
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

#install.packages("matrixcalc")
library(matrixcalc)
is.singular.matrix(Y,tol=1e-8) # precision
is.singular.matrix(Y,tol=1e-9) # here it woks
is.singular.matrix(L,tol=1e-12)

Leig <- eigen(L)$values
Yeig <- eigen(Y)$values
cbind(Leig,Yeig)
rads <- seq(0,2*pi,length=1000)

pdf("LYeigen.pdf",height=4.5,width=9)
par(mfrow=c(1,2))
plot(NULL, xlim=c(-1,1),ylim=c(-1,1),asp=1,axes=FALSE, xlab = "Real", ylab = "Complex",
        main = "eigenvalues on the unit circle")
lines(cos(rads),sin(rads))
points(Re(Leig),Im(Leig),col = "#0000FF60",pch=19,cex=.5)
points(Re(Yeig),Im(Yeig),col = "#FF000060",pch=19,cex=.5)

plot(0:110, sqrt(Re(Yeig)^2+Im(Yeig)^2), type='l',col="red",ylim=c(0,1),
        main = "radius")
lines(0:110, sqrt(Re(Leig)^2+Im(Leig)^2), col="blue")
dev.off()
getwd()


#install.packages("expm", repos="http://R-Forge.R-project.org")
library(expm)
all(Y%^%100 > 0) # yes, thano is regular
all(L%^%100 > 0) # Leslie not regular?


# mixing time

# http://stanford.edu/class/ee363/lectures/pf.pdf
1/log(1/rev(sort(abs(Re(eigen(Y)$values))))[2])
1/log(1/rev(sort(abs(Re(eigen(L)$values))))[2])

# decreasing magnitude (from JHJ 2006 standford workshop notes)
Ldamping <- Mod(eigen(L)$values[2])/Mod(eigen(L)$values[1])  
Ydamping <- Mod(eigen(Y)$values[2])/Mod(eigen(Y)$values[1])  
Yperiod       <- 2*pi/ Arg(eigen(Y)$values[2])
Lperiod       <- 2*pi/ Arg(eigen(L)$values[2])

# approx: number of steps over which deviation from equilibrium
#distribution decreases by factor e
#
# 5th decimal. Can it be the close-out?
#
## same first eigenvalue(roughly), but different first eigenvectors...
## these are the stable age structures
#plot(Re(eigen(L)$vectors[, 1]),type='l')
#lines(abs(Re(eigen(Y)$vectors[, 1])),col="green")

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

# tring to matricize redistribution. need to stagger differently...
#da2
#M <- diag(da2)
#col(M)
#M <- matrix(da2[col(M)],ncol=ncol(M))
#M[lower.tri(M)] <- 0
#
#
#M <- matrix(0,ncol=length(da2),nrow=length(da2))
#M[upper.tri(M,TRUE)] <- 1
#day <- diag(1/la2) %*% M %*% diag(da2)
#plot(rowSums(day*ca2)-ca2)
#
#matplot(day %*% diag(ca2),type='l')

# try igraph stuff
#install.packages("igraph")
library(igraph)

YG <- graph.adjacency(t(Y), mode="directed", weighted=TRUE)
LG <- graph.adjacency(t(L), mode="directed", weighted=TRUE)

is.dag(YG)
is.dag(LG)

is.connected(YG)
is.connected(LG)

# is.chordal(YG) # causes crash
# is.chordal(LG)

#largest.cliques(YG)
# takes a long time...
#plot(YG)
plot(LG)
average.path.length(YG, directed=TRUE, unconnected=FALSE)
average.path.length(LG, directed=TRUE, unconnected=FALSE)

betweenness(YG)
betweenness(LG)

closeness(YG)
closeness(LG)

#plot(YG)
get.all.shortest.paths(LG)
# vertex 
graph.cohesion(YG)
graph.cohesion(LG) # ?
# edge connectivity
graph.adhesion(YG)
graph.adhesion(LG)

girth(YG)
girth(LG)

diameter(YG,directed=TRUE)
diameter(LG,directed=TRUE)

barplot(eccentricity(YG))
plot(log(evcent(YG,directed=TRUE)$vector))


transitivity(YG)
transitivity(LG)

#install.packages("rgl")
#rglplot(YG)