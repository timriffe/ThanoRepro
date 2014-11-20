# for Tim, this will choke
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm", "tim-ThinkPad-L440")){
  # if I'm on the laptop
  setwd("/home/tim/git/ThanoRepro/ThanoRepro")
} else {
  # in that case I'm on Berkeley system, and other people in the dept can run this too
  setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ThanoRepro/ThanoRepro"))
}

source("R/Functions.R")
#install.packages("data.table")
library(data.table)
# read in data and select females
Data <- local(get(load("Data/DataAll.Rdata")))
Data <- Data[Data$Sex == "f", ]
Data$Fy[is.na(Data$Fy)] <- 0
Data$Fyf[is.nan(Data$Fyf)] <- 0
DATA <- data.table(Data)

# functions to optimize r quickish
#  thanatological age 
rThanoCoale <- compiler::cmpfun(function(Fy, da, .a = .5:110.5, maxit = 1e2, tol = 1e-15){  
    
    # Based on Coale (1957), modified.
    N    <- length(Fy)
    dxM  <- matrix(0, ncol = N, nrow = N)
    dxi  <- da
    for (i in 1:N){
        dxM[i, 1:length(dxi)  ] <- dxi 
        dxi <- dxi[2:length(dxi) ]
    }     
    R0      <- sum(dxM * Fy)
    T.guess <- sum(.a * dxM * Fy) / R0 # assuming r = 0
    r.i     <- log(R0) / T.guess
    # be careful to discount Fy by SRB appropriately for males / females
    # prior to specification
    for (i in 1:maxit){ # 15 is more than enough!
        #cat(r2,i,"\n")
        deltai <- 1 - sum(rowSums(t(t(dxM) * exp(-r.i * .a))) * Fy)
        # the mean generation time self-corrects 
        # according to the error produced by the Lotka equation
        r.i <- r.i - (deltai / (T.guess - (deltai / r.i))) 
        if (abs(deltai) <= tol){
            break
        }
    }
    return(r.i)  
})
# chronological age
Rmomentn <- compiler::cmpfun(function(fx,Lx,x,n=0){
            sum((x^n)*fx*Lx)
        })
rLotkaCoale <- compiler::cmpfun(function(fx,Lx,x=.5:110.5, maxit = 1e2, tol = 1e-15){
            # from Coale, Ansley J. (1957) A New Method for Calculating Lotka's r- the Intrinsic Rate of Growth in a Stable Population.
            # Population Studies, Vol. 11 no. 1, pp 92-94
            R0 <- Rmomentn(fx,Lx,x,0)
            # first assuming a mean generation time of 29
            ri <- log(R0)/29
            for (i in 1:maxit){ # 10 is more than enough!
                deltai <- sum(exp(-ri*x)*fx*Lx)-1
                # the mean generation time self-corrects 
                # according to the error produced by the Lotka equation
                ri <- ri+(deltai/(29-(deltai/ri)))
                if (abs(deltai) <= tol){
                    break
                }
            }
            return(ri)  
        })
# calculate out r Thanatos and r Lotka
rT <- DATA[, rThanoCoale(Fyf, dx,maxit=1e3), by = list(Code,Year)]
rL <- DATA[, rLotkaCoale(Fxf, Lx), by = list(Code,Year)]
library(reshape2)
rTM <- acast(rT,Year~Code,value.var = "V1")
rLM <- acast(rL,Year~Code,value.var = "V1")
rDiff <- rLM - rTM

pdf("Figures/rDist.pdf")
par(xaxs="i",yaxs="i",mai=c(1,.5,.5,.5))
plot(density(rL$V1),xlim=c(-.04,.03),ylim=c(0,55),col= "blue", main = "",axes = FALSE, ylab = "",
        panel.first = list(rect(-.04,0,.03,55,col=gray(.96), border = NA),
                segments(seq(-.04,.03,by=.01),0,seq(-.04,.03,by=.01),55,col="white",lwd=.8),
                segments(-.04,seq(0,55,by=10),.03,seq(0,55,by=10),col="white",lwd=.8),
                text(seq(-.04,.03,by=.01),-1, seq(-.04,.03,by=.01),cex=.7,pos=1,xpd = TRUE),
                text(-.041,seq(0,55,by=10), seq(0,55,by=10),cex=.7,pos=2,xpd = TRUE),
                segments(seq(-.04,.03,by=.01),0,seq(-.04,.03,by=.01),-1,xpd = TRUE),
                segments(-.04,seq(0,55,by=10),-.041,seq(0,55,by=10),xpd = TRUE)))
lines(density(rT$V1),col="red")
abline(v=mean(rL$V1),col="blue",lty=2)
abline(v=mean(rT$V1),col="red",lty=2)
text(-.042, 57, "Density", xpd = TRUE)
text(0,-5, "intrinsic growth rate, r", xpd = TRUE)
legend("topright", lty = 1, col = c("blue","red"), legend = c("Lotka","Thanatological"), bty = "n")
dev.off()


var(rT$V1)/var(rL$V1)

par(mfrow=c(1,2))
matplot(1891:2011, rTM, type = 'l', lty = 1, col = "#0000FF40")
abline(h=0)
matplot(1891:2011, rLM, type = 'l', lty = 1, col = "#0000FF40")
abline(h=0)

graphics.off()
for (i in 1:ncol(rTM)){
plot(NULL, type = "n",xlim=c(1891,2011),ylim=c(-.04,.02),main=colnames(rLM)[i])
abline(h=0)
lines(1891:2011, rLM[,i], col = "blue")
lines(1891:2011, rTM[,i], col = "red")
legend("bottomleft",lty=1,col=c("red","blue"),legend=c("Thano","Lotka"))
Sys.sleep(1)
}

matplot(1891:2011, rLM - rTM, type = 'l', lty = 1, col = "#0000FF40")
abline(h=0)

sum(abs(rDiff) < 1e-7, na.rm = TRUE) # only one case of virtually identical r values (tol is like 1e-12)
sum(rLM < rTM, na.rm = TRUE)
sum(rLM > 0, na.rm = TRUE)
sum(rTM > 0, na.rm = TRUE)
# when rL greater than 0
rLg0 <- sum(rLM > 0, na.rm = TRUE)
sum(rLM > 0 & rLM < rTM & abs(rDiff) > 1e-7, na.rm = TRUE) / rLg0
sum(rLM > 0 & rLM > rTM & abs(rDiff) > 1e-7, na.rm = TRUE) / rLg0

# when rL less than 0
rLl0 <- sum(rLM < 0, na.rm = TRUE)
sum(rLM < 0 & rLM > rTM & abs(rDiff) > 1e-7, na.rm = TRUE) / rLl0
sum(rLM < 0 & rLM < rTM & abs(rDiff) > 1e-7, na.rm = TRUE) / rLl0

sum(sign(rLM) != sign(rTM), na.rm = TRUE) # 138 times where sign was different
rLMd <- apply(rLM,2,diff)
rLTd <- apply(rTM,2,diff)



# 186 times where direction of change was different
sum(sign(rLMd) != sign(rLTd), na.rm = TRUE) 
# 13 cases where both sign and direction were different
sum((sign(rLMd) != sign(rLTd) ) & (sign(rLM) != sign(rTM))[-1,], na.rm = TRUE) 


matplot(1891:2011, t(t(rDiff) / colMeans(rDiff,na.rm=TRUE)), type = 'l', lty = 1, col = "#0000FF40")
stDiff <- t(t(rDiff) / colMeans(rDiff,na.rm=TRUE))
countries <- colnames(rDiff)

for (i in 1:length(countries)){
matplot(1891:2011, stDiff, type = 'l', lty = 1, col = "#00000040", lwd = 2, main = countries[i])
lines(1891:2011,stDiff[,countries[i]])
Sys.sleep(1)
}

matplot(1891:2011, stDiff, type = 'l', lty = 1, col = "#00000040", lwd = 2, main = countries[i])
lines(1891:2011,stDiff[,"DEUTW"])
scale(rDiff)
matplot(1891:2011, scale(rDiff), type = 'l', lty = 1, col = rainbow(31), lwd = 2)
ncol(rDiff)


# ---------------------------------------
ay      <- 0:110
dx2      <- c(dx, dx * 0) / sum(dx) # pad out with 0s
dM     <- matrix(dx2[col(matrix(nrow = 111, ncol = 111)) + ay], 
        nrow = 111, ncol = 111, dimnames = list(Ex = ay, Age = ay))

Fyx <- Minf0(Mna0(Thano(Dat$Bxf, dx) / Thano(Dat$Exposure, dx)))
r <- 0
sum(dM*exp(-r*.5:110.5)*Fyx)
sum(t(dM)*exp(-r*.5:110.5)*t(Fyx))
# i.e. if fertility is not collapsed to Fy only, then we get same results from age and thano

sum(dM*exp(-r*.5:110.5)*Fyx)
sum(t(dM)*exp(-r*.5:110.5)*t(Fyx))











