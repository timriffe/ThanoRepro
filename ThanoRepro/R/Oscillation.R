# for Tim, this will choke
if (system("hostname",intern=TRUE)=="triffe-N80Vm"){
  # if I'm on the laptop
  setwd("/home/tim/git/ThanoRepro/ThanoRepro")
} else {
  # in that case I'm on Berkeley system, and other people in the dept can run this too
  setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ThanoRepro/ThanoRepro"))
}
# seems like oscillation comparisons should be done on standard populations.

# could start with Px, get Py, and then derive Cohen's measure, but this adds an initial smoothing
# to Py. Instead I'll start with identical but randomly generated starting populations.
AllMatrices <- local(get(load("Data/AllMatrices.Rdata")))

CohenA <- compiler::cmpfun(function(P0, LorY, N=500){
    Pops <- matrix(nrow=length(P0), ncol = N)
    Pops[,1] <- P0
    for (i in 2:N){
        Pops[,i] <- c(LorY %*% Pops[, i - 1])
    }
    # make sum to 1
    Pops <- t(t(Pops) / colSums(Pops)) 
    sum(abs(Pops - Pops[,N]))/2
})

P0 <- runif(111)
CohenA(P0, AllMatrices[[1]][[2]])

library(parallel)
P0mat <- matrix(runif(111*1e3),ncol=1e3)


apply(P0mat,2,function(P0, AY){
       CohenA(P0, AY[["L"]]) / CohenA(P0, AY[["Y"]])
        }, AY = AllMatrices[[1]])

extremeCase1 <- c(1,rep(0,110))
CohenA(extremeCase1, AllMatrices[[1]][["L"]]) / CohenA(extremeCase1, AllMatrices[[1]][["Y"]])
CohenA(extremeCase1, AllMatrices[[1]][["L"]]) / CohenA(rev(extremeCase1), AllMatrices[[1]][["Y"]])

extremeCase2 <- c(rep(0,15),1,rep(0,95))
CohenA(extremeCase2, AllMatrices[[1]][["L"]]) / CohenA(extremeCase2, AllMatrices[[1]][["Y"]])
CohenA(extremeCase2, AllMatrices[[1]][["L"]]) / CohenA(rev(extremeCase2), AllMatrices[[1]][["Y"]])

extremeCase3 <- rep(1,111)
CohenA(extremeCase3, AllMatrices[[1]][["L"]]) / CohenA(extremeCase3, AllMatrices[[1]][["Y"]])


Y <- AllMatrices[[1]][["Y"]]
L <- AllMatrices[[1]][["L"]]
plot(Y %*% (Y %*% extremeCase3))
plot(L %*% (L %*% extremeCase3))






