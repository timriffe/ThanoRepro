
# note, Robert's neat dependency ratio parlor trick:

# P65+ / P18-64 = alpha
# how does alpha change with r?

Plotka <- function(Lx, r, a){
    sum(Lx * exp(-r*a))
}
xbarLotka <- function(Lx, r, a){
    sum(a * Lx * exp(-r*a)) /  sum(Lx * exp(-r*a))
}
Data <- local(get(load("/home/triffe/git/ThanoRepro/ThanoRepro/Data/DataAll.Rdata")))


Lx <- Data$Lx[with(Data, Code == "USA" & Year == 2000 & Sex == "m")]

age <- 0:110

r <- 0

Diffa <- function(Lx, r){
    xbarLotka(Lx[age >= 18 & age < 65], r, a=c(18.5:64.5)) -
            xbarLotka(Lx[age >= 65], r, a=c(65.5:110.5))
}
alphaLotka <- function(Lx,r){
    Plotka(Lx[age >= 65], r, a=c(65.5:110.5)) /
            Plotka(Lx[age >= 18 & age < 65], r, a=c(18.5:64.5))
    
}
r <- 0
DRLotka(Lx,r)

r <- .01
(alphaLotka(Lx,.01) - alphaLotka(Lx,0)) / alphaLotka(Lx,0) 
(alphaLotka(Lx,-.01) - alphaLotka(Lx,0)) / alphaLotka(Lx,0) 
Diffa(Lx,0.0) / 100
Diffa(Lx,.01) / 100
Diffa(Lx,-.01) / 100


