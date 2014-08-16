library(popbio)
setwd("/home/triffe/git/ThanoRepro/ThanoRepro")
# 
source("R/Functions.R")

Data <- local(get(load("Data/DataAll.Rdata")))
Data <- Data[Data$Sex == "f", ]

lambdaT <- local(get(load("Data/lambda.Rdata")))
lambdaL <- local(get(load("Data/lambdaLeslie.Rdata")))
# made in commented-out code below:
#AllMatrices <- local(get(load("Data/AllMatrices.Rdata")))
#print(object.size(AllMatrices),units = "Mb")
#AllMatrices <- lapply(split(Data, with(Data, paste(Code,Year))), function(Dat, lambdaL, lambdaT){
#           yr <- unique(Dat$Year)
#           Cd <- unique(Dat$Code)
#           N  <- nrow(Dat)
#           px <- Mna0(Dat$lx[2:N] / Dat$lx[1:(N-1)])
#           lL <- lambdaL$lambda.f[lambdaL$Code == Cd & lambdaL$Year == yr]
#           lT <- lambdaT$lambda.f[lambdaL$Code == Cd & lambdaL$Year == yr]
#         
#           list(Y = ThanoProjMatrix(Mna0(Dat$Fyf), Dat$dx, lT),
#               L = Leslie(Mna0(Dat$Fxf), px, lL))
#        },lambdaL = lambdaL, lambdaT = lambdaT)
# best to save this out, as it's quite big and take a while to calculate:
#save(AllMatrices, file="Data/AllMatrices.Rdata")
IDs <- names(AllMatrices)

# calculate Damping Ratios --- takes a while
# this is a recursive apply, because we have a nested list
Damping <- rapply(AllMatrices, function(X){
            eigen.analysis(X)$damping.ratio
        })
Damping <- matrix(Damping,ncol=2,byrow=TRUE,dimnames = list(IDs, c("DThano","DLeslie")))
# 100% of sample shows higher Damping ratio for thanatological matrix than for Leslie
sum(Damping[,1] > Damping[,2])

rMatrix <- rapply(AllMatrices, function(X){
            eigen.analysis(X)$lambda1
        })
rMatrix <- matrix(log(rMatrix),ncol=2,byrow=TRUE,dimnames = list(IDs, c("rThano","rLeslie")))
# or should it rMatrix - 1?

plot(density(Damping[,1]), xlim = c(1,1.1))
lines(density(Damping[,2]))

range(Damping[,1])
range(Damping[,2])
plot(density(Damping[,1]-Damping[,2]))

# what's the best way to show this?
par(mfrow=c(1,2))
plot(density((Damping[,1]-Damping[,2])/rowMeans(Damping)))
plot(density(Damping[,1]/Damping[,2]))

graphics.off()
plot(Damping[,1],Damping[,2])
# take a look at each country's trajectory
for (i in unique(lambdaT$Code)){
    plot(Damping[grepl(IDs,pattern=i),1],Damping[grepl(IDs,pattern=i),2],
            type='l',xlim=c(1,1.1),ylim=c(1,1.06),main=i)
    Sys.sleep(2)
}

