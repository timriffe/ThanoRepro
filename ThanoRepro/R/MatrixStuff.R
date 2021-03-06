
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm", "tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/ThanoRepro/ThanoRepro")
} else {
	if (system("hostname",intern=TRUE) == "PC-403478"){
		# on MPIDR PC
		setwd("U://git//ThanoRepro//ThanoRepro")
	} else {
		# in that case I'm on Berkeley system, and other people in the dept can run this too
		setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ThanoRepro/ThanoRepro"))
	}
}
getwd()
# 
# install.packages("popbio")
library(popbio)
devtools::load_all(file.path("R","RiffeFunctions"))

Data <- local(get(load("Data/DataAll.Rdata")))
Data <- Data[Data$Sex == "f", ]

lambdaT <- local(get(load("Data/lambda.Rdata")))
lambdaL <- local(get(load("Data/lambdaLeslie.Rdata")))

makeMatrices <- FALSE
if (makeMatrices){
# made in commented-out code below:
#AllMatrices <- local(get(load("Data/AllMatrices.Rdata")))
#print(object.size(AllMatrices),units = "Mb")
	
AllMatrices <- lapply(split(Data, with(Data, paste(Code,Year))), function(Dat, lambdaL, lambdaT){
           yr <- unique(Dat$Year)
           Cd <- unique(Dat$Code)
           N  <- nrow(Dat)
           px <- Mna0(Dat$lx[2:N] / Dat$lx[1:(N-1)])
           lL <- lambdaL$lambda.f[lambdaL$Code == Cd & lambdaL$Year == yr]
           lT <- lambdaT$lambda.f[lambdaL$Code == Cd & lambdaL$Year == yr]
         
           list(Y = ThanoProjMatrix(Mna0(Dat$Fyf), Dat$dx, lT),
               L = Leslie(Mna0(Dat$Fxf), px, lL))
        },lambdaL = lambdaL, lambdaT = lambdaT)
# best to save this out, as it's quite big and take a while to calculate:
	save(AllMatrices, file="Data/AllMatrices.Rdata")
    print(object.size(AllMatrices),units = "Mb")
}
AllMatrices <- local(get(load(file.path("Data","AllMatrices.Rdata"))))
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
plot(Damping[,1],Damping[,2],type='l')
dim(Damping)
colnames(cbind(Damping, Code = XYZ, Year = YYY))
XYZ <- unlist(lapply(strsplit(rownames(Damping), split = " "),"[[",1))
YYY <- unlist(lapply(strsplit(rownames(Damping), split = " "),"[[",2))
Damping2 <- data.frame(Damping, Code = XYZ, Year = YYY)
TD <- reshape2::acast(Damping2, Code~Year, value.var = "DThano")
CD <- reshape2::acast(Damping2, Code~Year, value.var = "DLeslie")

plot(NULL, type = "n", xlim=c(1.04,1.09),ylim=c(1.01,1.06))
for (i in 1:nrow(TD)){
	points(TD[i,],CD[i,], col = "#00000030", pch=19, cex= seq(.1,1.5,length.out=ncol(TD)))
}
