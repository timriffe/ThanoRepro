source("/home/tim/git/ThanoRepro/ThanoRepro/R/Functions.R")
Data <- local(get(load("/home/tim/git/ThanoRepro/ThanoRepro/Data/DataAll.Rdata")))
Data <- Data[Data$Sex == "f", ]

# 
Dat <- Data[Data$Year == 2000 & Data$Code == "SWE",]
Mna0(Minf0(rowSums(Thano(Dat$Births, Dat$dx)) / 
                        rowSums(Thano(Dat$Exposure, Dat$dx))))   -

plot(sapply(seq(-.015,-.005,length.out=100),ThanoMin,da=Dat$dx, Fy=Dat$Fyf), type = 'l')

