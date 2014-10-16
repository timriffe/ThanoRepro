# for Tim, this will choke
if (system("hostname",intern=TRUE)=="triffe-N80Vm"){
  # if I'm on the laptop
  setwd("/home/tim/git/ThanoRepro/ThanoRepro")
} else {
  # in that case I'm on Berkeley system, and other people in the dept can run this too
  setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ThanoRepro/ThanoRepro"))
}

source("R/Functions.R")
Data <- local(get(load("Data/DataAll.Rdata")))
Data <- Data[Data$Sex == "f", ]

# 
Dat <- Data[Data$Year == 2000 & Data$Code == "SWE",]
Mna0(Minf0(rowSums(Thano(Dat$Births, Dat$dx)) / 
                        rowSums(Thano(Dat$Exposure, Dat$dx))))   -

plot(sapply(seq(-.015,-.005,length.out=100),ThanoMin,da=Dat$dx, Fy=Dat$Fyf), type = 'l')

