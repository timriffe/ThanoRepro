
# for Tim, this will choke
if (system("hostname",intern=TRUE)=="triffe-N80Vm"){
  # if I'm on the laptop
  setwd("/home/tim/git/ThanoRepro/ThanoRepro")
} else {
  # in that case I'm on Berkeley system, and other people in the dept can run this too
  setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ThanoRepro/ThanoRepro"))
}

source("R/Functions.R")


HMDpath <- "/home/tim/DATA/HMD"
HFDpath <- "/home/tim/DATA/HFD"
#list.files(HFDpath)

# IRL and ESP are preliminary countries in HFD at time of this writing
HFDcountries <- unique(c(getHFDcountries(),c("IRL","ESP")))

# This will take a LONG time to generate. Inefficient web parsing. FYI
# 1)
if (!"pw" %in% ls()){
  cat("enter HFD password into console (no quotes) and press enter\n")
  pw <- userInput()
}
if (!"us" %in% ls()){
  cat("enter HFD username into console (no quotes) and press enter\n")
  us <- userInput()
}

# HFD data downloaded on:
Sys.Date()
# save it to a little file so you don't forget.
cat("HFD data downloaded on", as.character(Sys.Date()), file = "Data/HFDdate.txt")
Bx <- do.call(rbind, lapply(HFDcountries, function(XXX, us, pw){
      Dat <- readHFDweb(CNTRY = XXX, item = "birthsRR", username = us, password = pw)
      data.frame(Code = XXX, Dat, stringsAsFactors = FALSE)
    }, us = us, pw = pw))

colnames(Bx)[colnames(Bx) == "Total"] <-"Births"

# same thing for HFD exposures (can be slightly different from HMD exposures)
Ex <- do.call(rbind, lapply(HFDcountries, function(XXX, us, pw){
      Dat <- readHFDweb(CNTRY = XXX, item = "exposRR", username = us, password = pw)
      data.frame(Code = XXX, Dat, stringsAsFactors = FALSE)
    }, us = us, pw = pw))


# pad out HFD data to make conformable with HMD data (ages 0 - 110)
HFDpad <- function(DATAyr, colname = "Total"){
    young       <- DATAyr[1:12, ]
    young$Age   <- 0:11
    young[,colname] <- 0
    
    old         <- rbind(young[1:11, ], DATAyr)
    old$Age     <- 56:(55+nrow(old))
    old[,colname]   <- 0
    
    fyr <- rbind(young, DATAyr, old)  
}

# add ages from 0 and up to 110
library(data.table)
Bx          <- data.table(Bx)
Bx          <- Bx[,HFDpad(.SD,"Births"),by=list(Code,Year)]
Ex          <- data.table(Ex)
Ex          <- Ex[,HFDpad(.SD,"Exposure"),by=list(Code,Year)]

# merge, get ASFR
Bx$Exposure <- Ex$Exposure
HFD         <- Bx

HFD$Fx      <- HFD$Births /  HFD$Exposure
HFD$Fx[is.nan(HFD$Fx) | is.infinite(HFD$Fx)] <- 0

HFD$Sex     <- "f"
m           <- HFD
m$Sex       <- "m"
m$Births    <- m$Exposure <- m$Fx <- NA
HFD         <- rbind(HFD, m)

# -------------------------------------------------
# now prep HMD data

#list.files(HMDpath)
HMDcountries <- getHMDcountries()

Allcountries <- intersect(HFDcountries, HMDcountries)


# make similar long format for HMD data, Deaths, Exp, Pop Counts all smacked together. 
# Later we'll put on Fert data to make a single awesome object.

Data <- do.call(rbind, lapply(Allcountries, function(XXX,pw,us){
      
      flt <- readHMDweb(CNTRY = XXX, item = "fltper_1x1", username = us, password = pw)
      mlt <- readHMDweb(CNTRY = XXX, item = "mltper_1x1", username = us, password = pw)
      pop <- readHMDweb(CNTRY = XXX, item = "Population", username = us, password = pw)
      Exp <- readHMDweb(CNTRY = XXX, item = "Exposures_1x1", username = us, password = pw)
 
      flt$Sex       <- "f"
      mlt$Sex       <- "m"
      flt$Code      <- XXX
      mlt$Code      <- XXX
      flt$Pop1      <- pop$Female1
      flt$Pop2      <- pop$Female2
      mlt$Pop1      <- pop$Male1
      mlt$Pop2      <- pop$Male2
      mlt$Exposure  <- Exp$Male
      flt$Exposure  <- Exp$Female
      
      # stick sexes together
      rbind(flt, mlt)
    }, us = us, pw = pw))

###################################################### 
# left off cleanup here. Continue cleaning in future #
######################################################



# ------------------------------------
# work to combine
HMDxxyr <- paste(Data$Code, Data$Year)
HFDxxyr <- paste(HFD$Code, HFD$Year)
XXXyr   <- intersect(unique(HMDxxyr), unique(HFDxxyr))

HFD  <- HFD[HFDxxyr %in% XXXyr, ]
Data <- Data[HMDxxyr %in% XXXyr, ]

HFD  <- HFD[with(HFD, order(Code, Year, Sex, Age)),]
Data <- Data[with(Data, order(Code, Year, Sex, Age)),]

# universal Data object!
Data <- cbind(Data, HFD[, c("Births","Fx")])

# OK, one last step, make Fxf, assuming constant SRB over age of mother:

PFall <- do.call(rbind,lapply(Allcountries, function(XXX, HMDpath){
                    BMF <- read.table(file.path(HMDpath, "Births",paste0(XXX,".Births.txt")),
                            header = TRUE, as.is = TRUE, skip = 2)
                    BMF$PF <- BMF$Female / BMF$Total
                    BMF$Code <- XXX
                    BMF[, c("Code","Year","PF")]
                }, HMDpath = HMDpath))
PFvec <- PFall$PF
names(PFvec) <- paste(PFall$Code, PFall$Year)
DataID <- paste(Data$Code, Data$Year)
Data$Fxf <- Data$Fx * PFvec[DataID]
Data$Bxf <- Data$Births * PFvec[DataID]
# since we use dx, lets standardize it:
Data <- do.call(rbind,lapply(split(Data, list(Data$Code, Data$Sex, Data$Year)), function(Dat){
                                    Dat$dx <- Dat$dx / sum(Dat$dx)
                                    Dat$lx <- Dat$lx / 1e5
                                    Dat$Lx <- Dat$Lx / 1e5
                                    Dat$Tx <- Dat$Tx / 1e5
                                    Dat
                                }))
FyFun <- compiler::cmpfun(function(Exposure, Births,dx){
    rowSums(Thano(Births, dx)) /
            rowSums(Thano(Exposure, dx))
})
FyfFun <- compiler::cmpfun(function(Exposure, Bxf,dx){
    rowSums(Thano(Bxf, dx)) /
            rowSums(Thano(Exposure, dx))
})
library(data.table)
DATA <- data.table(Data)
DATA[,Fy := FyFun(Exposure, Births,dx),by = list(Code,Sex,Year)]
DATA[,Fyf := FyFun(Exposure, Bxf,dx),by = list(Code,Sex,Year)]
Data <- as.data.frame(DATA)



#-----------------------------------------
rownames(Data) <- NULL
#-----------------------------------------
# save it out for wider use
save(Data,file = "Data/DataAll.Rdata")

#unique(Data$Code)

# ------------------------------------------------------------
# now we need to read in deaths by Lexis triangle to calculate age 0 lambda for
# thanatological projection matrices.
dlpath <- "/home/tim/DATA/HMD/deaths/Deaths_lexis"

lambda <- do.call(rbind, lapply(Allcountries, function(XXX, dlpath){
                    Dat <- read.table(file.path(dlpath,paste0(XXX,".Deaths_lexis.txt")),
                            skip = 2, header = TRUE, as.is = TRUE)
                    out <- list()
                    for (sex in c("Female","Male")){ # sex = "Female"
                        D0 <- reshape2::acast(Dat[Dat$Age == "0",], Year ~ Cohort, value.var = sex)
                        D0[D0==0] <- NA
                        D0 <- do.call(cbind,apply(D0,2,function(x){
                                            if (sum(!is.na(x))==2){
                                                x[!is.na(x)]
                                            }
                                        }))
                        lambda        <- D0[1,] / colSums(D0)
                        years         <- as.integer(names(lambda))
                        years         <- c(years,years[length(years)]+1)
                        lambda        <- c(lambda,lambda[length(lambda)])
                        names(lambda) <- years
                        out[[sex]]  <- lambda
                    }
                    data.frame(Code = XXX, 
                            Year = years, 
                            lambda.m = out[["Male"]], 
                            lambda.f = out[["Female"]], 
                            stringsAsFactors = FALSE)
                }, dlpath = dlpath))
save(lambda, file = "Data/lambda.Rdata")

# something similar for use in Leslie matrices:

lambdaLeslie <- do.call(rbind,lapply(Allcountries, function(XXX, HMDpath){
                    TLTU <- read.table(file.path(HMDpath,"deaths","Deaths_lexis",paste0(XXX,".Deaths_lexis.txt")),
                            skip = 2, header = TRUE, as.is = TRUE)
                    yrs <- sort(unique(TLTU$Year))
                    DLm <- with(TLTU, Male[Age == "0" & Year == Cohort])
                    DLf <- with(TLTU, Female[Age == "0" & Year == Cohort])
                    names(DLm) <- names(DLf) <- yrs
                    
                    BMF <- read.table(file.path(HMDpath, "Births",paste0(XXX,".Births.txt")),
                            header = TRUE, as.is = TRUE, skip = 2)
                    Bf <- BMF$Female[BMF$Year %in% yrs]
                    Bm <- BMF$Male[BMF$Year %in% yrs]
                    
                    lambda.m <- 1 - DLm / Bm
                    lambda.f <- 1 - DLf / Bf
                    data.frame(Code = XXX,Year = yrs, lambda.m,lambda.f,stringsAsFactors = FALSE)
                }, HMDpath = HMDpath))
save(lambdaLeslie, file = "Data/lambdaLeslie.Rdata")

