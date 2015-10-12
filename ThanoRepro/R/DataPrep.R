
# for Tim, this will choke
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm", "tim-ThinkPad-L440")){
  # if I'm on the laptop
  setwd("/home/tim/git/ThanoRepro/ThanoRepro")
} else {
  # in that case I'm on Berkeley system, and other people in the dept can run this too
  setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ThanoRepro/ThanoRepro"))
}

source("R/Functions.R")
library(compiler)
library(data.table)
library(reshape2)
library(HMDHFDplus)
#list.files(HFDpath)

# IRL and ESP are preliminary countries in HFD at time of this writing
HFDcountries <- unique(c(getHFDcountries(),c("IRL","ESP")))
HMDcountries <- getHMDcountries()
Allcountries <- intersect(HFDcountries, HMDcountries)

# This will take a LONG time to generate. Inefficient web parsing. FYI
# 1)
if (!exists("pw")){
  cat("enter HFD password into console (no quotes) and press enter\n")
  pw <- userInput()
}

if (!exists("us")){
  cat("enter HFD username into console (no quotes) and press enter\n")
  us <- userInput()
}


# HFD data downloaded on:
Sys.Date()
# save it to a little file so you don't forget.
cat("HFD data downloaded on", as.character(Sys.Date()), file = "Data/HFDdate.txt")
Bx <- do.call(rbind, lapply(Allcountries, function(XXX, us, pw){
      Dat <- readHFDweb(CNTRY = XXX, item = "birthsRR", username = us, password = pw)
      data.frame(Code = XXX, Dat, stringsAsFactors = FALSE)
    }, us = us, pw = pw))

colnames(Bx)[colnames(Bx) == "Total"] <-"Births"

# same thing for HFD exposures (can be slightly different from HMD exposures)
Ex <- do.call(rbind, lapply(Allcountries, function(XXX, us, pw){
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

# make similar long format for HMD data, Deaths, Exp, Pop Counts all smacked together. 
# Later we'll put on Fert data to make a single awesome object.

# note when we download the data:
cat("HMD data downloaded on", as.character(Sys.Date()), file = "Data/HMDdate.txt")
HMD <- do.call(rbind, lapply(Allcountries, function(XXX,pw,us){
      
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

HMD <- as.data.frame(HMD)
HFD <- as.data.frame(HFD)

# ------------------------------------
# work to combine
# 1) find years/populations in common:
HMDxxyr <- paste(HMD$Code, HMD$Year)
HFDxxyr <- paste(HFD$Code, HFD$Year)
XXXyr   <- intersect(unique(HMDxxyr), unique(HFDxxyr))
#Bx[Bx$Year == 1956 & Bx$Code == "AUT" & Bx$Age > 10, ]
# cut down to only those pops with both fert and mort:
HFD     <- HFD[HFDxxyr %in% XXXyr, ]
HMD     <- HMD[HMDxxyr %in% XXXyr, ]

# reorder rows
HFD     <- HFD[with(HFD, order(Code, Year, Sex, Age)),]
HMD     <- HMD[with(HMD, order(Code, Year, Sex, Age)),]

# combine into universal Data object
Data    <- cbind(HMD, HFD[, c("Births","Fx")])

# OK, one last step, make Fxf, assuming constant SRB over age of mother:
PFall <- do.call(rbind,lapply(Allcountries, function(XXX, pw, us){
                    BMF <- readHMDweb(CNTRY = XXX, item = "Births", username = us, password = pw)
                    BMF$PF <- BMF$Female / BMF$Total
                    BMF$Code <- XXX
                    BMF[, c("Code","Year","PF")]
                }, pw = pw, us = us))
PFvec         <- PFall$PF
names(PFvec)  <- paste(PFall$Code, PFall$Year)

DataID        <- paste(Data$Code, Data$Year)
Data$Fxf      <- Data$Fx * PFvec[DataID]
Data$Bxf      <- Data$Births * PFvec[DataID]
# since we use dx, lets standardize it:

rescale <- function(x){
  x / sum(x, na.rm=TRUE)
}
reradix <- function(x){
  x / 1e5
}
Data <- data.table(Data)
Data$dx <- as.numeric(Data$dx) # data.table doesn't like converting integer to decimal...
Data$lx <- as.numeric(Data$lx)
Data$Lx <- as.numeric(Data$Lx)
Data$Tx <- as.numeric(Data$Tx)
Data[,dx := rescale(dx),by=list(Code,Sex,Year)]
Data[,lx := reradix(lx),by=list(Code,Sex,Year)]
Data[,Lx := reradix(Lx),by=list(Code,Sex,Year)]
Data[,Tx := reradix(Tx),by=list(Code,Sex,Year)]

# Now calculate thanatological fertility rates (asssuming fixed mort)                    
FyFun <- cmpfun(function(Exposure, Births,dx){
    rowSums(Thano(Births, dx)) /
            rowSums(Thano(Exposure, dx))
})
FyfFun <- cmpfun(function(Exposure, Bxf, dx){
    rowSums(Thano(Bxf, dx)) /
            rowSums(Thano(Exposure, dx))
})


Data[,Fy := FyFun(Exposure, Births, dx), by = list(Code,Sex,Year)]
Data[,Fyf := FyFun(Exposure, Bxf, dx), by = list(Code,Sex,Year)]
Data <- as.data.frame(Data)

head(Data[Data$Sex == "f" & Data$Age < 12,])
head(Data[Data$Sex == "f" & is.na(Data$Fyf),])
head(Data[Data$Sex == "f" & Data$Year == 1956 & Data$Code == "AUT",],15)
#-----------------------------------------
rownames(Data) <- NULL
#-----------------------------------------
# save it out for wider use
save(Data,file = "Data/DataAll.Rdata")

#unique(Data$Code)

# ------------------------------------------------------------
# now we need to read in deaths by Lexis triangle to calculate age 0 lambda for
# thanatological projection matrices. Lambda is a discount for those 
# births that die before making it to the end of the year
dlpath <- "/home/tim/DATA/HMD/deaths/Deaths_lexis"

lambda <- do.call(rbind, lapply(Allcountries, function(XXX, us, pw){
                    Dat <- readHMDweb(CNTRY = XXX, item = "Deaths_lexis", username = us, password = pw) 
                    years <- sort(unique(Dat$Year))
                    out <- list()
                    for (sex in c("Female","Male")){ # sex = "Female"
                        D0 <- acast(Dat[Dat$Age == "0",], Year ~ Cohort, value.var = sex)
                        
                        D0 <- rbind(D0[row(D0)+1==col(D0)],D0[row(D0)==col(D0)])
                       
                        lambda        <- D0[1,] / colSums(D0)
                        lambda[is.nan(lambda) | lambda == 0] <- mean(lambda[lambda!=0])
                        lambda        <- lambda
                        names(lambda) <- years
                        out[[sex]]    <- lambda
                    }
                    data.frame(Code = XXX, 
                            Year = years, 
                            lambda.m = out[["Male"]], 
                            lambda.f = out[["Female"]], 
                            stringsAsFactors = FALSE)
                }, pw = pw, us = us))
save(lambda, file = "Data/lambda.Rdata")

# something similar for use in Leslie matrices:
# for Leslie matrices lambda takes into account the loss of infants, but also the mortality of mothers.
lambdaLeslie <- do.call(rbind,lapply(Allcountries, function(XXX, us, pw){
                   
                    TLTU <- readHMDweb(CNTRY = XXX, item = "Deaths_lexis", username = us, password = pw)        
                          
                    yrs <- sort(unique(TLTU$Year))
                    DLm <- with(TLTU, Male[Age == "0" & Year == Cohort])
                    DLf <- with(TLTU, Female[Age == "0" & Year == Cohort])
                    names(DLm) <- names(DLf) <- yrs
                    
                    BMF <- readHMDweb(CNTRY = XXX, item = "Births", username = us, password = pw)
                    Bf  <- BMF$Female[BMF$Year %in% yrs]
                    Bm  <- BMF$Male[BMF$Year %in% yrs]
                    
                    lambda.m <- 1 - DLm / Bm
                    lambda.f <- 1 - DLf / Bf
                    data.frame(Code = XXX,Year = yrs, lambda.m,lambda.f,stringsAsFactors = FALSE)
                }, pw = pw, us = us))
save(lambdaLeslie, file = "Data/lambdaLeslie.Rdata")

