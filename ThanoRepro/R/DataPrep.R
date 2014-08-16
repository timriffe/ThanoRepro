source("/home/triffe/workspace/ThanoRepro/R/Functions.R")
HMDpath <- "/home/triffe/DATA/HMD"
HFDpath <- "/home/triffe/DATA/HFD"
#list.files(HFDpath)

Bx <- read.table(file.path(HFDpath,"birthsRR.txt"), 
        skip = 2, as.is = TRUE, header = TRUE)
others <- c("IRL","UKR","ESP","BLR")
XXX <- others[1]
BxOTHER <- do.call(rbind, lapply(others, function(XXX, HFDpath){
                    Dat <- read.table(file.path(HFDpath,"OTHER",paste0(XXX,"birthsRR.txt")),skip=2,as.is=TRUE,header=TRUE)
                    cbind(Code = XXX, Dat)
                },HFDpath=HFDpath))
Bx <- rbind(Bx, BxOTHER)

Bx$Age <- as.integer(gsub(pattern = "\\+", replacement = "",(gsub(pattern = "\\-", replacement = "", Bx$Age))))
colnames(Bx)[4] <-"Births"

Ex <- read.table(file.path(HFDpath,"exposRR.txt"), 
        skip = 2, as.is = TRUE, header = TRUE)

ExOTHER <- do.call(rbind, lapply(others, function(XXX, HFDpath){
                    Dat <- read.table(file.path(HFDpath,"OTHER",paste0(XXX,"exposRR.txt")),skip=2,as.is=TRUE,header=TRUE)
                    cbind(Code = XXX, Dat)
                },HFDpath=HFDpath))
Ex <- rbind(Ex, ExOTHER)

HFDcountries <- unique(Ex$Code)

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
# takes a min
Bx       <- do.call(rbind, lapply(split(Bx,list(Bx$Code, Bx$Year)), HFDpad, colname = "Births"))
Ex       <- do.call(rbind, lapply(split(Ex,list(Ex$Code, Ex$Year)), HFDpad, colname = "Exposure"))

HFD      <- cbind(Bx, Ex)
HFD$Fx   <- HFD$Births /  HFD$Exposure
HFD$Fx[is.nan(HFD$Fx) | is.infinite(HFD$Fx)] <- 0

HFD$Sex  <- "f"
m        <- HFD
m$Sex    <- "m"
m$Births <- m$Exposure <- m$Fx <- NA
HFD      <- rbind(HFD, m)

# -------------------------------------------------
# now prep HMD data

#list.files(HMDpath)
HMDcountries <- unlist(lapply(list.files(file.path(HMDpath, "Births")), function(x){
                    strsplit(x, split = "\\.")[[1]][[1]]
                }))

Allcountries <- intersect(HFDcountries, HMDcountries)


# make similar long format for HMD data, Deaths, Exp, Pop Counts all smacked together. 
# Later we'll put on Fert data to make a single awesome object.

Data <- do.call(rbind,lapply(Allcountries, function(XXX, HMDpath){
           
            flt <- read.table(file.path(HMDpath,"lt_female/fltper_1x1", paste0(XXX, ".fltper_1x1.txt")),
                    skip = 2, header = TRUE, as.is = TRUE, na.strings = ".")
            mlt <- read.table(file.path(HMDpath,"lt_male/mltper_1x1", paste0(XXX, ".mltper_1x1.txt")),
                    skip = 2, header = TRUE, as.is = TRUE, na.strings = ".")
            pop <- read.table(file.path(HMDpath,"population/Population", paste0(XXX, ".Population.txt")),
                    skip = 2, header = TRUE, as.is = TRUE, na.strings = ".")
            Exp <- read.table(file.path(HMDpath,"exposures/Exposures_1x1", paste0(XXX, ".Exposures_1x1.txt")),
                    skip = 2, header = TRUE, as.is = TRUE, na.strings = ".")
            
            flt$Sex <- "f"
            mlt$Sex <- "m"
            flt$Code <- XXX
            mlt$Code <- XXX
           
            # Jan 1 pops
            pop$Year    <- as.character(pop$Year)
            yrs         <- sort(unique(gsub("\\+","",gsub("\\-","",pop$Year))))
            Nyr         <- length(yrs)
    
            pop1        <- pop[!grepl("\\-", pop$Year), ]
            pop2        <- pop[!grepl("\\+", pop$Year), ]
            pop1$Year   <- gsub("\\+","",pop1$Year)
            pop2$Year   <- gsub("\\-","",pop2$Year)
            pop1        <- pop1[pop1$Year != yrs[Nyr], ]
            pop2        <- pop2[pop2$Year != yrs[1], ]

            # we can trust the ordering here, standard HMD data and dimensions
            flt$Pop1    <- pop1$Female
            flt$Pop2    <- pop2$Female
            mlt$Pop1    <- pop1$Male
            mlt$Pop2    <- pop2$Male
            mlt$Exposure <- Exp$Male
            flt$Exposure <- Exp$Female
            # return as a single obj
            rbind(flt, mlt)
        }, HMDpath = HMDpath))

Data$Age <- as.integer(gsub("\\+","",Data$Age))
Data$Year <- as.integer(Data$Year)

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
save(Data,file = "/home/triffe/workspace/ThanoRepro/Data/DataAll.Rdata")

#unique(Data$Code)

# ------------------------------------------------------------
# now we need to read in deaths by Lexis triangle to calculate age 0 lambda for
# thanatological projection matrices.
dlpath <- "/home/triffe/DATA/HMD/deaths/Deaths_lexis"

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
save(lambda, file = "/home/triffe/workspace/ThanoRepro/Data/lambda.Rdata")

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
save(lambdaLeslie, file = "/home/triffe/workspace/ThanoRepro/Data/lambdaLeslie.Rdata")

