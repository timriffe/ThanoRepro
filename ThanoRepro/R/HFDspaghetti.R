# for Tim, this will choke
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm", "tim-ThinkPad-L440")){
  # if I'm on the laptop
  setwd("/home/tim/git/ThanoRepro/ThanoRepro")
} else {
  # in that case I'm on Berkeley system, and other people in the dept can run this too
  setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ThanoRepro/ThanoRepro"))
}


source("R/Functions.R")
library(data.table)
# read in data and select females
Data <- local(get(load("Data/DataAll.Rdata")))
Data <- Data[Data$Sex == "f", ]
DT <- as.data.table(Data)


playWithData <- FALSE
if(playWithData){

# wow, the max age-specific rate is .29! WTF?
#max(Data$Fx)
#Data$Fx[Data$Fx==0] <- NA
plot(NULL, type="n",xlim = c(12,55),ylim=c(0,.3))
lines(Data$Age, Data$Fx, col = "#00000010")


# ----------------------------------------------

Data$Fx[Data$Fx == 0] <- NA
Data$Fy[Data$Fy == 0] <- NA

par(mfrow=c(1,2))
# -----------------------------------------------------
# all data
plot(NULL, type="n",xlim = c(0,95),ylim=c(0,.15))
lines(Data$Age, Data$Fy, col = "#00000005")

plot(NULL, type="n",xlim = c(12,55),ylim=c(0,.3))
lines(Data$Age, Data$Fx, col = "#FF000005")

# -----------------------------------------------------
# pre war
plot(NULL, type="n",xlim = c(12,55),ylim=c(0,.3))
lines(with(Data, Age[Year > 1920 & Year < 1939]), with(Data, Fx[Year > 1920 & Year < 1939]), col = "#00000050")

plot(NULL, type="n",xlim = c(0,95),ylim=c(0,.15))
lines(with(Data, Age[Year > 1920 & Year < 1939]), with(Data, Fy[Year > 1920 & Year < 1939]), col = "#00000050")

# -----------------------------------------------------
# post war to 1970
plot(NULL, type="n",xlim = c(12,55),ylim=c(0,.3))
lines(with(Data, Age[Year > 1945 & Year < 1970]), with(Data, Fx[Year > 1945 & Year < 1970]), col = "#00000010")

plot(NULL, type="n",xlim = c(0,95),ylim=c(0,.15))
lines(with(Data, Age[Year > 1945 & Year < 1970]), with(Data, Fy[Year > 1945 & Year < 1970]), col = "#00000010")

# -----------------------------------------------------
# after 1970
plot(NULL, type="n",xlim = c(12,55),ylim=c(0,.3))
lines(with(Data, Age[Year >= 1970]), with(Data, Fx[Year >= 1970]), col = "#00000008")

plot(NULL, type="n",xlim = c(0,95),ylim=c(0,.15))
lines(with(Data, Age[Year >= 1970]), with(Data, Fy[Year >= 1970]), col = "#00000008")


# -----------------------------------------------------
# after 1990
plot(NULL, type="n",xlim = c(12,55),ylim=c(0,.3))
lines(with(Data, Age[Year >= 1990]), with(Data, Fx[Year >= 1990]), col = "#00000008")

plot(NULL, type="n",xlim = c(0,95),ylim=c(0,.15))
lines(with(Data, Age[Year >= 1990]), with(Data, Fy[Year >= 1990]), col = "#00000008")

# -----------------------------------------------------
# after 2000
plot(NULL, type="n",xlim = c(12,55),ylim=c(0,.3))
lines(with(Data, Age[Year >= 2000]), with(Data, Fx[Year >= 2000]), col = "#00000008")

plot(NULL, type="n",xlim = c(0,95),ylim=c(0,.15))
lines(with(Data, Age[Year >= 2000]), with(Data, Fy[Year >= 2000]), col = "#00000008")


# 

st <- function(x){
    x/sum(x,na.rm=TRUE)
}

#update.packages()
DT[, FxSt := st(Fx), by = list(Code, Year)]
DT[, FySt := st(Fy), by = list(Code, Year)]
Data <- as.data.frame(DT)

par(mfrow=c(2,2))
plot(NULL, type="n",xlim = c(0,95),ylim=c(0,.15))
lines(Data$Age, Data$Fy, col = "#00000005")

plot(NULL, type="n",xlim = c(12,55),ylim=c(0,.3))
lines(Data$Age, Data$Fx, col = "#00000005")

plot(NULL, type="n",xlim = c(0,95),ylim=c(0,.11))
lines(Data$Age, Data$FySt, col = "#00000005")

plot(NULL, type="n",xlim = c(0,95),ylim=c(0,.11))
lines(Data$Age, Data$FxSt, col = "#00000005")
}

graphics.off()

# just for a cleaner figure
Data$Fy[Data$Age > 100] <- NA
# Fy 
png("Figures/FySpaghettiDraft.png")
#pdf("Figures/FySpaghetti.pdf") # high res vector
par(xaxs = "i", yaxs = "i", mai=c(.5,.5,.5,.5))
plot(NULL, type = "n",xlim = c(0,95),ylim=c(0,.13), axes = FALSE, xlab = "", ylab = "",
        panel.first = list(rect(0,0,95,.15,col = gray(.95),border = NA),
                segments(seq(0,95,by=10),0,seq(0,95,by=10),.13,col = "white",lwd=.5),
                segments(0,seq(0,.13,by=.02),95,seq(0,.13,by=.02),col = "white",lwd=.5),
                text(seq(0,95,by=10),0,seq(0,95,by=10),cex = .8,pos=1,xpd=TRUE),
                text(0,seq(0,.13,by=.02),seq(0,.13,by=.02),cex = .8,pos=2,xpd=TRUE)))
lines(Data$Age, Data$Fy, col = "#00000007")
text(-2,.1366,"Rate",xpd=TRUE)
text(mean(c(0,95)),-.01,"Thanatological age (years left)",xpd=TRUE)
#text(c(25,25,80,80),c(.03,.12,.03,.12),"DRAFT",cex=4,srt=45,col="#BBBBBB95",xpd=TRUE)
dev.off()

# just for a cleaner figure
Data$Fx[Data$Fx==0] <- NA
# Fx 
png("Figures/FxSpaghettiDraft.png")
#pdf("Figures/FxSpaghetti.pdf") # high res vector
par(xaxs = "i", yaxs = "i", mai=c(.5,.5,.5,.5))
plot(NULL, type = "n",xlim = c(12,55),ylim=c(0,.3), axes = FALSE, xlab = "", ylab = "",
        panel.first = list(rect(0,0,95,.3,col = gray(.95),border = NA),
                segments(seq(15,55,by=5),0,seq(15,55,by=5),.3,col = "white",lwd=.5),
                segments(12,seq(0,.3,by=.05),55,seq(0,.3,by=.05),col = "white",lwd=.5),
                text(seq(15,55,by=5),0,seq(15,55,by=5),cex = .8,pos=1,xpd=TRUE),
                text(12,seq(0,.3,by=.05),seq(0,.3,by=.05),cex = .8,pos=2,xpd=TRUE)))
lines(Data$Age, Data$Fx, col = "#00000007")
text(10,.315,"Rate",xpd=TRUE)
text(mean(c(12,55)),-.02,"Chronological age (years lived)",xpd=TRUE)
#text(c(25,25,45,45),c(.03,.25,.03,.25),"DRAFT",cex=4,srt=45,col="#BBBBBB95",xpd=TRUE)
dev.off()

length(unique(paste(Data$Code,Data$Year)))
#
## Fyst
#par(xaxs = "i", yaxs = "i", mai=c(1,1,1,1))
#plot(NULL, type = "n",xlim = c(0,95),ylim=c(0,.04), axes = FALSE, xlab = "", ylab = "",
#        panel.first = list(rect(0,0,95,.04,col = gray(.95),border = NA),
#                segments(seq(0,95,by=10),0,seq(0,95,by=10),.0,.04,col = "white",lwd=.5),
#                segments(0,seq(0,.04,by=.01),95,seq(0,.04,by=.01),col = "white",lwd=.5),
#                text(seq(0,95,by=10),0,seq(0,95,by=10),cex = .6,pos=1,xpd=TRUE),
#                text(0,seq(0,.01,by=.01),seq(0,.04,by=.01),cex = .6,pos=2,xpd=TRUE)))
#lines(Data$Age, Data$FySt, col = "#00000007")
#
#MRLB <- as.matrix(tapply(Data$Fy,list(Data$Code,Data$Year), function(fy){
#            sum(fy * (.5:110.5), na.rm = TRUE) / sum(fy, na.rm = TRUE)
#        }))
#
#ModalRLB <- as.matrix(tapply(Data$Fy,list(Data$Code,Data$Year), FindSmoothMode,x=0:110))
#
#
#matplot(t(ModalRLB), type= 'l', col = gray(.4), lty = 1, ylim = c(34,65))
#par(new=TRUE)
#matplot(t(MRLB), type= 'l', col = "#FF000050", lty = 1, axes =FALSE, ylab = "", xlab = "", ylim = c(34,70))
#
#plot(1891:2011, ModalRLB["SWE",], type = 'l', ylim = c(35,60))
#lines(1891:2011, MRLB["SWE",],col = "red")
#
#plot(1891:2011, ModalRLB["JPN",], type = 'l', ylim = c(35,65))
#lines(1891:2011, MRLB["JPN",],col = "red")
#
## difference between mean and mode
#matplot(t(ModalRLB)-t(MRLB), type= 'l', col = gray(.4), lty = 1)
#
#
## for caption
ranges <- tapply(Data$Year, Data$Code, function(yrs){
            paste0("$",paste(min(yrs),max(yrs),sep = "-"),"$")
        })
paste(paste(names(ranges),ranges, sep = ", "), collapse = "; ")
length(ranges)