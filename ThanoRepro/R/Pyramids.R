# for Tim, this will choke
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm", "tim-ThinkPad-L440")){
  # if I'm on the laptop
  setwd("/home/tim/git/ThanoRepro/ThanoRepro")
} else {
  # in that case I'm on Berkeley system, and other people in the dept can run this too
  setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ThanoRepro/ThanoRepro"))
}
source("R/Functions.R")
#source("R/PlotFunctions.R")
#
#
#devtools::load_all("/home/triffe/git/Leaves/PlosOne/R/RiffeetalFunctions")



colsFun1  <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"), space = "Lab")
colsFun2  <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "YlGn")), space = "Lab")


Data      <- local(get(load("Data/DataAll.Rdata")))
US2010    <- GetThanoMatrices(.Code = "USA", .Year = 2010, Data = Data)
USMat     <- prepare4plot(US2010,colsFun1,colsFun2)

percent.labs <- c("1.0%","0.8%","0.6%","0.4%","0.2%","","0.2%","0.4%","0.6%","0.8%","1.0%")

pdf("Figures/US2010Age.pdf",width=4.8,height=4.5)
par(xaxs = "i", yaxs = "i", mai = c(.4,.4,.4,.8), xpd = TRUE)
plot(NULL, type = 'n', xlim = c(-.01,.01),ylim = c(0,111), axes = FALSE, xlab = "", ylab = "",
        panel.first = list(rect(-.01,0,.01,111,col = gray(.94),border = gray(.8),lwd=.5,xpd=TRUE),
                segments(-.01,seq(0,110,by=10),.01,seq(0,110,by=10),col="white",lwd=.8),
                segments(-.01,seq(0,110,by=10),-.0102,seq(0,110,by=10),lwd=.8),
                segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),111,col="white",lwd=.8),
                segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),-1,lwd=.8),
                text(seq(-.01,.01,by=.002),-1,percent.labs,pos=1,cex=.7),
                text(-.0102,seq(0,110,by=10),seq(0,110,by=10),pos=2,cex=.7)
        ))
plotPyr(USMat,0,0,FALSE)
colorStrip(BrewerPal = "YlGnBu",w=.0015,h=80,lwd=.8,border=gray(.99),rising = TRUE,revcol=FALSE)
text(x=.01,91,"Years left",pos=4,cex=.8)
text(-.0125,117,"Years lived",pos=4,cex=.8)
dev.off()


pdf("Figures/US2010Thano.pdf",width=4.8,height=4.5)
par(xaxs = "i", yaxs = "i", mai = c(.4,.4,.4,.8),xpd=TRUE)
plot(NULL, type = 'n', xlim = c(-.01,.01),ylim = c(0,111), axes = FALSE, xlab = "", ylab = "",
        panel.first = list(rect(-.01,0,.01,111,col = gray(.94),border = gray(.8),lwd=.5,xpd=TRUE),
                segments(-.01,seq(0,110,by=10),.01,seq(0,110,by=10),col="white",lwd=.8),
                segments(-.01,seq(0,110,by=10),-.0102,seq(0,110,by=10),lwd=.8),
                segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),111,col="white",lwd=.8),
                segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),-1,lwd=.8),
                text(seq(-.01,.01,by=.002),-1,percent.labs,pos=1,cex=.7),
                text(-.0102,seq(0,110,by=10),seq(0,110,by=10),pos=2,cex=.7)
        ))
plotPyr(USMat,0,0,TRUE)
colorStrip(BrewerPal = "YlGn",w=.0017,h=80,lwd=.8,border=gray(.99), rising = FALSE)
text(x=.01,91,"Years lived",pos=4,cex=.8)
text(-.0125,117,"Years left",pos=4,cex=.8)
dev.off()



#US2010 <- GetThanoMatrices(.Code = "USA", .Year = 2010, Data = Data)
#USMat <- prepare4plot(US2010,colsFun1,colsFun2)
#
#postscript("/home/triffe/git/Leaves/PlosOne/Figures/Striking.eps",
#        width=4.3,height=4.5, paper = "special", onefile = FALSE, bg = "white",
#        horizontal = FALSE)
#par(xaxs = "i", yaxs = "i", mai = c(.2,.2,.2,.2),xpd=TRUE)
#plot(NULL, type = 'n', xlim = c(-.01,.01),ylim = c(0,111), axes = FALSE, xlab = "", ylab = "",
#        panel.first = list(rect(-.01,0,.01,111,col = NA,border = gray(.8),lwd=.5,xpd=TRUE),
#                segments(-.01,seq(0,110,by=10),.01,seq(0,110,by=10),col="white",lwd=.8),
#                segments(-.01,seq(0,110,by=10),-.0102,seq(0,110,by=10),lwd=.8),
#                segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),111,col="white",lwd=.8),
#                segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),-1,lwd=.8),
#                text(seq(-.01,.01,by=.002),-1,percent.labs,pos=1,cex=.7),
#                text(-.0102,seq(0,110,by=10),seq(0,110,by=10),pos=2,cex=.7)
#        )
#)
#plotPyr(USMat,0,0,TRUE)
#colorStrip(BrewerPal = "YlGn",w=.0017,h=80,lwd=.8,border=gray(.7), rising = FALSE)
#text(x=.01,91,"Years lived",pos=4,cex=.8)
#text(-.0125,117,"Years left",pos=4,cex=.8)
#dev.off()
#


#
#Data <- local(get(load("/home/triffe/workspace/ThanoRepro/Data/DataAll.Rdata")))
#source("/home/triffe/workspace/ThanoRepro/R/Functions.R")
#source("/home/triffe/workspace/ThanoRepro/R/PlotFunctions.R")
#
## Produces Figures 1a and 1b... and some other fun stuff.
#
#
#US2010 <- GetThanoMatrices(.Code = "USA", .Year = 2010, Data = Data)
#
#percent.labs <- c("1.0%","0.8%","0.6%","0.4%","0.2%","","0.2%","0.4%","0.6%","0.8%","1.0%")
#
#pdf("/home/triffe/workspace/ThanoRepro/Figures/US2010Age.pdf",width=4.8,height=4.5)
#par(xaxs = "i", yaxs = "i", mai = c(.4,.4,.4,.8), xpd = TRUE)
#plot(NULL, type = 'n', xlim = c(-.01,.01),ylim = c(0,111), axes = FALSE, xlab = "", ylab = "",
#        panel.first = list(rect(-.01,0,.01,111,col = gray(.97),border = NA),
#                           segments(-.01,seq(0,110,by=10),.01,seq(0,110,by=10),col="white",lwd=.8),
#                           segments(-.01,seq(0,110,by=10),-.0102,seq(0,110,by=10),lwd=.8),
#                           segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),111,col="white",lwd=.8),
#                           segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),-1,lwd=.8),
#                           text(seq(-.01,.01,by=.002),-1,percent.labs,pos=1,cex=.7),
#                           text(-.0102,seq(0,110,by=10),seq(0,110,by=10),pos=2,cex=.7)
#                           ))
#PyramidOrLeafWithHeterogeneity(Males = US2010$Males, Females = US2010$Females, revcol = TRUE)
#colorStrip(BrewerPal = "YlGnBu",w=.0015,h=80,lwd=.8,border=gray(.7),rising = TRUE,revcol=FALSE)
#text(x=.01,95,"Remaining\nyears\ngroups",pos=4,cex=.8)
#text(-.0125,117,"Age",pos=4,cex=.8)
#dev.off()
#
#
#pdf("/home/triffe/workspace/ThanoRepro/Figures/US2010Thano.pdf",width=4.8,height=4.5)
#par(xaxs = "i", yaxs = "i", mai = c(.4,.4,.4,.8),xpd=TRUE)
#plot(NULL, type = 'n', xlim = c(-.01,.01),ylim = c(0,111), axes = FALSE, xlab = "", ylab = "",
#        panel.first = list(rect(-.01,0,.01,111,col = gray(.96),border = NA),
#                segments(-.01,seq(0,110,by=10),.01,seq(0,110,by=10),col="white",lwd=.8),
#                segments(-.01,seq(0,110,by=10),-.0102,seq(0,110,by=10),lwd=.8),
#                segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),111,col="white",lwd=.8),
#                segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),-1,lwd=.8),
#                text(seq(-.01,.01,by=.002),-1,percent.labs,pos=1,cex=.7),
#                text(-.0102,seq(0,110,by=10),seq(0,110,by=10),pos=2,cex=.7)
#        ))
#
#PyramidOrLeafWithHeterogeneity(Males = US2010$Males, Females = US2010$Females, 
#        Pyramid =FALSE, revcol =TRUE, BrewerPal = "YlGn")
#colorStrip(BrewerPal = "YlGn",w=.0017,h=80,lwd=.8,border=gray(.7), rising = FALSE)
#text(x=.01,91,"Age groups",pos=4,cex=.8)
#text(-.0125,117,"Remaining\nyears (y)",pos=4,cex=.8)
#dev.off()
#
#
#
#
#
#
#
#
#
#
#
#
