
Data <- local(get(load("/home/triffe/workspace/ThanoRepro/Data/DataAll.Rdata")))
source("/home/triffe/workspace/ThanoRepro/R/Functions.R")
source("/home/triffe/workspace/ThanoRepro/R/PlotFunctions.R")

# Produces Figures 1a and 1b... and some other fun stuff.


US2010 <- GetThanoMatrices(.Code = "USA", .Year = 2010, Data = Data)

percent.labs <- c("1.0%","0.8%","0.6%","0.4%","0.2%","","0.2%","0.4%","0.6%","0.8%","1.0%")

pdf("/home/triffe/workspace/ThanoRepro/Figures/US2010Age.pdf",width=4.8,height=4.5)
par(xaxs = "i", yaxs = "i", mai = c(.4,.4,.4,.8), xpd = TRUE)
plot(NULL, type = 'n', xlim = c(-.01,.01),ylim = c(0,111), axes = FALSE, xlab = "", ylab = "",
        panel.first = list(rect(-.01,0,.01,111,col = gray(.97),border = NA),
                           segments(-.01,seq(0,110,by=10),.01,seq(0,110,by=10),col="white",lwd=.8),
                           segments(-.01,seq(0,110,by=10),-.0102,seq(0,110,by=10),lwd=.8),
                           segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),111,col="white",lwd=.8),
                           segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),-1,lwd=.8),
                           text(seq(-.01,.01,by=.002),-1,percent.labs,pos=1,cex=.7),
                           text(-.0102,seq(0,110,by=10),seq(0,110,by=10),pos=2,cex=.7)
                           ))
PyramidOrLeafWithHeterogeneity(Males = US2010$Males, Females = US2010$Females, revcol = TRUE)
colorStrip(BrewerPal = "YlGnBu",w=.0015,h=80,lwd=.8,border=gray(.7),rising = TRUE,revcol=FALSE)
text(x=.01,95,"Remaining\nyears\ngroups",pos=4,cex=.8)
text(-.0125,117,"Age",pos=4,cex=.8)
dev.off()


pdf("/home/triffe/workspace/ThanoRepro/Figures/US2010Thano.pdf",width=4.8,height=4.5)
par(xaxs = "i", yaxs = "i", mai = c(.4,.4,.4,.8),xpd=TRUE)
plot(NULL, type = 'n', xlim = c(-.01,.01),ylim = c(0,111), axes = FALSE, xlab = "", ylab = "",
        panel.first = list(rect(-.01,0,.01,111,col = gray(.96),border = NA),
                segments(-.01,seq(0,110,by=10),.01,seq(0,110,by=10),col="white",lwd=.8),
                segments(-.01,seq(0,110,by=10),-.0102,seq(0,110,by=10),lwd=.8),
                segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),111,col="white",lwd=.8),
                segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),-1,lwd=.8),
                text(seq(-.01,.01,by=.002),-1,percent.labs,pos=1,cex=.7),
                text(-.0102,seq(0,110,by=10),seq(0,110,by=10),pos=2,cex=.7)
        ))

PyramidOrLeafWithHeterogeneity(Males = US2010$Males, Females = US2010$Females, 
        Pyramid =FALSE, revcol =TRUE, BrewerPal = "YlGn")
colorStrip(BrewerPal = "YlGn",w=.0017,h=80,lwd=.8,border=gray(.7), rising = FALSE)
text(x=.01,91,"Age groups",pos=4,cex=.8)
text(-.0125,117,"Remaining\nyears (y)",pos=4,cex=.8)
dev.off()












