setwd("/home/triffe/git/ThanoRepro/ThanoRepro")
source("R/Functions.R")
source("R/PlotFunctions.R")
devtools::load_all("/home/triffe/git/Leaves/PlosOne/R/RiffeetalFunctions")

plotPyr <- function(XXX,x,y,Thano=TRUE){
    N <- ncol(XXX[["Cm"]])
    mgrab <- ifelse(Thano,"Tm","Cm")
    fgrab <- ifelse(Thano, "Tf", "Cf")
    cgrab <- ifelse(Thano, "Tcol", "Ccol")
    for (i in N:1){
        m <- rep(XXX[[mgrab]][, i], each = 2)
        f <- rep(XXX[[fgrab]][, i], each = 2)
        polygon(c(m, rev(f)) + x, XXX[["y"]] + y, border = NA, col = XXX[[cgrab]][i])
    }
    
    PyramidOutline(abs(XXX[[mgrab]][, N]), XXX[[fgrab]][, N], 
            x = x, y = y, scale = 1, border = gray(.1), lwd = .2)
}

# OK data object in order
# transform to plot-specific data
prepare4plot <- function(XXX, colsFun1, colsFun2){
    Mf          <- XXX$Females
    Mm          <- XXX$Males
    TOT         <- sum(Mf) + sum(Mm)
    
    Het1f       <- apply(Mf ,2, aggN, N = 10) / TOT             # age with ey het
    Het1m       <- apply(Mm ,2, aggN, N = 10) / TOT             # age with ey het
    Het2f       <- apply(Mf ,1, aggN, N = 10) / TOT             # age with ey het
    Het2m       <- apply(Mm ,1, aggN, N = 10) / TOT 
    
# cumulative sum for cleaner plotting
    cumHet1f    <- t(apply(Het1f, 2, cumsum))
    cumHet2f    <- t(apply(Het2f, 2, cumsum))
    cumHet1m    <- -t(apply(Het1m, 2, cumsum))
    cumHet2m    <- -t(apply(Het2m, 2, cumsum))
    
    cols1       <- colsFun1(ncol(cumHet1m))
    cols2       <- colsFun2(ncol(cumHet1m))
    y           <- c(0, rep(1:110, each = 2), 111)
    y           <- c(y, rev(y))
    list(Cf = cumHet1f, Tf = cumHet2f, Cm = cumHet1m, Tm = cumHet2m, Ccol = cols1, Tcol = cols2, y = y)
}

colsFun1        <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"), space = "Lab")
colsFun2        <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "YlGn")), space = "Lab")


Data <- local(get(load("Data/DataAll.Rdata")))
US2010 <- GetThanoMatrices(.Code = "USA", .Year = 2010, Data = Data)
USMat <- prepare4plot(US2010,colsFun1,colsFun2)

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



US2010 <- GetThanoMatrices(.Code = "USA", .Year = 2010, Data = Data)
USMat <- prepare4plot(US2010,colsFun1,colsFun2)

postscript("/home/triffe/git/Leaves/PlosOne/Figures/Striking.eps",
        width=4.3,height=4.5, paper = "special", onefile = FALSE, bg = "white",
        horizontal = FALSE)
par(xaxs = "i", yaxs = "i", mai = c(.2,.2,.2,.2),xpd=TRUE)
plot(NULL, type = 'n', xlim = c(-.01,.01),ylim = c(0,111), axes = FALSE, xlab = "", ylab = "",
        panel.first = list(rect(-.01,0,.01,111,col = NA,border = gray(.8),lwd=.5,xpd=TRUE),
                segments(-.01,seq(0,110,by=10),.01,seq(0,110,by=10),col="white",lwd=.8),
                segments(-.01,seq(0,110,by=10),-.0102,seq(0,110,by=10),lwd=.8),
                segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),111,col="white",lwd=.8),
                segments(seq(-.01,.01,by=.002),0,seq(-.01,.01,by=.002),-1,lwd=.8),
                text(seq(-.01,.01,by=.002),-1,percent.labs,pos=1,cex=.7),
                text(-.0102,seq(0,110,by=10),seq(0,110,by=10),pos=2,cex=.7)
        )
)
plotPyr(USMat,0,0,TRUE)
colorStrip(BrewerPal = "YlGn",w=.0017,h=80,lwd=.8,border=gray(.7), rising = FALSE)
text(x=.01,91,"Years lived",pos=4,cex=.8)
text(-.0125,117,"Years left",pos=4,cex=.8)
dev.off()




















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
