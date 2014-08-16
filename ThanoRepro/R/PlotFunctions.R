
makeVectorAgeGroups <- function(vec, ages = 0:110, N = 2){
    x.new  <- ages - ages %% N
    c(tapply(vec, x.new, sum))
}
# draw rect based on x,y of cetroid, plus width and height
makeRect <- function(x, y, w, h, ...){
    w2 <- w / 2
    h2 <- h / 2
    rect(x - w2, y - h2, x + w2, y + h2, ...)
}

PyramidOutline2 <- function(males, females, 
        scale = sum(c(males, females)), 
        x = 0, y = 0, ...){
    N       <- length(males)
    Total   <- sum(c(males, females), na.rm = TRUE)
    widths  <- rep(1, N)
    age     <- c(0,cumsum(widths)[-N])
    u.age   <- age[N] + widths[N]
    
    males   <- scale * males / Total
    females <- scale * females / Total
    
    polygon(x = c(0, rep(females, each = 2) + 0,0, rev(c(-0, rep(-males, each = 2) - 0, -0))) + x, 
            y =  c(rep(c(age, u.age), each = 2), rev(c(rep(c(age, u.age), each = 2)))) + y, ...)
}

# Males, Females should be the full thano x chrono population matrix, as returned by Thano()
# Pyramid for age-structure. FALSE for thanatological age structure instead
# N interval width for remaining-years heterogeneity (age pyramid) or age heterogeneity (thano leaf)
# BrewerPal the colo ramp from color brewer package.
# revcol you can decide direction of ramp explicitly
# x,y for positioning of pyramid
# ycex, for vertical scaling.
# ---------------------
# *only does pyramids / leaves scaled to N  <- 1
PyramidOrLeafWithHeterogeneity <- function(Males, Females, Pyramid = TRUE,
        N = 10, BrewerPal = "YlGnBu", revcol = FALSE,x=0,y=0){
    # define color functions
    if (revcol){
        colR        <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, BrewerPal)), space = "Lab")
    } else {
        colR        <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, BrewerPal), space = "Lab")
    }
    # if we're plotting a thanatological leaf
    if (!Pyramid){
        Males   <- t(Males)
        Females <- t(Females)
    }
   
    # reduces number of polygons we need to draw. 
    # *reduces file size if saving + renders better in pdf
    ages    <- 1:nrow(Males) - 1
    Males   <- apply(Males, 2, makeVectorAgeGroups, ages = ages, N = N)
    Females <- apply(Females, 2, makeVectorAgeGroups, ages = ages, N = N)
    
    # scale
    Tot     <- sum(Males) + sum(Females)
    Males   <- -Males / Tot
    Females <- Females / Tot
    
    # this if-statement is experimental
    # has to do with direction of shading gradient
    if (Pyramid){
        Males       <- Males[nrow(Males):1,]
        Females     <- Females[nrow(Females):1,]
    }
    
# first column = age 0
# take cumsum down the columns, step 1 in determining centroids
    MalesCa                 <- apply(Males, 2, cumsum)
    FemalesCa               <- apply(Females, 2, cumsum)
    
    Males[Males == 0]       <- NA
    Females[Females == 0]   <- NA
    MindNA                  <- !is.na(Males)
    FindNA                  <- !is.na(Females)
    
    MalesCa[!MindNA]        <- NA
    FemalesCa[!FindNA]      <- NA

# centroid x vals will be this minus 1/2 or orig values
    MalesCnta   <- MalesCa - Males / 2
    FemalesCnta <- FemalesCa - Females / 2
    
    # err, hard-coded to single-age data. Could be arg, would need to check rest
    heights     <- 1
# the prop values are now the widths
    
    x1m         <- MalesCnta[MindNA]
    y1m         <- col(Males)[MindNA] - .5 # midpoints
    
    x1f         <- FemalesCnta[FindNA]
    y1f         <- col(Females)[FindNA] - .5 # midpoints
    
    # work out colors
    colsm       <- colR(nrow(Males))[row(Males)] 
    dim(colsm)  <- dim(Males)
    colsm       <- colsm[MindNA]

    colsf       <- colR(nrow(Females))[row(Females)] 
    dim(colsf)  <- dim(Females)
    colsf       <- colsf[FindNA]
    
    # widths of rectangles
    wm          <- Males[MindNA]
    wf          <- Females[FindNA]
# % labels
    makeRect(x1m+x,y1m+y,wm,1,col = colsm, border = colsm, lwd = .5, xpd = TRUE)
    makeRect(x1f+x,y1f+y,wf,1,col = colsf, border = colsf, lwd = .5, xpd = TRUE)
    PyramidOutline2(-colSums(Males, na.rm=TRUE),colSums(Females, na.rm=TRUE),scale=1,border=gray(.1), x=x, y=y, lwd = .2)
    #segments(x,y,x,111+y,col="white",lwd=.2)
}

# where Data is the universal huge data object for this paper
GetThanoMatrices <- function(.Code = "SWE", .Year = 1891, Data){
    Dat     <- Data[with(Data, Code == .Code & Year == .Year), ]
    
    Males   <- Thano(Dat$Pop1[Dat$Sex == "m"], Dat$dx[Dat$Sex == "m"])
    Females <- Thano(Dat$Pop1[Dat$Sex == "f"], Dat$dx[Dat$Sex == "f"])
    list(Males = Males, Females = Females)
}

# for legends. function is quite rigid. only does vertical color strips
colorStrip <- function(x = .0105,y=5,w=.002,h = 100,BrewerPal,revcol = FALSE,rising = TRUE,...){
    # define color functions
    if (revcol){
        colR        <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, BrewerPal)), space = "Lab")
    } else {
        colR        <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, BrewerPal), space = "Lab")
    }
    
    y <- seq(0,1,length=13) * h + y
    
    rect(x, y[1:12], x + w, y[2:13], col = colR(12), ...)
    if (rising){
        text(x+w,y[1:12],seq(0,110,by=10),cex=.7,pos=4)
    } else {
        text(x+w,y[2:13],seq(110,0,by=-10),cex=.7,pos=4)
    }
    
}











