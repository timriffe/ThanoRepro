
# this script contains auximilary functions used in plotting in this paper. These are peculiar and not necessarily useful beyond their application in the paper scripts. Since I am lame and do all plotting in base, in invariably end up with myriad such custom figure functions...

#'
#' @title makeVectorAgeGroups move single age data into age groups
#' 
#' @param vec a vector of data to be grouped
#' @param ages vector of single ages. \code{0:110} by default to match HMD.
#' @param N the age interval for the new groups.
#' 
#' @return a vector of grouped counts
#' 
#' @export
#' 
#' 

makeVectorAgeGroups <- function(vec, ages = 0:110, N = 2){
    x.new  <- ages - ages %% N
    c(tapply(vec, x.new, sum))
}
# draw rect based on x,y of cetroid, plus width and height
#'
#' @title makeRect auxiliary function
#' 
#' @description draws a rectangle centered on the given coordinates and of specified width and height
#' 
#' @param x center x coord.
#' @param y center y coord.
#' @param w width.
#' @param h height.
#' @param ... optional arguments passed to \code{rect()}.
#' 
#' @return NULL function called for its plotting side effects
#' 
#' @export

makeRect <- function(x, y, w, h, ...){
    w2 <- w / 2
    h2 <- h / 2
    rect(x - w2, y - h2, x + w2, y + h2, ...)
}
#'
#' @title PyramidOrLeafWithHeterogeneity plots either a leaf or a pyramid
#' 
#' @description This function takes cross-classified age matrices (thanatological age in rows and chronological age in columns, and will plot it either as a pyramid with death-cohorts highlighted with a color ramp or a thanatological leaf, with birth cohorts highlighted with a color ramp. I might not be very useful outside of this paper. The figure is added to an already-open device of the appropriate dimension. Scaled to sum to 1, as used in the paper. That's fixed here, so you'l need to re-define the function in order to make it more flexible...
#' 
#' @param Males a matrix of counts cross-classified by single chronological and thanatological ages, as returned by \code{Thano()}. These plot on the left.
#' @param Females a matrix of counts cross-classified by single chronological and thanatological ages, as returned by \code{Thano()}. These plot on the right.
#' @param Pyramid logical. \code{TRUE} will make a chronological age pyramid and \code{FALSE} will make a thanatological leaf.
#' @param N the age groups for the cohort heterogeneity represented by color bands. default 10.
#' @param BrewerPal the palette from \code{RColorBrewer} to be used.
#' @param revcol logical. Should the color order be reversed?
#' @param x the center x coordinate (count = 0)
#' @param y the bottom y coordinate (age = 0).
#' 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' 
#' @export
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
        colR        <- colorRampPalette(rev(brewer.pal(9, BrewerPal)), space = "Lab")
    } else {
        colR        <- colorRampPalette(RColorBrewer::brewer.pal(9, BrewerPal), space = "Lab")
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
    PyramidOutline(-colSums(Males, na.rm=TRUE),colSums(Females, na.rm=TRUE),scale=1,border=gray(.1), x=x, y=y, lwd = .2)
    #segments(x,y,x,111+y,col="white",lwd=.2)
}

#'
#' @title GetThanoMatrices make male and female cross-classified age matrices
#' 
#' @description This function is just for the sake of modularity. We have the HMD-HFD data object created for this paper in the data prep script, and rather than pre-calculating several thanatological age matrices, we make them on the fly.
#' 
#' @param .Code the country abbreviation
#' @param .Year which year of data. These are period matrices
#' @param Data the large data.frame produced in the data prep script.
#' 
#' @return a list with male and female matrices. Thano age in rows and chrono age in columns.
#' 
#' @export
#' 

# where Data is the universal huge data object for this paper
GetThanoMatrices <- function(.Code = "SWE", .Year = 1891, Data){
    Dat     <- Data[with(Data, Code == .Code & Year == .Year), ]
    
    Males   <- Thano(Dat$Pop1[Dat$Sex == "m"], Dat$dx[Dat$Sex == "m"])
    Females <- Thano(Dat$Pop1[Dat$Sex == "f"], Dat$dx[Dat$Sex == "f"])
    list(Males = Males, Females = Females)
}

#'
#' @title colorStrip a legend function for pyramid plotting
#' 
#' @description This is a rigid function, not very applicable outside the code used in this paper. Plots the legend in 10-year age groups from age 0 to 110. Take it or leave it ;-)
#' 
#' @param x the lower left x coord of the vertical color strip
#' @param y the lower left y coord of the vertical color strip
#' @param x the width of the color strip in x units
#' @param h the height in ages
#' @param BrewerPal the \code{RColorBrewer} palette to be used.
#' @param revcol logical. reverse the color palette order?
#' @param rising logical. Should ages count up from the bottom or the other way around?
#' 
#' @return nothing is returned. Called for plotting side effects.
#' 
#' @export
#' 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' 
# for legends. function is quite rigid. only does vertical color strips
colorStrip <- function(x = .0105,y=5,w=.002,h = 100,BrewerPal,revcol = FALSE,rising = TRUE,...){
    # define color functions
    if (revcol){
        colR        <- colorRampPalette(rev(brewer.pal(9, BrewerPal)), space = "Lab")
    } else {
        colR        <- colorRampPalette(brewer.pal(9, BrewerPal), space = "Lab")
    }
    
    y <- seq(0,1,length=13) * h + y
    
    rect(x, y[1:12], x + w, y[2:13], col = colR(12), ...)
    if (rising){
        text(x+w,y[1:12],seq(0,110,by=10),cex=.7,pos=4)
    } else {
        text(x+w,y[2:13],seq(110,0,by=-10),cex=.7,pos=4)
    }
    
}

