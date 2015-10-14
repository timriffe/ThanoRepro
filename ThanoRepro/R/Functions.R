
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm", "tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/ThanoRepro/ThanoRepro")
} else {
	if (system("hostname",intern=TRUE) == "PC-403478"){
		# on MPIDR PC
		setwd("U://git//ThanoRepro//ThanoRepro")
	} else {
		# in that case I'm on Berkeley system, and other people in the dept can run this too
		setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ThanoRepro/ThanoRepro"))
	}
}
getwd()

PACKAGES <- installed.packages()[,"Package"]

# this will be tricky if you're not on a unix/linux system. MAC is more straightforward than Windows.
# see https://github.com/UCBdemography/DemogBerkeley for installation tips.
# requires some external programs be installed to run properly. If you try and fail, then email 
# the author. This is only used for reading HMD/HFD data into R.
#if (!"DemogBerkeley" %in% PACKAGES){
#  A <- try(devtools::install_github("UCBdemography/DemogBerkeley", subdir = "DemogBerkeley"))
#  if (class(A) == "try-error"){
#    cat("\nDemogBerkeley didn't install. \nsee https://github.com/UCBdemography/DemogBerkeley \n for installation tips")
#  }
#}

# TR: swapped in HMDHFDplus, much more stable!
require(HMDHFDplus)


devtools::load_all("R/RiffeFunctions")



#
#
##' @title Minf0 utility function to replace \code{Inf} with 0
##' 
##' @description This function is not the most rigorous computation practice, but it reduces headache when working with real data that may have 0 in the denominator.
##'  
##' @param M any vector of data
##' 
##' @importFrom compiler cmpfun
##' 
##' @return M the same vector or array, with 0s imputed.
##' 
##' @export 
##'
#
#Minf0 <- compiler::cmpfun(function(M){
#            M[is.infinite(M)]  <- 0
#            M
#        })
#
##' @title Mna0 utility function to replace \code{NA} or \code{NaN} with 0
##' 
##' @description This function is not the most rigorous computation practice, but it reduces headache when working with real data that may have 0 in the denominator.
##'  
##' @param M any vector of data
##' 
##' @importFrom compiler cmpfun
##' 
##' @return M the same vector or array, with 0s imputed.
##' 
##' @export 
##'
#Mna0 <- compiler::cmpfun(function(M){
#            M[is.na(M)]  <- 0
#            M[is.nan(M)] <- 0
#            M
#        })
#
#
#
##' @title Thano redistribute age an age-classified vector into an age-remaining years cross-classified matrix
##' 
##' @description A discrete approximation of the thanatological redistribution formula. If \code{length(Px) != length(dx)}, we pad out the shorter of the two with 0s.
##'  
##' @param Px any single-year age-classified vector of data. Age classes must align with those of \code{dx} starting at 0.
##' @param dx single-age density column of the lifetable. Age classes must align with those of \code{Px} starting at 0.
##' 
##' @importFrom compiler cmpfun
##' 
##' @export 
##' 
#Thano <- compiler::cmpfun(function(Px, dx, stagger = TRUE){
#            Np <- length(Px)
#            Nd <- length(dx)
#            if (Np != Nd){
#                N <- max(Np, Nd)
#                Px <- c(Px, rep(0, N - Np))
#                dx <- c(dx, rep(0, N - Nd))
#            } else {
#                N <- Np
#            }
#            
#            ay      <- 1:N - 1
#            
#            dx      <- Mna0(dx)   # remove NAs if any       
#            dx      <- c(dx, dx * 0) / sum(dx) # pad out with 0s
#            EDx     <- matrix(dx[col(matrix(nrow = N, 
#                                            ncol = N)) + ay], 
#                    nrow = N, 
#                    ncol = N, 
#                    dimnames = list(Ex = ay, 
#                            Age = ay)
#            )
#            if (stagger){
#                EDx <- (EDx + cbind(EDx[, 2:ncol(EDx)], 0)) / 2
#            }
#            t(Px * Minf0(Mna0(EDx / rowSums(EDx))))
#        })
#
#
#
## original
#ThanoSimple <- compiler::cmpfun(function(Px, dx){
#            ay      <- 0:110
#            dx      <- Mna0(dx)
#            dx      <- c(dx, dx * 0) / sum(dx) # pad out with 0s
#            EDx     <- matrix(dx[col(matrix(nrow = 111, ncol = 111)) + ay], 
#                    nrow = 111, ncol = 111, dimnames = list(Ex = ay, Age = ay))
#            t(Px * Minf0(Mna0(EDx / rowSums(EDx))))
#        })
#
#FindSmoothMode <- function(x, fy, delta = .05){
#    NA0ind <- fy == 0 | is.na(fy)
#    xy <- smooth.spline(x=x[!NA0ind],y = log(fy[!NA0ind]))
#    xy <- predict(xy, seq(min(x),max(x), by = delta))
#    xy$x[which.max(xy$y[xy$x<90])]
#}
#
#ThanoMin <- compiler::cmpfun(function(r, da, Fy, .a = .5:110.5){
#            # get the overlapped / staggered dx structure
#            N    <- length(Fy)
#            daM  <- matrix(0, ncol = N, nrow = N)
#            dai  <- da
#            for (i in 1:N){
#                daM[i, 1:length(dai)  ] <- dai 
#                dai <- dai[2:length(dai) ]
#            }   
#            1 - sum(rowSums(t(t(daM) * exp(-r * .a))) * Fy)
#        })
#
#ThanoProjMatrix <- compiler::cmpfun(function(Fy, da, lambda){
#    N       <- length(Fy)
#    # discount for part of infant mortality not surviving until end of year
#    da[1]   <- da[1] * (1-lambda)
#    
#    # NxN matrix
#    # fertility component
#    Y       <- outer(da, Fy, "*") + 
#            # add survival element-wise
#            rbind(cbind(0,diag(N - 1)), 0)
#    
#    # reduce e0 fertility by 1/2, as only exposed for part of year
#    Y[, 1]  <- Y[, 1] / 2
#    # do not allow for Inf or NA values: impute 0
#    Y       <- Mna0(Minf0(Y))
#    # return projection matrix
#    Y
#})
#
## 
#Leslie <- function(Fx, px, lambda){
#    Fx <- Fx * lambda # lambda = 1 - (DL / Births)
#    cbind(rbind(Fx[1:length(px)], diag(px)), 0)
#}
#
#mx2dxHMD <- compiler::cmpfun(function(mx){
#            mx                  <- Mna0(as.numeric(mx))
#            
#            # mean proportion of interval passed at death
#            ax                  <- mx * 0 + .5                      # ax = .5, pg 38 MPv5
#            
#            ax[1]   <- ((0.045 + 2.684 * mx[1]) + (0.053 + 2.800 * mx[1])) / 2 # hack, hard to pass in sex variable
#            
#            qx                  <- mx / (1 + (1 - ax) * mx)          # Eq 60 MPv5 (identity)
## ---------------------------------------------------------------------------------
## set open age qx to 1
#            i.openage           <- length(mx)+1 # removed argument OPENAGE
#            qx[i.openage]       <- ifelse(is.na(qx[i.openage]), NA, 1)
#            ax[i.openage]       <- 1 / mx[i.openage]                   
## ---------------------------------------------------------------------------------
## define remaining lifetable columns:
#            px                  <- 1 - qx                                                                                 # Eq 64 MPv5
#            px[is.nan(px)]      <- 0 # skips BEL NAs, as these are distinct from NaNs
## lx needs to be done columnwise over px, argument 2 refers to the margin.
#            lx                  <- c(1, cumprod(px[1:(i.openage-1)]))
#            # NA should only be possible if there was a death with no Exp below age 80- impossible, but just to be sure
#            # lx[is.na(lx)]   <- 0 # removed for BEL testing        
#            dx                  <- lx * qx                                                                                # Eq 66 MPv5
#            Mna0(dx)
#        })
#da2fya <- function(da, stagger = FALSE){
#    N       <- length(da)
#    ay      <- 1:N - 1
#    
#    da      <- Mna0(da)   # remove NAs if any       
#    da      <- c(da, da * 0) / sum(da) # pad out with 0s
#    fya     <- matrix(da[col(matrix(nrow = N, 
#                                    ncol = N)) + ay], 
#            nrow = N, 
#            ncol = N, 
#            dimnames = list(Ex = ay, 
#                    Age = ay)
#    )
#    if (stagger){
#        fya <- (fya + cbind(fya[, 2:ncol(fya)], 0)) / 2
#    }
#    fya <- Minf0(Mna0(fya / rowSums(fya)))
#    fya
#}
#
#getca <- function(La,r,a=.5:110.5){
#    La * exp(a * -r) / sum(La * exp(a * -r))
#}
#rThanoCoale <- compiler::cmpfun(function(Fy, da, .a = .5:110.5, maxit = 1e2, tol = 1e-15){  
#            
#            # Based on Coale (1957), modified.
#            N    <- length(Fy)
#            dxM  <- matrix(0, ncol = N, nrow = N)
#            dxi  <- da
#            for (i in 1:N){
#                dxM[i, 1:length(dxi)  ] <- dxi 
#                dxi <- dxi[2:length(dxi) ]
#            }     
#            R0      <- sum(dxM * Fy)
#            T.guess <- sum(.a * dxM * Fy) / R0 # assuming r = 0
#            r.i     <- log(R0) / T.guess
#            # be careful to discount Fy by SRB appropriately for males / females
#            # prior to specification
#            for (i in 1:maxit){ # 15 is more than enough!
#                #cat(r2,i,"\n")
#                deltai <- 1 - sum(sum(colSums(dxM * exp(-r.i * .a)) * Fy) * Fy)
#                # the mean generation time self-corrects 
#                # according to the error produced by the Lotka equation
#                r.i <- r.i - (deltai / (T.guess - (deltai / r.i))) 
#                if (abs(deltai) <= tol){
#                    break
#                }
#            }
#            return(r.i)  
#        })
## chronological age
#Rmomentn <- compiler::cmpfun(function(fx,Lx,x,n=0){
#            sum((x^n)*fx*Lx)
#        })
#rLotkaCoale <- compiler::cmpfun(function(fx,Lx,x=.5:110.5, T.guess = 29, maxit = 1e2, tol = 1e-15){
#            # from Coale, Ansley J. (1957) A New Method for Calculating Lotka's r- the Intrinsic Rate of Growth in a Stable Population.
#            # Population Studies, Vol. 11 no. 1, pp 92-94
#            R0 <- Rmomentn(fx,Lx,x,0)
#            # first assuming a mean generation time of 29
#            ri <- log(R0)/T.guess 
#            for (i in 1:maxit){ # 10 is more than enough!
#                deltai <- sum(exp(-ri*x)*fx*Lx)-1
#                # the mean generation time self-corrects 
#                # according to the error produced by the Lotka equation
#                ri <- ri+(deltai/(29-(deltai/ri)))
#                if (abs(deltai) <= tol){
#                    break
#                }
#            }
#            return(ri)  
#        })