#
# several functions used to optimize for r or create projection matrices.
#
#' @title getca get stable chronological age structure from La and r
#' 
#' @description This is a straightforward identity
#' 
#' @param .a age midpoints (same for thano/chrono), default \code{.5:110.5}.
#' @param La lifetable exposure (or survival, if you're lazy).
#' @param r intrinsic growth rate.
#' 
#' @return ca the stable age structure (proportion in each age class).
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export

getca <- function(La,r,a=.5:110.5){
  La * exp(a * -r) / sum(La * exp(a * -r))
}

#' @title rThanoCoale optimize r from thanatological fertility rates and da.
#' 
#' @description This is an adaptation of A Coale's 1957 fast-converging solution to find r. 
#' 
#' @param fy fertility rates classified by thanatological age.
#' @param da lifetable deaths distribution.
#' @param .a age midpoints (same for thano/chrono), default \code{.5:110.5}.
#' @param maxit maximum number of iterations to converge, default 100.
#' @param tol convergence tolerance, default 1e-15.
#' 
#' @return a single estimate of r.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export
#' 
#' @importFrom compiler cmpfun

rThanoCoale <- cmpfun(function(fy, da, .a = .5:110.5, maxit = 1e2, tol = 1e-15){  
    
    # Based on Coale (1957), modified.
    N    <- length(fy)
    dxM  <- matrix(0, ncol = N, nrow = N)
    dxi  <- da
    for (i in 1:N){
      dxM[i, 1:length(dxi)  ] <- dxi 
      dxi <- dxi[2:length(dxi) ]
    }     
    R0      <- sum(dxM * fy)
    T.guess <- sum(.a * dxM * fy) / R0 # assuming r = 0
    r.i     <- log(R0) / T.guess
    # be careful to discount Fy by SRB appropriately for males / females
    # prior to specification
    for (i in 1:maxit){ # 15 is more than enough!
      #cat(r2,i,"\n")
      deltai <- 1 - sum(colSums(dxM * exp(-r.i * .a)) * fy)
      # the mean generation time self-corrects 
      # according to the error produced by the Lotka equation
      r.i <- r.i - (deltai / (T.guess - (deltai / r.i))) 
      if (abs(deltai) <= tol){
        break
      }
    }
    return(r.i)  
  })

#' @title Rmomentn nth moment of Lotka equation.
#' 
#' @description This is a straightforward identity for chronological age.
#' 
#' @param fa ASFR, all ages, including zeros.
#' @param La lifetable exposure (or survival, if you're lazy).
#' @param a age (midpoints, usually).
#' @param n the moment. default is 0, which will give back net reproduction.
#' 
#' @return the value of the moment.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export
#' 
#' @importFrom compiler cmpfun

Rmomentn <- cmpfun(function(fa,La,a,n=0){
    sum((a^n)*fa*La)
  })

#' @title rLotkaCoale Coale's 1957 suggestion for optimizing r
#' 
#' @description This method works well. Check \code{1-sum(exp(-r*a)*Lx*fa) < tol}.
#' 
#' @param fa ASFR, all ages, including zeros.
#' @param La lifetable exposure (or survival, if you're lazy).
#' @param a age (midpoints, usually).
#' @param T.guess a guess at mean generation length. default 29.
#' @param maxit maximum number of iterations to converge, default 100.
#' @param tol convergence tolerance, default 1e-15.
#' 
#' @return r optimized
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export
#' 
#' @importFrom compiler cmpfun
rLotkaCoale <- cmpfun(function(fa,La,a=.5:110.5, T.guess = 29, maxit = 1e2, tol = 1e-15){
    # from Coale, Ansley J. (1957) A New Method for Calculating Lotka's r- the Intrinsic Rate of Growth in a Stable Population.
    # Population Studies, Vol. 11 no. 1, pp 92-94
    R0 <- Rmomentn(fa,La,a,0)
    # first assuming a mean generation time of 29
    ri <- log(R0)/T.guess 
    for (i in 1:maxit){ # 10 is more than enough!
      deltai <- sum(exp(-ri*a)*fa*La)-1
      # the mean generation time self-corrects 
      # according to the error produced by the Lotka equation
      ri <- ri+(deltai/(29-(deltai/ri)))
      if (abs(deltai) <= tol){
        break
      }
    }
    return(ri)  
  })

#' @title Leslie a super cheap Lesile matrix
#' 
#' @description This could be much more thoughtful.
#' 
#' @param fa ASFR, all ages, including zeros.
#' @param Sa vector of discrete survival probabilities, usually composed of L(a+1)/L(a)
#' @param lambda discount fertility for those births that don't make it to Dec 31, and mothers' between-year mortality.
#' 
#' @return the Leslie projection matrix
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export

Leslie <- function(fa, Sa, lambda){
  fa <- fa * lambda 
  cbind(rbind(fa[1:length(Sa)], diag(Sa)), 0)
}

#' @title ThanoMin residual functin for optimizing thanatological r
#' 
#' @description for use in \code{optim()} or similar
#' 
#' @param r intrinsic growth rate.
#' @param da lifetable deaths distribution.
#' @param fy fertility rates classified by thanatological age.
#' @param .y age (midpoints, usually). 
#' 
#' @return residual of given r guess.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export
#' 
#' @importFrom compiler cmpfun
ThanoMin <- cmpfun(function(r, da, fy, .y = .5:110.5){
    # get the overlapped / staggered dx structure
    N    <- length(fy)
    daM  <- matrix(0, ncol = N, nrow = N)
    dai  <- da
    for (i in 1:N){
      daM[i, 1:length(dai)  ] <- dai 
      dai <- dai[2:length(dai) ]
    }   
    1 - sum(rowSums(t(t(daM) * exp(-r * .y))) * fy)
  })
#' @title ThanoProjMatrix a super cheap Lesile matrix for thanatological model
#' 
#' @description This could be much more thoughtful.
#' 
#' @param fy fertility rates classified by thanatological age.
#' @param da lifetable deaths distribution.
#' @param lambda discount fertility the partof infant mortality that doesn't make it to the end of the year.
#' 
#' @return the projection matrix
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export
#' @importFrom compiler cmpfun
ThanoProjMatrix <- cmpfun(function(fy, da, lambda){
    N       <- length(fy)
    # discount for part of infant mortality not surviving until end of year
    da[1]   <- da[1] * (1-lambda)
    
    # NxN matrix
    # fertility component
    Y       <- outer(da, fy, "*") + 
               # add survival element-wise
               rbind(cbind(0,diag(N - 1)), 0)
    
    # reduce e0 fertility by 1/2, as only exposed for part of year
    Y[, 1]  <- Y[, 1] / 2
    # do not allow for Inf or NA values: impute 0
    Y       <- Mna0(Minf0(Y))
    # return projection matrix
    Y
  })

