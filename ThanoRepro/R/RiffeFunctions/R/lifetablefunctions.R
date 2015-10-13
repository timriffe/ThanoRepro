#'
#' @title AKm02q0 derive m0 from q0
#' 
#' @description Derive m0 from q0 according to the relevant segment of the Andreeve-Kingkade formula. This is elegant because it's an analytic solution, but ugly because, man, look at it. Carl Boe got this from Maple I think. This formula is only necessary because AK start with q0 whereas the HMD starts with m0, so we needed to adapt. This is just one expression of the AK formulas. In a forthcoming paper, the authors will simply offer a new set of piecewise coefficients to calculate a0 directly from m0, rather than this roundabout quadratic solution.
#' 
#' @param m0 a value of m0, per the period lifetable derivation
#' @param constant the intercept of the relevant Andreev-Kingkade segment
#' @param slope the slope of the relevant Andreev-Kingkade segment
#' 
#' @return q0 the estimate of q0 according to the identity between a0, m0, q0
#' 
#' @author Tim Riffe \email{tim.riffe@@gmail.com}
#' 
#' @export
#' 

AKm02q0 <- function(m0,constant,slope){
  -1 / slope / m0 * (-m0 +  (m0 * constant) - 0.1e1 + sqrt(((m0 ^ 2) - 2 * constant * (m0 ^ 2) + 2 * m0 + (constant ^ 2) * (m0 ^ 2) - 2 *  (m0 * constant) + 1 - 4 * slope * (m0 ^ 2)))) / 2
}

#' @title \code{AKm02a0} estimates a0 using the Andreev-Kinkade rule of thumb.
#'
#' @description
#' \code{AKm02a0} is an auxiliary function used by version 6 of the four HMD lifetable functions, \code{ltper_AxN()}, \code{ltcoh_AxN()}, \code{ltperBoth_AxN()}, \code{ltcohBoth_AxN()}. This function calls \code{AKm02q0()} to help get the work done, since the HMD needed to adapt the Andreev-Kingkade formulas to work with the period lifetable flow.
#'
#' @param m0 a value or vector of values of m0, the death rate for age 0 infants.
#' @param sex either "m" or "f"
#' 
#' @return a0, the estimated average age at death of those dying in the first year of life, either a single value or a vector of a_0 values.
#' 
#' @author Tim Riffe \email{tim.riffe@@gmail.com}
#' 
#' @export

AKm02a0 <- function(m0, sex = "m"){
  sex <- rep(sex,length(m0))
  ifelse(sex == "m",
    ifelse(m0 < 0.02306737, 0.1493 - 2.0367 * AKm02q0(m0, 0.1493, -2.0367),
      ifelse(m0 < 0.0830706, 0.0244 + 3.4994 * AKm02q0(m0, 0.0244, 3.4994), .2991)),
    ifelse(m0 < 0.01725977, 0.1490 - 2.0867 * AKm02q0(m0, 0.1490, -2.0867),
      ifelse(m0 < 0.06919348, 0.0438 + 4.1075 * AKm02q0(m0, 0.0438, 4.1075), 0.3141))
  )
}

#'
#' @title mx2dxHMD derive dx from mx using the HMD lifetable protocol
#' 
#' @description This is a typical lifetable function, but it produces only the dx column, rather than the full table. Sex is required as an argument for a0 estimation, where the new Andreev-Kingkade formula is used (to be introduced in version 6 of the HMD methods protocol).
#' 
#' @param mx a vector of event-exposure mortality rates.
#' @param sex either \code{"m"} or \code{"f"}.
#' 
#' @return dx a vector the same length as mx, radix 1.
#' 
#' @author Tim Riffe \email{tim.riffe@@gmail.com}
#' 
#' @export
#' 
#' @importFrom compiler cmpfun
#' 

mx2dxHMD <- cmpfun(function(mx, sex = "m"){
    mx                  <- Mna0(as.numeric(mx))
    # mean proportion of interval passed at death
    ax                  <- mx * 0 + .5                      # ax = .5, pg 38 MPv5
    
    ax[1]               <- AKm02a0(mx[1], sex)              # v6 a0 protocol
    
    qx                  <- mx / (1 + (1 - ax) * mx)          # Eq 60 MPv5 (identity)
# ---------------------------------------------------------------------------------
# set open age qx to 1
    i.openage           <- length(mx) # removed argument OPENAGE
    qx[i.openage]       <- 1
    ax[i.openage]       <- 1 / mx[i.openage]                   
# ---------------------------------------------------------------------------------
# define remaining lifetable columns:
    px                  <- 1 - qx       # Eq 64 MPv5
    px[is.nan(px)]      <- 0 # skips BEL NAs, as these are distinct from NaNs
    lx                  <- c(1, cumprod(px[1:(i.openage-1)]))
    # NA should only be possible if there was a death with no Exp below age 80- impossible, but just to be sure
    # lx[is.na(lx)]   <- 0 # removed for BEL testing        
    dx                  <- lx * qx # Eq 66 MPv5
# return result
    Mna0(dx)
  })

#'
#' @title derive DYA from d(a) 
#' 
#' @description This function makes a triangle out of d(a)
#' 
#' @param da the lifetable deaths distribution classified by chronological age
#' 
#' @return d(Y,A) the lifetable remaining lifetime distribution, not conditioned on surviva;
#' 
#' @author Tim Riffe \email{tim.riffe@@gmail.com}
#' 
#' @export
#' 

da2DYA <- function(da){
	N       <- length(da)
	ay      <- 1:N - 1
	
	da      <- Mna0(da)   # remove NAs if any       
	da      <- c(da, da * 0) / sum(da) # pad out with 0s
	fya     <- matrix(da[col(matrix(nrow = N, 
									ncol = N)) + ay], 
			nrow = N, 
			ncol = N, 
			dimnames = list(Ex = ay, 
					Age = ay)
	)
	fya
}

#'
#' @title derive f(y,a) from d(a) 
#' 
#' @description This function conditions d(a+y) on survival to age a.
#' 
#' @param da the lifetable deaths distribution classified by chronological age
#' @param stagger logical. default \code{TRUE} should remaining lifetime results be blended with neighboring ages?
#' 
#' @return fya the lifetable remaining lifetime distribution conditioned on survival.
#' 
#' @author Tim Riffe \email{tim.riffe@@gmail.com}
#' 
#' @export
#' 

da2fya <- function(da, stagger = FALSE){
  fya     <- da2DYA(da)
  if (stagger){
    fya <- (fya + cbind(fya[, 2:ncol(fya)], 0)) / 2
  }
  fya <- Minf0(Mna0(fya / rowSums(fya)))
  fya
}

