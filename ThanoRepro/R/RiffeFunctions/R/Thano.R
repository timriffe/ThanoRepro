#' @title Minf0 utility function to replace \code{Inf} with 0
#' 
#' @description This function is not the most rigorous computation practice, but it reduces headache when working with real data that may have 0 in the denominator.
#'  
#' @param M any vector of data
#' 
#' @importFrom compiler cmpfun
#' 
#' @return M the same vector or array, with 0s imputed.
#' 
#' @export 
#'

Minf0 <- cmpfun(function(M){
    M[is.infinite(M)]  <- 0
    M
  })

#' @title Mna0 utility function to replace \code{NA} or \code{NaN} with 0
#' 
#' @description This function is not the most rigorous computation practice, but it reduces headache when working with real data that may have 0 in the denominator.
#'  
#' @param M any vector of data
#' 
#' @importFrom compiler cmpfun
#' 
#' @return M the same vector or array, with 0s imputed.
#' 
#' @export 
#'
Mna0 <- compiler::cmpfun(function(M){
    M[is.na(M)]  <- 0
    M[is.nan(M)] <- 0
    M
  })



#' @title Thano redistribute age an age-classified vector into an age-remaining years cross-classified matrix
#' 
#' @description A discrete approximation of the thanatological redistribution formula. If \code{length(Px) != length(dx)}, we pad out the shorter of the two with 0s.
#'  
#' @param Px any single-year age-classified vector of data. Age classes must align with those of \code{dx} starting at 0.
#' @param dx single-age density column of the lifetable. Age classes must align with those of \code{Px} starting at 0.
#' 
#' @importFrom compiler cmpfun
#' 
#' @export 
#' 

Thano <- cmpfun(function(Px, dx, stagger = TRUE){
    Np <- length(Px)
    Nd <- length(dx)
    if (Np != Nd){
      N <- max(Np, Nd)
      Px <- c(Px, rep(0, N - Np))
      dx <- c(dx, rep(0, N - Nd))
    } else {
      N <- Np
    }
    
    ay      <- 1:N - 1
    
    dx      <- Mna0(dx)   # remove NAs if any       
    dx      <- c(dx, dx * 0) / sum(dx) # pad out with 0s
    EDx     <- matrix(dx[col(matrix(nrow = N, 
            ncol = N)) + ay], 
      nrow = N, 
      ncol = N, 
      dimnames = list(Ex = ay, 
        Age = ay)
    )
    if (stagger){
      EDx <- (EDx + cbind(EDx[, 2:ncol(EDx)], 0)) / 2
    }
    t(Px * Minf0(Mna0(EDx / rowSums(EDx))))
  })
