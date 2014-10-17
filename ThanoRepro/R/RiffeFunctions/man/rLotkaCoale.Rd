\name{rLotkaCoale}
\alias{rLotkaCoale}
\title{rLotkaCoale Coale's 1957 suggestion for optimizing r}
\usage{
  rLotkaCoale(fa, La, a = 0.5:110.5, T.guess = 29,
    maxit = 100, tol = 1e-15)
}
\arguments{
  \item{fa}{ASFR, all ages, including zeros.}

  \item{La}{lifetable exposure (or survival, if you're
  lazy).}

  \item{a}{age (midpoints, usually).}

  \item{T.guess}{a guess at mean generation length. default
  29.}

  \item{maxit}{maximum number of iterations to converge,
  default 100.}

  \item{tol}{convergence tolerance, default 1e-15.}
}
\value{
  r optimized
}
\description{
  This method works well. Check
  \code{1-sum(exp(-r*a)*Lx*fa) < tol}.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

