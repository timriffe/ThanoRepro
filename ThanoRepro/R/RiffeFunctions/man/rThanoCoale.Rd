\name{rThanoCoale}
\alias{rThanoCoale}
\title{rThanoCoale optimize r from thanatological fertility rates and da.}
\usage{
  rThanoCoale(fy, da, .a = 0.5:110.5, maxit = 100,
    tol = 1e-15)
}
\arguments{
  \item{fy}{fertility rates classified by thanatological
  age.}

  \item{da}{lifetable deaths distribution.}

  \item{.a}{age midpoints (same for thano/chrono), default
  \code{.5:110.5}.}

  \item{maxit}{maximum number of iterations to converge,
  default 100.}

  \item{tol}{convergence tolerance, default 1e-15.}
}
\value{
  a single estimate of r.
}
\description{
  This is an adaptation of A Coale's 1957 fast-converging
  solution to find r.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

