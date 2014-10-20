\name{Thano}
\alias{Thano}
\title{Thano redistribute age an age-classified vector into an age-remaining years cross-classified matrix}
\usage{
  Thano(Px, dx, stagger = TRUE)
}
\arguments{
  \item{Px}{any single-year age-classified vector of data.
  Age classes must align with those of \code{dx} starting
  at 0.}

  \item{dx}{single-age density column of the lifetable. Age
  classes must align with those of \code{Px} starting at
  0.}
}
\description{
  A discrete approximation of the thanatological
  redistribution formula. If \code{length(Px) !=
  length(dx)}, we pad out the shorter of the two with 0s.
}

