\name{AKm02q0}
\alias{AKm02q0}
\title{AKm02q0 derive m0 from q0}
\usage{
  AKm02q0(m0, constant, slope)
}
\arguments{
  \item{m0}{a value of m0, per the period lifetable
  derivation}

  \item{constant}{the intercept of the relevant
  Andreev-Kingkade segment}

  \item{slope}{the slope of the relevant Andreev-Kingkade
  segment}
}
\value{
  q0 the estimate of q0 according to the identity between
  a0, m0, q0
}
\description{
  Derive m0 from q0 according to the relevant segment of
  the Andreeve-Kingkade formula. This is elegant because
  it's an analytic solution, but ugly because, man, look at
  it. Carl Boe got this from Maple I think. This formula is
  only necessary because AK start with q0 whereas the HMD
  starts with m0, so we needed to adapt. This is just one
  expression of the AK formulas. In a forthcoming paper,
  the authors will simply offer a new set of piecewise
  coefficients to calculate a0 directly from m0, rather
  than this roundabout quadratic solution.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

