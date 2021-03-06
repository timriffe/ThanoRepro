\name{da2fya}
\alias{da2fya}
\title{mx2dxHMD derive dx from mx using the HMD lifetable protocol}
\usage{
  da2fya(da, stagger = FALSE)
}
\arguments{
  \item{da}{the lifetable deaths distribution classified by
  chronological age}

  \item{stagger}{logical. default \code{TRUE} should
  remaining lifetime results be blended with neighboring
  ages?}
}
\value{
  fya the lifetable remaining lifetime distribution
  conditioned on survival.
}
\description{
  This function conditions d(a+y) on survival to age a.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

