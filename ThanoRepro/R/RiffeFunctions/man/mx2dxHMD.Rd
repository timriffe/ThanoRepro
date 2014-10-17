\name{mx2dxHMD}
\alias{mx2dxHMD}
\title{mx2dxHMD derive dx from mx using the HMD lifetable protocol}
\usage{
  mx2dxHMD(mx, sex = "m")
}
\arguments{
  \item{mx}{a vector of event-exposure mortality rates.}

  \item{sex}{either \code{"m"} or \code{"f"}.}
}
\value{
  dx a vector the same length as mx, radix 1.
}
\description{
  This is a typical lifetable function, but it produces
  only the dx column, rather than the full table. Sex is
  required as an argument for a0 estimation, where the new
  Andreev-Kingkade formula is used (to be introduced in
  version 6 of the HMD methods protocol).
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

