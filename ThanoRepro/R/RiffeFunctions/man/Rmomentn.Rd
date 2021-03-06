\name{Rmomentn}
\alias{Rmomentn}
\title{Rmomentn nth moment of Lotka equation.}
\usage{
  Rmomentn(fa, La, a, n = 0)
}
\arguments{
  \item{fa}{ASFR, all ages, including zeros.}

  \item{La}{lifetable exposure (or survival, if you're
  lazy).}

  \item{a}{age (midpoints, usually).}

  \item{n}{the moment. default is 0, which will give back
  net reproduction.}
}
\value{
  the value of the moment.
}
\description{
  This is a straightforward identity for chronological age.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

