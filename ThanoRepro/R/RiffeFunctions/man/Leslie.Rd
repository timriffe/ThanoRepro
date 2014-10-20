\name{Leslie}
\alias{Leslie}
\title{Leslie a super cheap Lesile matrix}
\usage{
  Leslie(fa, Sa, lambda)
}
\arguments{
  \item{fa}{ASFR, all ages, including zeros.}

  \item{Sa}{vector of discrete survival probabilities,
  usually composed of L(a+1)/L(a)}

  \item{lambda}{discount fertility for those births that
  don't make it to Dec 31, and mothers' between-year
  mortality.}
}
\value{
  the Leslie projection matrix
}
\description{
  This could be much more thoughtful.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

