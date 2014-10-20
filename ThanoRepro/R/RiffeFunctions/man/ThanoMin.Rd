\name{ThanoMin}
\alias{ThanoMin}
\title{ThanoMin residual functin for optimizing thanatological r}
\usage{
  ThanoMin(r, da, fy, .y = 0.5:110.5)
}
\arguments{
  \item{r}{intrinsic growth rate.}

  \item{da}{lifetable deaths distribution.}

  \item{fy}{fertility rates classified by thanatological
  age.}

  \item{.y}{age (midpoints, usually).}
}
\value{
  residual of given r guess.
}
\description{
  for use in \code{optim()} or similar
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

