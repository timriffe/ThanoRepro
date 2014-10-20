\name{getca}
\alias{getca}
\title{getca get stable chronological age structure from La and r}
\usage{
  getca(La, r, a = 0.5:110.5)
}
\arguments{
  \item{.a}{age midpoints (same for thano/chrono), default
  \code{.5:110.5}.}

  \item{La}{lifetable exposure (or survival, if you're
  lazy).}

  \item{r}{intrinsic growth rate.}
}
\value{
  ca the stable age structure (proportion in each age
  class).
}
\description{
  This is a straightforward identity
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

