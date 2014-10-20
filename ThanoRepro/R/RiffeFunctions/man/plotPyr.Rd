\name{plotPyr}
\alias{plotPyr}
\title{plotPyr plots heterogenous pyramid using output from \code{prepare4plot()}}
\usage{
  plotPyr(XXX, x, y, Thano = TRUE)
}
\arguments{
  \item{XXX}{list returned by \code{GetThanoMatrices()}
  with male and female matrix elements}

  \item{x}{coord of center of figure}

  \item{y}{corrd of bottom of figure}

  \item{Thano}{logical do we draw a thanatological leaf or
  a chronological pyramid?}
}
\value{
  NULL. function called for its plotting side effects
}
\description{
  This is an internal function, very specific and rigid.
}

