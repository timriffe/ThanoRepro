\name{prepare4plot}
\alias{prepare4plot}
\title{prepare4plot gathers info needed for plotting and puts it in a list}
\usage{
  prepare4plot(XXX, colsFun1, colsFun2)
}
\arguments{
  \item{XXX}{list returned by \code{GetThanoMatrices()}
  with male and female matrix elements}

  \item{colsFun1}{color ramp function for chronological
  age}

  \item{colsFun2}{color ramp function for thanatological
  age}
}
\value{
  list with several components passed to special plotting
  function, \code{plotPyr()}
}
\description{
  This is an internal function, very specific and rigid.
}

