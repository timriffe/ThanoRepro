\name{Minf0}
\alias{Minf0}
\title{Minf0 utility function to replace \code{Inf} with 0}
\usage{
  Minf0(M)
}
\arguments{
  \item{M}{any vector of data}
}
\value{
  M the same vector or array, with 0s imputed.
}
\description{
  This function is not the most rigorous computation
  practice, but it reduces headache when working with real
  data that may have 0 in the denominator.
}

