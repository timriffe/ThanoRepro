\name{Mna0}
\alias{Mna0}
\title{Mna0 utility function to replace \code{NA} or \code{NaN} with 0}
\usage{
  Mna0(M)
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

