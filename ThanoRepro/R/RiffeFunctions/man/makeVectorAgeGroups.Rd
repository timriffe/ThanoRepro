\name{makeVectorAgeGroups}
\alias{makeVectorAgeGroups}
\title{makeVectorAgeGroups move single age data into age groups}
\usage{
  makeVectorAgeGroups(vec, ages = 0:110, N = 2)
}
\arguments{
  \item{vec}{a vector of data to be grouped}

  \item{ages}{vector of single ages. \code{0:110} by
  default to match HMD.}

  \item{N}{the age interval for the new groups.}
}
\value{
  a vector of grouped counts
}
\description{
  makeVectorAgeGroups move single age data into age groups
}

