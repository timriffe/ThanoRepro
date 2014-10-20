\name{GetThanoMatrices}
\alias{GetThanoMatrices}
\title{GetThanoMatrices make male and female cross-classified age matrices}
\usage{
  GetThanoMatrices(.Code = "SWE", .Year = 1891, Data)
}
\arguments{
  \item{.Code}{the country abbreviation}

  \item{.Year}{which year of data. These are period
  matrices}

  \item{Data}{the large data.frame produced in the data
  prep script.}
}
\value{
  a list with male and female matrices. Thano age in rows
  and chrono age in columns.
}
\description{
  This function is just for the sake of modularity. We have
  the HMD-HFD data object created for this paper in the
  data prep script, and rather than pre-calculating several
  thanatological age matrices, we make them on the fly.
}

