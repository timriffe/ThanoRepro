\name{PyramidOrLeafWithHeterogeneity}
\alias{PyramidOrLeafWithHeterogeneity}
\title{PyramidOrLeafWithHeterogeneity plots either a leaf or a pyramid}
\usage{
  PyramidOrLeafWithHeterogeneity(Males, Females,
    Pyramid = TRUE, N = 10, BrewerPal = "YlGnBu",
    revcol = FALSE, x = 0, y = 0)
}
\arguments{
  \item{Males}{a matrix of counts cross-classified by
  single chronological and thanatological ages, as returned
  by \code{Thano()}. These plot on the left.}

  \item{Females}{a matrix of counts cross-classified by
  single chronological and thanatological ages, as returned
  by \code{Thano()}. These plot on the right.}

  \item{Pyramid}{logical. \code{TRUE} will make a
  chronological age pyramid and \code{FALSE} will make a
  thanatological leaf.}

  \item{N}{the age groups for the cohort heterogeneity
  represented by color bands. default 10.}

  \item{BrewerPal}{the palette from \code{RColorBrewer} to
  be used.}

  \item{revcol}{logical. Should the color order be
  reversed?}

  \item{x}{the center x coordinate (count = 0)}

  \item{y}{the bottom y coordinate (age = 0).}
}
\description{
  Deprecated. Switched code to use \code{plotPyr()}. This
  function takes cross-classified age matrices
  (thanatological age in rows and chronological age in
  columns, and will plot it either as a pyramid with
  death-cohorts highlighted with a color ramp or a
  thanatological leaf, with birth cohorts highlighted with
  a color ramp. I might not be very useful outside of this
  paper. The figure is added to an already-open device of
  the appropriate dimension. Scaled to sum to 1, as used in
  the paper. That's fixed here, so you'l need to re-define
  the function in order to make it more flexible...
}

