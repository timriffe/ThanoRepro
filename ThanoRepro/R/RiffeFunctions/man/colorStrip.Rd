\name{colorStrip}
\alias{colorStrip}
\title{colorStrip a legend function for pyramid plotting}
\usage{
  colorStrip(x = 0.0105, y = 5, w = 0.002, h = 100,
    BrewerPal, revcol = FALSE, rising = TRUE, ...)
}
\arguments{
  \item{x}{the lower left x coord of the vertical color
  strip}

  \item{y}{the lower left y coord of the vertical color
  strip}

  \item{x}{the width of the color strip in x units}

  \item{h}{the height in ages}

  \item{BrewerPal}{the \code{RColorBrewer} palette to be
  used.}

  \item{revcol}{logical. reverse the color palette order?}

  \item{rising}{logical. Should ages count up from the
  bottom or the other way around?}
}
\value{
  nothing is returned. Called for plotting side effects.
}
\description{
  This is a rigid function, not very applicable outside the
  code used in this paper. Plots the legend in 10-year age
  groups from age 0 to 110. Take it or leave it ;-)
}

