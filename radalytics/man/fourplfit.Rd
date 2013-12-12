\name{fourplfit}
\alias{fourplfit}
\title{Radfit: fit the data in a frame}
\description{
Defines the function with the following variables
}
\usage{
fourplfit(data,ul,hlimL,hlimH,vlim)
}
\arguments{
  \item{data}{
    output datafile from the plate reader via SPS
  }
  \item{ul}{
   the upper limit of the assay as specified by the manufacturer
  }
  \item{hlimL}{
    horizonal limit low. This is the CV threshold for unknowns that have a concentration above the vertical threshold (vlim).
  }
  \item{hlimH}{
    horizontal limit high. This is the CV threshold for unknowns that have a concentration below the vertical threshold and above the limit of detection
  }
  \item{vlim}{
    vertical limit. This is selected based on a number of semi-arbitrary criteria. Examples; the second lowest calibrator, the upper limit divided by 10, or just by eyeballing the concentration vs. CV plot 
  }
}
