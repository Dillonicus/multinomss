% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{spatial_zones}
\alias{spatial_zones}
\title{Calculates the spatial zones, starting with each centroid and increasing size
from 1 to pre-specified k points. These represent circular zones in space.}
\usage{
spatial_zones(locs, spatial_prop)
}
\arguments{
\item{locs}{An n x 2 matrix of centroid point coordinates in latitude/longitude}

\item{spatial_prop}{A float/decimal specifying the maximum proportion of the total population to be included in the scanning windows}
}
\description{
Calculates the spatial zones, starting with each centroid and increasing size
from 1 to pre-specified k points. These represent circular zones in space.
}
