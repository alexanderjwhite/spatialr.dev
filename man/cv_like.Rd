% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv_like.R
\name{cv_like}
\alias{cv_like}
\title{Compute the Cross Validated Log Likelihood}
\usage{
cv_like(test_x, u, v, test_nn, index_map)
}
\arguments{
\item{test_x}{test matrix}

\item{u}{model computed u}

\item{v}{model computed v}

\item{test_nn}{nearest neighbor matrix}

\item{index_map}{index map}
}
\value{
double
}
\description{
Compute the Cross Validated Log Likelihood
}
