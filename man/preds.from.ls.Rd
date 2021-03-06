% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictor_generation.R
\name{preds.from.ls}
\alias{preds.from.ls}
\title{Predictor Generation}
\usage{
preds.from.ls(pred_choice, gridsize, code, n_preds_original, n_preds,
  base_ = NA, formu, effect.size = 1)
}
\arguments{
\item{pred_choice}{vector  of integers defining the set of predictors originally selected}

\item{gridsize}{integer determining size of the map}

\item{code}{code (e.g. "110") to be used in simData, check FReibier::simData documentation for further help}

\item{n_preds_original}{integer defining the number of predictors}

\item{n_preds}{integer defining total number of predictors, quadratic conditions and interactions}

\item{formu}{formula without an response such as ~ x1 + x2}

\item{effect.size}{max distance of the base coefficients to zero}

\item{base}{vector of numerics as the parameters coefficients}
}
\value{
matrix of predicotrs and a map suitable for plotting with ggplot2
}
\description{
Generate predictors_ from a landscape
}
\examples{
preds.from.ls()

}
\keyword{predictors}
