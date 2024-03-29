% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcLogLikC_prediction.R
\name{calcLogLikC_prediction}
\alias{calcLogLikC_prediction}
\title{calcLogLikC_prediction}
\usage{
calcLogLikC_prediction(par, c, returnOutcome = FALSE, returnPerMarker = FALSE)
}
\arguments{
\item{par}{parameters A list(mx,mu,omega,beta), The ordinary scale (ttheta)}

\item{c}{data object (returned from prepareC_calibration)}

\item{returnOutcome}{Whether to return all joint genotype outcome values (can be used for Deconvolution and marker-specific values)}

\item{returnPerMarker}{Whether to return logLik per marker}
}
\description{
Likelihood calculation of (evidence) prediction model (uses calibration info)
}
