% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FinalSize.R
\name{pExtremes}
\alias{pExtremes}
\title{Probability of extreme values given R}
\usage{
pExtremes(R, x, s0in, i0in, comp = `<`)
}
\arguments{
\item{R}{Value of R0}

\item{x}{Number of cases at the end of the outbreak}

\item{s0in}{Number of initial susceptibles}

\item{i0in}{Number of initial infectious}

\item{comp}{Direction of comparison i.e.(`>`or`<`)}

\item{r0in}{Number of initial recovered}
}
\value{
Exact probability that x cases or more extreme are found for this value of R
}
\description{
Probability of extreme values given R
}
\examples{
pExtremes(1, 1, 5,5, `<`)
pExtremes(1, 1, 5,5, `>`)

}
