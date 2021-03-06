% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example_data_generate.R
\name{example_data_generate}
\alias{example_data_generate}
\title{A function to simulate data for an example}
\usage{
example_data_generate(n, J, beta)
}
\arguments{
\item{n}{sample size.}

\item{J}{dimension for genotype matrix.}

\item{beta}{a vector of parameters.}
}
\value{
A list contains outcome vector Y and genotype matrix G.
}
\description{
This function returns a list contains the simulated outcome variable Y and genotype matrix G.
}
\examples{
\dontrun{
 library(MASS)
 J = 100
 n = 50
 beta=c(rep(2,3),rep(0,J-3))
 data=example_data_generate(n, J, beta)
}
}
