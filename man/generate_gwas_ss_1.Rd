% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theoretical_gwas.R
\name{generate_gwas_ss_1}
\alias{generate_gwas_ss_1}
\title{Create a GWAS summary dataset}
\usage{
generate_gwas_ss_1(params, nid, vy = 1, minmaf = 0.01, r = NULL)
}
\arguments{
\item{params}{Output from \code{generate_gwas_params}}

\item{nid}{sample size}

\item{vy}{Variance of trait}

\item{minmaf}{minimum allowed maf. default=0.01 to prevent instability}

\item{r}{LD correlation matrix. If NULL (default) then creates LD unaware sampling errors}
}
\value{
list of data frames
}
\description{
Determines SE and generates effect estimates given input parameters
}
