% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputeAges.R
\name{drawBaseline}
\alias{drawBaseline}
\title{Draw Ages Using the Inverse CDF Method from the baseline data}
\usage{
drawBaseline(baseline_data)
}
\arguments{
\item{baseline_data}{A data frame containing baseline data with columns 'cum_prob' and 'age'.}
}
\value{
A single age value drawn from the baseline data.
}
\description{
This function draws ages using the inverse CDF method from baseline data.
}
