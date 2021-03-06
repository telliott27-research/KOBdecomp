% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decomposition.R
\name{getExplained}
\alias{getExplained}
\title{Get Explained/Endowment Part of a Decomposition}
\usage{
getExplained(m, xes, ges, i = 1, j = 0)
}
\arguments{
\item{m}{a correctly specified model for an KOB decomposition}

\item{xes}{a character vector giving the names of the independent (or control) variables}

\item{ges}{a character vector giving the name(s) of the grouping variable and any polynomial transformations. The list should be ordered by increasing polynomial.}

\item{i, j}{values of the group variable to compare}
}
\value{
The amount of difference attributable to the explained part in a KOB decomposition
}
\description{
This function calculates the explained or endowment part of a Kitawaga-Oaxaca-Blinder (KOB) decomposition
generalized to a continuous grouping variable. Rarely used directly by the user.
}
