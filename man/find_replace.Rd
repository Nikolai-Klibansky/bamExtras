% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_replace.r
\name{find_replace}
\alias{find_replace}
\title{Find and replace single or multiple values in a vector. Usage is similar to gsub except that this function will not replace parts of character strings}
\usage{
find_replace(pattern, replacement, x)
}
\arguments{
\item{pattern}{vector of values to find}

\item{replacement}{vector of values to replace}

\item{x}{a character vector where matches are sought}
}
\description{
Find and replace single or multiple values in a vector. Usage is similar to gsub except that this function will not replace parts of character strings
}
\examples{
\dontrun{
# Find and replace names of multiple parameters
rdat <- rdat_BlackSeaBass
parms <- rdat$parms
parms_names <- names(parms)
parms_key <- c("R0"="BH.R0", "M.msst"="M.constant", "M.MSST"="M.constant")
find_replace(names(parms_key),parms_key,parms_names)
}
}
\author{
Nikolai Klibansky
}
\keyword{assessment}
\keyword{bam}
\keyword{fisheries}
\keyword{stock}
