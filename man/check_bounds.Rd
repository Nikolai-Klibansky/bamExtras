% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_bounds.r
\name{check_bounds}
\alias{check_bounds}
\title{Check to see if any estimated parameters are near bounds}
\usage{
check_bounds(x, p_cutoff = 0.01)
}
\arguments{
\item{x}{a BAM rdat object, i.e. x <- dget(rdat)}

\item{p_cutoff}{set minimum acceptable value for how close to a bound the estimated parameter can be before considering it too close to a bound, as a proportion of the range of bounds supplied to BAM}
}
\description{
Check to see if any estimated parameters are near bounds
}
\examples{
\dontrun{
check_bounds(rdat_BlackSeaBass)
}

}
\author{
Nikolai Klibansky
}
\keyword{assessment}
\keyword{bam}
\keyword{fisheries}
\keyword{stock}
