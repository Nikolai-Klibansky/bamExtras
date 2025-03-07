% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/M_Charnov.R
\name{M_Charnov}
\alias{M_Charnov}
\title{Compute age-varying natural mortality (M) based on Charnov et al (2013) and scale based on constant M}
\usage{
M_Charnov(a, Linf, K, t0, aP = 0, M_constant = NULL, aMin = 0, par_b = -1.5)
}
\arguments{
\item{a}{ages at which to compute M}

\item{Linf}{length infinity}

\item{K}{growth coefficient}

\item{t0}{age (time) at length zero}

\item{aP}{proportion of age (value between 0 and 1) at which to compute length (e.g. aP=0.5 to compute length at midyear)}

\item{M_constant}{constant M value used to scale age-varying M. Defaults to NULL. M will be scaled if a numeric value of M_constant is supplied.}

\item{aMin}{minimum age to include in computation of scaling factor when scaling M}

\item{par_b}{b parameter for Charnov equation ((L/Linf)^par_b)*K. Defaults to value from manuscript par_b=-1.5.}
}
\description{
Compute age-varying natural mortality (M) based on Charnov et al (2013) and scale based on constant M
}
\examples{
\dontrun{
# M for Black Sea Bass
M_Charnov(a=0:11,Linf = 502, K = 0.173, t0 = -0.97, M_constant=0.38, aMin=2)
}
}
\references{
Charnov, E. L., H. Gislason, and J. G. Pope. 2013. Evolutionary assembly rules for fish life histories. Fish and Fisheries 14:213-224.
}
\author{
Nikolai Klibansky
}
\keyword{assessment}
\keyword{bam}
\keyword{fisheries}
\keyword{stock}
