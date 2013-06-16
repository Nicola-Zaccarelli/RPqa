\name{Eps}
\alias{Eps}
\title{Threshold value - FINISH}
\description{
The procedure provides different ways of calculating a suitable threshold value $\epsilon$ for calculating a meaningful RP or CRP.
}
\usage{
Eps(dist.mat, type="perc.pd", param=0.1)
}
\arguments{
  \item{dist.mat}{A distance matrix from a time series or a vector of data.}
  \item{type}{The valid types of algorithm are: "max.pd" maximum phase state space diameter, "mean.pd" mean phase state space diameter, "median.pd" median phase state space diameter, "perc.pd" percentages of phase state space diameter, "perc.rr" percentages of recurrence rate in the final RP, "fix.nn" fixed numeb of neighbouring points, "noise" noise level; "range" to evaluate a range of possible value.}
  \item{param}{For the folowing type of threshold param means: for "perc.pd" the user specified percentage of phase state space diameter (i.e., a number between zero and one); for "perc.rr" the user specifed percentage of RR (i.e., a number between zero and one); for "fix.nn" the number of neighbouring points (by default it should be less than 1/10 of the total points number); for "range" a vector of three elements is necessary c(start, end, increment), so to explore the range start:end by increment.}
}

\details{
Computes the sample correlation integral on the provided time series for the specified length scale, and considering a time window \code{t} (see references). It  uses a naif algorithm: simply returns the fraction of points pairs nearer than \code{eps}. Normally, you would use \code{\link{d2}}, which takes roughly the same time, but computes the correlation sum for multiple length scales and embedding dimensions at once.
}
\value{
A list, a matrix or a vector along with the type of threshold selected.

FOR "fix.nn" is painfully slow as soon as n > 1000.
}

\references{
Hegger, R., Kantz, H., Schreiber, T., Practical implementation of nonlinear time series methods: The TISEAN package; CHAOS 9, 413-435 (1999)
}

\author{Nicola ZACCARELLI}
\keyword{RP}