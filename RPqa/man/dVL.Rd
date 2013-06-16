\name{dVL}
\alias{dVL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distribution of diagonal and vertical lines
}
\description{
This function is the workhorse of the RPqa package. It calculates the distribution of diagonal and vertical lines in a binary matrix. For the moment the Thieler window is not there, but I need to add it.
The underlying C function performes the work of evaluating the matrix and identifying the relevant elements.
Diagonal and vertical elements are calculated one after the other, but room to get this done in parallel is there.
}
\usage{
dVL(RPMatrix, Thieler = 0)
}
\arguments{
  \item{RPMatrix}{An input binary matrix}
  \item{Thieler}{The integer value of the Thieler window}
  }

\details{
The procedure is calculating the the value of diagonal and vertical lines of all length in the upper-left part of the matrix. Data must be bunary. The LOI has to go from [N,1] to [1,N].
}
\value{
  \item{distrLV }{A matrix of two rows containing the frequency of diagonal (first row) and vertical (second row) of a specific length. The length is in the name of the column.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicola ZACCARELLI
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (RPMatrix, Thieler = 0) 
{
    stopifnot(length(Thieler) == 1, Thieler >= 0, Thieler == 
        floor(Thieler), class(RPMatrix) \%in\% c("matrix"), dim(RMatrix)[1] == 
        dim(RMatrix)[2])
    if (is.numeric(RMatrix[1, 1]) == FALSE) 
        stop("The matrix has to be numeric!")
    cat("IMPORTANT: the LOI must go from [N,1] to [1, N] \n")
    num = dim(RPMatrix)[1]
    dLV = matrix(0, 2, num)
    distrLV <- .C("dVL", D = as.integer(RPMatrix), N = as.integer(num), 
        T = as.integer(Theiler), dLV = as.integer(dLV), PACKAGE = "RPqa")$dLV
    distrLV = as.matrix(distrLV, 2, num)
    colnames(distrLV) = 0:num
    delCol = colSums(distrLV) * (0:num)
    distrLV[, -delCol]
    distrLV[, -1]
    rownames(distrLV) = c("Diagonal", "Vertical")
    return(distrLV)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
