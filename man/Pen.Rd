\name{Pen}
\alias{Pen}
\title{Pen's Parade}
\description{
 plots Pen's Parade of a vector \code{x}}
\usage{
Pen(x, main="Pen's Parade", ylab=expression(x[(i)]/bar(x)), xlab=expression(i/n), col=4, lwd=2, las=1)
}
\arguments{
\item{x}{a vector containing non-negative elements}
 \item{...}{graphical parameters}
}
\details{
Pen's Parade is basically the inverse distribution function (standardized by
\code{mean(x)}.}
\references{F A Cowell: Measurement of Inequality, 2000, in A B Atkinson / F Bourguignon (Eds): Handbook of Income Distribution, Amsterdam,

F A Cowell: Measuring Inequality, 1995 Prentice Hall/Harvester Wheatshef,

J Pen: Income Distribution, 1971, Harmondsworth: Allen Lane.}
\author{Achim Zeileis \email{zeileis@ci.tuwien.ac.at}}

\seealso{\code{\link{Lc}}, \code{\link{Lc.plot}}}

\examples{
# generate income vector
x <- c(541, 1463, 2445, 3438, 4437, 5401, 6392, 8304, 11904, 22261)
# plot Pen's Parade of x
Pen(x)
}
\keyword{hplot}

