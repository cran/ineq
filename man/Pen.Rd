\name{Pen}
\alias{Pen}
\title{Pen's Parade}
\description{
 plots Pen's Parade of a vector \code{x}}
\usage{
Pen(x, n = rep(1, length(x)),
  scaled = TRUE, abline = TRUE, segments = FALSE,  
  main = "Pen's Parade", ylab = NULL, xlab = NULL, 
  col = 4, lwd = 2, las = 1, fill = NULL, \dots
  
}
\arguments{
 \item{x}{a vector containing non-negative elements.}
 \item{n}{a vector of frequencies or weights, must be same length as \code{x}.}
 \item{scaled}{logical. Should Pen's parade be divided by \code{mean(x)}?}
 \item{abline}{logical. Should a horizontal line for the mean be drawn?}
 \item{segments}{logical. Should histogram-like segments be drawn?}
 \item{fill}{a color for filling the area under the curve.}
 \item{xlab,ylab}{axis labels. Suitable defaults depending on
   \code{scaled} and \code{n} are chosen.}
 \item{main, col, lwd, las, \dots}{further high-level \code{\link{plot}} parameters.}
}
\details{
Pen's Parade is basically the inverse distribution function (standardized by
\code{mean(x)}).}
\references{F A Cowell: Measurement of Inequality, 2000, in A B Atkinson / F Bourguignon (Eds): Handbook of Income Distribution, Amsterdam,

F A Cowell: Measuring Inequality, 1995 Prentice Hall/Harvester Wheatshef,

J Pen: Income Distribution, 1971, Harmondsworth: Allen Lane.}
\author{Achim Zeileis \email{zeileis@ci.tuwien.ac.at}}

\seealso{\code{\link{Lc}}, \code{\link{plot.Lc}}}

\examples{
# load and attach Philippine income data
data(Ilocos)
attach(Ilocos)
# plot Pen's Parade of income
Pen(income)
}
\keyword{hplot}

