\name{Lc.plot}
\alias{Lc.plot}
\title{Plot Lorenz Curve}
\description{
 plots the Lorenz curve of a vector x (like computed by \code{Lc})}
\usage{
Lc.plot(Lc, general=F, lwd=2,xlab="p",ylab="L(p)",main="Lorenz curve", new=F, col=1, lty=1) 
}
\arguments{
 \item{Lc}{a list like provided by \code{Lc}}
 \item{general}{logical flag. if TRUE the generalized Lorenz curve will be
plotted}
 \item{new}{logical flag. if TRUE the curve will be plotted into an
 existing plot.}
 \item{...}{graphical parameters}
}
\references{B C Arnold: Majorization and the Lorenz Order: A Brief Introduction, 1987, Springer,

F A Cowell: Measurement of Inequality, 2000, in A B Atkinson / F Bourguignon (Eds): Handbook of Income Distribution, Amsterdam,

F A Cowell: Measuring Inequality, 1995 Prentice Hall/Harvester Wheatshef.}
\author{Achim Zeileis \email{zeileis@ci.tuwien.ac.at}}

\seealso{\code{\link{Lc}}, \code{\link{Lc.mehran}},
\code{\link{theor.Lc.plot}}}

\examples{
# income distribution of the USA in 1968 (in 10 classes)
# x vector of class means, n vector of class frequencies
x <- c(541, 1463, 2445, 3438, 4437, 5401, 6392, 8304, 11904, 22261)
n <- c(482, 825, 722, 690, 661, 760, 745, 2140, 1911, 1024)

# compute minimal Lorenz curve (= no inequality in each group)
Lc.min <- Lc(x, n=n)
# compute maximal Lorenz curve (limits of Mehran)
Lc.max <- Lc.mehran(x,n)
# plot both Lorenz curves in one plot                        
Lc.plot(Lc.min)
Lc.plot(Lc.max, new=T, col=4)

# add the theoretic Lorenz curve of a Lognormal-distribution with (var=0.78)
theor.Lc.plot(type="lognorm", parameter=0.78)
# add the theoretic Lorenz curve of a Dagum-distribution
theor.Lc.plot(type="Dagum", parameter=c(3.4,2.6))
}
\keyword{hplot}

