\name{theor.Lc.plot}
\alias{theor.Lc.plot}
\title{Plot Theoretical Lorenz Curves}
\description{
 add theoretical Lorenz curves of income distributions in a Lorenz curve plot}
\usage{
theor.Lc.plot(type=c("Singh-Maddala","Dagum","lognorm","Pareto","exponential"),parameter=0, col=2, lty=1)
}
\arguments{
\item{type}{character string giving the income distribution.
  must be one of the strings in the default argument (the first character
  is sufficient). defaults to "Singh-Maddala".}
 \item{parameter}{vector containing parameters of the distributions}
}
\references{C Dagum: Income Distribution Models, 1983, in: Johnson / Kotz (Eds): Encyclopedia of Statistical Sciences Vol.4, 27-34.

J B McDonald: Some generalized functions for the size distribution of income, 1984, Econometrica 52, 647-664.}
\author{Achim Zeileis \email{zeileis@statistik.uni-dortmund.de}}

\seealso{\code{\link{Lc}}, \code{\link{Lc.plot}}}

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
\keyword{aplot}

