\name{theor.Lc}
\alias{Lc.dagum}
\alias{Lc.singh}
\alias{Lc.pareto}
\alias{Lc.lognorm}
\alias{Lc.exp}
\title{Theoretical Lorenz Curves}
\description{
 theoretical Lorenz curves of income distributions}
\usage{
theor.Lc(p,parameter=0,type=c("Singh-Maddala","Dagum","lognorm","Pareto","exponential"))

Lc.dagum(p, parameter=c(2,2))
Lc.singh(p, parameter=c(2,2))
Lc.pareto(p, parameter=2)
Lc.lognorm(p, parameter=1)
Lc.exp(p)
}
\arguments{
\item{p}{vector with elements from [0,1]}
\item{parameter}{vector containing parameter(s) of the distributions}
\item{type}{character string giving the income distribution.
  must be one of the strings in the default argument (the first character
  is sufficient). defaults to "Singh-Maddala".}
}
\details{\code{theor.Lc} is just a wrapper for the theoretical Lorenz curves
         of income distributions \code{Lc.dagum}, \code{Lc.singh},
         \code{Lc.pareto}, \code{Lc.lognorm}, \code{Lc.exp}.

         \code{Lc.dagum} is the Lorenz curve of the Dagum distribution (2
         parameters), \code{Lc.singh} the one of the Singh-Maddala
         distribution (2 parameters), \code{Lc.pareto} the one of the Pareto
         distribution (1 parameter), \code{Lc.lognorm} the one of the Lognormal
         distribution (1 parameter) and \code{Lc.exp} the Lorenz curve of the
         exponential distribution (no parameter).
}
\value{value of the theoretical Lorenz curve at \code{p}}
\references{C Dagum: Income Distribution Models, 1983, in: Johnson / Kotz (Eds): Encyclopedia of Statistical Sciences Vol.4, 27-34.

J B McDonald: Some generalized functions for the size distribution of income, 1984, Econometrica 52, 647-664.}
\author{Achim Zeileis \email{zeileis@statistik.uni-dortmund.de}}

\seealso{\code{\link{Lc}}, \code{\link{Lc.plot}}, \code{\link{theor.Lc.plot}}}

\examples{
# vector of percentages
p <- (1:10)*0.1
# compute values of theoretic Lorenz curve of a Dagum-distribution
theor.Lc(p, parameter=c(3.4,2.6), type="Dagum")
# or
Lc.dagum(p, parameter=c(3.4,2.6))
}
\keyword{}

