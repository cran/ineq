\name{pov}
\alias{pov}
\alias{Watts}
\alias{Sen}
\alias{Foster}
\title{Poverty Measures}
\description{
 computes the poverty of an (income) vector according to the
specified poverty measure }
\usage{
pov(x,k,parameter=1,type=c("Watts", "Sen", "Foster")) 

Watts(x,k)
Sen(x,k)
Foster(x,k,parameter=1)
}
\arguments{
 \item{x}{a vector containing at least non-negative elements}
 \item{k}{a constant giving the absolute poverty boundary}
 \item{parameter}{parameter of the poverty measure}
 \item{type}{character string giving the measure used to compute poverty
    coefficient must be one of the strings in the default argument (the first
    character is sufficient). defaults to "Watts".}
}
\details{\code{pov} is just a wrapper for the poverty measures of
         \code{Watts}, \code{Sen} and \code{Foster}
         (Foster / Greer / Thorbecke).

         \code{Foster} gives for parameter 1 the head count ratio and for
         parameter 2 the poverty gap ratio.
}
\value{
  the value of the poverty measure
}
\references{Buhong Zheng: Aggregate Poverty Measures, 1997, Journal of Economic Surveys Vol.11 No.2, 123-162. 
 
J E Foster: On Economic Poverty: A Survey Of Aggregate Measures, 1984, Advances in Econometrics Vol. 3, 215-251. 
}
\author{Achim Zeileis \email{zeileis@ci.tuwien.ac.at}}

\seealso{\code{\link{ineq}}, \code{\link{conc}}}

\examples{
# generate vectors (of incomes)
x <- c(541, 1463, 2445, 3438, 4437, 5401, 6392, 8304, 11904, 22261)
y <- c(841, 2063, 2445, 3438, 4437, 5401, 6392, 8304, 11304, 21961)
# compute Watts index with poverty boundary 2000
pov(x, 2000)
pov(y, 2000)
# compute head count ratio with poverty boundary 2000
pov(x, 2000, parameter=1, type="Foster")
pov(y, 2000, parameter=1, type="Foster")}
\keyword{}

