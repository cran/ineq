ineq <- function(x, parameter=NULL, type=c("Gini", "RS", "Atkinson", "Theil",
    "Kolm", "var", "square.var", "entropy"))
{
  switch(match.arg(type),
  Gini = Gini(x),
  RS = RS(x),
  Atkinson = Atkinson(x, parameter=parameter),
  Theil = Theil(x, parameter=parameter),
  Kolm = Kolm(x, parameter=parameter),
  var = var.coeff(x),
  square.var = var.coeff(x, square=TRUE),
  entropy = entropy(x, parameter=parameter))
}

conc <- function(x, parameter=NULL, type=c("Herfindahl", "Rosenbluth"))
{
  switch(match.arg(type),
  Herfindahl = Herfindahl(x, parameter=parameter),
  Rosenbluth = Rosenbluth(x))
}

pov <- function(x, k, parameter=NULL, type=c("Watts", "Sen", "Foster"))
{
  switch(match.arg(type),
  Watts = Watts(x,k),
  Sen = Sen(x,k),
  Foster = Foster(x,k,parameter=parameter))
}

Lc <- function(x, n=rep(1,length(x)), plot=FALSE)
{
    k <- length(x)
    o <- order(x)
    x <- x[o]
    n <- n[o]
    x <- n*x
    p <- cumsum(n)/sum(n)
    L <- cumsum(x)/sum(x)
    p <- c(0,p)
    L <- c(0,L)
    L2 <- L * mean(x)
    Lc <- list(p,L,L2)
    names(Lc) <- c("p", "L", "L.general")
    class(Lc) <- "Lc"
    if(plot)
        plot(Lc)
    Lc
}

plot.Lc <- function(x, general=FALSE, lwd=2,xlab="p",ylab="L(p)",
  main="Lorenz curve", las=1, ...)
{
    if(!general)
      L <- x$L
    else
      L <- x$L.general
    plot(x$p, L, type="l", main=main, lwd=lwd, xlab=xlab, ylab=ylab, xaxs="i",
      yaxs="i", las=las, ...)
    abline(0,max(L))
}

lines.Lc <- function(x, general=FALSE, lwd=2, ...)
{
    if(!general)
      L <- x$L
    else
      L <- x$L.general
    lines(x$p, L, lwd=lwd, ...)
}

plot.theorLc <- function(x, parameter=NULL, xlab="p", ylab="L(p)", lwd=2, las=1, ...)
{
  dummy <- (0:1000)*0.001
  if(is.null(parameter))
    plot(dummy, x(dummy), type="l", xlab=xlab, ylab=ylab, xaxs="i", yaxs="i",
      lwd=lwd, las=las, ...)
  else
    plot(dummy, x(dummy,parameter=parameter), type="l", xlab=xlab,
      ylab=ylab, xaxs="i", yaxs="i", lwd=lwd, las=las, ...)
  abline(0,1)
}

lines.theorLc <- function(x, parameter=NULL, lwd=2, col=2, ...)
{
  dummy <- (0:1000)*0.001
  if(is.null(parameter))
    lines(dummy, x(dummy), lwd=lwd, col=col, ...)
  else
    lines(dummy, x(dummy,parameter=parameter), lwd=lwd, col=col, ...)
}

theorLc <- function(type=c("Singh-Maddala", "Dagum", "lognorm", "Pareto",
    "exponential"), parameter=0)
{
  switch(match.arg(type),
  "Singh-Maddala" = rval <- function(p) {Lc.singh(p,parameter=parameter)},
  Dagum = rval <- function(p) {Lc.dagum(p,parameter=parameter)},
  lognorm = rval <- function(p) {Lc.lognorm(p,parameter=parameter)},
  Pareto = rval <- function(p) {Lc.pareto(p,parameter=parameter)},
  exponential = rval <- function(p) {Lc.exp(p)})
  class(rval) <- "theorLc"
  return(rval)
}

Lc.exp <- function(p)
{
  elc <- 1/(1-p)
  elc <- (1-p)*log(elc)
  elc <- p - elc
  elc
}
class(Lc.exp) <- "theorLc"

Lc.lognorm <- function(p, parameter=1)
{
  if(parameter[1]>0)
    sigma <- parameter[1]
  else
  {
    warning("inadmissible parameter. default parameter=1 is used.")
    sigma <- 1
  }
  loglc <- p
  loglc[!loglc==0 & !loglc==1] <- pnorm(qnorm(loglc[!loglc==0 & !loglc==1]) - sigma)
  loglc
}
class(Lc.lognorm) <- "theorLc"

Lc.pareto <- function(p, parameter=2)
{
  if(parameter[1]>1) k<-(parameter[1]-1)/parameter[1]
  else
  {
    warning("inadmissible parameter. default parameter=2 is used.")
    k <- 0.5
  }
  parlc <- 1-((1-p)^k)
  parlc
}
class(Lc.pareto) <- "theorLc"

Lc.singh <- function(p, parameter=c(2,2))
{
  if(!(is.na(parameter[2]))&(parameter[1]>0)&(parameter[1]<2+parameter[2]))
  {
    b <- parameter[1]-1
    d <- 1/(parameter[2]+1)
  }
  else
  {
    warning("inadmissible parameter. default parameter=c(2,2) is used.")
    b <- 1
    d <- 1/3
  }
  smlc <- pbeta((1-(1-p)^b), (1+d), (b-d))
  smlc
}
class(Lc.singh) <- "theorLc"

Lc.dagum <- function(p, parameter=c(2,2))
{
  if(!(is.na(parameter[2]))&(parameter[1]>1))
  {
    a <- 1/parameter[1]
    b <- parameter[2]
  }
  else
  {
    warning("inadmissible parameter. default parameter=c(2,2) is used.")
    a <- 0.5
    b <- 2
  }
  daglc <- pbeta((p^b), (a+1/b), (1-a))
  daglc
}
class(Lc.dagum) <- "theorLc"

Lc.mehran <- function(x,n)
{
  Lc.min <- Lc(x,n=n)
  p <- Lc.min$p
  L <- Lc.min$L
  k <- length(p)
  q <- c(0,rep(1,k))
  K <- c(rep(0,k),1)
  for(i in k:2)
  {
    q[i] <- 2*p[i]-q[i+1]
  }
  for(i in 2:k)
  {
    K[i] <- 2*L[i-1] - K[i-1]
  }
  beta1 <- (L[2:k]-L[1:(k-1)])/(p[2:k]-p[1:(k-1)])
  beta2 <- (K[2:k]-K[1:(k-1)])/(q[2:k]-q[1:(k-1)])
  beta2 <- beta2[2:(k-1)]
  beta <- rep(0,(k-2))
  for(i in 1:(k-2))
  {
    if(beta1[i]>beta2[i])
      beta[i] <- beta1[i]
    else
    if(beta2[i]>beta1[i+1])
      beta[i] <- beta1[i+1]
    else
      beta[i] <- beta2[i]
  }
  d <- L[2:(k-1)] - beta*p[2:(k-1)]
  if(k==3)
    q <- NULL
  else
    q <- (d[2:(k-2)]-d[1:(k-3)])/(beta[1:(k-3)]-beta[2:(k-2)])
  q <- c(q,1)
  K <- beta*q + d
  L <- c(0,0,K,1)
  p <- c(0,-d[1]/beta[1],q,1)
  L <- L[is.finite(p)]
  p <- p[is.finite(p)]
  Lc.max <- list(p,L)
  names(Lc.max) <- c("p", "L")
  class(Lc.max) <- "Lc"
  Lc.max
}

major <- function(x,y)
{
    x <- sort(x)
    y <- sort(y)
    n <- length(x)
    if((length(y)==n)&(sum(x)==sum(y)))
        all((cumsum(x)-cumsum(y))<=0)
    else
        stop("incomparable arguments")
}

Pen <- function(x, n = rep(1, length(x)),
  scaled = TRUE, abline = TRUE, segments = FALSE,  
  main = "Pen's Parade", ylab = NULL, xlab = NULL, 
  col = 4, lwd = 2, las = 1, fill = NULL, ...) 
{
  o <- order(x)
  x <- x[o]
  n <- n[o]
    
  if(scaled) x <- x/mean(x)

  if(is.null(ylab)) {
    if(scaled) ylab <- expression(x[(i)]/bar(x))
      else ylab <- expression(x[(i)])
  }

  if(is.null(xlab)) {
    if(identical(all.equal(n, rep(1, length(x))), TRUE)) xlab <- expression(i/n)
      else xlab <- expression(Sigma[i](n[(i)]/n))
  }
  
  n <- cumsum(c(0, n))/sum(n)   

  plot(c(0, 1), c(0, max(x)), type = "n", main = main, ylab = ylab, xlab = xlab, 
      xaxs = "i", yaxs = "i", col = col, lwd = lwd, las = las,
      ...)

  ln <- length(n)
  n2 <- c(rep(n[-ln], rep(2, ln-1)), n[ln])
  x2 <- c(0, rep(x, rep(2, ln-1)))
      
  lines(n2, x2, col = col, lwd = lwd)
  
  if(!is.null(fill)) polygon(c(n2, 1, 0), c(x2, 0, 0), col = fill, border = col)

  if(abline) abline(h = mean(x), lty = 3)
  if (segments) segments(n, 0, n, x, col = col, lwd = lwd)
  box()
}


Gini <- function(x)
{
    n <- length(x)
    x <- sort(x)
    G <- sum(x * 1:n)
    G <- 2*G/(n*sum(x))
    G - 1 - (1/n)
}

RS <- function(x)
{
    d <- abs(x - mean(x))
    d <- mean(d)/(2*mean(x))
    d
}

Atkinson <- function(x, parameter=0.5)
{
    if(is.null(parameter)) parameter <- 0.5
    if(parameter==1)
        A <- 1 - (exp(mean(log(x)))/mean(x))
    else
        {
            x <- (x/mean(x))^(1-parameter)
            A <- 1 - mean(x)^(1/(1-parameter))
        }
    A
}

var.coeff <- function(x, square=FALSE)
{
    n <- length(x)
    V <- sqrt((n-1)*var(x)/n)/mean(x)
    if(square) V <- V^2
    V
}

Theil <- function(x, parameter=0)
{
  if(is.null(parameter)) parameter <- 0
  if(parameter==0)
  {
    x <- x[!(x==0)]
    Th <- x/mean(x)
    Th <- sum(x*log(Th))
    Th <- Th/sum(x)
  }
  else
  {
    Th <- exp(mean(log(x)))/mean(x)
    Th <- -log(Th)
  }
  Th
}

Herfindahl <- function(x, parameter=1)
{
  if(is.null(parameter))
    m <- 1
  else
    m <- parameter
  Herf <- x/sum(x)
  Herf <- Herf^(m+1)
  Herf <- sum(Herf)^(1/m)
  Herf
}

Kolm <- function(x, parameter=1)
{
  if(is.null(parameter)) parameter <- 1
  KM <- parameter * (mean(x)-x)
  KM <- mean(exp(KM))
  KM <- (1/parameter)*log(KM)
  KM
}

Rosenbluth <- function(x)
{
  n <- length(x)
  x <- sort(x)
  HT <- (n:1)*x
  HT <- 2*sum(HT/sum(x))
  HT <- 1/(HT-1)
  HT
}

entropy <- function(x, parameter=0.5)
{
  if(is.null(parameter)) parameter <- 0.5
  if(parameter==0)
    e <- Theil(x, parameter=1)
  else
  if(parameter==1)
    e <- Theil(x, parameter=0)
  else
  {
    k <- parameter
    e <- (x/mean(x))^k
    e <- mean(e-1)/(k*(k-1))
  }
  e
}

Sen <- function(x, k)
{
  x2 <- x[x<k]
  if(length(x2)<1)
    0
  else
  {
    H <- sum(x<k)/length(x)
    I <- sum((k-x2)/k)/length(x2)
    G <- Gini(x2)
    H*(I+(1-I)*G)
  }
}

Watts <- function(x, k)
{
  x2 <- x[x<k]
  if(length(x2)<1)
    0
  else
    sum(log(k/x2))/length(x)
}

Foster <- function(x, k, parameter=1)
{
  if(is.null(parameter)) parameter <- 1
  x2 <- x[x<k]
  if(length(x2)<1)
    0
  else
    sum(((k-x2)/k)^(parameter-1))/length(x)
}
