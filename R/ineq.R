ineq <-function(x,parameter=1,type=c("Gini","RS","Atkinson","Theil","Kolm","var","square.var","entropy"))
{
  switch(match.arg(type),
  Gini = Gini(x),
  RS = RS(x),
  Atkinson = Atkinson(x, parameter=parameter),
  Theil = Theil(x, parameter=parameter),
  Kolm = Kolm(x, parameter=parameter),
  var = var.coeff(x),
  square.var = var.coeff(x, square=T),
  entropy = entropy(x, parameter=parameter))
}

conc <- function(x, parameter=1, type=c("Herfindahl", "Rosenbluth"))
{
  switch(match.arg(type),
  Herfindahl = Herfindahl(x, parameter=parameter),
  Rosenbluth = Rosenbluth(x))
}

pov <- function(x,k,parameter=1,type=c("Watts", "Sen", "Foster"))
{
  switch(match.arg(type),
  Watts = Watts(x,k),
  Sen = Sen(x,k),
  Foster = Foster(x,k,parameter=parameter))
}

Lc <- function(x, n=rep(1,length(x)), plot=F)
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
    if(plot==T)
        Lc.plot(Lc)
    Lc
}

Lc.plot <- function(Lc, general=F, lwd=2,xlab="p",ylab="L(p)",main="Lorenz curve", new=F, col=1, lty=1)
{
    p <- Lc$p
    if(general==F)
      L <- Lc$L
    else
      L <- Lc$L.general
    if(new==F)
    {
     plot(p,L,type="l",main=main,lwd=lwd,xlab=xlab,ylab=ylab,xaxs="i",yaxs="i",las=1, col=col, lty=lty)
     abline(0,max(L))
    }
    else
    lines(p,L, type="l", lwd=lwd, col=col, lty=lty)
}

theor.Lc.plot<-function(type=c("Singh-Maddala","Dagum","lognorm","Pareto","exponential"),parameter=0, col=2, lty=1)
{
  dummy <- (0:1000)*0.001
  switch(match.arg(type),
  "Singh-Maddala"=lines(dummy,Lc.singh(dummy,parameter=parameter),col=col,type="l",lty=lty),
  Dagum=lines(dummy,Lc.dagum(dummy,parameter=parameter),col=col,type="l",lty=lty),
  lognorm=lines(dummy,Lc.lognorm(dummy,parameter=parameter),col=col,type="l",lty=lty),
  Pareto =lines(dummy,Lc.pareto(dummy,parameter), col=col,type="l",lty=lty),  
  exponential = lines(dummy, Lc.exp(dummy), col=col, type="l", lty=lty))
}

theor.Lc<-function(p,parameter=0,type=c("Singh-Maddala","Dagum","lognorm","Pareto","exponential"))
{ 
  switch(match.arg(type), 
  "Singh-Maddala" = Lc.singh(p,parameter=parameter),
  Dagum = Lc.dagum(p,parameter=parameter),
  lognorm = Lc.lognorm(p,parameter=parameter),
  Pareto = Lc.pareto(p,parameter=parameter),  
  exponential = Lc.exp(p))
}

Lc.exp <- function(p)
{
  elc <- 1/(1-p)
  elc <- (1-p)*log(elc)
  elc <- p - elc
  elc
}

Lc.lognorm <- function(p, parameter=1)
{
  if(parameter[1]>0)
    sigma <- parameter[1]
  else
  {
    cat("inadmissible parameter. default parameter=1 is used. \n")
    sigma <- 1
  }
  loglc <- p
  loglc[!loglc==0 & !loglc==1] <- pnorm(qnorm(loglc[!loglc==0 & !loglc==1]) - sigma)
  loglc
}

Lc.pareto <- function(p, parameter=2)
{
  if(parameter[1]>1) k<-(parameter[1]-1)/parameter[1]
  else
  {
    cat("inadmissible parameter. default parameter=2 is used. \n")
    k <- 0.5
  }
  parlc <- 1-((1-p)^k)
  parlc
}

Lc.singh <- function(p, parameter=c(2,2)) 
{ 
  if(!(is.na(parameter[2]))&(parameter[1]>0)&(parameter[1]<2+parameter[2]))
  {
    b <- parameter[1]-1
    d <- 1/(parameter[2]+1)
  }
  else 
  {
    cat("inadmissible parameter. default parameter=c(2,2) is used. \n")
    b <- 1 
    d <- 1/3
  }
  smlc <- pbeta((1-(1-p)^b), (1+d), (b-d))
  smlc 
} 

Lc.dagum <- function(p, parameter=c(2,2)) 
{ 
  if(!(is.na(parameter[2]))&(parameter[1]>1))
  {
    a <- 1/parameter[1]
    b <- parameter[2]
  }
  else 
  {
    cat("inadmissible parameter. default parameter=c(2,2) is used. \n")
    a <- 0.5 
    b <- 2
  }
  daglc <- pbeta((p^b), (a+1/b), (1-a))
  daglc 
} 

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
  c <- L[2:(k-1)] - beta*p[2:(k-1)]
  if(k==3)
    q <- NULL
  else
    q <- (c[2:(k-2)]-c[1:(k-3)])/(beta[1:(k-3)]-beta[2:(k-2)])
  q <- c(q,1)
  K <- beta*q + c
  L <- c(0,0,K,1)
  p <- c(0,-c[1]/beta[1],q,1)
  L <- L[is.finite(p)]
  p <- p[is.finite(p)]
  Lc.max <- list(p,L)
  names(Lc.max) <- c("p", "L")
  Lc.max
}

major <- function(x,y)
{
    x <- sort(x)
    y <- sort(y)
    n <- length(x)
    if((length(y)==n)&(sum(x)==sum(y)))
    {
        if(sum((cumsum(x)-cumsum(y))<=0)==n) result <- T
        else result <- F
        result
    }
    else
        cat("incomparable arguments \n")
}

Pen <- function(x, main="Pen Parade", ylab=expression(x[(i)]/bar(x)), 
xlab=expression(i/n), col=4, lwd=2, las=1)
{
  n <- length(x)
  a <- (0:n)/n
  b <- sort(c(0,x))/mean(x)
  plot(a,b, type="S", main=main, ylab=ylab,xlab=xlab, xaxs="i", yaxs="i",col=col, lwd=lwd, las=las)
  abline(1,0, lty=3)
}

Gini <- function(x)
{
    n <- length(x)
    x <- sort(x)
    i <- 1:n
    G <- t(x)%*%i
    G <- 2*G/(n*sum(x))
    G <- G - 1 - (1/n)
    as.vector(G)
}

RS <- function(x)
{
    d <- abs(x - mean(x))
    d <- mean(d)/(2*mean(x))
    d
}

Atkinson <- function(x, parameter=0.5)
{
    if(parameter==1)
        A <- 1 - ((prod(x)^(1/length(x)))/mean(x))
    else
        {
            x <- (x/mean(x))^(1-parameter)
            A <- 1 - mean(x)^(1/(1-parameter))
        }
    A
}

var.coeff <- function(x, square=F)
{
    n <- length(x)
    V <- sqrt((n-1)*var(x)/n)/mean(x)
    if(square==T) V <- V^2
    V
}

Theil <- function(x, parameter=0)
{
  if(parameter==0)
  {
    x <- x[!(x==0)]
    Th <- x/mean(x)
    Th <- sum(x*log(Th))
    Th <- Th/sum(x)
  }
  else
  {
    Th <- (prod(x)^(1/length(x)))/mean(x)
    Th <- -log(Th)
  }
  Th
}

Herfindahl <- function(x, parameter=1)
{
  m <- parameter
  Herf <- x/sum(x)
  Herf <- Herf^(m+1)
  Herf <- sum(Herf)^(1/m)
  Herf
}

Kolm <- function(x, parameter=1)
{
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
  if(parameter==0)
    e <- Theil(x, parameter=1)
  else
  if(parameter==1)
    e <- Theil(x, parameter=0)
  else
  {
    c <- parameter
    e <- (x/mean(x))^c
    e <- mean(e-1)/(c*(c-1))
  }
  e
}

Sen <- function(x,k)
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

Watts <- function(x,k)
{
  x2 <- x[x<k]
  if(length(x2)<1)
    0
  else
    sum(log(k/x2))/length(x)
}

Foster <- function(x,k,parameter=1)
{
  x2 <- x[x<k]
  if(length(x2)<1)
    0
  else
    sum(((k-x2)/k)^(parameter-1))/length(x)
}

