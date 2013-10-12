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

Gini <- function(x, corr = FALSE)
{
    x <- as.numeric(x)
    n <- length(x)
    x <- sort(x)
    G <- sum(x * 1:n)
    G <- 2*G/(n*sum(x))
    G <- G - 1 - (1/n)
    if(corr) G * n/(n-1) else G
}

RS <- function(x)
{
    x <- as.numeric(x)
    d <- abs(x - mean(x))
    d <- mean(d)/(2*mean(x))
    d
}

Atkinson <- function(x, parameter=0.5)
{
    x <- as.numeric(x)
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
    x <- as.numeric(x)
    n <- length(x)
    V <- sqrt((n-1)*var(x)/n)/mean(x)
    if(square) V <- V^2
    V
}

Theil <- function(x, parameter=0)
{
  x <- as.numeric(x)
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

Kolm <- function(x, parameter=1)
{
  x <- as.numeric(x)
  if(is.null(parameter)) parameter <- 1
  KM <- parameter * (mean(x)-x)
  KM <- mean(exp(KM))
  KM <- (1/parameter)*log(KM)
  KM
}

entropy <- function(x, parameter=0.5)
{
  x <- as.numeric(x)
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

