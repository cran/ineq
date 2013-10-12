pov <- function(x, k, parameter = NULL, type = c("Watts", "Sen", "SST", "Foster"))
{
  switch(match.arg(type),
  Watts = Watts(x, k),
  Sen = Sen(x, k),
  SST = SST(x, k),
  Foster = Foster(x, k, parameter = parameter))
}

Sen <- function(x, k)
{
  x <- as.numeric(x)
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

SST <- function(x, k)
{
  x <- as.numeric(x)
  x2 <- sort(x[x < k])
  n <- length(x)
  q <- length(x2)
  if(q < 1) 0 else {
    sum((2 * n - 2 * 1:q + 1) * (k - x2)/k)/n^2
  }
}

Watts <- function(x, k)
{
  x <- as.numeric(x)
  x2 <- x[x<k]
  if(length(x2)<1)
    0
  else
    sum(log(k/x2))/length(x)
}

Foster <- function(x, k, parameter=1)
{
  x <- as.numeric(x)
  if(is.null(parameter)) parameter <- 1
  x2 <- x[x<k]
  if(length(x2)<1)
    0
  else
    sum(((k-x2)/k)^(parameter-1))/length(x)
}
