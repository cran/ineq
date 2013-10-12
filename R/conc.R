conc <- function(x, parameter=NULL, type=c("Herfindahl", "Rosenbluth"))
{
  switch(match.arg(type),
  Herfindahl = Herfindahl(x, parameter=parameter),
  Rosenbluth = Rosenbluth(x))
}

Herfindahl <- function(x, parameter=1)
{
  x <- as.numeric(x)
  if(is.null(parameter))
    m <- 1
  else
    m <- parameter
  Herf <- x/sum(x)
  Herf <- Herf^(m+1)
  Herf <- sum(Herf)^(1/m)
  Herf
}

Rosenbluth <- function(x)
{
  x <- as.numeric(x)
  n <- length(x)
  x <- sort(x)
  HT <- (n:1)*x
  HT <- 2*sum(HT/sum(x))
  HT <- 1/(HT-1)
  HT
}

