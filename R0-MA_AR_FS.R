# This script calculates basic reproduction number (Ro).
# It is calculated using Martingale (Bouma et al., 1997),
# Attack rate (Obadia et al., 2012) and FinalSize algoritm
# (De Jong and Kimman, 1994).
# CI95% and CI68% are estimated for FinalSize algoritm.
# Use at own risk at extremal values (e.g. S1=0, I1>0, Ro>100).


# Input
# Io, I1 - Infected animals (start and end state)
# So, S1 - Susceptible animals (start and end state)
I0<-6
I1<-0
S0<-6
S1<-3

# Check-up of input
if (I0<1) stop('Wrong input (Io).')
if (I1<0) stop('Wrong input (I1).')
if (S0<1) stop('Wrong input (So).')
if (S1>S0) stop('Wrong input (S1).')

#FS (Final size algorithm)
FS<-function(S0, S1, I0, I1)
  {
  # Initial value for RR
  RR<-1
  # Total amount of animals
  N<-S0+I0
  # Function calculates probability of S1 at given R
  prob<-function(RR, S1)
    {
    # Creating matrix of probabilities
    # Matrix contains row 0 and column 0. So it is bigger than N x So
    # Besides this it has inverted notation for rows
    # Row 0, column 0 should be adressed as [N+1, 0]
    # Row N, column S0 should be addressed as [1, S0+1]
    prob.m<-matrix(0, nrow=N+1, ncol=S0+1, dimnames=list(c(N:0),c(0:S0)))
    # Filling with initial state (P=1)
    prob.m[N-I0+1, S0+1]<-1
    # beta is used as coefficient of probability. We can calculate it only once
    beta<-N/(RR*S0+N)
    # Filling the first column
    for (i in 1:I0) prob.m[N-I0+i+1, S0+1]<-prob.m[N-I0+i, S0+1]*beta
    # Diagonal filling
    for (S.temp in S0:1)
      {
      # We can calculate beta once per column
      beta<-(S.temp)/(RR*(S.temp)+N)
      # Filling first well
      prob.m[S.temp, S.temp]<-prob.m[S.temp+1, S.temp+1]*RR*beta
      # Filling all individual columns
      for (I.temp in 0:(N-S.temp-2))
        {
        beta2<-N/(RR*(S.temp-1)+N)
        prob.m[S.temp+I.temp+1, S.temp]<-prob.m[S.temp+I.temp+2, S.temp+1]*RR*beta+prob.m[S.temp+I.temp, S.temp]*beta2
        }
      # Filling lower two lines
      # It can be no execution of upper cycle. But beta should be calculated at least once.
      beta2<-N/(RR*(S.temp-1)+N)
      prob.m[N, S.temp]<-prob.m[N-1, S.temp]*beta2
      prob.m[N+1, S.temp]<-prob.m[N, S.temp]*beta2
    }
    # Internal check-up. P=1
    if ((round(as.numeric(sum(prob.m[N+1,])), 2))!=1) stop('Internal error (P<>1).')
    # End of creating matrix. We returning Prob(R) at given S1.
    return(prob.m[N+1,S1+1])
  }
  # Optimization of R to get Pmax for S1
  optimized.prob<-function(RR) prob(RR, S1)
  RR<-optimize(optimized.prob, interval=0:1000, maximum=TRUE)$maximum
  # Optimization of R to get lower CI95
  lowerCI.prob<-function(RR)
    {
    CI<-0
    for (i in 0:S1) CI<-CI+prob(RR, i)
    if (CI>0.025) CI<--CI
    return (CI)
  }
  lowerCI<-optimize(lowerCI.prob, interval=0:1000, maximum=TRUE)$maximum
  # Optimization of R to get higher CI95
  higherCI.prob<-function(RR)
  {
    CI<-0
    for (i in S1:S0) CI<-CI+prob(RR, i)
    if (CI>0.025) CI<--CI
    return (CI)
  }
  higherCI<-optimize(higherCI.prob, interval=0:1000, maximum=TRUE)$maximum
  # Optimization of R to get lower CI68
  lowerCI68.prob<-function(RR)
  {
    CI<-0
    for (i in 0:(S1)) CI<-CI+prob(RR, i)
    if (CI>0.16) CI<--CI
    return (CI)
  }
  lowerCI68<-optimize(lowerCI68.prob, interval=0:1000, maximum=TRUE)$maximum
  # Optimization of R to get higher CI68
  higherCI68.prob<-function(RR)
  {
    CI<-0
    for (i in S1:S0) CI<-CI+prob(RR, i)
    if (CI>0.16) CI<--CI
    return (CI)
  }
  higherCI68<-optimize(higherCI68.prob, interval=0:1000, maximum=TRUE)$maximum
  # Internal check-up of initial I0 and final R values
  if (RR<0) stop('Internal error (Ro<0).')
  if (I1==0) return(cat("Final-size: ", RR, " CI95% [", lowerCI, "-", higherCI, "]", " CI68% [", lowerCI68, "-", higherCI68, "]\n")) else return (cat("FS:         Cannot be estimated (I1>0).\n"))
}
#End FS (Final size algorithm)

#MA (Martingale algorithm)
MA<-function(S0, S1, I0, I1)
  {
  # Initial value for RR
  RR<-0
  for (i in (S1+1):S0) RR<-RR+1/i
  RR<-(S0+I0)/(I0+I1+S0-S1)*RR
  # At S1=0, R is not exactly determined
  if (S1==0) RR<-paste(">=", RR)
  # At S1=So, R should not be determined
  if (S1==S0) RR<-0 
  return(cat("Martingale :", RR, "\n"))
  }
#End MA (Martingale algorithm)

#AR (Attack rate algorithm)
AR<-function(S0, S1, I0)
{
  # AR - the percentage of the population eventually infected
  AR<-(I0+S0-S1)/(S0+I0)
  # SS - initial percentage of susceptible population
  SS<-S0/(S0+I0)
  RR<--log10((1-AR)/SS)/(AR-(1-SS))
  return(cat("Attack rate:", RR))
}
#End AR (Attack rate algorithm)

cat(FS(S0, S1, I0, I1), MA(S0, S1, I0, I1), AR(S0, S1, I0))
