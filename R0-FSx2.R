# This script calculates basic reproduction number (Ro).
# It is calculated using FinalSize algoritm (De Jong and Kimman, 1994)
# for two similar experiments with similar or different outcomes.
# CI95% and CI68% are estimated only at non-critical values.
# Use at own risk at extremal values (e.g. S1=0, S2=0, I1>0, Ro>100).

# Input
# Io, I1 - Infected animals (start and end state)
# So, S1, S2 - Susceptible animals (start and end state for the experiment 1 and 2)
I0<-6
I1<-0
S0<-6
S1<-1
S2<-1

# Check-up of input
if (I0<1) stop('Wrong input (Io).')
if (I1<0) stop('Wrong input (I1).')
if (S0<1) stop('Wrong input (So).')
if (S1>S0) stop('Wrong input (S1).')
if (S2>S0) stop('Wrong input (S2).')

#FS (Final size algorithm for 2 similar experiments)
FSx2<-function(S0, S1, S2, I0, I1)
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
    # Besides this it has inverted notation for raws
    # Row 0, column 0 should be adressed as [N+1, 0]
    # Row N, column S0 should be addressed as [1, S0+1]
    prob.m<-matrix(0, nrow=N+1, ncol=S0+1, dimnames=list(c(N:0),c(0:S0)))
    # Filling with initial state (P=1)
    prob.m[N-I0+1, S0+1]<-1
    # beta is used as coefficient of probability. We can calculate it only once
    beta<-N/(RR*S0+N)
    # Filling the first column
    for (i in 1:I0)
      {
        prob.m[N-I0+i+1, S0+1]<-prob.m[N-I0+i, S0+1]*beta
      }
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
        beta2<- N/(RR*(S.temp-1)+N)
        prob.m[S.temp+I.temp+1, S.temp]<-prob.m[S.temp+I.temp+2, S.temp+1]*RR*beta+prob.m[S.temp+I.temp, S.temp]*beta2
        }
      # Filling lower two lines
      # It can be no execution of upper cycle. But beta should be calculated at least once.
      beta2<- N/(RR*(S.temp-1)+N)
      prob.m[N, S.temp]<-prob.m[N-1, S.temp]*beta2
      prob.m[N+1, S.temp]<-prob.m[N, S.temp]*beta2
      }
    # Internal check-up. P=1
    if ((round(as.numeric(sum(prob.m[N+1,])), 2))!=1) stop('Internal error (P<>1).')
    # End of creating matrix. We returning Prob(R) at given S1.
    return(prob.m[N+1,])
  }
  # Making new matrix with probabilites of two experiments combined.
  prob2<-function(RR, S1, S2)
    {
    # Taking only lower line from previous table of probabilites for 1 experiment.
    prob.line<-prob(RR, S0)
    # Creating matrix of probabilities
    # Matrix contains row -1 and column -1 with initial values. So it is bigger than So+zero value x So+zero value.
    # Row 0, column 0 should be adressed as [2, 2]
    prob.m2<-matrix(0, nrow=S0+2, ncol=S0+2, dimnames=list(c(-1:S0),c(-1:S0)))
    # Filling first row and firts column.
    prob.m2[1, 2:(S0+2)]<-prob.line
    prob.m2[2:(S0+2), 1]<-prob.line
    # Filling all other wells.
    for (i in 2:(S0+2))
      {
        for (k in 2:(S0+2)) prob.m2[k, i]<-prob.m2[k, 1]*prob.m2[1, i]
      }
    # Internal check-up. P=1
    if ((round(as.numeric(sum(prob.m2[2:(S0+2), 2:(S0+2)])), 2))!=1) stop('Internal error (P for x2<>1).')
    # End of creating matrix. We returning Prob(R) at given S1 and S2.
    return(prob.m2[S1+2, S2+2])
    }
  # Optimization of R to get Pmax for S1 and S2
  optimized.prob<-function(RR) prob2(RR, S1, S2)
  RR<-optimize(optimized.prob, interval=0:1000, maximum=TRUE)$maximum
  # Optimization of R to get lower CI95
  lowerCI.prob<-function(RR)
  {
    CI<-0
    for (k in 0:(S1+S2))
    {
      for (i in 0:(S1+S2-k)) if (i<(S0+1)) CI<-CI+prob2(RR, i, k)
    }
    if (CI>0.025) CI<--CI
    return (CI)
  }
  lowerCI<-optimize(lowerCI.prob, interval=0:1000, maximum=TRUE)$maximum
  # Optimization of R to get higher CI95
  higherCI.prob<-function(RR)
  {
    CI<-0
    for (k in S0:0)
    {
      for (i in S0:0)
        {  
        if (S0-S1-S2+1<k) CI<-CI+prob2(RR, i, k)
        else if (S0-S1-S2+1<i+k) CI<-CI+prob2(RR, i, k)
        }
    }
    if (CI>0.025) CI<--CI
    return (CI)
  }
  higherCI<-optimize(higherCI.prob, interval=0:1000, maximum=TRUE)$maximum
  # Optimization of R to get lower CI68
  lowerCI68.prob<-function(RR)
  {
    CI<-0
    for (k in 0:(S1+S2))
    {
      for (i in 0:(S1+S2-k)) if (i<(S0+1)) CI<-CI+prob2(RR, i, k)
    }
    if (CI>0.16) CI<--CI
    return (CI)
  }
  lowerCI68<-optimize(lowerCI68.prob, interval=0:1000, maximum=TRUE)$maximum
  # Optimization of R to get higher CI68
  higherCI68.prob<-function(RR)
  {
    CI<-0
    for (k in S0:0)
    {
      for (i in S0:0)
      {  
        if (S0-S1-S2+1<k) CI<-CI+prob2(RR, i, k)
        else if (S0-S1-S2+1<i+k) CI<-CI+prob2(RR, i, k)
      }
    }
    if (CI>0.16) CI<--CI
    return (CI)
  }
  higherCI68<-optimize(higherCI68.prob, interval=0:1000, maximum=TRUE)$maximum
  # Internal check-up of initial I0, S1, S2 and final R values
  if (RR<0) stop('Internal error (Ro<0).')
  if (S1*S2==0) return(cat(RR))
  if (I1==0) return(cat((RR), " CI95% [", lowerCI, "-", higherCI, "]", " CI68% [", lowerCI68, "-", higherCI68, "]\n")) else return (cat("FS:         Cannot be estimated (I1>0).\n"))
}
#End FS (Final size algorithm)

cat(FSx2(S0, S1, S2, I0, I1))
