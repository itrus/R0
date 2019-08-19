# S - susceptible (negative) animals on each day. It includes "contact" animals as well.
S<-c(11, 9, 7, 4, 4, 4, 4, 4, 4, 4)
# C - "contact" animals. They turn to positives (I) on the next day. But today they are contact animals that are counted both as C and S.
C<-c(2, 2, 3, 0, 0, 0, 0, 0, 0, 0)
# I - infected (positive) animals.
I<-c(11, 13, 16, 19, 18, 18, 15, 10, 6, 2)
#N - total number of animals (N = S + I + R; R - recovered animals)
N<-c(22, 22, 23, 23, 23, 23, 23, 23, 23, 23)
# sampling interval
interval<-1

x<-glm(C/S~log(I*interval/N), family=quasibinomial(link=cloglog), weights=-log(I/N))
error95<-qnorm(0.975)*x$deviance[1]/sqrt(length(N))
error68<-qnorm(0.84)*x$deviance[1]/sqrt(length(N))

#printing out BETA and confidence intervals of BETA
print(exp(x$coefficient[1]))
cat("CI95% [", exp(x$coefficient[1]-error95), "-", exp(x$coefficient[1]+error95), "]\n")
cat("CI68% [", exp(x$coefficient[1]-error68), "-", exp(x$coefficient[1]+error68), "]")
