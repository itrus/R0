S<-c(11, 9, 7, 4, 4, 4, 4, 4, 4, 4)
I<-c(11, 13, 16, 19, 18, 18, 15, 10, 6, 2)
C<-c(2, 2, 3, 0, 0, 0, 0, 0, 0, 0)
N<-c(22, 22, 23, 23, 23, 23, 23, 23, 23, 23)

x<-glm(C/S~log(I/N), family=quasibinomial(link=cloglog), weights=-log(I/N))
error95<-qnorm(0.975)*x$deviance[1]/sqrt(length(N))
error68<-qnorm(0.84)*x$deviance[1]/sqrt(length(N))

print(exp(x$coefficient[1]))
cat("CI95% [", exp(x$coefficient[1]-error95), "-", exp(x$coefficient[1]+error95), "]")
cat("CI68% [", exp(x$coefficient[1]-error68), "-", exp(x$coefficient[1]+error68), "]")
