# R-GMM-bandwidth
this project provides an implementation of the bandwidth choice for robust GMM estimation as described in Wilhelm (2015)



The =R= function =bwWilhelm= computes the optimal bandwidth for HAC-robust GMM estimation as proposed in Wilhelm (2015, Econometric Theory, 31, 1054–1077, "Optimal Bandwidth Selection for Robust Generalized Method of Moments Estimation"). The arguments the function requires are identical to those of existing bandwidth selection methods such as the function =bwAndrews= in the =sandwich= package, except that the data matrix =x= (here the moment function evaluated at the data) must be an object of class =gmm=.

Here is an example of how to use the bandwidth selection procedure in two-step GMM estimation:

```R
rm(list = ls(all = TRUE))
set.seed(100)
library(gmm)

source("bwWilhelm.r")

# generate data
T <- 200		# sample size
p <- 2			# no. of parameters
l <- 5 		# no. of instruments (= no. of moments)
z <- matrix(rnorm(T*(l-1)), ncol=l-1)
x <- rowSums(z) + rnorm(T)
y <- x + rnorm(T)
dat <- cbind(y,1,x,1,z)

# moment functions
kron.IV <- function(A, B) {
	stopifnot(nrow(A)==nrow(B))
	if (is.vector(A)) A <- matrix(A, ncol=1)
	if (is.vector(B)) B <- matrix(B, ncol=1)
	B.enlarged <- B%x%t(rep(1,ncol(A)))
	A.enlarged <- t(rep(1,ncol(B)))%x%A
	return(A.enlarged*B.enlarged)
}
g <- function(b, dat) { 
	y <- dat[,1]
	x <- as.matrix(dat[,2:(p+1)])
	z <- dat[,(2+p):(p+l+1)]
	e <- y-x%*%b
	return(kron.IV(e,z))
}

# sample average of moment function derivatives
G <- function(b, dat) { 
	e <- -as.matrix(dat[,2:(p+1)])
	z <- dat[,(2+p):(p+l+1)]
	return(matrix(apply(kron.IV(e,z), 2, mean), byrow=TRUE, ncol=p))
}

# first step GMM estimator
print("first step estimation results:")
res1 <- gmm(g, dat, t0=c(0, 1), grad=G, type="twoStep", wmatrix="ident")
res1

# compute bandwidth for optimal HAC weighting matrix
print("Wilhelm (2015) bandwidth:")
optbw <- bwWilhelm(res1, order.by=NULL, kernel="Bartlett", approx="AR(1)", weights=1, prewhite=1, ar.method="ols")
optbw

# for comparison, compute Andrews (1991) bandwidth
print("Andrews (1991) bandwidth:")
bwAndrews(res1, order.by=NULL, kernel="Bartlett", approx="AR(1)", weights=1, prewhite=1, ar.method="ols")

# second step GMM estimator
print("second step estimation results:")
res2 <- gmm(g, dat, t0=c(0, 1), grad=G, type="twoStep", wmatrix="optimal", bw=optbw, kernel="Bartlett")
res2
```
