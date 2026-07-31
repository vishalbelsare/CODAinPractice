### chiPower function
chiPower <- function(X, close=TRUE, power=1, chi=TRUE, BoxCox=TRUE, CLR=TRUE)
{
# X: the compositional data matrix (it is closed in case)
# close: close the data after powering
# power: power of the transformation
# chi: apply chi-square standardization
# BoxCox: apply Box-Cox style of transformation
# CLR: translate columns so that convergence is to CLR transformation
foo <- as.matrix(X)
foo <- foo / rowSums(foo)
foo <- foo^power
if(close) foo <- foo / rowSums(foo)
if(chi) foo <- sweep(foo, 2, sqrt(colMeans(foo)), FUN="/")
if(BoxCox & !chi) foo <- (1/power)*(ncol(X)*foo - 1)
if(BoxCox & chi) foo <- (1/power)*(sqrt(ncol(X))*foo - 1)
if (BoxCox & chi & CLR) foo <- foo + rep(colMeans(foo), each=nrow(X))
X.chiPower <- foo
X.chiPower
}

