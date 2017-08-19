# The following code contains the equations for the allele frequency as 
# function of the selection coefficient, as provided by Boissinot et al. (2006)

# Load packages
library(Ryacas)
require("Ryacas")

# Set population size
n = 1000

# Set offset values for selection coefficient (not defined at s = 0) and 
# allele frequency (not defined at y = 1/(2n))
sOffSet <- 0.006
yOffSet <- 10^(-10)

# Generate a sequence of values for the selection coefficient
sVals <- c(seq(-0.1, -sOffSet, 0.001), seq(sOffSet, 0.1, 0.001))


#############
# I) y > 1/(2N) (where y denotes allele frequency)
#############

# 
# -(exp(-x)*(e^x-1)*(e^(2*n*x)-e^(2*n*y*x)))/(n*(y-1)*y*x*(e^(2*n*x)-1))
fa <- function(x, y) -(exp(-x)*(exp(x)-1)*(exp(2*n*x)-exp(2*n*y*x)))/(n*(y-1)*y*x*(exp(2*n*x)-1))
fa(0.1, 0.5)
x = -0.1
y = 0.5

# First derivative:
fa1 <- function(x, y) ((2*n*y-1)*e^(-x)*(((((2*n*y-2*n)*x-1)*e^x+(-2*n*y+2*n+1)*x+1)*e^(2*n*x)+(1-2*n*y*x)*e^x+(2*n*y-1)*x-1)*e^(2*n*y*x)+(e^x-x-1)*e^(4*n*x)+((2*n*x-1)*e^x+(1-2*n)*x+1)*e^(2*n*x)))/(2*n^2*(y-1)*y*x^2*(e^(2*n*x)-1)^2)

# Second derivative
fa2 <- function(x, y) {
  y *(exp(-x) * (((((8*n^3*y^3 + (-16*n^3 - 4*n^2)*y^2 + 
                (8*n^3 + 8*n^2)*y - 4*n^2)*x^2 + (-8*n^2*y^2 + (8*n^2 + 4*n)*y -
                 4*n)*x + 4*n*y - 2)*exp(x) + (-8*n^3*y^3 + (16*n^3 + 12*n^2)*y^2 +
                 (-8*n^3-16*n^2-6*n)*y+4*n^2+4*n+1)*x^2+(8*n^2*y^2+(-8*n^2-8*n)*y+4*n+2)*x-4*n*y+2)*exp(4*n*x)+(((-16*n^3*y^3+(16*n^3+8*n^2)*y^2+(8*n^3-8*n^2)*y-4*n^2)*x^2+(16*n^2*y^2+(-8*n^2-8*n)*y+4*n)*x-8*n*y+4)*exp(x)+(16*n^3*y^3+(-16*n^3-24*n^2)*y^2+(-8*n^3+16*n^2+12*n)*y+4*n^2-4*n-2)*x^2+(-16*n^2*y^2+(8*n^2+16*n)*y-4*n-4)*x+8*n*y-4)*exp(2*n*x)+((8*n^3*y^3-4*n^2*y^2)*x^2+(4*n*y-8*n^2*y^2)*x+4*n*y-2)*exp(x)+(-8*n^3*y^3+12*n^2*y^2-6*n*y+1)*x^2+(8*n^2*y^2-8*n*y+2)*x-4*n*y+2)*exp(2*n*y*x)+(-4*n*exp(x)+2*n*x^2+4*n*x+4*n)*exp(6*n*x)+((-8*n^3*x^2-8*n^2*x+8*n)*exp(x)+(8*n^3+8*n^2-4*n)*x^2+(8*n^2-8*n)*x-8*n)*exp(4*n*x)+((-8*n^3*x^2+8*n^2*x-4*n)*exp(x)+(8*n^3-8*n^2+2*n)*x^2+(4*n-8*n^2)*x+4*n)*exp(2*n*x)))/(2*n^2*(y-1)*y*x^3*(exp(2*n*x)-1)^3)
}
fa2(0.001, 0.5)

# Integrate second derivative numerically over y values from 1/(2*n) + yOffSet to 1
NIa <- sapply(sVals, function(s){
  integrate(function(y) fa2(s, y), lower = 1/(2*n) + yOffSet, upper = 1)$value
})

# Plot integral of second derivative for different selection coefficients
plot(sVals, NIa, type = "l")
plot(sVals, sapply(sVals, function(z) fa(z, 0.9)))

# Replace x by s and y by x in I) and multiply by x (this should provide the 
# mean frequency if integrated by x from 0 to 1)
#-x*(e^(-s)*(e^s-1)*(e^(2*n*s)-e^(2*n*x*s)))/(n*(x-1)*x*s*(e^(2*n*s)-1))

# Integral with respect to x
#-((e^s-1)*e^(2*n*s-s)*(ln(abs(x-1))+expintegral_e(1,-2*n*s*(x-1))))/(n*s*(e^(2*n*s)-1))

#############
# y < 1/(2N)
#############
# 
# -(e^(-x)*(e^x-1)*(e^(2*n*x)-e^(2*n*y*x)))/(n*(y-1)*y*x*(e^(2*n*x)-1))
fb <- function(x, y) -(exp(-x)*(exp(x)-exp(2*n*x))*(exp(2*n*y*x) - 1))/(n*(y-1)*y*x*(exp(2*n*x)-1))
fb(0.1, 0.5)
x = -0.1
y = 0.5

# First derivative:
fb1 <- function(x, y) (exp(-x)*((((2*n*y-1)*x-1)*exp(4*n*x)+(((2*n-2*n*y)*x+1)*exp(x)+(-2*n*y-2*n+1)*x+1)*exp(2*n*x)+(2*n*y*x-1)*exp(x))*exp(2*n*y*x)+(x+1)*exp(4*n*x)+((-2*n*x-1)*exp(x)+(2*n-1)*x-1)*exp(2*n*x)+exp(x)))/(n*(y-1)*y*x^2*(exp(2*n*x)-1)^2)

# Second derivative
fb2 <- function(x, y) {
  (exp(-x)*((4*n*((2*n*y-1)*x-1)*exp(4*n*x)+(2*n*y-1)*exp(4*n*x)+(((2*n-2*n*y)*x+1)*exp(x)+(2*n-2*n*y)*exp(x)-2*n*y-2*n+1)*exp(2*n*x)+2*n*(((2*n-2*n*y)*x+1)*exp(x)+(-2*n*y-2*n+1)*x+1)*exp(2*n*x)+(2*n*y*x-1)*exp(x)+2*n*y*exp(x))*exp(2*n*y*x)+2*n*y*(((2*n*y-1)*x-1)*exp(4*n*x)+(((2*n-2*n*y)*x+1)*exp(x)+(-2*n*y-2*n+1)*x+1)*exp(2*n*x)+(2*n*y*x-1)*exp(x))*exp(2*n*y*x)+4*n*(x+1)*exp(4*n*x)+exp(4*n*x)+((-2*n*x-1)*exp(x)-2*n*exp(x)+2*n-1)*exp(2*n*x)+2*n*((-2*n*x-1)*exp(x)+(2*n-1)*x-1)*exp(2*n*x)+exp(x)))/(n*(y-1)*y*x^2*(exp(2*n*x)-1)^2)-(4*exp(-x)*exp(2*n*x)*((((2*n*y-1)*x-1)*exp(4*n*x)+(((2*n-2*n*y)*x+1)*exp(x)+(-2*n*y-2*n+1)*x+1)*exp(2*n*x)+(2*n*y*x-1)*exp(x))*exp(2*n*y*x)+(x+1)*exp(4*n*x)+((-2*n*x-1)*exp(x)+(2*n-1)*x-1)*exp(2*n*x)+exp(x)))/((y-1)*y*x^2*(exp(2*n*x)-1)^3)-(exp(-x)*((((2*n*y-1)*x-1)*exp(4*n*x)+(((2*n-2*n*y)*x+1)*exp(x)+(-2*n*y-2*n+1)*x+1)*exp(2*n*x)+(2*n*y*x-1)*exp(x))*exp(2*n*y*x)+(x+1)*exp(4*n*x)+((-2*n*x-1)*exp(x)+(2*n-1)*x-1)*exp(2*n*x)+exp(x)))/(n*(y-1)*y*x^2*(exp(2*n*x)-1)^2)-(2*exp(-x)*((((2*n*y-1)*x-1)*exp(4*n*x)+(((2*n-2*n*y)*x+1)*exp(x)+(-2*n*y-2*n+1)*x+1)*exp(2*n*x)+(2*n*y*x-1)*exp(x))*exp(2*n*y*x)+(x+1)*exp(4*n*x)+((-2*n*x-1)*exp(x)+(2*n-1)*x-1)*exp(2*n*x)+exp(x)))/(n*(y-1)*y*x^3*(exp(2*n*x)-1)^2)}
fb2(-0.001, 0.5)

# Integrate second derivative numerically over y values from 1/(2*n) + yOffSet to 1
NIb <- sapply(sVals, function(s){
  integrate(function(y) fb2(s, y), lower = 0, upper = 1/(2*n) - yOffSet)$value
})

# Plot integral of second derivative for different selection coefficients
plot(sVals, NIb, type = "l")
plot(sVals, sapply(sVals, function(z) fb(z, 0.9)))
plot(sVals, NIa + NIb)
segments(0, -100, 0, 100, lty = 2)
segments(-1, 0, 1, 0)
