
rm(list = ls())
gc()
set.seed(1954)

# adjust to your working directory
setwd("~/Code/adaptive_malt/R")
source("pseudo_code.R")
source("step_size.R")


d=50
sigma=((d:1)/d)^(1/2)
init=rnorm(d)*sigma
# init=rep(0,d);init[2]=10
init = rnorm(d, mean = 1, sd = 1)
U=function(x){sum(0.5*x^2/sigma^2)}
grad=function(x){x/sigma^2}
n=10^4

output=adaptive_malt(n,init)

plot(output$eig_values,type="l")
plot(output$means,type="l")


# Test initial step size.
h0 = init_stepsize(U, grad, pos0 = init, tau0 = 1.8, w0 = init, h0 = 1)
output=adaptive_malt(n,init,h0)

plot(output$eig_values,type="l")
plot(output$means,type="l")
