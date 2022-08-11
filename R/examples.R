
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
n_warmup = 8000
delta = 0.8

output=adaptive_malt(n = n,init = init, n_warmup = n_warmup, delta = delta)

plot(output$eig_values,type="l")
plot(output$means,type="l")
plot(output$h, type = "l")
mean(output$alpha[n_warmup:n])  # if adaptation is good this should go to delta.

# Test initial step size.
h0 = init_stepsize(U, grad, pos0 = init, tau0 = 1.8, w0 = init, h0 = 1)
h0
output=adaptive_malt(n = n, init = init, h0 = h0, n_warmup = n_warmup, 
                     delta = delta)

plot(output$eig_values,type="l")
plot(output$means,type="l")
plot(output$h, type = "l")
mean(output$alpha[n_warmup:n])


# Check for divergent transitions (follow example in Inference Gym,
# https://github.com/tensorflow/probability/blob/main/spinoffs/inference_gym/notebooks/inference_gym_tutorial.ipynb)
divergence = output$Delta > 1000
sum(divergence)  # 3

# Check for divergences after warmup
sum(output$Delta[n_warmup:n] > 1000)  # 0
