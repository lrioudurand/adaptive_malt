#' osam
#'
#' @description Implements one step of the adaptive version of the sampling algorithm: Metropolis Adjusted Langevin Trajectories (MALT) as described in Riou-Durand and Vogrinc (2022).
#'
#' @param U A potential function to return the log-density of the distribution to be sampled from, up to an additive constant. It should input a real vector of the same length as \code{init} and output a scalar.
#' @param grad A function to return the gradient of the potential. It should input and output a real vector of the same length as \code{init}.
#' @param pos Real vector. current position.
#' @param tau The current length of a trajectory. Positive real number. The choice L=1 boils down to the Metropolis Adjusted Langevin Algorithm.
#' @param h The current time step. Positive real number.
#' @param w The current eigen vector.
#' @param m The current mean.
#' @param iter The iteration number.
#' @param H_bar Dual averaging statistics, used to adapt step size. Real number.
#' @param h_bar Proposed step size for sampling phase. Positive real number.
#' @param mu Reference initial step size for dual averaging. Positive real number.
#' @param n_warmup Iteration number after which adaptation should stop.
#' @param delta Target mean acceptance rate for proposal. Real number between 0 and 1.


# hyperparameters
eta_tau=0.05
eta_h=0.05
eta_w=3
kappa=8

# hyperparameters for dual averaging of step size
gamma = 0.05
t0 = 10
kappa_stepsize = 0.75

# test
# {
# pos = init
# tau = 1.8
# h = 0.2
# w = init
# m = init
# iter = 1
# H_bar = H_bar0
# h_bar = h_bar0
# h0=0.2
# mu = log(10 * h0)
# n_warmup = 1000
# delta = 0.8
# }


osam=function(U,grad,pos,tau,h,w,m,iter,
              H_bar, h_bar, mu, n_warmup, delta){
  eigen_max=sqrt(sum(w^2))
  eigen_vector=w/eigen_max
  g=2/sqrt(eigen_max) #look at the meads paper to improve
  eta=exp(-g*h)
  zeta=sqrt(1-eta^2)
  half=h/2
  small=h^2/8
  i=iter
  if(i<100){
    tau=h
  }
  L=ceiling(tau/h)
  steps=L+1
  #malt_x=pos
  x=pos
  grad_x=grad(x);grad_malt_x=grad_x
  gradsq_x=sum(grad_x^2);gradsq_malt_x=gradsq_x
  U_x=U(x);U_malt_x=U_x
  norm_draws=array(stats::rnorm(d*(L+1)),dim=c(d,L+1))
  expo_draw=stats::rexp(1)
  # malt update
  v=norm_draws[,steps]
  Delta=0
  for(j in 1:L){
    v=v-half*grad_x
    x=x+h*v
    grad_y=grad(x)
    Delta=Delta-half*sum(v*(grad_x+grad_y))
    v=v-half*grad_y
    v=eta*v+zeta*norm_draws[,j]
    grad_x=grad_y
  }
  prop=x
  gradsq_x=sum(grad_x^2)
  U_x=U(x)
  Delta=small*(gradsq_x-gradsq_malt_x)+U_x-U_malt_x+Delta
  if(expo_draw<Delta){x=pos;grad_x=grad_malt_x;gradsq_x=gradsq_malt_x;U_x=U_malt_x}else{grad_malt_x=grad_x;gradsq_malt_x=gradsq_x;U_malt_x=U_x}

  if (i <= n_warmup) {
    #update m
    eta_m=1/(ceiling(i/kappa)+1)
    m=(1-eta_m)*m+eta_m*x

    #update w
    x_c=x-m
    if(i>eta_w){
      w=w*(i-eta_w)/(i+1)+x_c*(sum(x_c*eigen_vector))*(eta_w+1)/(i+1)
    }else{w=x_c}

    #update tau
    pos_c=sum((pos-m)*eigen_vector)
    prop_c=sum((prop-m)*eigen_vector)
    v_c=sum(v*eigen_vector)
    diff_sq=((prop_c)^2-(pos_c)^2)
    noisy_grad=4*diff_sq*prop_c*v_c-(1/tau)*(diff_sq)^2
    tau=1.8  # Q: does this mean tau does not get updated?

    #update h
    H_bar = (1 - 1 / (i + t0)) * H_bar + (delta - min(1, exp(-Delta))) / (i + t0)
    log_h = mu - sqrt(i) / gamma * H_bar
    log_h_bar = i^(-kappa_stepsize) * log_h + (1 - i^(-kappa_stepsize)) * log(h_bar)
    h = exp(log_h)
    h_bar = exp(log_h_bar)
    if (i == n_warmup) h = h_bar
  }

  return(list(pos=x,tau=tau,h=h,w=w,m=m,iter=i+1,Delta=Delta,
              H_bar = H_bar, h_bar = h_bar))
}

# hyperparameters for adaptive step size
h_bar0 = 1
H_bar0 = 0

adaptive_malt=function(n,init,h0=0.2,n_warmup = 2000, delta = 0.8){
  update=list(pos=init,tau=1.8,h=h0,w=init,m=init,iter=1,
              H_bar = H_bar0, h_bar = h_bar0)
  mu = log(10 * h0)

  # store results for analysis
  eig_values=rep(NA,n)
  means=rep(NA,n)
  h = rep(NA, n)  # stepsize
  alpha = rep(NA, n)  # acceptance prob
  Delta = rep(NA, n)  # change in energy

  for(i in 1:n){
    update=osam(U,grad,update$pos,update$tau,update$h,update$w,update$m,update$iter,
                update$H_bar, update$h_bar, mu, n_warmup, delta)
    eig_values[i]=sqrt(sum(update$w^2))
    means[i]=sum(update$w*update$m)/sqrt(sum(update$w^2))
    h[i] = update$h
    alpha[i] = min(1, exp(-update$Delta))
    Delta[i] = update$Delta
  }

  return(list(eig_values=eig_values,means=means, h = h, alpha = alpha,
              Delta = Delta))
}

