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


# hyperparameters
eta_tau=0.05
eta_h=0.05
eta_w=3
kappa=8

osam=function(U,grad,pos,tau,h,w,m,iter){
  eigen_max=sqrt(sum(w^2))
  eigen_vector=w/eigen_max
  g=2/sqrt(eigen_max) #look at the meads paper to improve
  eta=exp(-g*h)
  zeta=sqrt(1-eta^2)
  half=h/2
  small=h^2/8
  L=ceiling(tau/h)
  steps=L+1
  i=iter
  if(i<100){
    L=1
  }
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

    #update m
    eta_m=1/(ceiling(i/kappa)+1)
    m=(1-eta_m)*m+eta_m*x
    #update w
    w=w*(i-eta_w)/(i+1)+x*(sum(x*eigen_vector))*(eta_w+1)/(i+1)
    #update tau
    pos_c=sum((pos-m)*eigen_vector)
    prop_c=sum((prop-m)*eigen_vector)
    v_c=sum(v*eigen_vector)
    diff_sq=((prop_c)^2-(pos_c)^2)
    noisy_grad=4*diff_sq*prop_c*v_c-(1/tau)(diff_sq)^2
    #update h

    return(list(pos=x,tau=tau,h=h,w=w,m=m,iter=i+1))
}

