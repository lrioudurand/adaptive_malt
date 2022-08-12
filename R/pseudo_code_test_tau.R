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
eta_noise=0.005

osam=function(U,grad,pos,tau,h,w,m,iter){
  eigen_max=sqrt(sum(w^2))
  eigen_vector=w/eigen_max
  g=2/sqrt(eigen_max) #look at the meads paper to improve
  eta=exp(-g*h)
  zeta=sqrt(1-eta^2)
  half=h/2
  small=h^2/8
  i=iter
  if(i<100){
    #tau=h
    tau=h
  }
  if(i==100){tau=1.2*sqrt(eigen_max)}
  if(tau<h){tau=h}
  if(tau>100){tau=100}
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
    v_0=v
    for(j in 1:L){
      v=eta*v+zeta*norm_draws[,j]
      v=v-half*grad_x
      x=x+h*v
      grad_y=grad(x)
      Delta=Delta-half*sum(v*(grad_x+grad_y))
      v=v-half*grad_y
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
    x_c=x-m
    if(i>eta_w){
    w=w*(i-eta_w)/(i+1)+x_c*(sum(x_c*eigen_vector))*(eta_w+1)/(i+1)
    }else{w=x_c}
    #update tau
    pos_c=sum((pos-m)*eigen_vector)
    prop_c=sum((prop-m)*eigen_vector)
    v_c=sum(v*eigen_vector)
    v_0=sum(v_0*eigen_vector)
    diff_sq=((prop_c)^2-(pos_c)^2)
    #noisy_grad=4*diff_sq*prop_c*v_c-(1/tau)*(diff_sq)^2
    #noisy_grad=4*diff_sq*(pos_c*v_0)-(1/tau)*(diff_sq)^2
    noisy_grad=2*diff_sq*(prop_c*v_c+pos_c*v_0)-(1/tau)*(diff_sq)^2
    #noisy_grad=-4*(pos_c)^2*prop_c*v_c-(1/tau)*(diff_sq)^2
    #noisy_grad=-4*((pos_c)^2-eigen_max)*prop_c*v_c-(1/tau)*(diff_sq)^2
    #noisy_grad=4*diff_sq*prop_c*v_c*(1/tau)-(1/tau^2)*(diff_sq)^2
    #log_tau=log(tau)+eta_noise*(noisy_grad)
    log_tau=log(tau)+eta_noise*exp(-2*i/n)*(noisy_grad)
    #log_tau=log(tau)+eta_noise*(noisy_grad)
    tau=exp(log_tau)
    #tau=tau+eta_noise*(noisy_grad)
    #update h
    h=0.05


    #storing noisy_grad
    return(list(pos=x,tau=tau,h=h,w=w,m=m,iter=i+1,noisy_grad=noisy_grad))
}




d=50
sigma=((d:1)/d)^(1/2)
init=rnorm(d)*sigma
#init=rep(0,d);init[2]=10
U=function(x){sum(0.5*x^2/sigma^2)}
grad=function(x){x/sigma^2}
n=10^4



adaptive_malt=function(n,init){
  update=list(pos=init,tau=1.8,h=0.2,w=init,m=init,iter=1)
  eig_values=rep(NA,n)
  means=rep(NA,n)
  tau_values=rep(NA,n)
  chain=rep(NA,n)
  noisyg=rep(NA,n)
  for(i in 1:n){
    update=osam(U,grad,update$pos,update$tau,update$h,update$w,update$m,update$iter)
    eig_values[i]=sqrt(sum(update$w^2))
    means[i]=sum(update$w*update$m)/sqrt(sum(update$w^2))
    tau_values[i]=update$tau
    chain[i]=update$pos[1]
    noisyg[i]=sqrt(sum(update$noisy_grad)^2)
  }
  return(list(eig_values=eig_values,means=means,tau_values=tau_values,chain=chain,noisyg=noisyg,last_eigen_vector=update$w))
}

output=adaptive_malt(n,init)

plot(output$eig_values,type="l")
plot(output$means,type="l")

plot(output$chain,type="l")
plot(output$noisyg,type="l")

plot(output$tau[1:150],type="l")
plot(output$tau[1:1000],type="l")
plot(output$tau,type="l")

output$last_eigen_vector
mean(output$tau[2000:n])
var(output$tau[2000:n])
var(output$noisyg[2000:n])