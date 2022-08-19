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
eta_noise=0.001
batch_size=2

alpha_adam=0.005
beta1_adam=0.9
beta2_adam=0.999
eps_adam=10^(-8)


osam=function(U,grad,pos,tau,h,w,m,m_sq,v_sq,cov_sq,deriv_sq,m_adam,v_adam,iter,noisy_grad){
  eigen_max=sqrt(sum(w^2))
  eigen_vector=w/eigen_max
  if(any(is.na(eigen_vector))){eigen_vector=w*0}
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
  if(i==100){tau=1.8*sqrt(eigen_max)}
  #if(is.na(tau)){tau=h}
  if(tau<h){tau=h}
  if(tau>10){tau=10}
  L=floor(tau/h)
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
  if(i>eta_w & eigen_max>0){
    w=w*(i-eta_w)/(i+1)+x_c*(sum(x_c*eigen_vector))*(eta_w+1)/(i+1)
  }else{w=x_c}
  #update tau
  pos_c=sum((pos-m)*eigen_vector)
  prop_c=sum((prop-m)*eigen_vector)
  v_c=sum(v*eigen_vector)
  v_0=sum(v_0*eigen_vector)
  diff_sq=((prop_c)^2-(pos_c)^2)
  sum_sq=((prop_c)^2+(pos_c)^2)
  #update m_sq, v_sq and rho
  m_sq=(1-eta_m)*m_sq+eta_m*(pos_c^2+prop_c^2)/2
  if(i>eta_w & eigen_max>0){
    v_sq=v_sq*(i-eta_w)/(i+1)+0.5*((pos_c^2-m_sq)^2+(pos_c^2-m_sq)^2)*(eta_w+1)/(i+1)
    cov_sq= cov_sq*(i-eta_w)/(i+1)+(pos_c^2-m_sq)*(prop_c^2-m_sq)*(eta_w+1)/(i+1)
    deriv_sq=deriv_sq*(i-eta_w)/(i+1)+diff_sq*(prop_c*v_c+pos_c*v_0)*(eta_w+1)/(i+1)
    rho=cov_sq/v_sq
    rho_prime=deriv_sq/v_sq
    #noisy_grad=rho_prime-((1-min(1,rho)^2))/(2*h*L)
    if(i%%batch_size!=0){
    noisy_grad=noisy_grad+2*diff_sq*(prop_c*v_c+pos_c*v_0)-(0.5*(1+min(1,rho))/(h*L))*(diff_sq)^2
    }else{noisy_grad=noisy_grad/batch_size}
    #if(i>20000){noisy_grad=2*v_sq*(rho_prime-((1-min(1,rho)^2))/(2*h*L))}
    #noisy_grad=2*diff_sq*(prop_c*v_c+pos_c*v_0)-(v_sq*(1-min(1,rho)^2)/(h*L))
  }else{noisy_grad=0}
  #noisy_grad=4*diff_sq*prop_c*v_c-(1/tau)*(diff_sq)^2
  #noisy_grad=4*diff_sq*(pos_c*v_0)-(1/tau)*(diff_sq)^2
  #=4*diff_sq*prop_c*v_c-(0.5/tau)*(diff_sq)^2
  #noisy_grad=2*diff_sq*(prop_c*v_c+pos_c*v_0)-(0.5*(1+min(1,rho))/(h*L))*(diff_sq)^2

  #noisy_grad=2*diff_sq*(prop_c*v_c+pos_c*v_0)-(v_sq/(h*L))*(1-min(1,rho)^2)
  #adam
  if(i%%batch_size==0){
  m_adam=beta1_adam*m_adam+(1-beta1_adam)*noisy_grad
  v_adam=beta2_adam*v_adam+(1-beta2_adam)*noisy_grad^2
  m_hat=m_adam/(1-beta1_adam^i)
  v_hat=v_adam/(1-beta2_adam^i)
  current_alpha=alpha_adam
  log_tau=log(tau)+current_alpha*m_hat/(sqrt(v_hat)+eps_adam)
  tau=exp(log_tau)
  noisy_grad=0
  }
  #tau=tau+alpha_adam*m_hat/(sqrt(v_hat)+eps_adam)
  #current_alpha=alpha_adam*5^(1-i/n)

  #noisy_grad=2*diff_sq*(prop_c*v_c+pos_c*v_0+sign(prop_c*v_c*pos_c*v_0)*sqrt(abs(prop_c*v_c*pos_c*v_0)))-(1/tau)*(diff_sq)^2
  #noisy_grad=2*diff_sq*(prop_c*v_c+pos_c*v_0)-(2/tau)*(m_sq-(pos_c*prop_c)^2)
  #noisy_grad=2*diff_sq*(prop_c*v_c+pos_c*v_0)-(0.75/tau)*(diff_sq)^2
  #noisy_grad=2*diff_sq*(prop_c*v_c+pos_c*v_0)-(1/tau)*(diff_sq)^2-(0.5/tau^2)
  #noisy_grad=-4*(pos_c)^2*prop_c*v_c-(1/tau)*(diff_sq)^2
  #noisy_grad=-4*((pos_c)^2-eigen_max)*prop_c*v_c-(1/tau)*(diff_sq)^2
  #noisy_grad=4*diff_sq*prop_c*v_c*(1/tau)-(1/tau^2)*(diff_sq)^2
  #log_tau=log(tau)+eta_noise*(noisy_grad)
  #log_tau=log(tau)+eta_noise*exp(-2*i/n)*(noisy_grad)
  #log_tau=log(tau)+eta_noise*(noisy_grad)

  #tau=tau+eta_noise*(noisy_grad)
  #update h
  h=0.05


  #storing noisy_grad
  return(list(pos=x,tau=tau,h=h,w=w,m=m,m_sq=m_sq,v_sq=v_sq,cov_sq=cov_sq,deriv_sq=deriv_sq,m_adam=m_adam,v_adam=v_adam,iter=i+1,noisy_grad=noisy_grad))
}




d=50
sigma=((d:1)/d)^(1/2)
init=rnorm(d)*sigma
#init=rep(0,d);init[2]=10
U=function(x){sum(0.5*x^2/sigma^2)}
grad=function(x){x/sigma^2}
n=10000



adaptive_malt=function(n,init){
  update=list(pos=init,tau=1.8,h=0.2,w=init,m=init,m_sq=sum(init^2),v_sq=0,cov_sq=0,deriv_sq=0,m_adam=0,v_adam=0,iter=1,noisy_grad=0)
  eig_values=rep(NA,n)
  means=rep(NA,n)
  tau_values=rep(NA,n)
  chain=rep(NA,n)
  noisyg=rep(NA,n)
  rho_values=rep(NA,n)
  rho_prime_values=rep(NA,n)
  for(i in 1:n){
    update=osam(U,grad,update$pos,update$tau,update$h,update$w,update$m,update$m_sq,update$v_sq,update$cov_sq,update$deriv_sq,update$m_adam,update$v_adam,update$iter,update$noisy_grad)
    eig_values[i]=sqrt(sum(update$w^2))
    means[i]=sum(update$w*update$m)/sqrt(sum(update$w^2))
    tau_values[i]=update$tau
    chain[i]=update$pos[1]
    noisyg[i]=sqrt(sum(update$noisy_grad)^2)
    rho_values[i]=update$cov_sq/update$v_sq
    if(is.na(rho_values[i])){rho_values[i]=0}
    rho_prime_values[i]=update$deriv_sq/update$v_sq
    if(is.na(rho_prime_values[i])){rho_prime_values[i]=0}
  }
  return(list(eig_values=eig_values,means=means,tau_values=tau_values,chain=chain,noisyg=noisyg,last_eigen_vector=update$w,rho_values=rho_values,rho_prime_values=rho_prime_values))
}

output=adaptive_malt(n,init)

plot(output$eig_values[1:150],type="l")
plot(output$eig_values,type="l")
plot(output$means,type="l")
plot(output$rho_values,type="l")
#plot(output$rho_values[2000:50000],type="l")
output$rho_values[n]

plot(output$chain,type="l")
plot(output$noisyg,type="l")

plot(output$tau[1:150],type="l")
plot(output$tau[1:1000],type="l")
plot(output$tau,type="l")



tau_bar=0
tau_average=rep(NA,n)
for(i in 1:n){
  tau_bar=beta2_adam*tau_bar+(1-beta2_adam)*output$tau[i]
  tau_average[i]=tau_bar/(1-beta2_adam^i)
}
plot(tau_average,type="l")
var(output$noisyg)
var(output$noisyg[(1:n/batch_size)*batch_size])
var(output$tau)
tau_average[n]
