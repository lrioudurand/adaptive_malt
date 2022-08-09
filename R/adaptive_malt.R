adaptive_malt=function(init,U,grad,n,g,h,L,warm=FALSE){
  d=length(init)
  if(d>=2){}else{stop("length(init) should be greater or equal than 2.")}
  if(n>=1){n=floor(n)}else{warning("The number of samples should be a positive integer, setting n=10000 instead.");n=10^4}
  if(g>=0){g=as.numeric(g)}else{warning("The friction/damping should be a non-negative real number, setting g=0 instead.");g=0}
  if(h>0){h=as.numeric(h)}else{stop("The time step should be a positive real number.")}
  if(L>=1){L=floor(L)}else{warning("The number of steps should be a positive integer, setting L=1 instead.");L=1}
  eta=exp(-g*h)
  zeta=sqrt(1-eta^2)
  half=h/2
  small=h^2/8
  steps=L+1
  chain=matrix(0,nrow=n,ncol=d)
  malt_x=init
  #warm up (doubles the computational time): draw unadjusted trajectories for n/2 iterations then adjusted trajectories for n/2 iterations.
  if(warm==T){
    x=malt_x
    grad_x=grad(x);grad_malt_x=grad_x
    gradsq_x=sum(grad_x^2);gradsq_malt_x=gradsq_x
    U_x=U(x);U_malt_x=U_x
    norm_draws=array(stats::rnorm(n*d*(L+1)),dim=c(n,d,L+1))
    expo_draws=stats::rexp(n)
    for(i in 1:(n-1)){
      v=norm_draws[i,,steps]
      Delta=0
      for(j in 1:L){
        v=v-half*grad_x
        x=x+h*v
        grad_y=grad(x)
        Delta=Delta-half*sum(v*(grad_x+grad_y))
        v=v-half*grad_y
        v=eta*v+zeta*norm_draws[i,,j]
        grad_x=grad_y
      }
      gradsq_x=sum(grad_x^2)
      U_x=U(x)
      Delta=small*(gradsq_x-gradsq_malt_x)+U_x-U_malt_x+Delta
      if(i<n/2){Delta=0}
      if(expo_draws[i]<Delta){x=chain[i,];grad_x=grad_malt_x;gradsq_x=gradsq_malt_x;U_x=U_malt_x}else{grad_malt_x=grad_x;gradsq_malt_x=gradsq_x;U_malt_x=U_x}
      chain[i+1,]=x
    }
    chain[1,]=chain[n,]
  }

  x=chain[1,]
  grad_x=grad(x);grad_malt_x=grad_x
  gradsq_x=sum(grad_x^2);gradsq_malt_x=gradsq_x
  U_x=U(x);U_malt_x=U_x
  norm_draws=array(stats::rnorm(n*d*(L+1)),dim=c(n,d,L+1))
  expo_draws=stats::rexp(n)
  for(i in 1:(n-1)){
    v=norm_draws[i,,steps]
    Delta=0
    for(j in 1:L){
      v=v-half*grad_x
      x=x+h*v
      grad_y=grad(x)
      Delta=Delta-half*sum(v*(grad_x+grad_y))
      v=v-half*grad_y
      v=eta*v+zeta*norm_draws[i,,j]
      grad_x=grad_y
    }
    gradsq_x=sum(grad_x^2)
    U_x=U(x)
    Delta=small*(gradsq_x-gradsq_malt_x)+U_x-U_malt_x+Delta
    if(expo_draws[i]<Delta){x=chain[i,];grad_x=grad_malt_x;gradsq_x=gradsq_malt_x;U_x=U_malt_x}else{grad_malt_x=grad_x;gradsq_malt_x=gradsq_x;U_malt_x=U_x}
    chain[i+1,]=x
  }
  output=list(samples=chain, draw=chain[n,], accept=mean(chain[2:n,1]!=chain[1:(n-1),1]),param=list(g=g,h=h,L=L))
  class(output)=c("malt")
  return(output)
}
