# Rejection Sampler - COM-Poisson
# Algorithm 2 - Benson & Friel (2021)
# Input: theta = (mu, nu)
# Return sample y

com_sampler = function(mu, nu){
  if(nu >= 1){
    #Enveloping Bound
    #B_fg = (mu^floor(mu)/factorial(floor(mu)))^(nu-1)
    log_B_fg = (nu-1)*(floor(mu)*log(mu)-lfactorial(floor(mu)))
    repeat{
      y_sam = rpois(1,mu) #Proposal
      #alpha = (mu^y_sam/factorial(y_sam))^nu/(B_fg*(mu^(y_sam)/factorial(y_sam)))
      log_alpha = nu*(y_sam*log(mu)-lfactorial(y_sam))-(log_B_fg + y_sam*log(mu) - lfactorial(y_sam))
      if(runif(1) <= exp(log_alpha)) break
    }
  } else{
    #Enveloping Bound
    p = (2*nu)/(2*nu*mu + 1 + nu)
    ratio = floor(mu/((1-p)^(1/nu)))
    #B_fg = (1/p)*(mu^(nu*ratio))/((1-p)^ratio*(factorial(ratio))^nu)
    log_B_fg = -log(p) + nu*ratio*log(mu) - (ratio*log(1-p) + nu*lfactorial(ratio))
    repeat{
      u0 = runif(1)
      y_sam = floor(log(u0)/log(1-p)) #Proposal
      #alpha = ((mu^y_sam/factorial(y_sam))^nu)/(B_fg*(1-p)^(y_sam)*p)
      log_alpha = nu*(y_sam*log(mu)-lfactorial(y_sam)) - (log_B_fg + y_sam*log(1-p) + log(p))
      if(runif(1) <= exp(log_alpha)) break
    }
  }
  return(y_sam)
}
