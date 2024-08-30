weibull <- function(time, alpha = 3.4, beta = 39.34) {
  pow <- -(time/beta)^alpha
  y <- exp(pow)
  return(y)
}
get_smc_profile<-function(start_sim,end_sim,SMC_start,nround,gap){
  prop_prof<-c(rep(1-weibull(1:gap),nround-1),1-weibull(1:100))
  return(list(
    SMC_times=c(start_sim,seq(SMC_start,SMC_start+length(prop_prof)-1,by=1),end_sim),
    SMC_vals=c(1,prop_prof,1)
  ))
}

smc_prof<-get_smc_profile(0,1800,900,4,30)
plot(smc_prof$SMC_times,smc_prof$SMC_vals)

smc_prof
