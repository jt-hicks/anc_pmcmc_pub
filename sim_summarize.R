library("coda")
library(binom)
library(ggplot2)
library(bayesplot)
library(reshape2)
library(ggpubr)
library(zoo)
library(patchwork)
library(RColorBrewer)
source('utils.R')
sim10_runs$results()[[1]]$mcmc
##First 10 sim runs
sims_results_1to10 <- lapply(1:length(sim10_runs$results()),function(i){
  combine_runs(sim10_runs$results()[[i]],n_steps=1000,n_chains=2,burnin=200)
  })
sims_results_11to30 <- lapply(1:length(sim20_runs$results()),function(i){
  combine_runs(sim20_runs$results()[[i]],n_steps=1000,n_chains=2,burnin=200)
})
sims_results_1to30 <- append(sims_results_1to10,sims_results_11to30)
11/2*5/4
mean(sapply(1:30, function(i){sims_results_1to30[[i]]$run_time}))/3600

#Diagnostics by volatility and EIR
ar_df <- bind_rows(lapply(1:length(sims_results_1to30),function(i){
  ar <- as.data.frame(t(1 - coda::rejectionRate(as.mcmc(sims_results_1to30[[i]]$mcmc))))
  ar$true_volatility <- sim_dataset_100$variables[[i,'volatility']]
  ar$true_init_EIR <- sim_dataset_100$variables[[i,'init_EIR']]
  return(ar)
}))
ar_df_long <- melt(ar_df,id=c('true_volatility','true_init_EIR'))
ar_by_vol <- ggplot(ar_df_long)+
  geom_point(aes(x=true_volatility,y=value))+
  facet_wrap(.~variable)+
  coord_cartesian(ylim=c(0,NA))+
  labs(x='True Volatility',y='Acceptance Rate')
ar_by_eir <- ggplot(ar_df_long)+
  geom_point(aes(x=true_init_EIR,y=value))+
  facet_wrap(.~variable)+
  scale_x_log10()+
  coord_cartesian(ylim=c(0,NA))+
  labs(x='True initial EIR',y='Acceptance Rate')

ess_df <- bind_rows(lapply(1:length(sims_results_1to30),function(i){
  ess <- as.data.frame(t(coda::effectiveSize(as.mcmc(sims_results_1to30[[i]]$mcmc))))
  ess$true_volatility <- sim_dataset_100$variables[[i,'volatility']]
  ess$true_init_EIR <- sim_dataset_100$variables[[i,'init_EIR']]
  return(ess)
}))
ess_df_long <- melt(ess_df,id=c('true_volatility','true_init_EIR'))
ess_by_vol <- ggplot(ess_df_long)+
  geom_point(aes(x=true_volatility,y=value))+
  facet_wrap(.~variable)+
  coord_cartesian(ylim=c(0,NA))+
  labs(x='True Volatility',y='Effective Sample Size')
ess_by_eir <- ggplot(ess_df_long)+
  geom_point(aes(x=true_init_EIR,y=value))+
  facet_wrap(.~variable)+
  scale_x_log10()+
  coord_cartesian(ylim=c(0,NA))+
  labs(x='True Initial EIR',y='Effective Sample Size')

rhat_df <- bind_rows(lapply(1:length(sims_results_1to30),function(i){
  chain1 <- as.matrix(sims_results_1to30[[i]]$mcmc[1:800,])
  chain2 <- as.matrix(sims_results_1to30[[i]]$mcmc[801:1600,])
  chain_array <- array(c(chain1,chain2),dim=c(800,5,2),
                       dimnames = list(rownames(as.matrix(sims_results_1to30[[i]]$mcmc[1:800,])),
                                       colnames(as.matrix(sims_results_1to30[[i]]$mcmc[1:800,])),
                                       c('chain 1','chain 2'))
                       )

  rhat <- data.frame(log_prior = posterior::rhat(chain_array[,'log_prior',]),
                     log_posterior = posterior::rhat(chain_array[,'log_posterior',]),
                     log_likelihood = posterior::rhat(chain_array[,'log_likelihood',]),
                     log_init_EIR = posterior::rhat(chain_array[,'log_init_EIR',]),
                     volatility = posterior::rhat(chain_array[,'volatility',]))
  rhat$true_volatility <- sim_dataset_100$variables[[i,'volatility']]
  rhat$true_init_EIR <- sim_dataset_100$variables[[i,'init_EIR']]
  return(rhat)
}))
rhat_df_long <- melt(rhat_df,id=c('true_volatility','true_init_EIR'))
rhat_by_vol <- ggplot(rhat_df_long)+
  geom_hline(yintercept = 1.05,color='darkgrey',size=1,linetype='dashed')+
  geom_point(aes(x=true_volatility,y=value))+
  facet_wrap(.~variable)+
  coord_cartesian(ylim=c(1,NA))+
  labs(x='True Volatility',y='R-hat')
rhat_by_eir <- ggplot(rhat_df_long)+
  geom_hline(yintercept = 1.05,color='darkgrey',size=1,linetype='dashed')+
  geom_point(aes(x=true_init_EIR,y=value))+
  facet_wrap(.~variable)+
  scale_x_log10()+
  coord_cartesian(ylim=c(1,NA))+
  labs(x='True Initial EIR',y='R-hat')

time_df <- bind_rows(lapply(1:length(sims_results_1to30),function(i){
  ess <- as.data.frame(t(coda::effectiveSize(as.mcmc(sims_results_1to30[[i]]$mcmc))))
  ess$run_time <- sims_results_1to30[[i]]$run_time
  ess$true_volatility <- sim_dataset_100$variables[[i,'volatility']]
  ess$true_init_EIR <- sim_dataset_100$variables[[i,'init_EIR']]
  df <- ess %>%
    mutate(ess_per_sec = min(ess[1,1:5])/as.numeric(run_time),
           run_rate = as.numeric(run_time)/200/2000)
  return(df)
}))
mean(sapply(1:30, function(i){sims_results_1to30[[i]]$run_time}))/200/2000

eff_by_vol <- ggplot(time_df)+
  geom_point(aes(x=true_volatility,y=ess_per_sec*3600))+
  labs(x='True volatility',y='ESS per hour')+
  coord_cartesian(ylim=c(0,NA))
eff_by_eir <- ggplot(time_df)+
  geom_point(aes(x=true_init_EIR,y=ess_per_sec*3600))+
  scale_x_log10()+
  labs(x='True init_EIR',y='ESS per hour')+
  coord_cartesian(ylim=c(0,NA))
rate_by_vol <- ggplot(time_df)+
  geom_point(aes(x=true_volatility,y=run_rate))+
  labs(x='True volatility',y='Seconds per particle per step')+
  coord_cartesian(ylim=c(0,NA))
ggsave(filename='comp_time_vol.tiff',dpi=300,width=4,height=4,plot=rate_by_vol)
rate_by_eir <- ggplot(time_df)+
  geom_point(aes(x=true_init_EIR,y=run_rate))+
  scale_x_log10()+
  labs(x='True init_EIR',y='Seconds per particle per step')+
  coord_cartesian(ylim=c(0,NA))
ggsave(filename='comp_time_eir.tiff',dpi=300,width=4,height=4,plot=rate_by_eir)

#Comparison of real versus fitted variables
fit_df <- bind_rows(lapply(1:length(sims_results_1to30),function(i){
  ci_range_vol <- quantile(sims_results_1to30[[i]]$mcmc$volatility,probs = c(0.025,0.5,0.975))
  ci_range_eir <- quantile(sims_results_1to30[[i]]$mcmc$log_init_EIR,probs = c(0.025,0.5,0.975))
  df <- data.frame(volatility_median=ci_range_vol[[2]],
                   volatility_lower = ci_range_vol[[1]],
                   volatility_upper = ci_range_vol[[3]],
                   eir_median=ci_range_eir[[2]],
                   eir_lower = ci_range_eir[[1]],
                   eir_upper = ci_range_eir[[3]],
                   true_volatility = sim_dataset_100$variables[[i,'volatility']],
                   true_init_EIR = sim_dataset_100$variables[[i,'init_EIR']],
                   draw=factor(i))
  return(df)
}))

vol_fit <- ggplot(fit_df,aes(y=fct_reorder(draw,true_volatility)))+
  geom_point(aes(x=volatility_median))+
  geom_linerange(aes(xmin=volatility_lower,xmax=volatility_upper))+
  geom_point(aes(x=true_volatility),color='red')+
  labs(x='Volatility',y='Simulation')
eir_fit <- ggplot(fit_df,aes(y=fct_reorder(draw,true_init_EIR)))+
  geom_point(aes(x=exp(eir_median)))+
  geom_linerange(aes(xmin=exp(eir_lower),xmax=exp(eir_upper)))+
  geom_point(aes(x=true_init_EIR),color='red')+
  scale_x_log10()+
  labs(x='EIR',y='Simulation')

#traces
trace_df <- bind_rows(lapply(1:length(sims_results_1to30),function(i){
  df <- sims_results_1to30[[i]]$mcmc
  df$true_volatility <- sim_dataset_100$variables[[i,'volatility']]
  df$true_init_EIR <- sim_dataset_100$variables[[i,'init_EIR']]
  df$index <- c(1:nrow(df))
  df$draw <- factor(i)
  return(df)
}))
vol_trace <- ggplot(trace_df)+
  geom_line(aes(x=index,y=volatility))+
  facet_wrap(.~draw,scale='free_y')
eir_trace <- ggplot(trace_df)+
  geom_line(aes(x=index,y=log_init_EIR))+
  facet_wrap(.~fct_reorder(draw,true_init_EIR),scale='free_y')

##Informed runs
sims_results_1to10_informed <- lapply(1:length(sim10_informed_runs$results()),function(i){
  combine_runs(sim10_informed_runs$results()[[i]],n_steps=1000,n_chains=2,burnin=200)
})
sims_results_1to10_informed[[1]]$mcmc$
trace_df <- bind_rows(lapply(1:length(sims_results_1to10_informed),function(i){
  df <- sims_results_1to10_informed[[i]]$mcmc
  df$true_volatility <- sim_dataset_100$variables[[i,'volatility']]
  df$init_betaa <- sim_dataset_100$variables[[i,'init_betaa']]
  df$index <- c(1:nrow(df))
  df$draw <- factor(i)
  return(df)
}))
vol_trace <- ggplot(trace_df)+
  geom_line(aes(x=index,y=volatility))+
  facet_wrap(.~draw,scale='free_y')
eir_trace <- ggplot(trace_df)+
  geom_line(aes(x=index,y=init_betaa))+
  facet_wrap(.~fct_reorder(draw,true_init_EIR),scale='free_y')
ar_df <- bind_rows(lapply(1:length(sims_results_1to10_informed),function(i){
  ar <- as.data.frame(t(1 - coda::rejectionRate(as.mcmc(sims_results_1to10_informed[[i]]$mcmc))))
  ar$true_volatility <- sim_dataset_100$variables[[i,'volatility']]
  ar$true_init_EIR <- sim_dataset_100$variables[[i,'init_EIR']]
  return(ar)
}))
ar_df_long <- melt(ar_df,id=c('true_volatility','true_init_EIR'))
ar_by_vol <- ggplot(ar_df_long)+
  geom_point(aes(x=true_volatility,y=value))+
  facet_wrap(.~variable)+
  coord_cartesian(ylim=c(0,NA))+
  labs(x='True Volatility',y='Acceptance Rate')
ar_by_eir <- ggplot(ar_df_long)+
  geom_point(aes(x=true_init_EIR,y=value))+
  facet_wrap(.~variable)+
  scale_x_log10()+
  coord_cartesian(ylim=c(0,NA))+
  labs(x='True initial EIR',y='Acceptance Rate')

time_df <- bind_rows(lapply(1:length(sims_results_1to10_informed),function(i){
  ess <- as.data.frame(t(coda::effectiveSize(as.mcmc(sims_results_1to10_informed[[i]]$mcmc))))
  ess$run_time <- sims_results_1to10_informed[[i]]$run_time
  ess$true_volatility <- sim_dataset_100$variables[[i,'volatility']]
  ess$true_init_EIR <- sim_dataset_100$variables[[i,'init_EIR']]
  df <- ess %>%
    mutate(ess_per_sec = min(ess[1,1:5])/as.numeric(run_time),
           run_rate = as.numeric(run_time)/200/2000)
  return(df)
}))
eff_by_vol <- ggplot(time_df)+
  geom_point(aes(x=true_volatility,y=ess_per_sec*3600))+
  labs(x='True volatility',y='ESS per hour')+
  coord_cartesian(ylim=c(0,NA))
eff_by_eir <- ggplot(time_df)+
  geom_point(aes(x=true_init_EIR,y=ess_per_sec*3600))+
  scale_x_log10()+
  labs(x='True init_EIR',y='ESS per hour')+
  coord_cartesian(ylim=c(0,NA))
rate_by_vol <- ggplot(time_df)+
  geom_point(aes(x=true_volatility,y=run_rate))+
  labs(x='True volatility',y='Seconds per particle per step')+
  coord_cartesian(ylim=c(0,NA))
rate_by_eir <- ggplot(time_df)+
  geom_point(aes(x=true_init_EIR,y=run_rate))+
  scale_x_log10()+
  labs(x='True init_EIR',y='Seconds per particle per step')+
  coord_cartesian(ylim=c(0,NA))

fit_df <- bind_rows(lapply(1:length(sims_results_1to10_informed),function(i){
  ci_range_vol <- quantile(sims_results_1to10_informed[[i]]$mcmc$volatility,probs = c(0.025,0.5,0.975))
  ci_range_eir <- quantile(sims_results_1to10_informed[[i]]$mcmc$init_betaa,probs = c(0.025,0.5,0.975))
  df <- data.frame(volatility_median=ci_range_vol[[2]],
                   volatility_lower = ci_range_vol[[1]],
                   volatility_upper = ci_range_vol[[3]],
                   eir_median=ci_range_eir[[2]],
                   eir_lower = ci_range_eir[[1]],
                   eir_upper = ci_range_eir[[3]],
                   true_volatility = sim_dataset_100$variables[[i,'volatility']],
                   true_init_EIR = sim_dataset_100$variables[[i,'init_EIR']],
                   draw=factor(i))
  return(df)
}))

vol_fit <- ggplot(fit_df,aes(y=fct_reorder(draw,true_volatility)))+
  geom_point(aes(x=volatility_median))+
  geom_linerange(aes(xmin=volatility_lower,xmax=volatility_upper))+
  geom_point(aes(x=true_volatility),color='red')+
  labs(x='Volatility',y='Simulation')
eir_fit <- ggplot(fit_df,aes(y=fct_reorder(draw,true_init_EIR)))+
  geom_point(aes(x=exp(eir_median)))+
  geom_linerange(aes(xmin=exp(eir_lower),xmax=exp(eir_upper)))+
  geom_point(aes(x=true_init_EIR),color='red')+
  scale_x_log10()+
  labs(x='EIR',y='Simulation')

sim10_informed_longerflex_runs
sims_results_1to10_informed_longerflex <- lapply(1:length(sim10_informed_longerflex_runs$results()),function(i){
  combine_runs(sim10_informed_longerflex_runs$results()[[i]],n_steps=1000,n_chains=2,burnin=200)
})
trace_df_longerflex <- bind_rows(lapply(1:length(sims_results_1to10_informed_longerflex),function(i){
  df <- sims_results_1to10_informed_longerflex[[i]]$mcmc
  df$true_volatility <- sim_dataset_100$variables[[i,'volatility']]
  df$true_init_EIR <- sim_dataset_100$variables[[i,'init_EIR']]
  df$index <- c(1:nrow(df))
  df$draw <- factor(i)
  return(df)
}))
vol_trace_lf <- ggplot(trace_df_longerflex)+
  geom_line(aes(x=index,y=volatility))+
  facet_wrap(.~draw,scale='free_y')
eir_trace_vf <- ggplot(trace_df_longerflex)+
  geom_line(aes(x=index,y=init_betaa))+
  facet_wrap(.~fct_reorder(draw,true_init_EIR),scale='free_y')
