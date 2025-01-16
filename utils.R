cluster_setup <- function(context_name,template,cores){
  path <- paste0('T:/jth/',context_name)
  sources <- c("submit_sbc.R")
  ctx <- context::context_save(path,
                               sources = sources,
                               packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde','RecordLinkage','parallel','future'),
                               package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate','jt-hicks/mamasante','hyunjimoon/SBC')))
  config <- didehpc::didehpc_config(template = template,cores=cores, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
  obj_sbc <- didehpc::queue_didehpc(ctx,config = config)
  return(obj)
}
addCIs<-function(df,Ys,Ns){
  df$mean<-NA
  df$upper<-NA
  df$lower<-NA
  CIs<-binom::binom.confint(Ys[!is.na(Ys)],Ns[!is.na(Ns)],method="exact")
  df[which(!is.na(Ns)),]$mean<-CIs$mean
  df[which(!is.na(Ns)),]$upper<-CIs$upper
  df[which(!is.na(Ns)),]$lower<-CIs$lower
  return(df)
}
## fn to return prevalence from log_odds
get_prev_from_log_odds<-function(log_odds){
  return(exp(log_odds)/(1+exp(log_odds)))
}

## fn to return odds from prevalence
get_odds_from_prev<-function(prev){
  return(prev/(1-prev))
}

addCIs_anc<-function(df,Ys.cs,Ns.cs,Ys.anc,Ns.anc){
  df$mean.cs<-NA
  df$upper.cs<-NA
  df$lower.cs<-NA
  CIs.cs<-binom.confint(Ys.cs,Ns.cs,method="exact")
  df$mean.cs[Ns.cs>0]<-CIs.cs$mean[Ns.cs>0]
  df$upper.cs[Ns.cs>0]<-CIs.cs$upper[Ns.cs>0]
  df$lower.cs[Ns.cs>0]<-CIs.cs$lower[Ns.cs>0]
  df$mean.anc<-NA
  df$upper.anc<-NA
  df$lower.anc<-NA
  CIs.anc<-binom.confint(Ys.anc,Ns.anc,method="exact")
  df$mean.anc[Ns.anc>0]<-CIs.anc$mean[Ns.anc>0]
  df$upper.anc[Ns.anc>0]<-CIs.anc$upper[Ns.anc>0]
  df$lower.anc[Ns.anc>0]<-CIs.anc$lower[Ns.anc>0]
  return(df)
}
################## generate the data ######################
# generate random walk of betaa (recursive fn)
genRandWalk <- function(x,vol,randWalk) {
  if (x == 0)    return (randWalk)
  else return(genRandWalk(x-1,vol,c(randWalk,min(randWalk[length(randWalk)]+rnorm(1)*vol,log(max_param)))))
}
get_t <- function(date,origin=as.Date('2015-01-01')){
  return(date-origin)
}
datasim4pmcmc <- function(N = 24,
                          max_param=125,
                          model_file="MiP_odin_model_nodelay.R"){
  volatility <- rgamma(1, shape = 3.4, rate = 3.1)
  log_init_EIR <- rnorm(1, mean = 4, sd = 3)
  init_EIR <- exp(log_init_EIR)

  cat('volatility = ',volatility,' init_EIR = ',init_EIR,'\n')
  init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

  prop_treated <- 0.4
  het_brackets <- 5

  ##Convert init_EIR into an init_betaa
  mpl_initial <- mamasante::model_param_list_create(init_EIR = init_EIR,
                                                 init_ft = prop_treated
  )

  pars_initial <- mamasante::equilibrium_init_create_stripped(age_vector = init_age,
                                                           init_EIR = init_EIR,
                                                           ft = prop_treated,
                                                           model_param_list = mpl_initial,
                                                           het_brackets = het_brackets)
  init_betaa <- pars_initial$betaa_eq

  betaa_dates <- seq.Date(from = as.Date('2015-01-01'),by = 'month',length.out = N)
  betaa_times <- sapply(betaa_dates,get_t)

  ### just a random walk on logscale
  log_betaa_vals <- genRandWalk(length(betaa_times)-1,volatility,log(init_betaa))

  betaa_vals <- exp(log_betaa_vals)
  ##set up the simulation for the simualted data
  time<- max(betaa_times) + 30
  out_step=1

  mpl <- mamasante::model_param_list_create(init_EIR = init_EIR,
                                         init_ft = prop_treated,
                                         betaa_times=betaa_times,
                                         betaa_vals=betaa_vals,
                                         lag_rates = 10
  )

  pars <- mamasante::equilibrium_init_create_stripped(age_vector = init_age,
                                                   init_EIR = init_EIR,
                                                   ft = prop_treated,
                                                   model_param_list = mpl,
                                                   het_brackets = het_brackets)

  ##The malaria model
  generator <- odin(model_file)
  state_use <- pars[names(pars) %in% coef(generator)$name]

  # create model with initial values
  mod <- generator(user = state_use, use_dde = TRUE)
  tt <- seq(0, time, out_step)

  # run the simulation to base the data
  mod_run <- mod$run(tt, step_max_n = 1e7,
                     atol = 1e-5,
                     rtol = 1e-5)

  # shape output
  out <- mod$transform_variables(mod_run)

  # plot data and generate data
  # plot(out$t,out$betaa_out,type='l',col="red",ylim=c(0,125))
  # lines(out$t,out$prev*100,col="blue",lwd=4)
  out_df <- data.frame(t=out$t,
                       prev05=out$prev05,
                       date=as.character(seq.Date(from = as.Date('2015-01-01'),by = 'day',length.out = length(tt))),
                       EIR_true=out$EIR_out,
                       betaa_true=out$betaa_out,
                       inc05_true=out$inc05,
                       inc_all_true = out$inc,
                       prev_all_true = out$prev_all,
                       mv = out$mv_out,
                       spz_rate = out$spz_rate)
  months <- unique(as.yearmon(out_df$date))
  midmonth_dates <- data.frame(date=as.character(as.Date(months,frac=0.5)))
  monthly_data <- left_join(midmonth_dates,out_df,by='date')

  monthly_data$tested<-round(rnorm(nrow(monthly_data),500,50))
  monthly_data$positive<-rbinom(nrow(monthly_data),monthly_data$tested,monthly_data$prev05)
  monthly_data$month <- zoo::as.yearmon(monthly_data$date)
  if(any(monthly_data$positive>monthly_data$tested)) {stop('Number of positive is greater than number tested')}

  return(list(variables = list(volatility=volatility,
              init_EIR = init_EIR),
              generated=monthly_data))
}
combine_runs <- function(results,n_chains,n_steps,burnin){
  first_chain <- c((burnin+1):n_steps)
  steps2keep <- first_chain
  if(n_chains > 1){
    for(i in 2:n_chains){
      more_steps <- first_chain + (i-1)*n_steps
      steps2keep <- append(steps2keep,more_steps)
    }
  }
  mcmc <- results$mcmc[steps2keep,]
  threads <- results$threads
  run_time <- results$run_time
  history <- results$history[, steps2keep, -1]
  return(list(n_chains=n_chains,
              n_steps=n_steps,
              burnin=burnin,
              mcmc=mcmc,
              history=history,
              run_time=run_time,
              threads=threads))
}
sim_history4plots <- function(results,sim_data,burnin){
  length <- length(results$history[1,,1])
  start_chain <- round(burnin*length)+1
  history <- results$history[,start_chain:length, -c(1,2)]
  # history <- results$history[,start_chain:length, -1]
  sim_data$date <- as.Date(sim_data$date)
  dates <- as.Date(sim_data$date)
  if(length>100){
    index_sample <- sample(length(history[1,,1]), 100)
  }else{
    index_sample <- c(1:length)
  }

  history.df.prev <- as.data.frame(t(history['prev_05',,]))
  prev_history <- history.df.prev%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(prev.median=median(value),
                     prev.mean=mean(value),
                     prev.upper=quantile(value,probs=0.975),
                     prev.lower=quantile(value,probs=0.025))%>%
    mutate(date=as.Date(date))
  prev_sample <- history.df.prev[,index_sample] %>%
    mutate(t=c(1:nrow(prev_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100))

  prev_plot <- ggplot()+
    geom_point(data=sim_data,aes(x=date,y=prev05_true))+
    geom_line(data = prev_history, aes(x=date,y=prev.median),color='darkgrey',linewidth=1)+
    geom_ribbon(data=prev_history,aes(x=date,ymin=prev.lower,ymax=prev.upper),fill='lightgrey',alpha=0.5)+
    coord_cartesian(ylim = c(0,NA))
  prev_plot_sample <- ggplot()+
    geom_line(data=prev_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
    geom_line(data = prev_history, aes(x=date,y=prev.median),color='darkgrey',linewidth=1)+
    geom_point(data=sim_data,aes(x=date,y=prev05))+
    coord_cartesian(ylim = c(0,NA))

  history.df.incall <- as.data.frame(t(history['clininc_all',,]))
  inc_history <- history.df.incall%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(incall.median=median(value),
                     incall.mean=mean(value),
                     incall.upper=quantile(value,probs=0.975),
                     incall.lower=quantile(value,probs=0.025))
  inc_sample <- history.df.incall[,index_sample] %>%
    mutate(t=c(1:nrow(inc_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100))
  inc_plot <- ggplot()+
    geom_point(data=sim_data,aes(x=date,y=inc_all_true))+
    geom_line(data = inc_history, aes(x=date,y=incall.median),color='darkgrey')+
    geom_ribbon(data=inc_history,aes(x=date,ymin=incall.lower,ymax=incall.upper),fill='lightgrey',alpha=0.5)+
    coord_cartesian(ylim = c(0,NA))
  inc_plot_sample <- ggplot()+
    geom_line(data=inc_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
    geom_line(data = inc_history, aes(x=date,y=incall.median),color='darkgrey',linewidth=1)+
    geom_point(data=sim_data,aes(x=date,y=inc_all_true))+
    coord_cartesian(ylim = c(0,NA))

  history.df.EIR <- as.data.frame(t(history['EIR',,]))
  EIR_history <- history.df.EIR%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(EIR.median=median(value),
                     EIR.mean=mean(value),
                     EIR.upper=quantile(value,probs=0.975),
                     EIR.lower=quantile(value,probs=0.025))
  EIR_sample <- history.df.EIR[,index_sample] %>%
    mutate(t=c(1:nrow(EIR_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100))
  EIR_plot <- ggplot()+
    geom_point(data=sim_data,aes(x=date,y=EIR_true))+
    geom_line(data = EIR_history, aes(x=date,y=EIR.median),color='darkgrey')+
    geom_ribbon(data=EIR_history,aes(x=date,ymin=EIR.lower,ymax=EIR.upper),fill='lightgrey',alpha=0.5)+
    coord_cartesian(ylim = c(0,NA))
  EIR_plot_sample <- ggplot()+
    geom_line(data=EIR_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
    geom_line(data = EIR_history, aes(x=date,y=EIR.median),color='darkgrey',linewidth=1)+
    geom_point(data=sim_data,aes(x=date,y=EIR_true))+
    coord_cartesian(ylim = c(0,NA))
#
#   history.df.mv <- as.data.frame(t(history['eff_moz_pop',,]))
#   mv_history <- history.df.mv%>%
#     dplyr::mutate(date=dates)%>%
#     melt(id='date')%>%
#     group_by(date)%>%
#     dplyr::summarise(mv.median=median(value),
#                      mv.mean=mean(value),
#                      mv.upper=quantile(value,probs=0.975),
#                      mv.lower=quantile(value,probs=0.025))
#
#   mv_sample <- history.df.mv[,index_sample] %>%
#     mutate(t=c(1:nrow(mv_history)))%>%
#     melt(id='t')%>%
#     dplyr::rename(time=t)%>%
#     mutate(date = rep(dates,100))
#   mv_plot <- ggplot()+
#     geom_point(data=sim_data,aes(x=date,y=mv))+
#     geom_line(data = mv_history, aes(x=date,y=mv.median),color='darkgrey')+
#     geom_ribbon(data=mv_history,aes(x=date,ymin=mv.lower,ymax=mv.upper),fill='lightgrey',alpha=0.5)+
#     coord_cartesian(ylim = c(0,NA))
#   mv_plot_sample <- ggplot()+
#     geom_line(data=mv_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
#     geom_line(data = mv_history, aes(x=date,y=mv.median),color='darkgrey',linewidth=1)+
#     geom_point(data=sim_data,aes(x=date,y=mv))+
#     coord_cartesian(ylim = c(0,NA))
#
#   history.df.spzrate <- as.data.frame(t(history['spz_rate',,]))
#   spz_history <- history.df.spzrate%>%
#     dplyr::mutate(date=dates)%>%
#     melt(id='date')%>%
#     group_by(date)%>%
#     dplyr::summarise(spzrate.median=median(value),
#                      spzrate.mean=mean(value),
#                      spzrate.upper=quantile(value,probs=0.975),
#                      spzrate.lower=quantile(value,probs=0.025))
#   spz_sample <- history.df.spzrate[,index_sample] %>%
#     mutate(t=c(1:nrow(mv_history)))%>%
#     melt(id='t')%>%
#     dplyr::rename(time=t)%>%
#     mutate(date = rep(dates,100))
#   spz_plot <- ggplot()+
#     geom_point(data=sim_data,aes(x=date,y=spz_rate))+
#     geom_line(data = spz_history, aes(x=date,y=spzrate.median),color='darkgrey')+
#     geom_ribbon(data=spz_history,aes(x=date,ymin=spzrate.lower,ymax=spzrate.upper),fill='lightgrey',alpha=0.5)+
#     coord_cartesian(ylim = c(0,NA))
#   spz_plot_sample <- ggplot()+
#     geom_line(data=spz_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
#     geom_line(data = spz_history, aes(x=date,y=spzrate.median),color='darkgrey',linewidth=1)+
#     geom_point(data=sim_data,aes(x=date,y=spz_rate))+
#     coord_cartesian(ylim = c(0,NA))

  history.df.betaa <- as.data.frame(t(history['betaa',,]))
  betaa_history <- history.df.betaa%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(betaa.median=median(value),
                     betaa.mean=mean(value),
                     betaa.upper=quantile(value,probs=0.975),
                     betaa.lower=quantile(value,probs=0.025))
  betaa_sample <- history.df.betaa[,index_sample] %>%
    mutate(t=c(1:nrow(betaa_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100))
  betaa_plot <- ggplot()+
    geom_point(data=sim_data,aes(x=date,y=betaa_true))+
    geom_line(data = betaa_history, aes(x=date,y=betaa.median),color='darkgrey')+
    geom_ribbon(data=betaa_history,aes(x=date,ymin=betaa.lower,ymax=betaa.upper),fill='lightgrey',alpha=0.5)+
    coord_cartesian(ylim = c(0,NA))
  betaa_plot_sample <- ggplot()+
    geom_line(data=betaa_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
    geom_line(data = betaa_history, aes(x=date,y=betaa.median),color='darkgrey',linewidth=1)+
    geom_point(data=sim_data,aes(x=date,y=betaa_true))+
    coord_cartesian(ylim = c(0,NA))

  return(list(prev_history=prev_history,
              prev_sample=prev_sample,
              inc_history=inc_history,
              inc_sample=inc_sample,
              EIR_history=EIR_history,
              EIR_sample=EIR_sample,
              # mv_history=mv_history,
              # mv_sample=mv_sample,
              # spz_history=spz_history,
              # spz_sample=spz_sample,
              betaa_history=betaa_history,
              betaa_sample=betaa_sample))
}
gen_seasonal_sim <- function(init_EIR=100,
                             max_param=125,
                             model_file= "init/odin_model_stripped_seasonal.R",
                             country = 'Burkina Faso',
                             admin_unit = 'Cascades',
                             sim_length = 2){
  #Provide age categories, proportion treated, and number of heterogeneity brackets
  #Create model parameter list. Also loads seasonality profile data file to match to desired admin_unit and country
  cat('Country = ',country,' Admin = ',admin_unit,' init_EIR = ',init_EIR,'\n')
  init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

  prop_treated <- 0.4
  het_brackets <- 5

  ##set up the simulation for the simualted data
  time<- 10*365
  out_step=1

  mpl <- mamasante::model_param_list_create(init_EIR = init_EIR,
                                         init_ft = prop_treated,
                                         country=country,
                                         admin_unit = admin_unit,
                                         lag_rates = 10,
                                         max_param = max_param,
                                         volatility = 1,
                                         state_check = 0,
                                         comparison = 'u5'
  )

  pars <- mamasante::equilibrium_init_create_stripped(age_vector = init_age,
                                                   init_EIR = init_EIR,
                                                   ft = prop_treated,
                                                   model_param_list = mpl,
                                                   het_brackets = het_brackets)

  ##The malaria model
  generator <- odin(model_file)
  state_use <- pars[names(pars) %in% coef(generator)$name]

  # create model with initial values
  mod <- generator(user = state_use, use_dde = TRUE)
  tt <- seq(0, time, out_step)

  # run the simulation to base the data
  mod_run <- mod$run(tt, step_max_n = 1e7,
                     atol = 1e-5,
                     rtol = 1e-5)

  # shape output
  out <- mod$transform_variables(mod_run)

  # plot data and generate data
  # plot(out$t,out$betaa_out,type='l',col="red",ylim=c(0,125))
  # lines(out$t,out$prev*100,col="blue",lwd=4)
  out_df <- data.frame(t=out$t,
                       prev=out$prev,
                       date=as.character(seq.Date(from = as.Date('2015-01-01'),by = 'day',length.out = length(tt))))
  months <- unique(as.yearmon(out_df$date))
  midmonth_dates <- data.frame(date=as.character(as.Date(months,frac=0.5)))
  monthly_data <- left_join(midmonth_dates,out_df,by='date')
  ##Use only last 2 years for simulated data
  sim_months <- sim_length*12-1
  monthly_data <- monthly_data[(nrow(monthly_data)-sim_months):nrow(monthly_data),]

  monthly_data$tested<-round(rnorm(nrow(monthly_data),500,50))
  monthly_data$positive<-rbinom(nrow(monthly_data),monthly_data$tested,monthly_data$prev)
  monthly_data$month <- zoo::as.yearmon(monthly_data$date)
  if(any(monthly_data$positive>monthly_data$tested)) {stop('Number of positive is greater than number tested')}

  true_data <- data.frame(t=out$t,
                          date=seq.Date(from = as.Date('2015-01-01'),by = 'day',length.out = length(tt)),
                          prev05_true=out$prev,
                          EIR_true=out$EIR_out,
                          betaa_true=out$betaa_out,
                          inc05_true=out$inc05,
                          inc_all_true = out$inc,
                          prev_all_true = out$prev_all)

  return(list(init_EIR = init_EIR,
              true_data=true_data,
              data_raw=monthly_data))
}
create_diag_figs <- function(result,folderpath,name,country,district,burnin=0){
  length <- length(result$mcmc[,1])
  start_chain <- round(burnin*length)+1
  mcmc <- result$mcmc[start_chain:length,]
  print('acceptance rate')
  ar <- 1 - coda::rejectionRate(as.mcmc(mcmc))
  print(ar)
  print('effective size')
  ess <- coda::effectiveSize(as.mcmc(mcmc))
  print(ess)

  title <- paste0('Diagnostic plots for seasonal model - ',district,', ',country)
  pars_list <- names(mcmc)
  trace_plots <- lapply(pars_list, function(x){
    bayesplot::mcmc_trace(mcmc,pars = x) +
      ggtitle(paste0(x,' / AR: ',round(ar[x],3),' / ESS: ',round(ess[x],1)))+
      theme(title = element_text(size=6),
            axis.title.y = element_blank())
  })
  dense_plots <- lapply(pars_list, function(x){
    bayesplot::mcmc_dens(mcmc,pars = x) +
      ggtitle(paste0(x,' / AR: ',round(ar[x],2),' / ESS: ',round(ess[x],1)))+
      theme(title = element_text(size=6),
            axis.title.x = element_blank())
  })
  diag <- (trace_plots[[1]]+dense_plots[[1]])/
    (trace_plots[[2]]+dense_plots[[2]])/
    (trace_plots[[3]]+dense_plots[[3]])/
    (trace_plots[[4]]+dense_plots[[4]]) +
    plot_layout(guides = "collect") + plot_annotation(title = title)

  ggsave(filename=paste0(folderpath,'/',name,'-',district,'-',country,'.tiff'),plot=diag,dpi=300,height = 11,width=8,units = 'in')
  return(diag)
}
run_history4plots <- function(results,data,burnin){
  length <- length(results$history[1,,1])
  start_chain <- round(burnin*length)+1
  history <- results$history[,start_chain:length, -1]
  dates <- data$date
  data$date <- as.Date(data$date)
  data$prev05 <- data$positive/data$tested
  if(length>100){
    index_sample <- sample(length(history[1,,1]), 100)
  }else{
    index_sample <- c(1:length)
  }

  history.df.prev <- as.data.frame(t(history['prev_05',,]))
  prev_history <- history.df.prev%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(prev.median=median(value),
                     prev.mean=mean(value),
                     prev.upper=quantile(value,probs=0.975),
                     prev.lower=quantile(value,probs=0.025))%>%
    mutate(date=as.Date(date))
  prev_sample <- history.df.prev[,index_sample] %>%
    mutate(t=c(1:nrow(prev_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100))

  prev_plot <- ggplot()+
    geom_point(data=data,aes(x=date,y=prev05))+
    geom_line(data = prev_history, aes(x=date,y=prev.median),color='darkgrey',linewidth=1)+
    geom_ribbon(data=prev_history,aes(x=date,ymin=prev.lower,ymax=prev.upper),fill='lightgrey',alpha=0.5)+
    coord_cartesian(ylim = c(0,NA))
  prev_plot_sample <- ggplot()+
    geom_line(data=prev_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
    geom_line(data = prev_history, aes(x=date,y=prev.median),color='darkgrey',linewidth=1)+
    geom_point(data=data,aes(x=date,y=prev05))+
    coord_cartesian(ylim = c(0,NA))

  history.df.incall <- as.data.frame(t(history['clininc_all',,]))
  inc_history <- history.df.incall%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(incall.median=median(value),
                     incall.mean=mean(value),
                     incall.upper=quantile(value,probs=0.975),
                     incall.lower=quantile(value,probs=0.025))
  inc_sample <- history.df.incall[,index_sample] %>%
    mutate(t=c(1:nrow(inc_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100))
  inc_plot <- ggplot()+
    geom_point(data=data,aes(x=date,y=inc_all_true))+
    geom_line(data = inc_history, aes(x=date,y=incall.median),color='darkgrey')+
    geom_ribbon(data=inc_history,aes(x=date,ymin=incall.lower,ymax=incall.upper),fill='lightgrey',alpha=0.5)+
    coord_cartesian(ylim = c(0,NA))
  inc_plot_sample <- ggplot()+
    geom_line(data=inc_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
    geom_line(data = inc_history, aes(x=date,y=incall.median),color='darkgrey',linewidth=1)+
    geom_point(data=data,aes(x=date,y=inc_all_true))+
    coord_cartesian(ylim = c(0,NA))

  history.df.EIR <- as.data.frame(t(history['EIR',,]))
  EIR_history <- history.df.EIR%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(EIR.median=median(value),
                     EIR.mean=mean(value),
                     EIR.upper=quantile(value,probs=0.975),
                     EIR.lower=quantile(value,probs=0.025))
  EIR_sample <- history.df.EIR[,index_sample] %>%
    mutate(t=c(1:nrow(mv_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100))
  EIR_plot <- ggplot()+
    geom_point(data=data,aes(x=date,y=EIR_true))+
    geom_line(data = EIR_history, aes(x=date,y=EIR.median),color='darkgrey')+
    geom_ribbon(data=EIR_history,aes(x=date,ymin=EIR.lower,ymax=EIR.upper),fill='lightgrey',alpha=0.5)+
    coord_cartesian(ylim = c(0,NA))
  EIR_plot_sample <- ggplot()+
    geom_line(data=EIR_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
    geom_line(data = EIR_history, aes(x=date,y=EIR.median),color='darkgrey',linewidth=1)+
    geom_point(data=data,aes(x=date,y=EIR_true))+
    coord_cartesian(ylim = c(0,NA))

  history.df.mv <- as.data.frame(t(history['eff_moz_pop',,]))
  mv_history <- history.df.mv%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(mv.median=median(value),
                     mv.mean=mean(value),
                     mv.upper=quantile(value,probs=0.975),
                     mv.lower=quantile(value,probs=0.025))

  mv_sample <- history.df.mv[,index_sample] %>%
    mutate(t=c(1:nrow(mv_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100))
  mv_plot <- ggplot()+
    geom_point(data=data,aes(x=date,y=mv))+
    geom_line(data = mv_history, aes(x=date,y=mv.median),color='darkgrey')+
    geom_ribbon(data=mv_history,aes(x=date,ymin=mv.lower,ymax=mv.upper),fill='lightgrey',alpha=0.5)+
    coord_cartesian(ylim = c(0,NA))
  mv_plot_sample <- ggplot()+
    geom_line(data=mv_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
    geom_line(data = mv_history, aes(x=date,y=mv.median),color='darkgrey',linewidth=1)+
    geom_point(data=data,aes(x=date,y=mv))+
    coord_cartesian(ylim = c(0,NA))

  history.df.spzrate <- as.data.frame(t(history['spz_rate',,]))
  spz_history <- history.df.spzrate%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(spzrate.median=median(value),
                     spzrate.mean=mean(value),
                     spzrate.upper=quantile(value,probs=0.975),
                     spzrate.lower=quantile(value,probs=0.025))
  spz_sample <- history.df.spzrate[,index_sample] %>%
    mutate(t=c(1:nrow(mv_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100))
  spz_plot <- ggplot()+
    geom_point(data=data,aes(x=date,y=spz_rate))+
    geom_line(data = spz_history, aes(x=date,y=spzrate.median),color='darkgrey')+
    geom_ribbon(data=spz_history,aes(x=date,ymin=spzrate.lower,ymax=spzrate.upper),fill='lightgrey',alpha=0.5)+
    coord_cartesian(ylim = c(0,NA))
  spz_plot_sample <- ggplot()+
    geom_line(data=spz_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
    geom_line(data = spz_history, aes(x=date,y=spzrate.median),color='darkgrey',linewidth=1)+
    geom_point(data=data,aes(x=date,y=spz_rate))+
    coord_cartesian(ylim = c(0,NA))

  history.df.betaa <- as.data.frame(t(history['betaa',,]))
  betaa_history <- history.df.betaa%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(betaa.median=median(value),
                     betaa.mean=mean(value),
                     betaa.upper=quantile(value,probs=0.975),
                     betaa.lower=quantile(value,probs=0.025))
  betaa_sample <- history.df.betaa[,index_sample] %>%
    mutate(t=c(1:nrow(mv_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100))
  betaa_plot <- ggplot()+
    geom_point(data=data,aes(x=date,y=betaa_true))+
    geom_line(data = betaa_history, aes(x=date,y=betaa.median),color='darkgrey')+
    geom_ribbon(data=betaa_history,aes(x=date,ymin=betaa.lower,ymax=betaa.upper),fill='lightgrey',alpha=0.5)+
    coord_cartesian(ylim = c(0,NA))
  betaa_plot_sample <- ggplot()+
    geom_line(data=betaa_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
    geom_line(data = betaa_history, aes(x=date,y=betaa.median),color='darkgrey',linewidth=1)+
    geom_point(data=data,aes(x=date,y=betaa_true))+
    coord_cartesian(ylim = c(0,NA))

  return(list(prev_history=prev_history,
              prev_sample=prev_sample,
              inc_history=inc_history,
              inc_sample=inc_sample,
              EIR_history=EIR_history,
              EIR_sample=EIR_sample,
              mv_history=mv_history,
              mv_sample=mv_sample,
              spz_history=spz_history,
              spz_sample=spz_sample,
              betaa_history=betaa_history,
              betaa_sample=betaa_sample))
}
daily2monthly <- function(out_df,sim_length=2){

  out_df$date <- as.character(out_df$date)
  months <- unique(as.yearmon(as.Date(out_df$date)))
  midmonth_dates <- data.frame(date=as.character(as.Date(months,frac=0.5)))
  monthly_data <- left_join(midmonth_dates,out_df,by='date')
  monthly_data$month <- zoo::as.yearmon(monthly_data$date)
  sim_months <- sim_length*12-1
  monthly_data <- monthly_data[(nrow(monthly_data)-sim_months):nrow(monthly_data),]

  return(monthly_data)
}
get_chain_info <- function(n_chains, length, burnin){
  chain_length <- length/n_chains
  start_chain <- round(burnin*chain_length)+1
  ind <- c()
  for(i in 1:n_chains){
    shift <- chain_length*(i-1)
    ind <- append(ind,c((start_chain+shift):(chain_length+shift)))
  }
  return(ind)
}
prep_results <- function(results,sim_data,burnin,n_chains=1,country,district,timelength=NA,anc=FALSE){
  ind <- get_chain_info(n_chains=n_chains,
                        length=length(results$history[1,,1]),
                        burnin=burnin)
  num_months <- length(results$history[1,1,])
  # history <- results$history[,start_chain:length, -c(1,2)]
  if(!('date'%in%names(sim_data))){
    sim_data$date <- as.Date(sim_data$month,frac=0.5)
  }
  sim_data$date <- as.Date(sim_data$date)
  dates <- as.Date(sim_data$date)
  if(is.na(timelength)){
    timelength=length(dates)
  }
  history <- results$history[,ind, (num_months-timelength+1):num_months]
  times <- results$times[(num_months-timelength+1):num_months]

  if(length(ind)>100){
    index_sample <- sample(length(history[1,,1]), 100)
  }else{
    index_sample <- c(1:length)
  }

  history.df.prev <- as.data.frame(t(history['prev_05',,]))
  prev_history <- history.df.prev%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(median=median(value),
                     mean=mean(value),
                     upper=quantile(value,probs=0.975),
                     lower=quantile(value,probs=0.025))%>%
    mutate(date=as.Date(date),
           measure = 'prev_05')
  prev_means <- history.df.prev%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(variable)%>%
    summarise(mean = mean(value))
  prev_summary <- prev_means %>%
    summarise(median = median(mean))
  prev_hpd <- HPDinterval(as.mcmc(prev_means$mean))
  prev_summary$lower_hpd <- prev_hpd[1]
  prev_summary$upper_hpd <- prev_hpd[2]
  prev_summary$measure <- 'prev_05'
  prev_sample <- history.df.prev[,index_sample] %>%
    mutate(t=c(1:nrow(prev_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100),
           measure = 'prev_05')

  history.df.inc05 <- as.data.frame(t(history['clininc_05',,]))
  inc_history <- history.df.inc05%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(median=median(value),
                     mean=mean(value),
                     upper=quantile(value,probs=0.975),
                     lower=quantile(value,probs=0.025))%>%
    mutate(date=as.Date(date),
           measure = 'inc05')
  inc_sample <- history.df.inc05[,index_sample] %>%
    mutate(t=c(1:nrow(inc_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100),
           measure = 'inc05')
  inc_means <- history.df.inc05%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(variable)%>%
    summarise(mean = mean(value))
  inc_summary <- inc_means %>%
    summarise(median = median(mean))
  inc_hpd <- HPDinterval(as.mcmc(inc_means$mean))
  inc_summary$lower_hpd <- inc_hpd[1]
  inc_summary$upper_hpd <- inc_hpd[2]
  inc_summary$measure <- 'inc05'


  history.df.EIR <- as.data.frame(t(history['EIR',,]))
  EIR_history <- history.df.EIR%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(median=median(value),
                     mean=mean(value),
                     upper=quantile(value,probs=0.975),
                     lower=quantile(value,probs=0.025))%>%
    mutate(date=as.Date(date),
           measure = 'EIR')
  EIR_sample <- history.df.EIR[,index_sample] %>%
    mutate(t=c(1:nrow(EIR_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100),
           measure = 'EIR')
  EIR_means <- history.df.EIR%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(variable)%>%
    summarise(mean = mean(value))
  EIR_summary <- EIR_means %>%
    summarise(median = median(mean))
  EIR_hpd <- HPDinterval(as.mcmc(EIR_means$mean))
  EIR_summary$lower_hpd <- EIR_hpd[1]
  EIR_summary$upper_hpd <- EIR_hpd[2]
  EIR_summary$measure <- 'EIR'

  history.df.betaa <- as.data.frame(t(history['betaa',,]))
  betaa_history <- history.df.betaa%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(median=median(value),
                     mean=mean(value),
                     upper=quantile(value,probs=0.975),
                     lower=quantile(value,probs=0.025))%>%
    mutate(date=as.Date(date),
           measure = 'betaa')
  betaa_means <- history.df.betaa%>%
    dplyr::mutate(date=dates)%>%
    melt(id='date')%>%
    group_by(variable)%>%
    summarise(mean = mean(value))
  betaa_summary <- betaa_means %>%
    summarise(median = median(mean))
  betaa_hpd <- HPDinterval(as.mcmc(betaa_means$mean))
  betaa_summary$lower_hpd <- betaa_hpd[1]
  betaa_summary$upper_hpd <- betaa_hpd[2]
  betaa_summary$measure <- 'betaa'
  betaa_sample <- history.df.betaa[,index_sample] %>%
    mutate(t=c(1:nrow(betaa_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100),
           measure = 'betaa')
  history_summary <- bind_rows(prev_history,inc_history,EIR_history,betaa_history)
  sample <- bind_rows(prev_sample,inc_sample,EIR_sample,betaa_sample)
  point_est_summary <- bind_rows(prev_summary,inc_summary,EIR_summary,betaa_summary)

  if(anc){
    history.df.prev.pg <- as.data.frame(t(history['prev_pg',,]))
    prev_history_pg <- history.df.prev.pg%>%
      dplyr::mutate(date=dates)%>%
      melt(id='date')%>%
      group_by(date)%>%
      dplyr::summarise(median=median(value),
                       mean=mean(value),
                       upper=quantile(value,probs=0.975),
                       lower=quantile(value,probs=0.025))%>%
      mutate(date=as.Date(date),
             measure = 'prev_pg')
    prev_pg_means <- history.df.prev.pg%>%
      dplyr::mutate(date=dates)%>%
      melt(id='date')%>%
      group_by(variable)%>%
      summarise(mean = mean(value))
    prev_pg_summary <- prev_pg_means %>%
      summarise(median = median(mean))
    prev_pg_hpd <- HPDinterval(as.mcmc(prev_pg_means$mean))
    prev_pg_summary$lower_hpd <- prev_pg_hpd[1]
    prev_pg_summary$upper_hpd <- prev_pg_hpd[2]
    prev_pg_summary$measure <- 'prev_pg'
    prev_sample_pg <- history.df.prev.pg[,index_sample] %>%
      mutate(t=c(1:nrow(prev_history_pg)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      mutate(date = rep(dates,100),
             measure = 'prev_pg')
    sgmg <- ifelse('prev_sg' %in% names(history[,1,1]),'prev_sg',ifelse('prev_mg' %in% names(history[,1,1]),'prev_mg',NA))
    history.df.prev.mg <- as.data.frame(t(history[sgmg,,]))
    prev_history_mg <- history.df.prev.mg%>%
      dplyr::mutate(date=dates)%>%
      melt(id='date')%>%
      group_by(date)%>%
      dplyr::summarise(median=median(value),
                       mean=mean(value),
                       upper=quantile(value,probs=0.975),
                       lower=quantile(value,probs=0.025))%>%
      mutate(date=as.Date(date),
             measure = 'prev_mg')
    prev_mg_means <- history.df.prev.mg%>%
      dplyr::mutate(date=dates)%>%
      melt(id='date')%>%
      group_by(variable)%>%
      summarise(mean = mean(value))
    prev_mg_summary <- prev_mg_means %>%
      summarise(median = median(mean))
    prev_mg_hpd <- HPDinterval(as.mcmc(prev_mg_means$mean))
    prev_mg_summary$lower_hpd <- prev_mg_hpd[1]
    prev_mg_summary$upper_hpd <- prev_mg_hpd[2]
    prev_mg_summary$measure <- 'prev_mg'
    prev_sample_mg <- history.df.prev.mg[,index_sample] %>%
      mutate(t=c(1:nrow(prev_history_pg)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      mutate(date = rep(dates,100),
             measure = 'prev_mg')

    history.df.prev.all <- as.data.frame(t(history['prev_anc_all',,]))
    prev_history_all <- history.df.prev.all%>%
      dplyr::mutate(date=dates)%>%
      melt(id='date')%>%
      group_by(date)%>%
      dplyr::summarise(median=median(value),
                       mean=mean(value),
                       upper=quantile(value,probs=0.975),
                       lower=quantile(value,probs=0.025))%>%
      mutate(date=as.Date(date),
             measure = 'prev_anc_all')
    prev_all_means <- history.df.prev.all%>%
      dplyr::mutate(date=dates)%>%
      melt(id='date')%>%
      group_by(variable)%>%
      summarise(mean = mean(value))
    prev_all_summary <- prev_all_means %>%
      summarise(median = median(mean))
    prev_all_hpd <- HPDinterval(as.mcmc(prev_all_means$mean))
    prev_all_summary$lower_hpd <- prev_all_hpd[1]
    prev_all_summary$upper_hpd <- prev_all_hpd[2]
    prev_all_summary$measure <- 'prev_anc_all'
    prev_sample_all <- history.df.prev.all[,index_sample] %>%
      mutate(t=c(1:nrow(prev_history_pg)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      mutate(date = rep(dates,100),
             measure = 'prev_anc_all')
    history_summary <- bind_rows(history_summary,prev_history_pg,prev_history_mg,prev_history_all)
    sample <- bind_rows(sample,prev_sample_pg,prev_sample_mg,prev_sample_all)
    point_est_summary <- bind_rows(point_est_summary,prev_pg_summary,prev_mg_summary,prev_all_summary)

  }
  history_summary$country <- country
  history_summary$district <- district
  sample$country <- country
  sample$district <- district
  point_est_summary$country <- country
  point_est_summary$district <- district

  return(list(summary=history_summary,
              sample = sample,
              times = times,
              point_est = point_est_summary))
}
get_u5prev_fromanc <- function(avg_prev,comparison=NULL){
  ## Gravidity prevalence conversion coefficients:
  coefs_pgsgmg_df <- apply(mamasante::load_file('pgsgmg_corr_sample.RDS'),2,median)
  coefs_all_df <- apply(mamasante::load_file('all_corr_sample.RDS'),2,median)

  av_lo_child <- coefs_pgsgmg_df[['av_lo_child']]
  gradient_pg <- coefs_pgsgmg_df[['gradient_pg']]
  intercept_pg <- coefs_pgsgmg_df[['intercept_pg']]
  gradient_sg <- coefs_pgsgmg_df[['gradient_sg']]
  intercept_sg <- coefs_pgsgmg_df[['intercept_sg']]
  gradient_mg <- coefs_pgsgmg_df[['gradient_mg']]
  intercept_mg <- coefs_pgsgmg_df[['intercept_mg']]
  av_lo_child_all <- coefs_all_df[['av_lo_child']]
  gradient_all <- coefs_all_df[['gradient']]
  intercept_all <- coefs_all_df[['intercept']]

  #Determine average log-odds of childhood prevalence depending on first year of available data
  if(comparison=='ancall'){
    log_odds_pall <- log(mamasante::get_odds_from_prev(avg_prev))
    log_odds_child <- ((log_odds_pall - intercept_all) + av_lo_child_all*gradient_all)/(gradient_all + 1)
    prev_u5 <- mamasante::get_prev_from_log_odds(log_odds_child)
  }else {
    log_odds_pg <- log(mamasante::get_odds_from_prev(avg_prev[1]))
    log_odds_child <- ((log_odds_pg - intercept_pg) + av_lo_child*gradient_pg)/(gradient_pg + 1)
    prev_u5 <- mamasante::get_prev_from_log_odds(log_odds_child)
  }

  return(prev_u5)
}
get_param_sum <- function(results,param_name,posterior_name,n_chains,burnin,country){
  ind <- get_chain_info(n_chains=n_chains,length=nrow(results$mcmc),burnin=burnin)
  dist <- results$mcmc[ind,param_name]
  max_post <- which.max(results$mcmc[ind,posterior_name])
  MAP <- dist[max_post]
  hpd <- HPDinterval(as.mcmc(data.frame(param=dist)))
  return(data.frame(measure=param_name,
                    country=country,
                    mean=mean(dist),
                    median=median(dist),
                    map=MAP,
                    lower_ci = quantile(dist,0.025,names=FALSE),
                    upper_ci = quantile(dist,0.975,names=FALSE),
                    lower_hpd = hpd[1],
                    upper_hpd = hpd[2]))
}
datasim4smc <- function(betaa_vals,
                        betaa_times,
                        init_EIR,
                        model_file="MiP_odin_model_nodelay.R",
                        comparison_group=NULL){
  init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
  prop_treated <- 0.4
  het_brackets <- 5
  ##Convert init_EIR into an init_betaa
  mpl_initial <- mamasante::model_param_list_create(init_EIR = init_EIR,
                                                 init_ft = prop_treated,
                                                 comparison = comparison_group
  )

  pars_initial <- mamasante::equilibrium_init_create_stripped(age_vector = init_age,
                                                           init_EIR = init_EIR,
                                                           ft = prop_treated,
                                                           model_param_list = mpl_initial,
                                                           het_brackets = het_brackets)
  init_betaa <- pars_initial$betaa_eq

  betaa_times <- append(c(0),betaa_times)

  betaa_vals <- append(c(init_betaa),betaa_vals)

  ##set up the simulation for the simualted data
  time<- max(betaa_times) + 30
  out_step=1

  mpl <- mamasante::model_param_list_create(init_EIR = init_EIR,
                                         init_ft = prop_treated,
                                         betaa_times=betaa_times,
                                         betaa_vals=betaa_vals,
                                         lag_rates = 10,
                                         comparison = comparison_group
  )

  pars <- mamasante::equilibrium_init_create_stripped(age_vector = init_age,
                                                   init_EIR = init_EIR,
                                                   ft = prop_treated,
                                                   model_param_list = mpl,
                                                   het_brackets = het_brackets)

  ##The malaria modelMiP
  generator <- odin(model_file)
  state_use <- pars[names(pars) %in% coef(generator)$name]

  # create model with initial values
  mod <- generator(user = state_use, use_dde = TRUE)
  tt <- seq(0, time, out_step)

  # run the simulation to base the data
  mod_run <- mod$run(tt, step_max_n = 1e7,
                     atol = 1e-5,
                     rtol = 1e-5)

  # shape output
  out <- mod$transform_variables(mod_run)

  # plot data and generate data
  # plot(out$t,out$betaa_out,type='l',col="red",ylim=c(0,125))
  # lines(out$t,out$prev*100,col="blue",lwd=4)
  out_df <- data.frame(t=out$t,
                       prev05=out$prev05,
                       EIR_true=out$EIR_out,
                       betaa_true=out$betaa_out,
                       inc05_true=out$inc05,
                       inc_all_true = out$inc,
                       prev_all_true = out$prev_all,
                       mv = out$mv_out,
                       spz_rate = out$spz_rate)

  return(out_df)
}
return_msi <- function(dataframe,single_year=FALSE){
  dataframe$year <- year(dataframe$date)
  dataframe$month <- month(dataframe$date)
  if(single_year){
    if(nrow(dataframe)!=12){
      stop('There need to be 12 months of data to calculate MSI for a single year.')
    }
    dataframe$year = 1
  }
  output <- dataframe %>%
    left_join(month_arc, by=join_by(month))%>%
    mutate(sin = value*sin(radians),
           cos = value*cos(radians))%>%
    group_by(year)%>%
    summarise(rk = sqrt(sum(sin)^2 + sum(cos)^2),
              theta_k = atan2(sum(sin),sum(cos)),
              theta_k_deg = (theta_k*180/pi + 360) %% 360,
              msi = rk/sum(value),
              month_num=n()) %>%
    filter(month_num==12)

  return(output)
}
return_vectors <- function(dataframe,single_year=FALSE){
  # dataframe$year <- year(dataframe$date)
  # dataframe$month <- month(dataframe$date)
  if(single_year){
    if(nrow(dataframe)!=12){
      stop('There need to be 12 months of data to calculate MSI for a single year.')
    }
    dataframe$year = 1
  }
  output <- dataframe %>%
    left_join(month_arc, by=join_by(month))%>%
    mutate(sin = value*sin(radians),
           cos = value*cos(radians))%>%
    group_by(year)%>%
    mutate(month_num = n(),
           x_coord = NA,
           y_coord = NA) %>%
    filter(month_num==12)
  output <- bind_rows(output,
                      data.frame(year=unique(output$year),
                                 month=0,
                                 x_coord=0,
                                 y_coord=0)) %>%
    arrange(year,month)
  for(x in unique(output$year)){
    for(y in 1:12){
      if(y==1){
        output[output$year==x & output$month==y,'x_coord']=0 + output[output$year==x & output$month==y,'value'] * cos(output[output$year==x & output$month==y,'radians'])
        output[output$year==x & output$month==y,'y_coord']=0 + output[output$year==x & output$month==y,'value'] * sin(output[output$year==x & output$month==y,'radians'])
      } else{
        output[output$year==x & output$month==y,'x_coord']=output[output$year==x & output$month==y-1,'x_coord'] + output[output$year==x & output$month==y,'value'] * cos(output[output$year==x & output$month==y,'radians'])
        output[output$year==x & output$month==y,'y_coord']=output[output$year==x & output$month==y-1,'y_coord'] + output[output$year==x & output$month==y,'value'] * sin(output[output$year==x & output$month==y,'radians'])

      }
    }
  }

  return(output)
}

calc_case_conc <- function(dataframe,window_size){
  dataframe$year <- year(dataframe$date)
  dataframe$month <- month(dataframe$date)

  if(nrow(dataframe)<12){
    stop('There need to be at least 12 months of data to calculate seasonality.')
  } else if(nrow(dataframe)>12){
    dataframe <- dataframe %>%
      group_by(month)%>%
      summarise(value = mean(value))
  }
  annual_total <- sum(dataframe$value)
  period_props <- bind_rows(lapply(1:12,function(start,window_size){
    index <- get_period_index(x=start,window_size=window_size)
    data.frame(months = paste0(index[1],'-',index[window_size]),
               prop = sum(dataframe[dataframe$month %in% index,]$value)/annual_total)
  },window_size=window_size))
  if(anyDuplicated(period_props$prop)>0){
    warning('Mulitple periods with equal value')
  }
  return(period_props[period_props$prop==max(period_props$prop),])
}
calc_annual_proportion <- function(dataframe){
  dataframe$year <- year(dataframe$date)
  dataframe$month <- month(dataframe$date)

  if(nrow(dataframe)<12){
    stop('There need to be at least 12 months of data to calculate seasonality.')
  } else if(nrow(dataframe)>12){
    dataframe <- dataframe %>%
      group_by(month)%>%
      summarise(value = mean(value))
  }
  annual_total <- sum(dataframe$value)
  period_props <- bind_rows(lapply(1:12,function(x){
    data.frame(prop = sum(dataframe[dataframe$month == x,]$value)/annual_total,
               month=x)
  }))
  if(anyDuplicated(period_props$prop)>0){
    warning('Mulitple periods with equal value')
  }
  return(period_props)
}

get_period_index <- function(x,window_size){
  index_1 <- c(x:(x+window_size-1))
  index_2 <- index_1 %% 12
  return(replace(index_2,index_2==0,12))
}
# Function to calculate the mean of a 12-row block
calc_summary <- function(start_row, data, window_size = 12) {
  end_row <- start_row + window_size - 1
  if (end_row > nrow(data)) {
    return(NULL)  # Handle cases where the window exceeds the data frame
  }
  return(return_msi(data[start_row:end_row,],single_year = TRUE))
}
process_group <- function(group_data) {
  start_rows <- 1:(nrow(group_data) - 12 + 1)
  results <- lapply(start_rows, calc_summary, data = group_data)
  results <- do.call(rbind, results)  # Combine results into a single dataframe
  results_out <- results %>%
    summarise(msi = mean(msi),
              theta_k_deg = mean(theta_k_deg))
  results_out$country <- unique(group_data$country)
  results_out$variable <- unique(group_data$variable)
  return(results_out)
}
process_group_peak <- function(group_data) {
  group_data$month <- month(as.Date(group_data$date))
  results_3 <- calc_case_conc(group_data,3)
  results_3$period = '3 months'
  results_4 <- calc_case_conc(group_data,4)
  results_4$period = '4 months'
  results_out <- bind_rows(results_3,results_4)
  results_out$country <- unique(group_data$country)
  results_out$variable <- unique(group_data$variable)
  return(results_out)
}
return_rainfall_df <- function(shape,country){
  as.data.frame(t(terra::extract(my_raster, shape, fun = mean, na.rm=TRUE)))%>%
    rownames_to_column()%>%
    filter(rowname!='ID')%>%
    mutate(month=as.yearmon(unique_months$date),
           country=country)%>%
    rename(rainfall='V1')%>%
    select(!rowname)
}
align_dates <- function(data,yearshift,site_name){
  data[data$country==site_name,]$date_aligned <- data[data$country==site_name,]$date_aligned %m-% years(yearshift)
}
