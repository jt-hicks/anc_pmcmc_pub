prep_results_sim <- function(results,sim_data_true,sim_data_raw,burnin,n_chains=1,site,timelength=NA,anc=FALSE){
  ind <- get_chain_info(n_chains=n_chains,
                        length=length(results$history[1,,1]),
                        burnin=burnin)
  num_months <- length(results$history[1,1,])
  # history <- results$history[,start_chain:length, -c(1,2)]

  if(!('date'%in%names(sim_data_raw))){
    sim_data_raw$date <- zoo::as.Date.yearmon(sim_data_raw$month,frac=0.5)
  }
  sim_data_raw$date <- as.character(as.Date(sim_data_raw$date))
  dates <- as.Date(sim_data_raw$date)

  months <- unique(zoo::as.yearmon(sim_data_true$date))
  midmonth_dates <- data.frame(date=as.character(as.Date(months,frac=0.5)))
  sim_data_true$date <- as.character(sim_data_true$date)
  monthly_data_true <- left_join(midmonth_dates,sim_data_true,by='date')

  sim_data <- left_join(monthly_data_true,sim_data_raw,by=c('date','admin','country','site'),suffix=c('.true','.raw'))%>%
    select(-prev_05)%>%
    rename(prev_05 = prev05_true,
           clininc_05=inc05_true,
           clininc_all=inc_all_true,
           EIR=EIR_true,
           betaa=betaa_true)%>%
    select(date,prev_05,clininc_05,clininc_all,EIR,betaa,country,admin,site,init_EIR)%>%
    reshape2::melt(id=c('date','site','country','admin','init_EIR'))%>%
    rename(true_value = value,
           measure = variable)%>%
    mutate(date=as.Date(date),
           month=as.yearmon(date))

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
    reshape2::melt(id='date')%>%
    group_by(date)%>%
    dplyr::summarise(median=median(value),
                     mean=mean(value),
                     upper=quantile(value,probs=0.975),
                     lower=quantile(value,probs=0.025))%>%
    mutate(date=as.Date(date),
           measure = 'prev_05')
  prev_means <- history.df.prev%>%
    dplyr::mutate(date=dates)%>%
    reshape2::melt(id='date')%>%
    group_by(variable)%>%
    summarise(mean = mean(value))
  prev_summary <- prev_means %>%
    summarise(median = median(mean))
  prev_hpd <- coda::HPDinterval(coda::as.mcmc(prev_means$mean))
  prev_summary$lower_hpd <- prev_hpd[1]
  prev_summary$upper_hpd <- prev_hpd[2]
  prev_summary$measure <- 'prev_05'
  prev_sample <- history.df.prev[,index_sample] %>%
    mutate(t=c(1:nrow(prev_history)))%>%
    reshape2::melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100),
           measure = 'prev_05')


  # history.df.inc05 <- as.data.frame(t(history['clininc_05',,]))
  # inc_history <- history.df.inc05%>%
  #   dplyr::mutate(date=dates,
  #                 month=as.yearmon(date))%>%
  #   reshape2::melt(id='month')%>%
  #   left_join(sim_data%>%filter(measure=='clininc_05'),by='month')%>%
  #   group_by(month,true_value)%>%
  #   dplyr::summarise(median=median(value),
  #                    mean=mean(value),
  #                    upper=quantile(value,probs=0.975),
  #                    lower=quantile(value,probs=0.025),
  #                    upper_50=quantile(value,probs=0.75),
  #                    lower_50=quantile(value,probs=0.25),
  #                    med_rel_error = median((value-true_value)/true_value),
  #                    med_abs_error = median(abs(value-true_value)),
  #                    rmse = sqrt(mean((value-true_value)^2)),
  #                    nrmse_range = rmse/(max(value)-min(value)),
  #                    nrmse_mean = rmse/mean(value),
  #                    mad = median(abs(value-median(value))))%>%
  #   mutate(date=as.Date(month,frac=0.5),
  #          measure = 'inc05',
  #          truth_in_quant = ifelse(lower<=true_value&upper>=true_value,1,0),
  #          truth_in_quant_50 = ifelse(lower_50<=true_value&upper_50>=true_value,1,0))
  # inc_sample <- history.df.inc05[,index_sample] %>%
  #   mutate(t=c(1:nrow(inc_history)))%>%
  #   reshape2::melt(id='t')%>%
  #   dplyr::rename(time=t)%>%
  #   mutate(date = rep(dates,100),
  #          measure = 'inc05',
  #          month = as.yearmon(date))%>%
  #   left_join(sim_data%>%filter(measure=='clininc_05'),by='month')%>%
  #   mutate(abs_error = abs(value-true_value),
  #          rel_error = (value-true_value)/true_value,
  #          rel_abs_error = abs(value-true_value)/true_value)
  #
  # inc_means <- history.df.inc05%>%
  #   dplyr::mutate(date=dates,
  #                 month=as.yearmon(date))%>%
  #   reshape2::melt(id='month')%>%
  #   left_join(sim_data%>%filter(measure=='clininc_05'),by='month')%>%
  #   mutate(error = value-true_value,
  #          abs_error = abs(value-true_value),
  #          rel_error = (value-true_value)/true_value,
  #          rel_abs_error = abs(value-true_value)/true_value)%>%
  #   group_by(variable)%>%
  #   summarise(mean = mean(value),
  #             median = median(value),
  #             median_abs_error = median(abs_error),
  #             median_rel_error = median(rel_error),
  #             median_rel_abs_error = median(rel_abs_error),
  #             rmse = sqrt(mean((error)^2)),
  #             nrmse_mean = rmse/mean
  #             )
  # inc_summary <- inc_means %>%
  #   summarise(median_mean = median(mean),
  #             median_median = median(median),
  #             median_mae = median(median_abs_error),
  #             median_mre = median(median_rel_error),
  #             median_mrae = median(median_rel_abs_error),
  #             median_rmse = median(rmse),
  #             median_nrmse = median(nrmse_mean))
  # inc_hpd <- HPDinterval(as.mcmc(inc_means$mean))
  # inc_summary <- get_hpd(df_means=inc_means,df_summary=inc_summary)
  # inc_summary$truth_captured <- mean(inc_history$truth_in_quant)
  # inc_summary$truth_captured_50 <- mean(inc_history$truth_in_quant_50)
  # inc_summary$measure <- 'inc05'
  #
  # history.df.incall <- as.data.frame(t(history['clininc_all',,]))
  # incall_history <- history.df.incall%>%
  #   dplyr::mutate(date=dates,
  #                 month=as.yearmon(date))%>%
  #   reshape2::melt(id='month')%>%
  #   left_join(sim_data%>%filter(measure=='clininc_all'),by='month')%>%
  #   group_by(month,true_value)%>%
  #   dplyr::summarise(median=median(value),
  #                    mean=mean(value),
  #                    upper=quantile(value,probs=0.975),
  #                    lower=quantile(value,probs=0.025),
  #                    med_rel_error = median((value-true_value)/true_value),
  #                    med_abs_error = median(abs(value-true_value)),
  #                    rmse = sqrt(mean((value-true_value)^2)),
  #                    nrmse_range = rmse/(max(value)-min(value)),
  #                    nrmse_mean = rmse/mean(value),
  #                    mad = median(abs(value-median(value))))%>%
  #   mutate(date=as.Date(month,frac=0.5),
  #          measure = 'incall',
  #          truth_in_quant = ifelse(lower<=true_value&upper>=true_value,1,0))
  # incall_sample <- history.df.incall[,index_sample] %>%
  #   mutate(t=c(1:nrow(incall_history)))%>%
  #   reshape2::melt(id='t')%>%
  #   dplyr::rename(time=t)%>%
  #   mutate(date = rep(dates,100),
  #          measure = 'incall')
  # incall_means <- history.df.incall%>%
  #   dplyr::mutate(date=dates)%>%
  #   reshape2::melt(id='date')%>%
  #   group_by(variable)%>%
  #   summarise(mean = mean(value))
  # incall_summary <- inc_means %>%
  #   summarise(median = median(mean))
  # incall_hpd <- HPDinterval(as.mcmc(inc_means$mean))
  # incall_summary$lower_hpd <- incall_hpd[1]
  # incall_summary$upper_hpd <- incall_hpd[2]
  # incall_summary$measure <- 'incall'
  # incall_summary$truth_captured <- mean(incall_history$truth_in_quant)
  #
  # history.df.EIR <- as.data.frame(t(history['EIR',,]))
  # EIR_history <- history.df.EIR%>%
  #   dplyr::mutate(date=dates,
  #                 month=as.yearmon(date))%>%
  #   reshape2::melt(id='month')%>%
  #   left_join(sim_data%>%filter(measure=='EIR'),by='month')%>%
  #   group_by(month,true_value)%>%
  #   dplyr::summarise(median=median(value),
  #                    mean=mean(value),
  #                    upper=quantile(value,probs=0.975),
  #                    lower=quantile(value,probs=0.025),
  #                    med_rel_error = median((value-true_value)/true_value),
  #                    med_abs_error = median(abs(value-true_value)),
  #                    rmse = sqrt(mean((value-true_value)^2)),
  #                    nrmse_range = rmse/(max(value)-min(value)),
  #                    nrmse_mean = rmse/mean(value),
  #                    mad = median(abs(value-median(value))))%>%
  #   mutate(date=as.Date(month,frac=0.5),
  #          measure = 'EIR',
  #          truth_in_quant = ifelse(lower<=true_value&upper>=true_value,1,0))
  # EIR_sample <- history.df.EIR[,index_sample] %>%
  #   mutate(t=c(1:nrow(EIR_history)))%>%
  #   reshape2::melt(id='t')%>%
  #   dplyr::rename(time=t)%>%
  #   mutate(date = rep(dates,100),
  #          measure = 'EIR')
  # EIR_means <- history.df.EIR%>%
  #   dplyr::mutate(date=dates)%>%
  #   reshape2::melt(id='date')%>%
  #   group_by(variable)%>%
  #   summarise(mean = mean(value))
  # EIR_summary <- EIR_means %>%
  #   summarise(median = median(mean))
  # EIR_hpd <- HPDinterval(as.mcmc(EIR_means$mean))
  # EIR_summary$lower_hpd <- EIR_hpd[1]
  # EIR_summary$upper_hpd <- EIR_hpd[2]
  # EIR_summary$measure <- 'EIR'
  # EIR_summary$truth_captured <- mean(EIR_history$truth_in_quant)
  #
  # history.df.betaa <- as.data.frame(t(history['betaa',,]))
  # betaa_history <- history.df.betaa%>%
  #   dplyr::mutate(date=dates,
  #                 month=as.yearmon(date))%>%
  #   reshape2::melt(id='month')%>%
  #   left_join(sim_data%>%filter(measure=='betaa'),by='month')%>%
  #   group_by(month,true_value)%>%
  #   dplyr::summarise(median=median(value),
  #                    mean=mean(value),
  #                    upper=quantile(value,probs=0.975),
  #                    lower=quantile(value,probs=0.025),
  #                    med_rel_error = median((value-true_value)/true_value),
  #                    med_abs_error = median(abs(value-true_value)),
  #                    rmse = sqrt(mean((value-true_value)^2)),
  #                    nrmse_range = rmse/(max(value)-min(value)),
  #                    nrmse_mean = rmse/mean(value),
  #                    mad = median(abs(value-median(value))))%>%
  #   mutate(date=as.Date(month,frac=0.5),
  #          measure = 'betaa',
  #          truth_in_quant = ifelse(lower<=true_value&upper>=true_value,1,0))
  # betaa_means <- history.df.betaa%>%
  #   dplyr::mutate(date=dates)%>%
  #   reshape2::melt(id='date')%>%
  #   group_by(variable)%>%
  #   summarise(mean = mean(value))
  # betaa_summary <- betaa_means %>%
  #   summarise(median = median(mean))
  # betaa_hpd <- HPDinterval(as.mcmc(betaa_means$mean))
  # betaa_summary$lower_hpd <- betaa_hpd[1]
  # betaa_summary$upper_hpd <- betaa_hpd[2]
  # betaa_summary$measure <- 'betaa'
  # betaa_summary$truth_captured <- mean(betaa_history$truth_in_quant)
  # betaa_sample <- history.df.betaa[,index_sample] %>%
  #   mutate(t=c(1:nrow(betaa_history)))%>%
  #   reshape2::melt(id='t')%>%
  #   dplyr::rename(time=t)%>%
  #   mutate(date = rep(dates,100),
  #          measure = 'betaa')
  #
  inc05_summaries <- get_pmcmc_summaries(history=history,measure_name = 'clininc_05',index_sample = index_sample,sim_data = sim_data, dates = dates)
  incall_summaries <- get_pmcmc_summaries(history=history,measure_name = 'clininc_all',index_sample = index_sample,sim_data = sim_data, dates = dates)
  EIR_summaries <- get_pmcmc_summaries(history=history,measure_name = 'EIR',index_sample = index_sample,sim_data = sim_data, dates = dates)
  betaa_summaries <- get_pmcmc_summaries(history=history,measure_name = 'betaa',index_sample = index_sample,sim_data = sim_data, dates = dates)

  history_summary <- bind_rows(prev_history,
                               inc05_summaries[['output_history']],
                               incall_summaries[['output_history']],
                               EIR_summaries[['output_history']],
                               betaa_summaries[['output_history']])
  sample <- bind_rows(prev_sample,
                      inc05_summaries[['output_sample']],
                      incall_summaries[['output_sample']],
                      EIR_summaries[['output_sample']],
                      betaa_summaries[['output_sample']])
  point_est_summary <- bind_rows(prev_summary,
                                 inc05_summaries[['output_summary']],
                                 incall_summaries[['output_summary']],
                                 EIR_summaries[['output_summary']],
                                 betaa_summaries[['output_summary']])
  means_summary <- bind_rows(prev_means,
                                 inc05_summaries[['output_means']],
                                 incall_summaries[['output_means']],
                                 EIR_summaries[['output_means']],
                                 betaa_summaries[['output_means']])

  history_summary$site <- site
  sample$site <- site
  point_est_summary$site <- site
  means_summary$site <- site

  return(list(summary=history_summary,
              sample = sample,
              times = times,
              point_est = point_est_summary,
              means = means_summary))
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
get_hpd <- function(df_means,df_summary){
  hpd_all <- bind_cols(lapply(names(df_means),function(name){
    if(name == 'variable'){
      return()
    }
    hpd <- coda::HPDinterval(coda::as.mcmc(df_means[,name]))
    hpd_df <- data.frame(hpd[1],hpd[2])
    names(hpd_df) <- c(paste0(name,'_hpd_lower'),paste0(name,'_hpd_upper'))
    return(hpd_df)
  }))

  return(bind_cols(df_summary,hpd_all))

}
get_pmcmc_summaries <- function(history,measure_name,index_sample,sim_data,dates){
  history.df <- as.data.frame(t(history[measure_name,,]))
  output_history <- history.df%>%
    dplyr::mutate(date=dates,
                  month=as.yearmon(date))%>%
    reshape2::melt(id='month')%>%
    left_join(sim_data%>%filter(measure==measure_name),by='month')%>%
    group_by(month,true_value)%>%
    dplyr::summarise(median=median(value),
                     mean=mean(value),
                     upper=quantile(value,probs=0.975),
                     lower=quantile(value,probs=0.025),
                     upper_50=quantile(value,probs=0.75),
                     lower_50=quantile(value,probs=0.25),
                     med_rel_error = median((value-true_value)/true_value),
                     med_abs_error = median(abs(value-true_value)),
                     rmse = sqrt(mean((value-true_value)^2)),
                     range = max(value)-min(value),
                     nrmse_range = rmse/(max(value)-min(value)),
                     nrmse_mean = rmse/mean(value),
                     mad = median(abs(value-median(value))))%>%
    ungroup()%>%
    mutate(date=as.Date(month,frac=0.5),
           measure = measure_name,
           truth_in_quant = ifelse(lower<=true_value&upper>=true_value,1,0),
           truth_in_quant_50 = ifelse(lower_50<=true_value&upper_50>=true_value,1,0))
  output_sample <- history.df[,index_sample] %>%
    mutate(t=c(1:nrow(output_history)))%>%
    reshape2::melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(date = rep(dates,100),
           measure = measure_name,
           month = as.yearmon(date))%>%
    left_join(sim_data%>%filter(measure==measure_name),by=c('month','measure'))%>%
    mutate(abs_error = abs(value-true_value),
           rel_error = (value-true_value)/true_value,
           rel_abs_error = abs(value-true_value)/true_value)

  output_means <- history.df%>%
    dplyr::mutate(date=dates,
                  month=as.yearmon(date))%>%
    reshape2::melt(id='month')%>%
    left_join(sim_data%>%filter(measure==measure_name),by='month')%>%
    mutate(error = value-true_value,
           abs_error = abs(value-true_value),
           rel_error = (value-true_value)/true_value,
           rel_abs_error = abs(value-true_value)/true_value)%>%
    group_by(variable)%>%
    summarise(mean = mean(value),
              median = median(value),
              median_abs_error = median(abs_error),
              median_rel_error = median(rel_error),
              median_rel_abs_error = median(rel_abs_error),
              rmse = sqrt(mean((error)^2)),
              nrmse_mean = rmse/mean
    )
  output_summary <- output_means %>%
    summarise(median_mean = median(mean),
              median_median = median(median),
              median_mae = median(median_abs_error),
              median_mre = median(median_rel_error),
              median_mrae = median(median_rel_abs_error),
              median_rmse = median(rmse),
              median_nrmse = median(nrmse_mean))
  output_summary <- get_hpd(df_means=output_means,df_summary=output_summary)
  output_summary$truth_captured <- mean(output_history$truth_in_quant)
  output_summary$truth_captured_50 <- mean(output_history$truth_in_quant_50)
  output_summary$measure <- measure_name
  return(list(output_history=output_history,output_means=output_means,output_summary=output_summary,output_sample=output_sample))

}
