weibull <- function(time, alpha = 3.4, beta = 39.34) {
  pow <- -(time/beta)^alpha
  y <- exp(pow)
  return(y)
}
start_sim = 1
end_sim=900
smc_times = c(561,592,622)
get_smc_profile<-function(start_sim,end_sim,smc_times){
  gaps <- diff(c(smc_times, end_sim))
  prop_prof <- unlist(lapply(gaps, function(gap_length) {
    1 - weibull(1:gap_length)
  }))
  SMC_vals=c(rep(1,times=smc_times[1]),prop_prof)

  return(SMC_vals)
}

smc_prof<-get_smc_profile(0,900,smc_times)
plot(smc_prof)

smc_prof


##Get betaa values from runs
##Run simulation to get daily incidence rates

sim_from_betaa <- function(n,init_EIR,sample){
  sims <- lapply(1:n,function(x){
    df <- data_gen_piecewise(init_EIR =init_EIR,
                             betaa_sample = sample[,c(1,x+1)])
    return(df[,c('t','inc05')])
  })

  combined_df <- sims[[1]][, "t", drop = FALSE]

  # Iteratively join each dataframe in the list to the combined dataframe
  for (i in seq_along(sims)) {
    combined_df <- combined_df %>%
      left_join(sims[[i]], by = "t") %>%
      rename(!!paste0("inc05.", i) := inc05)
  }
  return(combined_df)
}
mcmc_sample <- sample(c(101:1000), 100)
ghana_times <- trials_short_informed_results[[3]]$times
trials_short_informed_prepped[[3]]$times
ghana_betaa_sample <- trials_short_informed_prepped[[3]]$sample%>%
  filter(measure=='betaa')%>%
  mutate(t=rep(trials_short_informed_prepped[[3]]$times,100))
trials_short_informed_prepped[[3]]$summary
ghana_betaa_sample <- bind_cols(c(data.frame(t=trials_short_informed_results[[3]]$times[-1]),
                             data.frame(t(trials_short_informed_results[[3]]$history['betaa',mcmc_sample,-1]))))
ghana_init_eir <- trials_short_informed_results[[3]]$history['EIR',1,1]

ghana_sim <- sim_from_betaa(n=100,init_EIR=ghana_init_eir,sample=ghana_betaa_sample)

burkinafaso_betaa_sample <- bind_cols(c(data.frame(t=trials_short_informed_results[[1]]$times[-1]),
                                  data.frame(t(trials_short_informed_results[[1]]$history['betaa',mcmc_sample,-1]))))
burkinafaso_init_eir <- trials_short_informed_results[[1]]$history['EIR',1,1]

burkinafaso_sim <- sim_from_betaa(n=100,init_EIR=burkinafaso_init_eir,sample=burkinafaso_betaa_sample)

mali_betaa_sample <- bind_cols(c(data.frame(t=trials_short_informed_results[[4]]$times[-1]),
                                        data.frame(t(trials_short_informed_results[[4]]$history['betaa',mcmc_sample,-1]))))
mali_init_eir <- trials_short_informed_results[[4]]$history['EIR',1,1]

mali_sim <- sim_from_betaa(n=100,init_EIR=mali_init_eir,sample=mali_betaa_sample)

gambia_betaa_sample <- bind_cols(c(data.frame(t=trials_short_informed_results[[2]]$times[-1]),
                                 data.frame(t(trials_short_informed_results[[2]]$history['betaa',mcmc_sample,-1]))))
gambia_init_eir <- trials_short_informed_results[[2]]$history['EIR',1,1]

gambia_sim <- sim_from_betaa(n=100,init_EIR=gambia_init_eir,sample=gambia_betaa_sample)


matplot(combined_df[,2:101],type='l')
##Get SMC dates from lit
get_smc_times <- function(months,data){
  start_obs <- min(zoo::as.Date(zoo::as.yearmon(data$month)))#Month of first observation (in Date format)
  time_origin <- zoo::as.Date(paste0(lubridate::year(start_obs)-1,'-01-01')) #January 1 of year before observation (in Date format)
  data$date <- zoo::as.Date(zoo::as.yearmon(data$month), frac = 0.5) #Convert dates to middle of month
  data$t <- as.integer(difftime(data$date,time_origin,units="days")) #Calculate date as number of days since January 1 of year before observation
  smc_times <- data[month(data$date)%in%months,]$t
  names(smc_times) <- data[month(data$date)%in%months,]$date
  return(smc_times)
}
bf_smc_times <- get_smc_times(months=c(7,8,9),data=trials_pg_data_all_list$`Burkina Faso`)
mali_smc_times <- get_smc_times(months=c(8,9,10),data=trials_pg_data_all_list$Mali)
gambia_smc_times <- get_smc_times(months=c(9,10,11),data=trials_pg_data_all_list$Gambia)

##Calculate SMC profiles
smc_prof_bf<-get_smc_profile(1,nrow(burkinafaso_sim),bf_smc_times)
smc_prof_mali<-get_smc_profile(1,nrow(mali_sim),mali_smc_times)
smc_prof_gambia<-get_smc_profile(1,nrow(gambia_sim),gambia_smc_times)

##Apply SMC profile to daily incidence rates
create_smc_sim <- function(sim_df,smc_profile,exclude_col=c('t','date')){
  df_excluded <- sim_df[,exclude_col,drop=FALSE]
  df_to_multiply <- sim_df[,!names(sim_df)%in%exclude_col]
  df_multiplied <- df_to_multiply * smc_profile
  df_result <- cbind(df_excluded,df_multiplied)
  return(df_result)
}
burkinafaso_sim_smc <- create_smc_sim(sim_df=burkinafaso_sim,smc_profile = smc_prof_bf)
mali_sim_smc <- create_smc_sim(sim_df=mali_sim,smc_profile = smc_prof_mali)
gambia_sim_smc <- create_smc_sim(sim_df=gambia_sim,smc_profile = smc_prof_gambia)

##Plot Incidence time series
burkinafaso_sim_long <- melt(burkinafaso_sim,id='t')%>%
  mutate(country='Burkina Faso',
         smc = 'Control')
mali_sim_long <- melt(mali_sim,id='t')%>%
  mutate(country='Mali',
         smc = 'Control')
gambia_sim_long <- melt(gambia_sim,id='t')%>%
  mutate(country='The Gambia',
         smc = 'Control')

burkinafaso_sim_smc_long <- melt(burkinafaso_sim_smc,id='t')%>%
  mutate(country='Burkina Faso',
         smc = 'SP+AQ')
mali_sim_smc_long <- melt(mali_sim_smc,id='t')%>%
  mutate(country='Mali',
         smc = 'SP+AQ')
gambia_sim_smc_long <- melt(gambia_sim_smc,id='t')%>%
  mutate(country='The Gambia',
         smc = 'SP+AQ')

smc_sim_all <- bind_rows(burkinafaso_sim_long,mali_sim_long,gambia_sim_long,
                         burkinafaso_sim_smc_long,mali_sim_smc_long,gambia_sim_smc_long)%>%
  group_by(country,smc,t)%>%
  summarise(median=median(value),
            lower=quantile(value,0.025),
            upper=quantile(value,0.975))
time_breaks <- get_smc_times(months=c(1:12),data=trials_pg_data_all_list$`Burkina Faso`)[1:13]
time_labels <- c('Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun')
smc_sim_plots <- ggplot(smc_sim_all)+
  geom_line(aes(x=t,y=median*1000,color=smc),linewidth=1)+
  facet_wrap(.~country)+
  scale_x_continuous(breaks=time_breaks,labels = time_labels)+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(xlim=c(500,900),ylim=c(0,18))+
  labs(y='Daily clinical malaria case\nper 1000 children under 5')+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())
ggsave('./trial_figs/output/smc_sim_plots.tiff',plot=smc_sim_plots,width=7,height=3.5,units='in')

##Calculate SMC effectiveness
gambia_smc_times
smc_sim_control <- bind_rows(burkinafaso_sim_long,mali_sim_long,gambia_sim_long)
smc_sim_spaq <- bind_rows(burkinafaso_sim_smc_long,mali_sim_smc_long,gambia_sim_smc_long)
smc_sim_rr <- left_join(smc_sim_control,smc_sim_spaq,by=c('t','variable','country'),suffix=c('.control','.spaq'))%>%
  mutate(fu_period =case_when(
    country == 'Burkina Faso' & t >= 561 & t < 622+45 ~ "Yes",
    country == 'Mali' & t >= 592 & t < 653+45 ~ "Yes",
    country == 'The Gambia' & t >= 622 & t < 683+45 ~ "Yes",
    TRUE ~ 'No'))%>%
  filter(fu_period=='Yes')%>%
  group_by(variable,country)%>%
  summarise(inc_control = sum(value.control),
            inc_spaq = sum(value.spaq))%>%
  mutate(pe = 1-(inc_spaq/inc_control))%>%
  group_by(country)%>%
  summarise(median.pe = median(pe),
            lower.pe=quantile(pe,0.025),
            upper.pe=quantile(pe,0.975)
            )
ggplot(smc_sim_rr)+
  geom_violin(aes(x=country,y=pe))
##Plot SMC effectiveness
smc_trial_pe <- data.frame(country = c('Burkina Faso','Mali','The Gambia'),
                           median.obs = c(0.71,0.83,0.93),
                           lower.obs = c(0.68,0.80,0.80),
                           upper.obs = c(0.74,0.86,0.98))
smc_pe_both <- left_join(smc_sim_rr,smc_trial_pe,by='country')

pe_comparison_smc <- ggplot(smc_pe_both)+
  geom_point(aes(x=median.obs,y=median.pe))+
  geom_errorbar(aes(x=median.obs,ymin=lower.pe,ymax=upper.pe),width=0)+
  geom_errorbarh(aes(y=median.pe,xmin=lower.obs,xmax=upper.obs),height=0)+
  geom_abline(linetype='dashed')+
  coord_cartesian(xlim=c(0.5,1),ylim=c(0.5,1))+
  labs(x='Protective efficacy in trial',
       y='Protective efficacy in simulation')

##Malawi SMC sim
##3 Rounds
malawi_smc_times_3 <- get_smc_times(months=c(1,2,12),data=trials_pg_data_all_list$Malawi)
amoah_smc_times_3 <- get_smc_times(months=c(1,2,12),data=amoah_prev_u5)
rogerson_smc_times_3 <- get_smc_times(months=c(1,2,12),data=rogerson_prev)

smc_months <- as.yearmon(as.Date(c(names(malawi_smc_times_3),names(amoah_smc_times_3),names(rogerson_smc_times_3))))
malawi_betaa_sample <- bind_cols(c(data.frame(t=trials_short_informed_results[[6]]$times[-1]),
                                   data.frame(t(trials_short_informed_results[[6]]$history['betaa',mcmc_sample,-1]))))
malawi_init_eir <- trials_short_informed_results[[6]]$history['EIR',1,1]

malawi_sim <- sim_from_betaa(n=100,init_EIR=malawi_init_eir,sample=malawi_betaa_sample)
malawi_sim$date <- seq.Date(from=as.Date(names(malawi_smc_times_3)[1])-malawi_smc_times_3[[1]],length.out = nrow(malawi_sim),by='day')
malawi_sim_long <- malawi_sim%>%
  melt(id=c('t','date'))%>%
  mutate(source=malawi_sources[1],
         smc = 'Control')

amoah_betaa_sample <- bind_cols(c(data.frame(t=amoah_u5_short_2_results$times[-1]),
                                   data.frame(t(amoah_u5_short_2_results$history['betaa',mcmc_sample,-1]))))
amoah_init_eir <- amoah_u5_short_2_results$history['EIR',1,1]

amoah_sim <- sim_from_betaa(n=100,init_EIR=amoah_init_eir,sample=amoah_betaa_sample)
amoah_sim$date <- seq.Date(from=as.Date(names(amoah_smc_times_3)[1])-amoah_smc_times_3[[1]],length.out = nrow(amoah_sim),by='day')
amoah_sim_long <- amoah_sim%>%
  melt(id=c('date','t'))%>%
  mutate(source=malawi_sources[2],
         smc = 'Control')

rogerson_betaa_sample <- bind_cols(c(data.frame(t=rogerson_prev_short_2_results$times[-1]),
                                  data.frame(t(rogerson_prev_short_2_results$history['betaa',mcmc_sample,-1]))))
rogerson_init_eir <- rogerson_prev_short_2_results$history['EIR',1,1]

rogerson_sim <- sim_from_betaa(n=100,init_EIR=rogerson_init_eir,sample=rogerson_betaa_sample)
rogerson_sim$date <- seq.Date(from=as.Date(names(rogerson_smc_times_3)[1])-rogerson_smc_times_3[[1]],length.out = nrow(rogerson_sim),by='day')
rogerson_sim_long <- rogerson_sim%>%
  melt(id=c('date','t'))%>%
  mutate(source=malawi_sources[3],
         smc = 'Control')


smc_prof_malawi_3<-get_smc_profile(1,nrow(malawi_sim),malawi_smc_times_3)
smc_prof_amoah_3<-get_smc_profile(1,nrow(amoah_sim),amoah_smc_times_3)
smc_prof_rogerson_3<-get_smc_profile(1,nrow(rogerson_sim),rogerson_smc_times_3)

malawi_sim_smc_3 <- create_smc_sim(sim_df=malawi_sim,smc_profile = smc_prof_malawi_3)%>%
  melt(id=c('date','t'))%>%
  mutate(source=malawi_sources[1],
         smc = 'SP+AQ')
amoah_sim_smc_3 <- create_smc_sim(sim_df=amoah_sim,smc_profile = smc_prof_amoah_3)%>%
  melt(id=c('date','t'))%>%
  mutate(source=malawi_sources[2],
         smc = 'SP+AQ')
rogerson_sim_smc_3 <- create_smc_sim(sim_df=rogerson_sim,smc_profile = smc_prof_rogerson_3)%>%
  melt(id=c('date','t'))%>%
  mutate(source=malawi_sources[3],
         smc = 'SP+AQ')

malawi_smc_sim_all <- bind_rows(malawi_sim_long,amoah_sim_long,rogerson_sim_long,
                         malawi_sim_smc_3,amoah_sim_smc_3,rogerson_sim_smc_3)%>%
  group_by(source,smc,t,date)%>%
  summarise(median=median(value),
            lower=quantile(value,0.025),
            upper=quantile(value,0.975))%>%
  mutate(source_smc = factor(paste0(source,'-',smc),levels=c( "Rogerson, et al. 2000-SP+AQ","Madanitsa, et al. 2016-SP+AQ","Amoah, et al. 2021-SP+AQ",
                                                              "Rogerson, et al. 2000-Control","Madanitsa, et al. 2016-Control","Amoah, et al. 2021-Control"))
  )

colors_sources <- c(`Madanitsa, et al. 2016`=tri_vir[2],`Amoah, et al. 2021`=tri_vir[3],`Rogerson, et al. 2000`=tri_vir[1])

smc_source_pal <- c(tri_vir,lighten(tri_vir,0.4))
names(smc_source_pal) <- c("Rogerson, et al. 2000-Control","Madanitsa, et al. 2016-Control","Amoah, et al. 2021-Control",
                           "Rogerson, et al. 2000-SP+AQ","Madanitsa, et al. 2016-SP+AQ","Amoah, et al. 2021-SP+AQ")
str(malawi_smc_sim_all)
multiplier_inc <- 1000
malawi_rainfall_smc <- malawi_rainfall%>%
  mutate(rainfall_rel = rainfall_rel*15/500)
malawi_rainfall_smc_only <- malawi_rainfall_smc%>%
  filter(month %in% smc_months)

plot_malawi_smc <- ggplot()+
  geom_col(data=malawi_rainfall_smc,aes(x=date,y=rainfall_rel),alpha = 1,fill = 'darkgrey',just=0)+
  geom_col(data=malawi_rainfall_smc_only,aes(x=date,y=rainfall_rel),alpha = 1,fill = '#C3A995FF',just=0)+
  scale_y_continuous(expand=c(0,0),limits=c(0,12.5),sec.axis = sec_axis(~ . /(rainfall_multiplier*15/max(rainfall$rainfall)), name = "Monthly rainfall (cm)"))+
  theme(strip.placement = "outside")+
  geom_line(data=malawi_smc_sim_all,aes(x=as.Date(date),y=median*multiplier_inc,color=source_smc,group=source_smc),linewidth=0.8)+
  scale_color_manual(values=smc_source_pal)+
  # facet_grid(.~source,scales = 'free_x')+
  scale_x_break(c(as.Date('1998-07-15'),as.Date('2011-08-15')), space = 0.2) +
  scale_x_break(c(as.Date('2012-11-15'),as.Date('2016-04-15')), space = 0.2) +
  scale_x_date(date_labels = "%b '%y", date_breaks = '3 months',limits =c(min(results_eir_summary$date),as.Date('2017-11-15')) )+
  labs(y = y_axis_label)+
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.length = unit(3, "pt"),
        legend.position = 'none',
        panel.spacing.y = unit(3, "mm"),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()
  )
ggsave('./trial_figs/output/plot_malawi_smc.tiff',plot=plot_malawi_smc,width=7,height=2,units='in')
