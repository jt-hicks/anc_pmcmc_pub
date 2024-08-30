library(mamasante)
library(SBC)
library('tidyverse')
library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(didehpc)
library(pkgdepends)
library("coda")
library(binom)
library(ggplot2)
library(bayesplot)
library(reshape2)
library(ggpubr)
library(zoo)
library(patchwork)
library(RColorBrewer)
library(viridis)
source('utils.R')

theme_set(theme_minimal(base_size = 7)+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))


admin_units_seasonal <- mamasante::load_file("admin_units_seasonal.rds")
get_seasonal_profile <- function(admin,country,years=1){
  admin_matches <- mamasante::admin_match(admin_unit = admin, country = country,
                               admin_units_seasonal = admin_units_seasonal)

  ssa0 <- admin_units_seasonal$a0[admin_matches]
  ssa1 <- admin_units_seasonal$a1[admin_matches]
  ssa2 <- admin_units_seasonal$a2[admin_matches]
  ssa3 <- admin_units_seasonal$a3[admin_matches]
  ssb1 <- admin_units_seasonal$b1[admin_matches]
  ssb2 <- admin_units_seasonal$b2[admin_matches]
  ssb3 <- admin_units_seasonal$b3[admin_matches]
  theta_c <- admin_units_seasonal$theta_c[admin_matches]

  t <- c(1:365*years)
  return(data.frame(rainfall=pmax((ssa0+ssa1*cos(2*pi*-t/365)+ssa2*cos(2*2*pi*-t/365)+ssa3*cos(3*2*pi*-t/365)+ssb1*sin(2*pi*-t/365)+ssb2*sin(2*2*pi*-t/365)+ ssb3*sin(3*2*pi*-t/365))/theta_c,0.001),
                    t=t,
                    admin=admin,
                    country=country))
}
admin_list <- list('Tanga','Upper East','Fatick','Equateur')
country_list <- list('Tanzania','Ghana','Senegal','Democratic Republic of Congo')

seasonal_profiles <- bind_rows(lapply(1:4,function(x){
  get_seasonal_profile(admin=admin_list[[x]],country=country_list[[x]])
}))

ggplot(seasonal_profiles)+
  geom_line(aes(x=t,y=rainfall,group=admin,color=admin),size=1)+
  scale_color_viridis_d()

tanga_sim <- gen_seasonal_sim(init_EIR=50,
                              max_param=125,
                              model_file= "init/odin_model_stripped_seasonal.R",
                              country = 'Tanzania',
                              admin_unit = 'Tanga',
                              sim_length = 8)
tanga_sim$true_data$country <- 'Tanzania'
tanga_sim$true_data$admin <- 'Tanga'
tanga_sim$true_data$site <- 'Tanga, Tanzania'
tanga_sim$data_raw$country <- 'Tanzania'
tanga_sim$data_raw$admin <- 'Tanga'
tanga_sim$data_raw$site <- 'Tanga, Tanzania'
tanga_sim$true_monthly <- daily2monthly(out_df=tanga_sim$true_data)

plot(tanga_sim$true_data$t,tanga_sim$true_data$prev05_true)
plot(tanga_sim$data_raw$t,tanga_sim$data_raw$prev)
tanga_sim$data_raw$date
ue_ghana_sim <- gen_seasonal_sim(init_EIR=50,
                                  max_param=125,
                                  model_file= "init/odin_model_stripped_seasonal.R",
                                  country = 'Ghana',
                                  admin_unit = 'Upper East',
                                  sim_length = 8)
ue_ghana_sim$true_data$country <- 'Ghana'
ue_ghana_sim$true_data$admin <- 'Upper East'
ue_ghana_sim$true_data$site <- 'Upper East, The Ghana'
ue_ghana_sim$data_raw$country <- 'Ghana'
ue_ghana_sim$data_raw$admin <- 'Upper East'
ue_ghana_sim$data_raw$site <- 'Upper East, The Ghana'
ue_ghana_sim$true_monthly <- daily2monthly(out_df=ue_ghana_sim$true_data)

fatick_sen_sim <- gen_seasonal_sim(init_EIR=50,
                                   max_param=125,
                                   model_file= "init/odin_model_stripped_seasonal.R",
                                   country = 'Senegal',
                                   admin_unit = 'Fatick',
                                   sim_length = 8)
fatick_sen_sim$true_data$country <- 'Senegal'
fatick_sen_sim$true_data$admin <- 'Fatick'
fatick_sen_sim$true_data$site <- 'Fatick, Senegal'
fatick_sen_sim$data_raw$country <- 'Senegal'
fatick_sen_sim$data_raw$admin <- 'Fatick'
fatick_sen_sim$data_raw$site <- 'Fatick, Senegal'
fatick_sen_sim$true_monthly <- daily2monthly(out_df=fatick_sen_sim$true_data)

equateur_drc_sim <- gen_seasonal_sim(init_EIR=50,
                                     max_param=125,
                                     model_file= "init/odin_model_stripped_seasonal.R",
                                     country = 'Democratic Republic of Congo',
                                     admin_unit = 'Equateur',
                                     sim_length = 8)
equateur_drc_sim$true_data$country <- 'Democratic Republic of Congo'
equateur_drc_sim$true_data$admin <- 'Equateur'
equateur_drc_sim$true_data$site <- 'Equateur, DRC'
equateur_drc_sim$data_raw$country <- 'Democratic Republic of Congo'
equateur_drc_sim$data_raw$admin <- 'Equateur'
equateur_drc_sim$data_raw$site <- 'Equateur, DRC'
equateur_drc_sim$true_monthly <- daily2monthly(out_df=equateur_drc_sim$true_data)


season_true_sim_all <- bind_rows(tanga_sim$true_data,
                                 ue_ghana_sim$true_data,
                                 fatick_sen_sim$true_data,
                                 equateur_drc_sim$true_data)
season_truemonthly_sim_all <- bind_rows(tanga_sim$true_monthly,
                                        ue_ghana_sim$true_monthly,
                                        fatick_sen_sim$true_monthly,
                                        equateur_drc_sim$true_monthly)%>%
  rename(prev_05=prev05_true,
         incall=inc_all_true,
         EIR=EIR_true,
         betaa=betaa_true)%>%
  select(t,date,country,admin,site,month,prev_05,incall,EIR,betaa)%>%
  melt(id=c('t','date','country','admin','site','month'))%>%
  rename(measure=variable)%>%
  mutate(date=as.Date(date))

season_truemonthly_sim_list <- list(tanga_sim$true_monthly,
                                    ue_ghana_sim$true_monthly,
                                    fatick_sen_sim$true_monthly,
                                    equateur_drc_sim$true_monthly)
season_data_raw_sim_all <- bind_rows(tanga_sim$data_raw,
                                     ue_ghana_sim$data_raw,
                                     fatick_sen_sim$data_raw,
                                     equateur_drc_sim$data_raw)
season_data_raw_sim_list <- list(tanga_sim$data_raw,
                                 ue_ghana_sim$data_raw,
                                 fatick_sen_sim$data_raw,
                                 equateur_drc_sim$data_raw)
admin_list <- list('Tanga','Upper East','Fatick','Equateur')
country_list <- list('Tanzania','Ghana','Senegal','Democratic Republic of Congo')
season_data_raw_sim_all <- addCIs(season_data_raw_sim_all,Ys=season_data_raw_sim_all$positive,Ns=season_data_raw_sim_all$tested)
saveRDS(season_data_raw_sim_all,'./sim_seasonal_data/season_data_raw_sim_all.rds')
saveRDS(season_data_raw_sim_list,'./sim_seasonal_data/season_data_raw_sim_list.rds')
saveRDS(season_truemonthly_sim_list,'./sim_seasonal_data/season_truemonthly_sim_list.rds')
saveRDS(season_truemonthly_sim_all,'./sim_seasonal_data/season_truemonthly_sim_all.rds')
saveRDS(season_true_sim_all,'./sim_seasonal_data/season_true_sim_all.rds')

season_data_raw_sim_all <- readRDS('./sim_seasonal_data/season_data_raw_sim_all.rds')
season_data_raw_sim_list <- readRDS('./sim_seasonal_data/season_data_raw_sim_list.rds')
season_truemonthly_sim_list <- readRDS('./sim_seasonal_data/season_truemonthly_sim_list.rds')
season_truemonthly_sim_all <- readRDS('./sim_seasonal_data/season_truemonthly_sim_all.rds')
season_true_sim_all <- readRDS('./sim_seasonal_data/season_true_sim_all.rds')

season_data_raw_sim_all_summary <- season_data_raw_sim_all%>%
  group_by(site)%>%
  summarise(sum_total = sum(tested),
            sum_positive = sum(positive))
season_data_raw_sim_all_summary <- addCIs(season_data_raw_sim_all_summary,season_data_raw_sim_all_summary$sum_positive,season_data_raw_sim_all_summary$sum_total)

compare_prev <- ggplot(season_data_raw_sim_all)+
  geom_line(aes(x=date,y=mean,color=site,group=site),linewidth=1)+
  # geom_point(aes(x=date,y=mean,color=site))+
  # geom_errorbar(aes(x=date,ymin=lower,ymax=upper,color=site))+
  scale_color_viridis_d()+
  scale_y_continuous(limits=c(0,.7))
compare_eir <- ggplot(season_true_sim_all)+
  geom_line(aes(x=date,y=EIR_true,color=site,group=site))+
  # geom_point(aes(x=date,y=EIR_true,color=site))+
  # geom_errorbar(aes(x=date,ymin=lower,ymax=upper,color=site))+
  scale_color_viridis_d()


###Set up cluster
obj_sim_seas <- cluster_setup(context_name = 'sim_seas',template='20Core', cores=10)
obj_sim_seas$login()

windows_authenticate()
lapply(season_data_raw_sim_list,function(data) data$t)
season_data_raw_sim_list[[1]]$
sim_seas_short_run <- task_create_expr(parallel::clusterApply(NULL,season_data_raw_sim_list,function(data){
  mamasante::run_pmcmc(data_raw=data,
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = diag(0.5,2),
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 10,
                    n_chains = 1,
                    n_workers = 1,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'fitted')
}),
parallel = hipercow_parallel('parallel',cores_per_process = 8),
resources=resources)
sim_seas_short_result <- task_result(sim_seas_short_run)
sim_seas_short_result[[1]]$run_time/3600
task_info(sim_seas_short_run)
task_log_show(sim_seas_short_run)

sim_seas_short_run_plots <- lapply(1:4,function(i) create_diag_figs(sim_seas_short_result[[i]],
                                                                         country = country_list[[i]],
                                                                         district = admin_list[[i]],
                                                                         folderpath = './sim_seasonal_figs/diag',
                                                                         name = 'sim_seas_short_fitted_update'))
sim_seas_short_informed_props <- lapply(1:5, function(id){
  cov(sim_seas_short_informed_results[[id]]$pars[500:1000,])
})

task_result(id)

sim_seas_short_run$status() #'candied_seabird'
length(season_data_raw_sim_list)

sim_seas_short_run$tasks[[1]]$log()
sim_seas_short_run$tasks[[2]]$log()
obj_sim_seas$unsubmit(sim_seas_short_run$tasks[[3]]$id)
sim_seas_short_run$tasks[[4]]$log()
obj_sim_seas$unsubmit(sim_seas_short_run$tasks[[5]]$id)

tanga_prelim_results <- sim_seas_short_run$tasks[[1]]$result()
tanga_prelim_results$mcmc
mcmc_trace(tanga_prelim_results$mcmc)
drc_prelim_results <- sim_seas_short_run$tasks[[4]]$result()
mcmc_trace(drc_prelim_results$mcmc)

cov(tanga_prelim_results$pars)
cov(tanga_prelim_results$pars[500:1000,])
cov(drc_prelim_results$pars[500:1000,])

annual_prev <- season_data_raw_sim_all%>%
  mutate(year = year(date))%>%
  group_by(site,year)%>%
  summarise(annual_prev = sum(positive)/sum(tested))

sim_seas_short_informed_run <- obj_sim_seas$enqueue_bulk(1:5,function(x,data){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = matrix(c(5,0.001,0.001,0.1),nrow=2),
                    target_prev = 0.4,
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 10,
                    n_chains = 1,
                    n_workers = 1,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed')
},data=season_data_raw_sim_list)
sim_seas_short_informed_run$status() #uncapitalistic_malamute
sim_seas_short_informed_run <- obj_sim_seas$task_bundle_get('uncapitalistic_malamute')
sim_seas_short_informed_run$tasks[[1]]$log()
sim_seas_short_informed_run$tasks[[4]]$log()

ghana_results <- sim_seas_short_informed_run$tasks[[5]]$result()
mcmc_trace(ghana_results$mcmc)
names(ghana_results$mcmc)
cov(ghana_results$pars)
drc_results <- sim_seas_short_informed_run$tasks[[4]]$result()
mcmc_trace(drc_results$mcmc)
cov(drc_results$pars)

sim_seas_short_informed_results <- lapply(1:5, function(id){
  sim_seas_short_informed_run$tasks[[id]]$result()
})
sim_seas_short_informed_plots <- lapply(1:4,function(i) create_diag_figs(sim_seas_short_informed_results[[i]],
                                                                         country = country_list[[i]],
                                                                         district = admin_list[[i]],
                                                                         folderpath = './sim_seasonal_figs/diag',
                                                                         name = 'sim_seas_short_informed_up'))
sim_seas_short_informed_props <- lapply(1:5, function(id){
  cov(sim_seas_short_informed_results[[id]]$pars[500:1000,])
})

sim_seas_short_informed_prepped <- lapply(1:5, function(x) {
  prep_results(results=sim_seas_short_informed_results[[x]],
               sim_data=season_truemonthly_sim_list[[x]],
               burnin=0.5,
               country=country_list[[x]],
               district=admin_list[[x]])
})
sim_seas_short_informed_summary <- bind_rows(lapply(1:5, function(x){
  sim_seas_short_informed_prepped[[x]]$summary
}))
sim_seas_short_informed_sample <- bind_rows(lapply(1:5, function(x){
  sim_seas_short_informed_prepped[[x]]$sample
}))
measure_levels <- c('prev_05','incall','EIR','betaa')
ggplot()+
  geom_line(data=sim_seas_short_informed_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
  geom_line(data = sim_seas_short_informed_summary, aes(x=date,y=median),color='darkgrey',linewidth=1)+
  geom_point(data=season_truemonthly_sim_all,aes(x=date,y=value))+
  facet_wrap(country~factor(measure,levels=measure_levels),scales = 'free_y',nrow=5,ncol=4)+
  scale_y_continuous(limits = c(0,NA))
unique(season_truemonthly_sim_all$measure)
unique(season_truemonthly_sim_all$country)

first_annual_prev <- season_data_raw_sim_all%>%
  mutate(year = year(date))%>%
  group_by(country,year)%>%
  summarise(annual_prev = sum(positive)/sum(tested))%>%
  filter(year==2023)%>%
  select(country,annual_prev)
targetprevs <- first_annual_prev$annual_prev
names(targetprevs) <- first_annual_prev$country
country_list
sim_seas_short_informed_run_2 <- obj_sim_seas$enqueue_bulk(1:5,function(x,data,annual_prev,country){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = matrix(c(5,0.001,0.001,0.1),nrow=2),
                    target_prev = annual_prev[[country[[x]]]],
                    target_prev_group='u5',
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 10,
                    n_chains = 1,
                    n_workers = 1,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed',
                    check_flexibility=TRUE
  )
},data=season_data_raw_sim_list,annual_prev=targetprevs,country=country_list)
sim_seas_short_informed_run_2$status() #plump_skimmer
sim_seas_short_informed_run_2$tasks[[5]]$log()
country_list[5]
sim_seas_short_informed_run_2 <- obj_sim_seas$task_bundle_get('plump_skimmer')

sim_seas_short_2_informed_results <- lapply(1:5, function(id){
  sim_seas_short_informed_run_2$tasks[[id]]$result()
})
sim_seas_short_2_informed_plots <- lapply(1:5,function(i) create_diag_figs(sim_seas_short_2_informed_results[[i]],
                                                                           country = country_list[[i]],
                                                                           district = admin_list[[i]],
                                                                           folderpath = './sim_seasonal_figs/diag',
                                                                           name = 'sim_seas_short_2_informed'))
sim_seas_short_2_informed_props <- lapply(1:5, function(id){
  cov(sim_seas_short_2_informed_results[[id]]$pars[500:1000,])
})

sim_seas_short_2_informed_prepped <- lapply(1:5, function(x) {
  prep_results(results=sim_seas_short_2_informed_results[[x]],
               sim_data=season_truemonthly_sim_list[[x]],
               burnin=0.5,
               country=country_list[[x]],
               district=admin_list[[x]])
})
sim_seas_short_2_informed_summary <- bind_rows(lapply(1:5, function(x){
  sim_seas_short_2_informed_prepped[[x]]$summary
}))
sim_seas_short_2_informed_sample <- bind_rows(lapply(1:5, function(x){
  sim_seas_short_2_informed_prepped[[x]]$sample
}))
measure_levels <- c('prev_05','incall','EIR','betaa')
ggplot()+
  geom_line(data=sim_seas_short_2_informed_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
  geom_line(data = sim_seas_short_2_informed_summary, aes(x=date,y=median),color='darkgrey',linewidth=1)+
  geom_point(data=season_truemonthly_sim_all,aes(x=date,y=value))+
  facet_wrap(country~factor(measure,levels=measure_levels),scales = 'free_y',nrow=5,ncol=4)+
  scale_y_continuous(limits = c(0,NA))

sim_seas_short_informed_run_3 <- obj_sim_seas$enqueue_bulk(1:5,function(x,data,annual_prev,country){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = matrix(1),
                    target_prev = annual_prev[[country[[x]]]],
                    target_prev_group='u5',
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 10,
                    n_chains = 1,
                    n_workers = 1,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed',
                    check_flexibility=TRUE
  )
},data=season_data_raw_sim_list,annual_prev=targetprevs,country=country_list)
sim_seas_short_informed_run_3$status() #chivalrous_longhornbeetle
sim_seas_short_informed_run_3$tasks[[4]]$log()

windows_authenticate()
hipercow_environment_create(name ='sim_seas',
                            globals = c('season_data_raw_sim_list','targetprevs','country_list','x'))
hipercow_environment_show(name='sim_seas')
x <- 100
str(season_data_raw_sim_list)
names(season_data_raw_sim_list) <- unlist(country_list)
sim_seas_short_informed_run_3 <- task_create_expr({
  parallel::clusterApply(cl=NULL,season_data_raw_sim_list,function(x){
    first_annual_prev <- sum(x[1:12,'positive'])/sum(x[1:12,'tested'])
    mamasante::run_pmcmc(data_raw=x,
                         init_EIR = 100,
                         n_particles=200,
                         proposal_matrix = matrix(1),
                         target_prev = first_annual_prev,
                         target_prev_group='u5',
                         max_param=125,
                         prop_treated = 0.4,
                         n_steps = 1000,
                         n_threads = 10,
                         n_chains = 1,
                         n_workers = 1,
                         state_check = 0,## Run equilibrium checks
                         seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                         seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                         seed = 1L,
                         start_pf_time = 30*12,
                         particle_tune = FALSE,
                         comparison = 'u5',
                         initial = 'informed',
                         check_flexibility=TRUE)})},
  parallel = hipercow_parallel('parallel',cores_per_process = 8),
resources=resources)
# task_cancel(sim_seas_short_informed_run_3)
id_guess<- '0cdf45e47788b1f0a963e8c56953de53'

id_guess_3<- 'c50382768a90584a4bfbad4046453872'
task_status( '0cdf45e47788b1f0a963e8c56953de53')
sim_seas_short_informed_run_3 <- id_guess_3
task_status(sim_seas_short_informed_run_3)
task_info(sim_seas_short_informed_run_3)
task_log_show(sim_seas_short_informed_run_3)
task_log_watch(sim_seas_short_informed_run_3)
task_info(id_guess)

task_info(id_guess_3)
result_guess <- task_result(id_guess)
result_guess_2 <- task_result(id_guess)
result_guess_3 <- task_result(id_guess)

sim_seas_short_3_informed_results <- task_result(sim_seas_short_informed_run_3)
sim_seas_short_3_informed_plots <- lapply(1:4,function(i) create_diag_figs(sim_seas_short_3_informed_results[[i]],
                                                                           country = country_list[[i]],
                                                                           district = admin_list[[i]],
                                                                           folderpath = './sim_seasonal_figs/diag',
                                                                           name = 'sim_seas_short_3_informed_updated'))

sim_seas_short_3_informed_props <- lapply(1:4, function(id){
  var(sim_seas_short_3_informed_results[[id]]$pars[500:1000,])
})

sim_seas_short_3_informed_prepped <- lapply(1:4, function(x) {
  prep_results(results=sim_seas_short_3_informed_results[[x]],
               sim_data=season_truemonthly_sim_list[[x]],
               burnin=0.5,
               country=country_list[[x]],
               district=admin_list[[x]])
})
sim_seas_short_3_informed_summary <- bind_rows(lapply(1:4, function(x){
  sim_seas_short_3_informed_prepped[[x]]$summary
}))
sim_seas_short_3_informed_sample <- bind_rows(lapply(1:4, function(x){
  sim_seas_short_3_informed_prepped[[x]]$sample
}))
measure_levels <- c('prev_05','incall','EIR','betaa')
ggplot()+
  geom_line(data=sim_seas_short_3_informed_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
  geom_line(data = sim_seas_short_3_informed_summary, aes(x=date,y=median),color='darkgrey',linewidth=1)+
  geom_point(data=season_truemonthly_sim_all,aes(x=date,y=value))+
  facet_wrap(country~factor(measure,levels=measure_levels),scales = 'free_y',nrow=5,ncol=4)+
  scale_y_continuous(limits = c(0,NA))

sim_short_hc_prev <- ggplot()+
  geom_line(data=sim_seas_short_3_informed_sample[sim_seas_short_3_informed_sample$measure=='prev_05',],aes(x=date,y=value,group=variable),color="#FDBF6F",alpha=0.1,linewidth=0.2)+
  geom_line(data = sim_seas_short_3_informed_summary[sim_seas_short_3_informed_summary$measure=='prev_05',], aes(x=date,y=median),color='#FF7F00',linewidth=0.8)+
  geom_line(data=season_truemonthly_sim_all_8y[season_truemonthly_sim_all_8y$measure=='prev_05',],aes(x=date,y=value),size=0.5)+
  facet_wrap(country~.,nrow=5,ncol=1)+
  scale_y_continuous(limits = c(0,1), expand = c(0, 0))+
  labs(y='Prevalence <5yo')+
  theme(axis.title.x = element_blank())+
  coord_cartesian(ylim = c(0,NA))

sim_short_hc_inc <- ggplot()+
  geom_line(data=sim_seas_short_3_informed_sample[sim_seas_short_3_informed_sample$measure=='inc05',],aes(x=date,y=value*1000,group=variable),color="#B2DF8A",alpha=0.1,linewidth=0.2)+
  geom_line(data = sim_seas_short_3_informed_summary[sim_seas_short_3_informed_summary$measure=='inc05',], aes(x=date,y=median*1000),color="#33A02C",linewidth=0.8)+
  geom_line(data=season_truemonthly_sim_all_8y[season_truemonthly_sim_all_8y$measure=='inc05',],aes(x=date,y=value*1000),size=0.5)+
  facet_wrap(country~.,nrow=5,ncol=1)+
  # scale_y_continuous(limits = c(0,1))+
  scale_y_continuous(expand = c(0, 0))+
  labs(y='Clinical incidence (under 5 years) per 1000 people per day')+
  theme(axis.title.x = element_blank())+
  coord_cartesian(ylim = c(0,20))
sim_short_hc_betaa <- ggplot()+
  geom_line(data=sim_seas_short_3_informed_sample[sim_seas_short_3_informed_sample$measure=='betaa',],aes(x=date,y=value,group=variable),color="#FB9A99",alpha=0.1,linewidth=0.2)+
  geom_line(data = sim_seas_short_3_informed_summary[sim_seas_short_3_informed_summary$measure=='betaa',], aes(x=date,y=median),color="#E31A1C",linewidth=0.8)+
  geom_line(data=season_truemonthly_sim_all_8y[season_truemonthly_sim_all_8y$measure=='betaa',],aes(x=date,y=value),size=0.5)+
  facet_wrap(country~.,nrow=5,ncol=1)+
  # scale_y_continuous(limits = c(0,1))+
  labs(y='Mosquito emergence per day')+
  # scale_y_log10()+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0, 0))+
  coord_cartesian(ylim = c(0,30))
sim_short_hc_plots <- sim_short_hc_prev + sim_short_hc_inc +sim_short_hc_betaa + plot_layout(ncol=3)
ggsave('./sim_seasonal_figs/output/sim_short_hc_plots.tiff',plot=sim_short_hc_plots,width=7,height=5,units='in')


obj_sim_seas$login()
sim_seas_long_informed_run <- obj_sim_seas$enqueue_bulk(1:5,function(x,data,props){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = props[[x]],
                    target_prev = 0.4,
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 10000,
                    n_threads = 10,
                    n_chains = 1,
                    n_workers = 1,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed')
},data=season_data_raw_sim_list,props=sim_seas_short_informed_props)
sim_seas_long_informed_run$status() #'lavender_mountainlion'
sim_seas_long_informed_run <- obj_sim_seas$task_bundle_get('lavender_mountainlion')
sim_seas_long_informed_run$status()

sim_seas_long_informed_results <- lapply(1:5, function(id){
  sim_seas_long_informed_run$tasks[[id]]$result()
})
sim_seas_long_informed_plots <- lapply(1:5,function(i) create_diag_figs(sim_seas_long_informed_results[[i]],
                                                                        country = country_list[[i]],
                                                                        district = admin_list[[i]],
                                                                        folderpath = './sim_seasonal_figs/diag',
                                                                        name = 'sim_seas_long_informed_burnin0.1',
                                                                        burnin=0.1))
sim_seas_long_informed_prepped <- lapply(1:5, function(x) {
  prep_results(results=sim_seas_long_informed_results[[x]],
               sim_data=season_truemonthly_sim_list[[x]],
               burnin=0.1,
               country=country_list[[x]],
               district=admin_list[[x]])
})
sim_seas_long_informed_summary <- bind_rows(lapply(1:5, function(x){
  sim_seas_long_informed_prepped[[x]]$summary
}))
sim_seas_long_informed_sample <- bind_rows(lapply(1:5, function(x){
  sim_seas_long_informed_prepped[[x]]$sample
}))
measure_levels <- c('prev_05','incall','EIR','betaa')
ggplot()+
  geom_line(data=sim_seas_long_informed_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
  geom_line(data = sim_seas_long_informed_summary, aes(x=date,y=median),color='darkgrey',linewidth=1)+
  geom_point(data=season_truemonthly_sim_all,aes(x=date,y=value))+
  facet_wrap(country~factor(measure,levels=measure_levels),scales = 'free_y',nrow=5,ncol=4)+
  scale_y_continuous(limits = c(0,NA))

sim_seas_long_2_informed_run <- obj_sim_seas$enqueue_bulk(1:5,function(x,data,props,annual_prev,country){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = props[[x]],
                    target_prev = annual_prev[[country[[x]]]],
                    target_prev_group='u5',
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 10000,
                    n_threads = 10,
                    n_chains = 1,
                    n_workers = 1,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed',
                    check_flexibility=TRUE)
},data=season_data_raw_sim_list,props=sim_seas_short_2_informed_props,annual_prev=targetprevs,country=country_list)
sim_seas_long_2_informed_run$status() #'diving_americangoldfinch'
sim_seas_long_2_informed_run$tasks[[4]]$log()
sim_seas_long_2_informed_run$tasks[[2]] <- sim_seas_long_2_informed_run_2$tasks[[1]]

sim_seas_long_2_informed_run <- obj_sim_seas$task_bundle_get('diving_americangoldfinch')


sim_seas_long_2_informed_run_2 <- obj_sim_seas$enqueue_bulk(2,function(x,data,props,annual_prev,country){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = props[[x]],
                    target_prev = annual_prev[[country[[x]]]],
                    target_prev_group='u5',
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 10000,
                    n_threads = 10,
                    n_chains = 1,
                    n_workers = 1,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed',
                    check_flexibility=TRUE)
},data=season_data_raw_sim_list,props=sim_seas_short_2_informed_props,annual_prev=targetprevs,country=country_list)
sim_seas_long_2_informed_run_2$status() #'unpronounceable_amurratsnake'
sim_seas_long_2_informed_run_2 <- obj_sim_seas$task_bundle_get('unpronounceable_amurratsnake')

sim_seas_long_2_informed_results <- lapply(1:5, function(id){
  sim_seas_long_2_informed_run$tasks[[id]]$result()
})
sim_seas_long_2_informed_plots <- lapply(1:5,function(i) create_diag_figs(sim_seas_long_2_informed_results[[i]],
                                                                          country = country_list[[i]],
                                                                          district = admin_list[[i]],
                                                                          folderpath = './sim_seasonal_figs/diag',
                                                                          name = 'sim_seas_long_2_informed'))
sim_seas_long_2_informed_props <- lapply(1:5, function(id){
  cov(sim_seas_long_2_informed_results[[id]]$pars[1000:10000,])
})

sim_seas_long_2_informed_prepped <- lapply(1:5, function(x) {
  prep_results(results=sim_seas_long_2_informed_results[[x]],
               sim_data=season_truemonthly_sim_list[[x]],
               burnin=0.5,
               country=country_list[[x]],
               district=admin_list[[x]])
})
sim_seas_long_2_informed_summary <- bind_rows(lapply(1:5, function(x){
  sim_seas_long_2_informed_prepped[[x]]$summary
}))
sim_seas_long_2_informed_sample <- bind_rows(lapply(1:5, function(x){
  sim_seas_long_2_informed_prepped[[x]]$sample
}))
measure_levels <- c('prev_05','incall','EIR','betaa')
ggplot()+
  geom_line(data=sim_seas_long_2_informed_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
  geom_line(data = sim_seas_long_2_informed_summary, aes(x=date,y=median),color='darkgrey',linewidth=1)+
  geom_point(data=season_truemonthly_sim_all,aes(x=date,y=value))+
  facet_wrap(country~factor(measure,levels=measure_levels),scales = 'free_y',nrow=5,ncol=4)+
  scale_y_continuous(limits = c(0,NA))


sim_seas_short_parttune_run <- obj_sim_seas$enqueue_bulk(1:5,function(x,data){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = matrix(c(5,0.001,0.001,0.1),nrow=2),
                    target_prev = 0.4,
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 10,
                    n_chains = 1,
                    n_workers = 1,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*12,
                    particle_tune = TRUE,
                    comparison = 'u5',
                    initial = 'informed')
},data=season_data_raw_sim_list)
sim_seas_short_parttune_run$status() #irongray_ballpython
sim_seas_short_parttune_run$tasks[[1]]$log()


##Mis-specification examples
tanga_sim <- gen_seasonal_sim(init_EIR=50,
                              max_param=125,
                              model_file= "init/odin_model_stripped_seasonal.R",
                              country = 'Tanzania',
                              admin_unit = 'Tanga',
                              sim_length = 8)
tanga_sim$true_data$country <- 'Tanzania'
tanga_sim$true_data$admin <- 'Tanga'
tanga_sim$true_data$site <- 'Tanga, Tanzania'
tanga_sim$data_raw$country <- 'Tanzania'
tanga_sim$data_raw$admin <- 'Tanga'
tanga_sim$data_raw$site <- 'Tanga, Tanzania'
tanga_sim$true_monthly <- daily2monthly(out_df=tanga_sim$true_data,sim_length=8)


fatick_sen_sim <- gen_seasonal_sim(init_EIR=50,
                                   max_param=125,
                                   model_file= "init/odin_model_stripped_seasonal.R",
                                   country = 'Senegal',
                                   admin_unit = 'Fatick',
                                   sim_length = 8)
fatick_sen_sim$true_data$country <- 'Senegal'
fatick_sen_sim$true_data$admin <- 'Fatick'
fatick_sen_sim$true_data$site <- 'Fatick, Senegal'
fatick_sen_sim$data_raw$country <- 'Senegal'
fatick_sen_sim$data_raw$admin <- 'Fatick'
fatick_sen_sim$data_raw$site <- 'Fatick, Senegal'
fatick_sen_sim$true_monthly <- daily2monthly(out_df=fatick_sen_sim$true_data,sim_length=8)

equateur_drc_sim <- gen_seasonal_sim(init_EIR=50,
                                     max_param=125,
                                     model_file= "init/odin_model_stripped_seasonal.R",
                                     country = 'Democratic Republic of Congo',
                                     admin_unit = 'Equateur',
                                     sim_length = 8)
equateur_drc_sim$true_data$country <- 'Democratic Republic of Congo'
equateur_drc_sim$true_data$admin <- 'Equateur'
equateur_drc_sim$true_data$site <- 'Equateur, DRC'
equateur_drc_sim$data_raw$country <- 'Democratic Republic of Congo'
equateur_drc_sim$data_raw$admin <- 'Equateur'
equateur_drc_sim$data_raw$site <- 'Equateur, DRC'
equateur_drc_sim$true_monthly <- daily2monthly(out_df=equateur_drc_sim$true_data,sim_length=8)

ue_ghana_sim <- gen_seasonal_sim(init_EIR=50,
                              max_param=125,
                              model_file= "init/odin_model_stripped_seasonal.R",
                              country = 'Ghana',
                              admin_unit = 'Upper East',
                              sim_length = 8)
ue_ghana_sim$true_data$country <- 'Ghana'
ue_ghana_sim$true_data$admin <- 'Upper East'
ue_ghana_sim$true_data$site <- 'Upper East, Ghana'
ue_ghana_sim$data_raw$country <- 'Ghana'
ue_ghana_sim$data_raw$admin <- 'Upper East'
ue_ghana_sim$data_raw$site <- 'Upper East, Ghana'
ue_ghana_sim$true_monthly <- daily2monthly(out_df=ue_ghana_sim$true_data,sim_length=8)


season_true_sim_all_8y <- bind_rows(tanga_sim$true_data,
                                    fatick_sen_sim$true_data,
                                    equateur_drc_sim$true_data,
                                    ue_ghana_sim$true_data)
season_truemonthly_sim_all_8y <- bind_rows(tanga_sim$true_monthly,
                                           fatick_sen_sim$true_monthly,
                                           equateur_drc_sim$true_monthly,
                                           ue_ghana_sim$true_monthly)%>%
  rename(prev_05=prev05_true,
         incall=inc_all_true,
         inc05=inc05_true,
         EIR=EIR_true,
         betaa=betaa_true)%>%
  select(t,date,country,admin,site,month,prev_05,incall,inc05,EIR,betaa)%>%
  melt(id=c('t','date','country','admin','site','month'))%>%
  rename(measure=variable)%>%
  mutate(date=as.Date(date))
table(season_truemonthly_sim_all_8y$country)
truth_summary <- bind_rows(tanga_sim$true_monthly,
                           fatick_sen_sim$true_monthly,
                           equateur_drc_sim$true_monthly,
                           ue_ghana_sim$true_monthly)%>%
  rename(prev_05=prev05_true,
         incall=inc_all_true,
         inc05=inc05_true,
         EIR=EIR_true,
         betaa=betaa_true)%>%
  select(t,date,country,admin,site,month,prev_05,incall,inc05,EIR,betaa) %>%
  group_by(country,admin,site)%>%
  summarise(mean_prev_05=mean(prev_05),
            mean_incall = mean(incall),
            mean_inc05 = mean(inc05),
            mean_EIR = mean(EIR),
            mean_betaa = mean(betaa))

season_truemonthly_sim_list_8y <- list(tanga_sim$true_monthly,
                                       fatick_sen_sim$true_monthly,
                                       equateur_drc_sim$true_monthly,
                                       ue_ghana_sim$true_monthly)
season_data_raw_sim_all_8y <- bind_rows(tanga_sim$data_raw,
                                        fatick_sen_sim$data_raw,
                                        equateur_drc_sim$data_raw,
                                        ue_ghana_sim$data_raw)
season_data_raw_sim_list_8y <- list(tanga_sim$data_raw,
                                    fatick_sen_sim$data_raw,
                                    equateur_drc_sim$data_raw,
                                    ue_ghana_sim$data_raw)
names(season_data_raw_sim_list_8y) <- c('Tanga, TZ','Fatick, SN','Equateur, CD','Upper East, GH')
saveRDS(season_data_raw_sim_all_8y,'./sim_seasonal_data/season_data_raw_sim_all_8y.rds')
saveRDS(season_data_raw_sim_list_8y,'./sim_seasonal_data/season_data_raw_sim_list_8y.rds')
saveRDS(season_truemonthly_sim_list_8y,'./sim_seasonal_data/season_truemonthly_sim_list_8y.rds')
saveRDS(season_truemonthly_sim_all_8y,'./sim_seasonal_data/season_truemonthly_sim_all_8y.rds')
saveRDS(season_true_sim_all_8y,'./sim_seasonal_data/season_true_sim_all_8y.rds')
season_truemonthly_sim_list_8y <- readRDS('./sim_seasonal_data/season_truemonthly_sim_list_8y.rds')
season_truemonthly_sim_all_8y <- readRDS('./sim_seasonal_data/season_truemonthly_sim_all_8y.rds')

season_data_raw_sim_all_8y_cis <- addCIs(season_data_raw_sim_all_8y,season_data_raw_sim_all_8y$positive,season_data_raw_sim_all_8y$tested)

test <- readRDS('./sim_seasonal_data/season_truemonthly_sim_all_8y.rds')
first_annual_prev <- season_data_raw_sim_all_8y%>%
  mutate(year = year(date))%>%
  group_by(country,year)%>%
  summarise(annual_prev = sum(positive)/sum(tested))%>%
  filter(year==2017)%>%
  select(country,annual_prev)
targetprevs <- first_annual_prev$annual_prev
names(targetprevs) <- first_annual_prev$country
targetprevs_toohigh <- targetprevs+targetprevs*.25
targetprevs_toolow <- targetprevs-targetprevs*.25

sim_seas_short_informed_run_toohigh <- obj_sim_seas$enqueue_bulk(1:5,function(x,data,annual_prev,country){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = matrix(0.05),
                    target_prev = annual_prev[[country[[x]]]],
                    target_prev_group='u5',
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 10,
                    n_chains = 1,
                    n_workers = 1,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed',
                    check_flexibility=TRUE
  )
},data=season_data_raw_sim_list_8y,annual_prev=targetprevs_toohigh,country=country_list)
sim_seas_short_informed_run_toohigh$status() #contractile_afghanhound
sim_seas_short_informed_run_toohigh <- obj_sim_seas$task_bundle_get('contractile_afghanhound')
sim_seas_short_informed_run_toohigh$tasks[[1]]$log()
obj_sim_seas$unsubmit(sim_seas_short_informed_run_toohigh$ids)

sim_seas_short_informed_toohigh_results <- lapply(1:5, function(id){
  sim_seas_short_informed_run_toohigh$tasks[[id]]$result()
})
sim_seas_short_informed_toohigh_plots <- lapply(1:5,function(i) create_diag_figs(sim_seas_short_informed_toohigh_results[[i]],
                                                                                country = country_list[[i]],
                                                                                district = admin_list[[i]],
                                                                                folderpath = './sim_seasonal_figs/diag',
                                                                                name = 'sim_seas_short_toohigh'))
sim_seas_short_informed_toohigh_props <- lapply(1:5, function(id){
  var(sim_seas_short_informed_toohigh_results[[id]]$pars[500:1000,])
})

sim_seas_short_informed_toohigh_prepped <- lapply(1:5, function(x) {
  prep_results(results=sim_seas_short_informed_toohigh_results[[x]],
               sim_data=season_truemonthly_sim_list_8y[[x]],
               burnin=0.5,
               country=country_list[[x]],
               district=admin_list[[x]],
               timelength = 12*8)
})
sim_seas_short_informed_toohigh_summary <- bind_rows(lapply(1:5, function(x){
  sim_seas_short_informed_toohigh_prepped[[x]]$summary
}))
sim_seas_short_informed_toohigh_sample <- bind_rows(lapply(1:5, function(x){
  sim_seas_short_informed_toohigh_prepped[[x]]$sample
}))
measure_levels <- c('prev_05','incall','EIR','betaa')
ggplot()+
  geom_line(data=sim_seas_short_informed_toohigh_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
  geom_line(data = sim_seas_short_informed_toohigh_summary, aes(x=date,y=median),color='darkgrey',linewidth=1)+
  geom_point(data=season_truemonthly_sim_all_8y,aes(x=date,y=value))+
  facet_wrap(country~factor(measure,levels=measure_levels),scales = 'free_y',nrow=5,ncol=4)+
  scale_y_continuous(limits = c(0,NA))

sim_seas_short_informed_run_toolow <- obj_sim_seas$enqueue_bulk(1:5,function(x,data,annual_prev,country){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = matrix(0.05),
                    target_prev = annual_prev[[country[[x]]]],
                    target_prev_group='u5',
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 10,
                    n_chains = 1,
                    n_workers = 1,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed',
                    check_flexibility=TRUE
  )
},data=season_data_raw_sim_list_8y,annual_prev=targetprevs_toolow,country=country_list)
sim_seas_short_informed_run_toolow$status() #fulltime_gonolek
sim_seas_short_informed_run_toolow <- obj_sim_seas$task_bundle_get('fulltime_gonolek')
obj_sim_seas$unsubmit(sim_seas_short_informed_run_toolow$ids)

sim_seas_short_informed_toolow_results <- lapply(1:5, function(id){
  sim_seas_short_informed_run_toolow$tasks[[id]]$result()
})
sim_seas_short_informed_toolow_plots <- lapply(1:5,function(i) create_diag_figs(sim_seas_short_informed_toolow_results[[i]],
                                                                                 country = country_list[[i]],
                                                                                 district = admin_list[[i]],
                                                                                 folderpath = './sim_seasonal_figs/diag',
                                                                                 name = 'sim_seas_short_toolow'))
sim_seas_short_informed_toolow_props <- lapply(1:5, function(id){
  var(sim_seas_short_informed_toolow_results[[id]]$pars[500:1000,])
})

sim_seas_short_informed_toolow_prepped <- lapply(1:5, function(x) {
  prep_results(results=sim_seas_short_informed_toolow_results[[x]],
               sim_data=season_truemonthly_sim_list_8y[[x]],
               burnin=0.5,
               country=country_list[[x]],
               district=admin_list[[x]],
               timelength = 12*8)
})
sim_seas_short_informed_toolow_summary <- bind_rows(lapply(1:5, function(x){
  sim_seas_short_informed_toolow_prepped[[x]]$summary
}))
sim_seas_short_informed_toolow_sample <- bind_rows(lapply(1:5, function(x){
  sim_seas_short_informed_toolow_prepped[[x]]$sample
}))
measure_levels <- c('prev_05','incall','EIR','betaa')
ggplot()+
  geom_line(data=sim_seas_short_informed_toolow_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
  geom_line(data = sim_seas_short_informed_toolow_summary, aes(x=date,y=median),color='darkgrey',linewidth=1)+
  geom_point(data=season_truemonthly_sim_all_8y,aes(x=date,y=value))+
  facet_wrap(country~factor(measure,levels=measure_levels),scales = 'free_y',nrow=5,ncol=4)+
  scale_y_continuous(limits = c(0,NA))

sim_seas_short_informed_run_bangon <- obj_sim_seas$enqueue_bulk(1:5,function(x,data,annual_prev,country){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = matrix(0.05),
                    target_prev = annual_prev[[country[[x]]]],
                    target_prev_group='u5',
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 10,
                    n_chains = 1,
                    n_workers = 1,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed',
                    check_flexibility=TRUE
  )
},data=season_data_raw_sim_list_8y,annual_prev=targetprevs,country=country_list)
names_df <- data.frame(name=names(season_data_raw_sim_list_8y))

test <- task_create_expr({
  print('this worked')
},resources = hipercow::hipercow_resources(cores=1))

hipercow_configuration()
hipercow
task_status(test)
task_info(test)
windows_authenticate()
sim_seas_short_informed_run_bangon <- task_create_bulk_expr({
  print(name)
  x <- season_data_raw_sim_list_8y[[name]]
  first_annual_prev <- sum(x[1:12,'positive'])/sum(x[1:12,'tested'])
  mamasante::run_pmcmc(data_raw=x,
                       init_EIR = 100,
                       n_particles=200,
                       proposal_matrix = matrix(0.05),
                       target_prev = first_annual_prev,
                       target_prev_group='u5',
                       max_param=125,
                       prop_treated = 0.4,
                       n_steps = 1000,
                       n_threads = 32,
                       n_chains = 1,
                       n_workers = 1,
                       state_check = 0,## Run equilibrium checks
                       seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                       seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                       seed = 1L,
                       start_pf_time = 30*12,
                       particle_tune = FALSE,
                       comparison = 'u5',
                       initial = 'informed',
                       check_flexibility=TRUE)
},resources=hipercow::hipercow_resources(cores=32),data=names_df)
#dynamic_comet
task_status(sim_seas_short_informed_run_bangon$ids)
task_info(sim_seas_short_informed_run_bangon$ids[4])
task_log_show(sim_seas_short_informed_run_bangon$ids[4])
task_result(sim_seas_short_informed_run_bangon$ids[1])$trace
sim_seas_short_informed_run_bangon$status() #vermivorous_asiaticmouflon
sim_seas_short_informed_run_bangon <- obj_sim_seas$task_bundle_get('vermivorous_asiaticmouflon')
obj_sim_seas$unsubmit(sim_seas_short_informed_run_bangon$ids)
country_list_sim <- c('Tanzania','Senegal','Democratic Republic of Congo','Ghana')
admin_list_sim <- c('Tanga','Fatick','Equateur','Upper East')

sim_seas_short_informed_bangon_results <- lapply(1:4, function(id){
  task_result(sim_seas_short_informed_run_bangon$ids[id])
  })
sim_seas_short_informed_bangon_plots <- lapply(1:4,function(i) create_diag_figs(sim_seas_short_informed_bangon_results[[i]],
                                                                          country = country_list_sim[[i]],
                                                                          district = admin_list_sim[[i]],
                                                                          folderpath = './sim_seasonal_figs/diag',
                                                                          name = 'sim_seas_short_bangon'))
sim_seas_short_informed_bangon_props <- lapply(1:4, function(id){
  var(sim_seas_short_informed_bangon_results[[id]]$pars[500:1000,])
})

sim_seas_short_informed_bangon_prepped <- lapply(1:4, function(x) {
  prep_results(results=sim_seas_short_informed_bangon_results[[x]],
               sim_data=season_truemonthly_sim_list_8y[[x]],
               burnin=0.5,
               country=country_list_sim[[x]],
               district=admin_list_sim[[x]],
               timelength = 12*8)
})
sim_seas_short_informed_bangon_summary <- bind_rows(lapply(1:4, function(x){
  sim_seas_short_informed_bangon_prepped[[x]]$summary
}))
sim_seas_short_informed_bangon_sample <- bind_rows(lapply(1:4, function(x){
  sim_seas_short_informed_bangon_prepped[[x]]$sample
}))

measure_levels <- c('prev_05','inc05','EIR','betaa')
ggplot()+
  geom_line(data=sim_seas_short_informed_bangon_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
  geom_line(data = sim_seas_short_informed_bangon_summary, aes(x=date,y=median),color='darkgrey',linewidth=1)+
  geom_point(data=season_truemonthly_sim_all_8y[season_truemonthly_sim_all_8y$measure!='incall',],aes(x=date,y=value))+
  facet_wrap(country~factor(measure,levels=measure_levels),scales = 'free_y')+
  scale_y_continuous(limits = c(0,NA))
sim_seas_short_informed_bangon_4diff <- left_join(sim_seas_short_informed_bangon_summary,season_truemonthly_sim_all_8y,by=c('date','country','measure'))
sim_seas_short_informed_bangon_4diff$difference <- (sim_seas_short_informed_bangon_4diff$median-sim_seas_short_informed_bangon_4diff$value)
sim_seas_short_informed_bangon_4diff_explore <- sim_seas_short_informed_bangon_4diff %>%
  filter(difference>5)
daily_dates <- seq.Date(from=as.Date(min(sim_seas_short_informed_bangon_4diff$month)),to=as.Date(max(sim_seas_short_informed_bangon_4diff$month),frac=1),by='day')
daily_df <- data.frame(date=daily_dates,
                       month=as.yearmon(daily_dates))
sim_seas_short_informed_bangon_annualcases <- left_join(daily_df,sim_seas_short_informed_bangon_4diff[sim_seas_short_informed_bangon_4diff$measure=='incall',],by = join_by(month))%>%
  mutate(year=year(date.x))%>%
  group_by(year,country)%>%
  summarise(est_incidence = sum(median),
            est_incidence_lower = sum(lower),
            est_incidence_upper = sum(upper),
            true_incidence = sum(value))%>%
  mutate(percent_diff = (est_incidence - true_incidence)/true_incidence,
         percent_diff_lower = (est_incidence_lower - true_incidence)/true_incidence,
         percent_diff_upper = (est_incidence_upper - true_incidence)/true_incidence)
ggplot(sim_seas_short_informed_bangon_annualcases)+
  geom_point(aes(x=true_incidence,y=est_incidence,color=country))+
  # geom_errorbar(aes(x=true_incidence,ymin=est_incidence_lower,ymax=est_incidence_upper),width=0)+
  geom_abline()+
  scale_x_continuous(limits=c(0.6,1),expand=c(0,0))+
  scale_y_continuous(limits=c(0.6,1),expand=c(0,0))

ggplot(sim_seas_short_informed_bangon_annualcases)+
  geom_point(aes(x=year,y=percent_diff,color=country))+
  geom_errorbar(aes(x=year,ymin=percent_diff_lower,ymax=percent_diff_upper),width=0)+
  facet_wrap(.~country)



ggplot(sim_seas_short_informed_bangon_4diff)+
  geom_point(aes(x=value,y=median))+
  geom_abline()+
  facet_wrap(.~factor(measure,levels=measure_levels),scales = 'free')+
  # scale_y_continuous(limits = c(0,NA))+
  labs(x='True Value',y='Estimated Value')
ggplot(sim_seas_short_informed_bangon_4diff)+
  geom_point(aes(x=date,y=difference))+
  facet_wrap(country~factor(measure,levels=measure_levels),scales = 'free_y',nrow=5,ncol=4)+
  # scale_y_continuous(limits = c(0,NA))+
  labs(x='Date',y='% Difference')
ggplot()+
  geom_line(data=sim_seas_short_informed_bangon_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
  geom_line(data = sim_seas_short_informed_bangon_summary, aes(x=date,y=median),color='darkgrey',linewidth=1)+
  geom_point(data=season_truemonthly_sim_all_8y,aes(x=date,y=value))+
  facet_wrap(country~factor(measure,levels=measure_levels),scales = 'free_y',nrow=5,ncol=4)+
  scale_y_continuous(limits = c(0,NA))

colors_three <- c(viridis(3,begin=0.5,end=0.9))
"#21908CFF""#43BF71FF""#BBDF27FF"
"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00""#CAB2D6" "#6A3D9A"
bangon_prev <- ggplot()+
  geom_line(data=sim_seas_short_informed_bangon_sample[sim_seas_short_informed_bangon_sample$measure=='prev_05',],aes(x=date,y=value,group=variable),color="#FDBF6F",alpha=0.1,linewidth=0.2)+
  geom_line(data = sim_seas_short_informed_bangon_summary[sim_seas_short_informed_bangon_summary$measure=='prev_05',], aes(x=date,y=median),color='#FF7F00',linewidth=1)+
  geom_line(data=season_truemonthly_sim_all_8y[season_truemonthly_sim_all_8y$measure=='prev_05',],aes(x=date,y=value),size=0.5)+
  facet_wrap(country~.,nrow=5,ncol=1)+
  scale_y_continuous(limits = c(0,1), expand = c(0, 0))+
  labs(y='Prevalence <5yo')+
  theme(axis.title.x = element_blank())+
  coord_cartesian(ylim = c(0,NA))

bangon_inc <- ggplot()+
  geom_line(data=sim_seas_short_informed_bangon_sample[sim_seas_short_informed_bangon_sample$measure=='inc05',],aes(x=date,y=value*1000,group=variable),color="#B2DF8A",alpha=0.1,linewidth=0.2)+
  geom_line(data = sim_seas_short_informed_bangon_summary[sim_seas_short_informed_bangon_summary$measure=='inc05',], aes(x=date,y=median*1000),color="#33A02C",linewidth=1)+
  geom_line(data=season_truemonthly_sim_all_8y[season_truemonthly_sim_all_8y$measure=='inc05',],aes(x=date,y=value*1000),size=0.5)+
  facet_wrap(country~.,nrow=5,ncol=1)+
  # scale_y_continuous(limits = c(0,1))+
  scale_y_continuous(expand = c(0, 0))+
  labs(y='Clinical incidence (under 5 years) per 1000 people per day')+
  theme(axis.title.x = element_blank())+
  coord_cartesian(ylim = c(0,30))
bangon_betaa <- ggplot()+
  geom_line(data=sim_seas_short_informed_bangon_sample[sim_seas_short_informed_bangon_sample$measure=='betaa',],aes(x=date,y=value,group=variable),color="#FB9A99",alpha=0.1,linewidth=0.2)+
  geom_line(data = sim_seas_short_informed_bangon_summary[sim_seas_short_informed_bangon_summary$measure=='betaa',], aes(x=date,y=median),color="#E31A1C",linewidth=1)+
  geom_line(data=season_truemonthly_sim_all_8y[season_truemonthly_sim_all_8y$measure=='betaa',],aes(x=date,y=value),size=0.5)+
  facet_wrap(country~.,nrow=5,ncol=1)+
  # scale_y_continuous(limits = c(0,1))+
  labs(y='Mosquito emergence per day')+
  # scale_y_log10()+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0, 0))+
  coord_cartesian(ylim = c(0,30))
bangon_plots <- bangon_prev + bangon_inc +bangon_betaa + plot_layout(ncol=3)
ggsave('./sim_seasonal_figs/output/bangon_plots_4.tiff',plot=bangon_plots,width=7,height=5,units='in')

sim_seas_long_informed_run_bangon <- obj_sim_seas$enqueue_bulk(1:5,function(x,data,annual_prev,country,props){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = matrix(props[[x]]),
                    target_prev = annual_prev[[country[[x]]]],
                    target_prev_group='u5',
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 10000,
                    n_threads = 10,
                    n_chains = 2,
                    n_workers = 2,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed',
                    check_flexibility=TRUE
  )
},data=season_data_raw_sim_list_8y,annual_prev=targetprevs,country=country_list,props=sim_seas_short_informed_bangon_props)
sim_seas_long_informed_run_bangon$status() #wearisome_dairycow
sim_seas_long_informed_run_bangon <- obj_sim_seas$task_bundle_get('wearisome_dairycow')
sim_seas_long_informed_run_bangon$tasks[[1]]$log()

create_diag_figs(sim_seas_long_informed_run_bangon$tasks[[4]]$result(),
                 country = country_list[[4]],
                 district = admin_list[[4]],
                 folderpath = './sim_seasonal_figs/diag',
                 name = 'sim_seas_long_bangon')

obj_sim_seas$login()
sim_seas_long_informed_run_bangon_ghana <- obj_sim_seas$enqueue_bulk(5,function(x,data,annual_prev,country,props){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = matrix(props[[x]]),
                    target_prev = annual_prev[[country[[x]]]],
                    target_prev_group='u5',
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 10000,
                    n_threads = 10,
                    n_chains = 2,
                    n_workers = 2,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 2L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed',
                    check_flexibility=TRUE
  )
},data=season_data_raw_sim_list_8y,annual_prev=targetprevs,country=country_list,props=sim_seas_short_informed_bangon_props)
sim_seas_long_informed_run_bangon_ghana$status() #splendid_alaskanmalamute
sim_seas_long_informed_run_bangon_ghana <- obj_sim_seas$task_bundle_get('splendid_alaskanmalamute')
sim_seas_long_informed_run_bangon_ghana$tasks[[1]]$log()
create_diag_figs(sim_seas_long_informed_run_bangon_ghana$tasks[[1]]$result(),
                 country = country_list[[5]],
                 district = admin_list[[5]],
                 folderpath = './sim_seasonal_figs/diag',
                 name = 'sim_seas_long_bangon')
obj_sim_seas$login()
sim_seas_long_informed_run_bangon_senegal <- obj_sim_seas$enqueue_bulk(3,function(x,data,annual_prev,country,props){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = matrix(props[[x]]),
                    target_prev = annual_prev[[country[[x]]]],
                    target_prev_group='u5',
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 10000,
                    n_threads = 10,
                    n_chains = 2,
                    n_workers = 2,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 3L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed',
                    check_flexibility=TRUE
  )
},data=season_data_raw_sim_list_8y,annual_prev=targetprevs,country=country_list,props=sim_seas_short_informed_bangon_props)
sim_seas_long_informed_run_bangon_senegal <- obj_sim_seas$task_bundle_get('egalitarian_meadowlark')
sim_seas_long_informed_run_bangon_senegal$status() #egalitarian_meadowlark
create_diag_figs(sim_seas_long_informed_run_bangon_senegal$tasks[[1]]$result(),
                 country = country_list[[3]],
                 district = admin_list[[3]],
                 folderpath = './sim_seasonal_figs/diag',
                 name = 'sim_seas_long_bangon')

sim_seas_long_informed_bangon_results <- lapply(1:5, function(id){
  sim_seas_long_informed_run_bangon$tasks[[id]]$result()
})
sim_seas_long_informed_bangon_results[[3]] <- sim_seas_long_informed_run_bangon_senegal$tasks[[1]]$result()
sim_seas_long_informed_bangon_results[[5]] <- sim_seas_long_informed_run_bangon_ghana$tasks[[1]]$result()

sim_seas_long_informed_bangon_prepped <- lapply(1:5, function(x) {
  prep_results(results=sim_seas_long_informed_bangon_results[[x]],
               sim_data=season_truemonthly_sim_list_8y[[x]],
               n_chains=2,
               burnin=0.2,
               country=country_list[[x]],
               district=admin_list[[x]],
               timelength = 12*8)
})
sim_seas_long_informed_bangon_summary <- bind_rows(lapply(1:5, function(x){
  sim_seas_long_informed_bangon_prepped[[x]]$summary
}))
sim_seas_long_informed_bangon_sample <- bind_rows(lapply(1:5, function(x){
  sim_seas_long_informed_bangon_prepped[[x]]$sample
}))
sim_seas_long_informed_bangon_point_est <- bind_rows(lapply(1:5, function(x){
  sim_seas_long_informed_bangon_prepped[[x]]$point_est
}))
sim_seas_long_informed_bangon_point_est %>%
  filter(measure=='prev_05')
sim_seas_long_informed_bangon_point_est %>%
  filter(measure=='inc05')
sim_seas_long_informed_bangon_point_est %>%
  filter(measure=='betaa')
volatility_sum <- bind_rows(lapply(1:5, function(x){
  ind <- get_chain_info(n_chains=2,length=nrow(sim_seas_long_informed_bangon_results[[x]]$mcmc),burnin=0.2)
  dist <- sim_seas_long_informed_bangon_results[[x]]$mcmc[ind,'volatility']
  max_post <- which.max(sim_seas_long_informed_bangon_results[[x]]$mcmc[ind,'log_posterior'])
  MAP <- dist[max_post]
  hpd <- HPDinterval(as.mcmc(data.frame(volatility=dist)))
  return(data.frame(country=country_list[[x]],
                    mean=mean(dist),
                    median=median(dist),
                    map=MAP,
                    lower_ci = quantile(dist,0.025),
                    upper_ci = quantile(dist,0.975),
                    lower_hpd = hpd[1],
                    upper_hpd = hpd[2]))
}))

volatility_sum <- bind_rows(lapply(1:5, function(x){
  get_param_sum(results=sim_seas_long_informed_bangon_results[[x]],
                param_name='volatility',
                posterior_name='log_posterior',
                n_chains=2,
                burnin=0.2,
                country=country_list[[x]])
}))
bangon_long_prev <- ggplot()+
  geom_line(data=sim_seas_long_informed_bangon_sample[sim_seas_long_informed_bangon_sample$measure=='prev_05',],aes(x=date,y=value,group=variable),color="#FDBF6F",alpha=0.1,linewidth=0.2)+
  geom_line(data = sim_seas_long_informed_bangon_summary[sim_seas_long_informed_bangon_summary$measure=='prev_05',], aes(x=date,y=median),color='#FF7F00',linewidth=0.8)+
  geom_line(data=season_truemonthly_sim_all_8y[season_truemonthly_sim_all_8y$measure=='prev_05',],aes(x=date,y=value),size=0.5)+
  facet_wrap(country~.,nrow=5,ncol=1)+
  scale_y_continuous(limits = c(0,1), expand = c(0, 0))+
  labs(y='Prevalence <5yo')+
  theme(axis.title.x = element_blank())+
  coord_cartesian(ylim = c(0,NA))

bangon_long_inc <- ggplot()+
  geom_line(data=sim_seas_long_informed_bangon_sample[sim_seas_long_informed_bangon_sample$measure=='incall',],aes(x=date,y=value*1000,group=variable),color="#B2DF8A",alpha=0.1,linewidth=0.2)+
  geom_line(data = sim_seas_long_informed_bangon_summary[sim_seas_long_informed_bangon_summary$measure=='incall',], aes(x=date,y=median*1000),color="#33A02C",linewidth=0.8)+
  geom_line(data=season_truemonthly_sim_all_8y[season_truemonthly_sim_all_8y$measure=='incall',],aes(x=date,y=value*1000),size=0.5)+
  facet_wrap(country~.,nrow=5,ncol=1)+
  # scale_y_continuous(limits = c(0,1))+
  scale_y_continuous(expand = c(0, 0))+
  labs(y='Clinical incidence (all ages) per 1000 people per day')+
  theme(axis.title.x = element_blank())+
  coord_cartesian(ylim = c(0,20))
bangon_long_betaa <- ggplot()+
  geom_line(data=sim_seas_long_informed_bangon_sample[sim_seas_long_informed_bangon_sample$measure=='betaa',],aes(x=date,y=value,group=variable),color="#FB9A99",alpha=0.1,linewidth=0.2)+
  geom_line(data = sim_seas_long_informed_bangon_summary[sim_seas_long_informed_bangon_summary$measure=='betaa',], aes(x=date,y=median),color="#E31A1C",linewidth=0.8)+
  geom_line(data=season_truemonthly_sim_all_8y[season_truemonthly_sim_all_8y$measure=='betaa',],aes(x=date,y=value),size=0.5)+
  facet_wrap(country~.,nrow=5,ncol=1)+
  # scale_y_continuous(limits = c(0,1))+
  labs(y='Mosquito emergence per day')+
  # scale_y_log10()+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0, 0))+
  coord_cartesian(ylim = c(0,30))
bangon_long_plots <- bangon_long_prev + bangon_long_inc +bangon_long_betaa + plot_layout(ncol=3)
ggsave('./sim_seasonal_figs/output/bangon_long_plots.tiff',plot=bangon_long_plots,width=7,height=5,units='in')

ghana_bangon_prepped <- prep_results(results=sim_seas_long_informed_run_bangon_ghana$tasks[[1]]$result(),
               sim_data=season_data_raw_sim_list_8y[[5]],
               burnin=0.1,
               country='Ghana',
               district=NA,
               timelength = 12*8)

ghana_betaa <- ghana_bangon_prepped$sample%>%
  filter(measure=='betaa')%>%
  group_by(country,variable)%>%
  mutate(last_value=lag(value),
         step=log(value/last_value))%>%
  filter(value!=125)
ggplot(ghana_betaa)+
  geom_histogram(aes(x=step))

bangon_prev <- ggplot()+
  geom_line(data=sim_seas_short_informed_bangon_sample[sim_seas_short_informed_bangon_sample$measure=='prev_05',],aes(x=date,y=value,group=variable),color=colors_three[1],alpha=0.1,linewidth=0.2)+
  geom_line(data = sim_seas_short_informed_bangon_summary[sim_seas_short_informed_bangon_summary$measure=='prev_05',], aes(x=date,y=median),color=colors_three[1],linewidth=1)+
  geom_line(data=season_truemonthly_sim_all_8y[season_truemonthly_sim_all_8y$measure=='prev_05',],aes(x=date,y=value),size=0.5)+
  facet_wrap(country~.,nrow=5,ncol=1)+
  scale_y_continuous(limits = c(0,1), expand = c(0, 0))+
  labs(y='Prevalence <5yo')+
  theme(axis.title.x = element_blank())+
  coord_cartesian(ylim = c(0,NA))

bangon_inc <- ggplot()+
  geom_line(data=sim_seas_short_informed_bangon_sample[sim_seas_short_informed_bangon_sample$measure=='incall',],aes(x=date,y=value*1000,group=variable),color=colors_three[2],alpha=0.1,linewidth=0.2)+
  geom_line(data = sim_seas_short_informed_bangon_summary[sim_seas_short_informed_bangon_summary$measure=='incall',], aes(x=date,y=median*1000),color=colors_three[2],linewidth=1)+
  geom_line(data=season_truemonthly_sim_all_8y[season_truemonthly_sim_all_8y$measure=='incall',],aes(x=date,y=value*1000),size=0.5)+
  facet_wrap(country~.,nrow=5,ncol=1)+
  # scale_y_continuous(limits = c(0,1))+
  scale_y_continuous(expand = c(0, 0))+
  labs(y='Clinical incidence (all ages) per 1000 people per day')+
  theme(axis.title.x = element_blank())+
  coord_cartesian(ylim = c(0,NA))
bangon_betaa <- ggplot()+
  geom_line(data=sim_seas_short_informed_bangon_sample[sim_seas_short_informed_bangon_sample$measure=='betaa',],aes(x=date,y=value,group=variable),color=colors_three[3],alpha=0.1,linewidth=0.2)+
  geom_line(data = sim_seas_short_informed_bangon_summary[sim_seas_short_informed_bangon_summary$measure=='betaa',], aes(x=date,y=median),color=colors_three[3],linewidth=1)+
  geom_line(data=season_truemonthly_sim_all_8y[season_truemonthly_sim_all_8y$measure=='betaa',],aes(x=date,y=value),size=0.5)+
  facet_wrap(country~.,nrow=5,ncol=1)+
  # scale_y_continuous(limits = c(0,1))+
  labs(y='Mosquito emergence per day')+
  # scale_y_log10()+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0, 0))+
  coord_cartesian(ylim = c(0,40))

#####Variety of EIR levels
ue_gambia_sim_highburden <- gen_seasonal_sim(init_EIR=200,
                                             max_param=125,
                                             model_file= "init/odin_model_stripped_seasonal.R",
                                             country = 'Gambia',
                                             admin_unit = 'Upper East',
                                             sim_length = 8)
ue_gambia_sim_highburden$true_data$country <- 'Gambia'
ue_gambia_sim_highburden$true_data$admin <- 'Upper East'
ue_gambia_sim_highburden$true_data$site <- 'Upper East, The Gambia'
ue_gambia_sim_highburden$true_data$burden <- 'High'
ue_gambia_sim_highburden$data_raw$country <- 'Gambia'
ue_gambia_sim_highburden$data_raw$admin <- 'Upper East'
ue_gambia_sim_highburden$data_raw$site <- 'Upper East, The Gambia'
ue_gambia_sim_highburden$data_raw$burden <- 'High'
ue_gambia_sim_highburden$true_monthly <- daily2monthly(out_df=ue_gambia_sim_highburden$true_data,sim_length = 8)

ue_gambia_sim_lowburden <- gen_seasonal_sim(init_EIR=10,
                                            max_param=125,
                                            model_file= "init/odin_model_stripped_seasonal.R",
                                            country = 'Gambia',
                                            admin_unit = 'Upper East',
                                            sim_length = 8)
ue_gambia_sim_lowburden$true_data$country <- 'Gambia'
ue_gambia_sim_lowburden$true_data$admin <- 'Upper East'
ue_gambia_sim_lowburden$true_data$site <- 'Upper East, The Gambia'
ue_gambia_sim_lowburden$true_data$burden <- 'Low'
ue_gambia_sim_lowburden$data_raw$country <- 'Gambia'
ue_gambia_sim_lowburden$data_raw$admin <- 'Upper East'
ue_gambia_sim_lowburden$data_raw$site <- 'Upper East, The Gambia'
ue_gambia_sim_lowburden$data_raw$burden <- 'Low'
ue_gambia_sim_lowburden$true_monthly <- daily2monthly(out_df=ue_gambia_sim_lowburden$true_data,sim_length = 8)
ue_gambia_sim$data_raw$burden <- 'Medium'
ue_gambia_sim$true_data$burden <- 'Medium'
ue_gambia_sim$true_monthly$burden <- 'Medium'

equateur_drc_sim_highburden <- gen_seasonal_sim(init_EIR=200,
                                                max_param=125,
                                                model_file= "init/odin_model_stripped_seasonal.R",
                                                country = 'Democratic Republic of Congo',
                                                admin_unit = 'Equateur',
                                                sim_length = 8)
equateur_drc_sim_highburden$true_data$country <- 'Democratic Republic of Congo'
equateur_drc_sim_highburden$true_data$admin <- 'Equateur'
equateur_drc_sim_highburden$true_data$site <- 'Equateur, DRC'
equateur_drc_sim_highburden$true_data$burden <- 'High'
equateur_drc_sim_highburden$data_raw$country <- 'Democratic Republic of Congo'
equateur_drc_sim_highburden$data_raw$admin <- 'Equateur'
equateur_drc_sim_highburden$data_raw$site <- 'Equateur, DRC'
equateur_drc_sim_highburden$data_raw$burden <- 'High'
equateur_drc_sim_highburden$true_monthly <- daily2monthly(out_df=equateur_drc_sim_highburden$true_data,sim_length = 8)

equateur_drc_sim_lowburden <- gen_seasonal_sim(init_EIR=10,
                                               max_param=125,
                                               model_file= "init/odin_model_stripped_seasonal.R",
                                               country = 'Democratic Republic of Congo',
                                               admin_unit = 'Equateur',
                                               sim_length = 8)
equateur_drc_sim_lowburden$true_data$country <- 'Democratic Republic of Congo'
equateur_drc_sim_lowburden$true_data$admin <- 'Equateur'
equateur_drc_sim_lowburden$true_data$site <- 'Equateur, DRC'
equateur_drc_sim_lowburden$true_data$burden <- 'Low'
equateur_drc_sim_lowburden$data_raw$country <- 'Democratic Republic of Congo'
equateur_drc_sim_lowburden$data_raw$admin <- 'Equateur'
equateur_drc_sim_lowburden$data_raw$site <- 'Equateur, DRC'
equateur_drc_sim_lowburden$data_raw$burden <- 'Low'
equateur_drc_sim_lowburden$true_monthly <- daily2monthly(out_df=equateur_drc_sim_lowburden$true_data,sim_length = 8)
equateur_drc_sim$true_data$burden <- 'Medium'
equateur_drc_sim$data_raw$burden <- 'Medium'
equateur_drc_sim$true_monthly$burden <- 'Medium'

season_burden_true_sim_all <- bind_rows(ue_gambia_sim_lowburden$true_data,
                                        ue_gambia_sim$true_data,
                                        ue_gambia_sim_highburden$true_data,
                                        equateur_drc_sim_lowburden$true_data,
                                        equateur_drc_sim$true_data,
                                        equateur_drc_sim_highburden$true_data,)
season_burden_truemonthly_sim_all <- bind_rows(ue_gambia_sim_lowburden$true_monthly,
                                               ue_gambia_sim$true_monthly,
                                               ue_gambia_sim_highburden$true_monthly,
                                               equateur_drc_sim_lowburden$true_monthly,
                                               equateur_drc_sim$true_monthly,
                                               equateur_drc_sim_highburden$true_monthly,)%>%
  rename(prev_05=prev05_true,
         incall=inc_all_true,
         EIR=EIR_true,
         betaa=betaa_true)%>%
  select(t,date,country,admin,site,month,prev_05,incall,EIR,betaa,burden)%>%
  melt(id=c('t','date','country','admin','site','month','burden'))%>%
  rename(measure=variable)%>%
  mutate(date=as.Date(date))

season_burden_truemonthly_sim_list <- list(`Gambia-Low`=ue_gambia_sim_lowburden$true_monthly,
                                           `Gambia-Medium`=ue_gambia_sim$true_monthly,
                                           `Gambia-High`=ue_gambia_sim_highburden$true_monthly,
                                           `Democratic Republic of Congo-Low`=equateur_drc_sim_lowburden$true_monthly,
                                           `Democratic Republic of Congo-Medium`=equateur_drc_sim$true_monthly,
                                           `Democratic Republic of Congo-High`=equateur_drc_sim_highburden$true_monthly)
site_names <- names(season_burden_truemonthly_sim_list)
season_burden_data_raw_sim_all <- bind_rows(ue_gambia_sim_lowburden$data_raw,
                                            ue_gambia_sim$data_raw,
                                            ue_gambia_sim_highburden$data_raw,
                                            equateur_drc_sim_lowburden$data_raw,
                                            equateur_drc_sim$data_raw,
                                            equateur_drc_sim_highburden$data_raw)
season_burden_data_raw_sim_list <- list(`Gambia-Low`=ue_gambia_sim_lowburden$data_raw,
                                        `Gambia-Medium`=ue_gambia_sim$data_raw,
                                        `Gambia-High`=ue_gambia_sim_highburden$data_raw,
                                        `Democratic Republic of Congo-Low`=equateur_drc_sim_lowburden$data_raw,
                                        `Democratic Republic of Congo-Medium`=equateur_drc_sim$data_raw,
                                        `Democratic Republic of Congo-High`=equateur_drc_sim_highburden$data_raw)
saveRDS(season_burden_data_raw_sim_all,'./sim_seasonal_data/season_burden_data_raw_sim_all.rds')
saveRDS(season_burden_data_raw_sim_list,'./sim_seasonal_data/season_burden_data_raw_sim_list.rds')
saveRDS(season_burden_truemonthly_sim_list,'./sim_seasonal_data/season_burden_truemonthly_sim_list.rds')
saveRDS(season_burden_truemonthly_sim_all,'./sim_seasonal_data/season_burden_truemonthly_sim_all.rds')
saveRDS(season_burden_true_sim_all,'./sim_seasonal_data/season_burden_true_sim_all.rds')

season_burden_data_raw_sim_all <- readRDS('./sim_seasonal_data/season_burden_data_raw_sim_all.rds')
season_burden_data_raw_sim_list <- readRDS('./sim_seasonal_data/season_burden_data_raw_sim_list.rds')
season_burden_truemonthly_sim_list <- readRDS('./sim_seasonal_data/season_burden_truemonthly_sim_list.rds')
season_burden_truemonthly_sim_all <- readRDS('./sim_seasonal_data/season_burden_truemonthly_sim_all.rds')
season_burden_true_sim_all <- readRDS('./sim_seasonal_data/season_burden_true_sim_all.rds')

compare_prev <- ggplot(season_burden_data_raw_sim_all)+
  geom_line(aes(x=date,y=positive/tested,color=site,group=site),linewidth=1)+
  # geom_point(aes(x=date,y=mean,color=site))+
  # geom_errorbar(aes(x=date,ymin=lower,ymax=upper,color=site))+
  scale_color_viridis_d()+
  scale_y_continuous(limits=c(0,1))+
  facet_wrap(.~burden)
compare_eir <- ggplot(season_burden_true_sim_all)+
  geom_line(aes(x=date,y=EIR_true,color=site,group=site))+
  scale_color_viridis_d()+
  facet_wrap(.~burden)

first_annual_prev_burden <- season_burden_data_raw_sim_all%>%
  mutate(year = year(date))%>%
  group_by(country,burden,year)%>%
  summarise(annual_prev = sum(positive)/sum(tested))%>%
  filter(year==2017)%>%
  select(country,burden,annual_prev)
targetprevs_burden <- first_annual_prev_burden$annual_prev
names(targetprevs_burden) <- paste0(first_annual_prev_burden$country,'-',first_annual_prev_burden$burden)

obj_sim_seas$login()
sim_seas_short_informed_run_burden <- obj_sim_seas$enqueue_bulk(1:6,function(x,data,annual_prev,country){
  mamasante::run_pmcmc(data_raw=data[[x]],
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = matrix(0.05),
                    target_prev = annual_prev[[country[[x]]]],
                    target_prev_group='u5',
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 10,
                    n_chains = 1,
                    n_workers = 1,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*12,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed',
                    check_flexibility=TRUE
  )
},data=season_burden_data_raw_sim_list,annual_prev=targetprevs_burden,country=site_names)
sim_seas_short_informed_run_burden <- obj_sim_seas$task_bundle_get('fuchsia_firebelliedtoad')
sim_seas_short_informed_run_burden$status() #fuchsia_firebelliedtoad
sim_seas_short_informed_run_burden$tasks[[1]]$result()$history[1,1,]

sim_seas_short_informed_burden_results <- lapply(1:6, function(id){
  sim_seas_short_informed_run_burden$tasks[[id]]$result()
})
sim_seas_short_informed_burden_plots <- lapply(1:6,function(i) create_diag_figs(sim_seas_short_informed_burden_results[[i]],
                                                                                 country = site_names[[i]],
                                                                                district = 1,
                                                                                 folderpath = './sim_seasonal_figs/diag',
                                                                                 name = 'sim_seas_short_burden'))
sim_seas_short_informed_burden_props <- lapply(1:6, function(id){
  var(sim_seas_short_informed_burden_results[[id]]$pars[500:1000,])
})

sim_seas_short_informed_burden_prepped <- lapply(1:6, function(x) {
  prep_results(results=sim_seas_short_informed_burden_results[[x]],
               sim_data=season_burden_truemonthly_sim_list[[x]],
               burnin=0.5,
               country = site_names[[x]],
               district = NA,
               timelength = 12*8)
})
sim_seas_short_informed_burden_summary <- bind_rows(lapply(1:6, function(x){
  sim_seas_short_informed_burden_prepped[[x]]$summary
}))
sim_seas_short_informed_burden_sample <- bind_rows(lapply(1:6, function(x){
  sim_seas_short_informed_burden_prepped[[x]]$sample
}))
measure_levels <- c('prev_05','incall','EIR','betaa')
country_levels <- unique(season_burden_truemonthly_sim_all$country)
burden_comp <- ggplot()+
  geom_line(data=sim_seas_short_informed_burden_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.1,linewidth=0.2)+
  geom_line(data = sim_seas_short_informed_burden_summary, aes(x=date,y=median),color='darkgrey',linewidth=0.8)+
  geom_line(data=season_burden_truemonthly_sim_all,aes(x=date,y=value),linewidth=0.5)+
  facet_wrap(factor(country,levels=country_levels)~factor(measure,levels=measure_levels),scales = 'free_y',nrow=6,ncol=4)+
  scale_y_continuous(limits = c(0,NA))
ggsave('./sim_seasonal_figs/output/diff_burdens_short.tiff',plot=burden_comp,width=7,height=5,units='in')

table(season_burden_truemonthly_sim_all$country)
table(is.na(season_burden_truemonthly_sim_all$measure))
unique(season_burden_truemonthly_sim_all$country)
