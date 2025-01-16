#Supplemental Figures

library(lubridate)
library(dplyr)
library(tidyr)
library("colorspace")
library(scales)
library(patchwork)
library(ggplot2)
library(orderly2)
library(odin)
library(zoo)
library(hipercow)
library(RColorBrewer)
library(ggpubr)

devtools::install_github("jt-hicks/mamasante")

source('utils.R')
theme_set(theme_minimal(base_size = 7)+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))

id_data_gen <- orderly2::orderly_run('create_sim_data')

hipercow_init()
hipercow_configure('windows')

hipercow_provision()
hipercow_environment_create(packages=c('dplyr','ggplot2'))
windows_authenticate()
resources <- hipercow_resources(cores=32)

short_seas_id <- task_create_bulk_expr(orderly2::orderly_run('run_pmcmc',parameters=list(name=name,length=1000)),data=data.frame(name=c(1:16)),
                                     resources=resources)
task_status(short_seas_id$ids)
task_log_show(short_seas_id$ids)

short_diag_id <- orderly2::orderly_run('run_diagnostics',parameters=list(length=1000,n_datasets=16))

short_figs_id <- orderly2::orderly_run('run_seasonal_figures')



long_seas_id <- task_create_bulk_expr(orderly2::orderly_run('run_pmcmc',parameters=list(name=name,length=10000,proposal_matrix='from1000')),data=data.frame(name=c(1:16)),
                                      resources=resources)


long_seas_id <- hipercow_bundle_load('ungeometric_bighornedsheep')
task_status(long_seas_id$ids)#'ungeometric_bighornedsheep'
task_log_show(long_seas_id$ids[3])

##Mis-specification
misspecification <- c(-0.2,0.2)
name <- c(5:8) #All four sites at EIR=50
start_pf_time <- c(30,90,180,360)

parameter_df <- expand.grid(misspecification=misspecification,name=name,start_pf_time=start_pf_time)

short_misp_id <- task_create_bulk_expr(orderly2::orderly_run('run_pmcmc',parameters=list(name=name,
                                                                                         length=1000,
                                                                                         start_pf_time=start_pf_time,
                                                                                         misspecification=misspecification)),
                                       data=parameter_df,
                                       resources=resources)
parameter_df_redo <- parameter_df[31,]
parameter_df_redo$seed = 1712
short_misp_id_redo <- task_create_bulk_expr(orderly2::orderly_run('run_pmcmc',parameters=list(name=name,length=1000,start_pf_time=start_pf_time,misspecification=misspecification,seed=seed)),data=parameter_df_redo,
                                       resources=resources)
task_status(short_misp_id_redo$ids) #semicontinuous_duckbillcat

short_misp_diag_id <- orderly2::orderly_run('run_diagnostics_misp',parameters=list(length=1000,n_datasets=32))

task_status(short_misp_id$ids) #apolitical_americanquarterhorse

misspecification <- c(0)
name <- c(5:8) #All four sites at EIR=50
start_pf_time <- c(30,90,180,360)

parameter_control_df <- expand.grid(misspecification=misspecification,name=name,start_pf_time=start_pf_time)

short_misp_control_id <- task_create_bulk_expr(orderly2::orderly_run('run_pmcmc',parameters=list(name=name,
                                                                                         length=1000,
                                                                                         start_pf_time=start_pf_time,
                                                                                         misspecification=misspecification)),
                                       data=parameter_control_df,
                                       resources=resources)
short_misp_control_id<- hipercow::hipercow_bundle_load('semimonarchical_dipper')
task_status(short_misp_control_id$ids)#semimonarchical_dipper

short_misp_control_diag_id <- orderly2::orderly_run('run_diagnostics_misp',parameters=list(length=1000,n_datasets=16))

short_misp_figures_id <- orderly2::orderly_run('run_seasonal_figures_misp')


##Data for vignette
sim_data_true <- sim_data_EIR10_Tanga_Tanzania
sim_data_raw <- sim_seasonal_dataraw_list[[1]]
months <- unique(zoo::as.yearmon(sim_data_true$date))
midmonth_dates <- data.frame(date=as.character(zoo::as.Date(months,frac=0.5)))
sim_data_true$date <- as.character(sim_data_true$date)
monthly_data_true <- dplyr::left_join(midmonth_dates,sim_data_true,by='date')
sim_data_raw$init_EIR <- 10

library('dplyr')
library('zoo')
sim_data <- left_join(monthly_data_true,sim_data_raw,by=c('date','admin','country','site'),suffix=c('.true','.raw'))%>%
  rename(prev_05_true = prev05_true,
         clininc_05_true=inc05_true,
         clininc_all_true=inc_all_true,
         EIR_true=EIR_true,
         betaa_true=betaa_true)%>%
  select(-(c(t.true,t.raw)))%>%
  filter(!is.na(tested))%>%
  mutate(date=as.Date(date),
         month=as.yearmon(date))
saveRDS(sim_data, file='Q:/mamasante/data/sim_data_tanzania.rds')
?saveRDS
