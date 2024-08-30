library(lubridate)
library(dplyr)
library(tidyr)
library(sifter)
library(ggplot2)

theme_set(theme_minimal(base_size = 7)+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))

years <- 2
N <- 12*years
init_EIR <- 2000

init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

prop_treated <- 0.4
het_brackets <- 5

##Convert init_EIR into an init_betaa
mpl_initial <- sifter::model_param_list_create(init_EIR = init_EIR,
                                               init_ft = prop_treated
)

pars_initial <- sifter::equilibrium_init_create_stripped(age_vector = init_age,
                                                         init_EIR = init_EIR,
                                                         ft = prop_treated,
                                                         model_param_list = mpl_initial,
                                                         het_brackets = het_brackets)
pars_initial$init_EIR
################## generate the data ######################
betaa_dates <- seq.Date(from = as.Date('2015-01-01'),by = 'month',length.out = N)
betaa_times <- sapply(betaa_dates,get_t)

### just a random walk on logscale
betaa_vals <- rep(c(50, 300, 300, 300, 4, 9, 22, 30, 56,
                       71, 23, 8),years)
##set up the simulation for the simualted data
time<- max(betaa_times) + 30
out_step=1

mpl <- sifter::model_param_list_create(betaa_times=betaa_times,
                                       betaa_vals=betaa_vals,
                                       lag_rates = 10
)

pars <- sifter::equilibrium_init_create_stripped(age_vector = init_age,
                                                 init_EIR = init_EIR,
                                                 ft = prop_treated,
                                                 model_param_list = mpl,
                                                 het_brackets = het_brackets)

##The malaria model
model_file="MiP_odin_model_nodelay.R"
generator <- odin(model_file)
state_use <- pars[names(pars) %in% coef(generator)$name]

# create model with initial values
mod <- generator(user = pars, use_dde = TRUE)
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

ggplot(monthly_data)+
  geom_line(aes(x=month,y=inc_all_true))+
  coord_cartesian(ylim = c(0,NA))
ggplot(monthly_data)+
  geom_line(aes(x=month,y=prev05))+
  coord_cartesian(ylim = c(0,NA))
ggplot(monthly_data)+
  geom_line(aes(x=month,y=mv))+
  coord_cartesian(ylim = c(0,NA))

four_years <- monthly_data[c((16*12+1):(20*12)),]%>%
  mutate(date_shifted = ymd(date) - years(10),
         month_shifted = month - 10)
all_dates <- four_years %>%
  group_by(month_shifted) %>%
  tidyr::expand(date = seq(floor_date(date_shifted, unit = "month"),
                    ceiling_date(date_shifted, unit="month")-days(1), by="day"))
four_years_byday <- four_years %>%
  select(date_shifted,month_shifted,betaa_true)%>%
  right_join(all_dates,by='month_shifted')

incidence_plot <- ggplot(four_years)+
  geom_line(data=sim_seasonal_hist4plots$inc_sample,aes(x=date,y=value*1000,group=variable),alpha=0.5,color='lightgrey')+
  geom_line(aes(x=date_shifted,y=inc_all_true*1000))+
  geom_line(data=sim_seasonal_hist4plots$inc_history,aes(x=date,y=incall.median*1000),linetype='dashed')+
  labs(x='Year',y='Clinical cases\nper 1000 persons')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
incidence_truth_plot <- ggplot(four_years)+
  geom_line(data=sim_seasonal_hist4plots$inc_sample,aes(x=date,y=value*1000,group=variable),alpha=0,color='lightgrey')+
  geom_line(aes(x=date_shifted,y=inc_all_true*1000))+
  labs(x='Year',y='Clinical cases\nper 1000 persons')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
four_years_cis <- addCIs(four_years,Ys=four_years$positive,Ns=four_years$tested)
prev05_plot <- ggplot(four_years)+
  geom_line(data=sim_seasonal_hist4plots$prev_sample,aes(x=date,y=value,group=variable),alpha=0.5,color='lightgrey')+
  geom_line(aes(x=date_shifted,y=prev05))+
  geom_line(data=sim_seasonal_hist4plots$prev_history,aes(x=date,y=prev.median),linetype='dashed')+
  labs(x='Year',y='Prevalence (<5 years old)')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
prev05_truth_plot <- ggplot(four_years_cis)+
  geom_line(data=sim_seasonal_hist4plots$prev_sample,aes(x=date,y=value,group=variable),alpha=0,color='lightgrey')+
  geom_line(aes(x=date_shifted,y=mean))+
  geom_errorbar(aes(x=date_shifted,ymin=lower,ymax=upper),width=0,alpha=0)+
  labs(x='Year',y='Prevalence (<5 years old)')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
prev05_data_only_plot <- ggplot(four_years_cis)+
  geom_line(data=sim_seasonal_hist4plots$prev_sample,aes(x=date,y=value,group=variable),alpha=0,color='lightgrey')+
  geom_point(aes(x=date_shifted,y=mean))+
  geom_errorbar(aes(x=date_shifted,ymin=lower,ymax=upper),width=0)+
  labs(x='Year',y='Prevalence (<5 years old)')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
prev05_data_plot <- ggplot(four_years_cis)+
  geom_line(data=sim_seasonal_hist4plots$prev_sample,aes(x=date,y=value,group=variable),alpha=0.5,color='lightgrey')+
  geom_line(data=sim_seasonal_hist4plots$prev_history,aes(x=date,y=prev.median),linetype='dashed')+
  geom_point(aes(x=date_shifted,y=mean))+
  geom_errorbar(aes(x=date_shifted,ymin=lower,ymax=upper),width=0)+
  labs(x='Year',y='Prevalence (<5 years old)')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
eir_plot <- ggplot(four_years)+
  geom_line(data=sim_seasonal_hist4plots$EIR_sample,aes(x=date,y=value,group=variable),alpha=0.5,color='lightgrey')+
  geom_line(aes(x=date_shifted,y=EIR_true))+
  geom_line(data=sim_seasonal_hist4plots$EIR_history,aes(x=date,y=EIR.median),linetype='dashed')+
  labs(x='Year',y='Infectious mosquito bites\nper person per year')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
eir_truth_plot <- ggplot(four_years)+
  geom_line(data=sim_seasonal_hist4plots$EIR_sample,aes(x=date,y=value,group=variable),alpha=0,color='lightgrey')+
  geom_line(aes(x=date_shifted,y=EIR_true))+
  labs(x='Year',y='Infectious mosquito bites\nper person per year')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

moz_plot <- ggplot(four_years)+
  geom_line(data=sim_seasonal_hist4plots$mv_sample,aes(x=date,y=value,group=variable),alpha=0.5,color='lightgrey')+
  geom_line(aes(x=date_shifted,y=mv))+
  geom_line(data=sim_seasonal_hist4plots$mv_history,aes(x=date,y=mv.median),linetype='dashed')+
  coord_cartesian(ylim = c(0,NA))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

betaa_plot <- ggplot(four_years_byday)+
  geom_line(data=sim_seasonal_hist4plots$betaa_sample,aes(x=date,y=value,group=variable),alpha=0.5,color='lightgrey')+
  geom_line(aes(x=date,y=betaa_true))+
  geom_line(data=sim_seasonal_hist4plots$betaa_history,aes(x=date,y=betaa.median),linetype='dashed')+
  labs(x='Year',y='# adult mosquitoes emerging\nfrom pupae per person per day')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

betaa_truth_plot <- ggplot(four_years_byday)+
  geom_line(data=sim_seasonal_hist4plots$betaa_sample,aes(x=date,y=value,group=variable),alpha=0,color='lightgrey')+
  geom_line(aes(x=date,y=betaa_true))+
  labs(x='Year',y='# adult mosquitoes emerging\nfrom pupae per person per day')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

windows(15,6.5)
truth <- betaa_truth_plot + eir_truth_plot + incidence_truth_plot + prev05_truth_plot + plot_layout(nrow=1)
observed <- betaa_truth_plot + eir_truth_plot + incidence_truth_plot + prev05_data_only_plot + plot_layout(nrow=1,ncol=4)
estimates <- betaa_plot + eir_plot + incidence_plot + prev05_data_plot + plot_layout(nrow=1)
ggsave('true_sim_demo.tiff',plot = truth,units = 'cm',width = 6,height=12,dpi=300)
ggsave('true_sim_demo_wide.tiff',plot = truth,units = 'cm',width = 17,height=4,dpi=300)
ggsave('estimate_sim_demo_wide.tiff',plot = estimates,units = 'cm',width = 17,height=4,dpi=300)
ggsave('observed_sim_demo_wide.tiff',plot = observed,units = 'cm',width = 17,height=4,dpi=300)
regr_test_df <- cars
str(cars)
str(datasets::iris)
lm(Species~.,data=iris)

names(sim_seasonal_run_result$history[,1,1])

results <- sim_seasonal_run_result
sim_data <- four_years

sim_seasonal_hist4plots <- sim_history4plots(sim_seasonal_run_result,sim_data=four_years,burnin = 0.2)

###What does a simulation based on real ento data look like?
site_list <- c('bf_banfora','bf_gaoua','bf_orodara',
               'mz_changara','mz_chemba','mz_guro',
               'ng_asa','ng_ejigbo','ng_ifenorth','ng_moro')
names(site_list) <- c('Banfora','Gaoua','Orodara','Changara','Chemba','Guro','Asa','Ejigbo','Ife North','Moro')
ento_data <- bind_rows(lapply(1:length(site_list),function(i,site_list){
  df <- read.table(file=paste0('Q:/anc_pmcmc/nnp/data/',site_list[[i]],'_ento.txt'),skip=1,col.names=c('month','value','plus','minus'),sep=',',header=FALSE)%>%
    mutate(month = as.yearmon(month,'%b-%Y'),
           minus = ifelse(value<=0,0,minus),
           plus = ifelse(value<=0,0,plus),
           value = ifelse(value<=0,0,value),
           lower = value-minus,
           upper = value+plus,
           district = names(site_list[i]))
  return(df)
},site_list=site_list))
ento_data_asa <- ento_data %>%
  filter(district=='Asa')

N <- 21
init_EIR <- 500

init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

prop_treated <- 0.4
het_brackets <- 5

##Convert init_EIR into an init_betaa
mpl_initial <- sifter::model_param_list_create(init_EIR = init_EIR,
                                               init_ft = prop_treated
)

pars_initial <- sifter::equilibrium_init_create_stripped(age_vector = init_age,
                                                         init_EIR = init_EIR,
                                                         ft = prop_treated,
                                                         model_param_list = mpl_initial,
                                                         het_brackets = het_brackets)
pars_initial$init_EIR
################## generate the data ######################
betaa_dates <- seq.Date(from = as.Date('2020-11-01'),by = 'month',length.out = N)
betaa_times <- sapply(betaa_dates,get_t)

### just a random walk on logscale
betaa_vals <- ento_data_asa$value
betaa_vals[1] <- 200
##set up the simulation for the simualted data
time<- max(betaa_times) + 30
out_step=1

mpl <- sifter::model_param_list_create(betaa_times=betaa_times,
                                       betaa_vals=betaa_vals,
                                       lag_rates = 10
)

pars <- sifter::equilibrium_init_create_stripped(age_vector = init_age,
                                                 init_EIR = init_EIR,
                                                 ft = prop_treated,
                                                 model_param_list = mpl,
                                                 het_brackets = het_brackets)

##The malaria model
model_file="MiP_odin_model_nodelay.R"
generator <- odin(model_file)
state_use <- pars[names(pars) %in% coef(generator)$name]

# create model with initial values
mod <- generator(user = pars, use_dde = TRUE)
tt <- seq(min(betaa_times), max(betaa_times), out_step)
rm(mod_run)
# run the simulation to base the data
mod_run <- mod$run(tt, step_max_n = 1e7,
                   atol = 1e-5,
                   rtol = 1e-5)

# shape output
out <- mod$transform_variables(mod_run)
out$betaa_out[1]
# plot data and generate data
# plot(out$t,out$betaa_out,type='l',col="red",ylim=c(0,125))
# lines(out$t,out$prev*100,col="blue",lwd=4)
out_df <- data.frame(t=out$t,
                     prev05=out$prev05,
                     date=as.character(seq.Date(from = as.Date('2020-11-01'),by = 'day',length.out = length(tt))),
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

data_raw_ng_pg_asa <- readRDS('Q:/anc_pmcmc/nnp/data/data_raw_ng_pg_asa.RDS')
ggplot(out_df)+
  geom_line(aes(x=as.Date(date),y=prev05),size=1)+
  geom_line(data=data_raw_ng_pg_asa,aes(x=as.Date(month,frac=0.5),y=positive/tested),linetype='dashed',size=1)+
  scale_y_continuous(limits=c(0,1))
ggplot(out_df)+
  geom_line(aes(x=as.Date(date),y=betaa_true))
ggplot(out_df)+
  geom_line(aes(x=as.Date(date),y=EIR_true))
ggplot(out_df)+
  geom_line(aes(x=as.Date(date),y=inc_all_true))

ggplot(out_df)+
  geom_line(aes(x=as.Date(date),y=mv))+
  geom_col(data=ento_data_asa,aes(x=as.Date(month),y=value))
