library(lubridate)
library(dplyr)
library(tidyr)
devtools::install_github("jt-hicks/mamasante")

theme_set(theme_minimal(base_size = 7)+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))

years <- 20
N <- 12*years
init_EIR <- 80

init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

prop_treated <- 0.4
het_brackets <- 5

##Convert init_EIR into an init_betaa
mpl_initial <- mamasante::model_param_list_create(init_EIR = init_EIR,
                                               init_ft = prop_treated,
                                               comparison='u5'
)

pars_initial <- mamasante::equilibrium_init_create_stripped(age_vector = init_age,
                                                         init_EIR = init_EIR,
                                                         ft = prop_treated,
                                                         model_param_list = mpl_initial,
                                                         het_brackets = het_brackets)

################## generate the data ######################
betaa_dates <- seq.Date(from = as.Date('2015-01-01'),by = 'month',length.out = N)
betaa_times <- sapply(betaa_dates,get_t)

### just a random walk on logscale
betaa_vals <- rep(c(2, 2, 2, 2, 4, 9, 22, 30, 56,
                       71, 23, 8),years)
##set up the simulation for the simualted data
time<- max(betaa_times) + 30
out_step=1

mpl <- mamasante::model_param_list_create(init_EIR = init_EIR,
                                       init_ft = prop_treated,
                                       betaa_times=betaa_times,
                                       betaa_vals=betaa_vals,
                                       lag_rates = 10,
                                       comparison = 'u5'
)

pars <- mamasante::equilibrium_init_create_stripped(age_vector = init_age,
                                                 init_EIR = init_EIR,
                                                 ft = prop_treated,
                                                 model_param_list = mpl,
                                                 het_brackets = het_brackets)
##The malaria model
model_file<-"MiP_odin_model_nodelay.R"
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
  mutate(date = ymd(date) - years(10),
         month = month - 10)
all_dates <- four_years %>%
  group_by(month) %>%
  tidyr::expand(date = seq(floor_date(date, unit = "month"),
                    ceiling_date(date, unit="month")-days(1), by="day"))
four_years_byday <- four_years %>%
  select(date,month,betaa_true)%>%
  right_join(all_dates,by='month')
sim_seasonal_hist4plots <- sim_history4plots(sim_seasonal_run_result,sim_data=four_years,burnin = 0.2)
four_years_cis <- addCIs(four_years,Ys=four_years$positive,Ns=four_years$tested)

library("colorspace")
test_palette <- viridis::viridis(4,begin=0,end=0.4)
show_col(test_palette)
show_col(c(test_palette,lighten(test_palette,0.3)),ncol=4)
show_col(viridis::mako(12,begin=0,end=1))

true_palette <- viridis::viridis(4,begin=0.1,end=0.35)
part_palette <- lighten(true_palette,0.4)
show_col(c(true_palette,part_palette),ncol=4)
alpha_var <-0.1
betaa_plot <- ggplot(four_years_byday)+
  geom_line(data=sim_seasonal_hist4plots$betaa_sample,aes(x=date,y=value,group=variable),alpha=alpha_var,color=part_palette[1])+
  geom_line(aes(x=date.y,y=betaa_true),color=true_palette[1])+
  geom_line(data=sim_seasonal_hist4plots$betaa_history,aes(x=date,y=betaa.median),linetype='dashed',color=true_palette[1])+
  labs(x='Year',y='# adult mosquitoes emerging\nfrom pupae per person per day')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

betaa_truth_plot <- ggplot(four_years_byday)+
  geom_line(data=sim_seasonal_hist4plots$betaa_sample,aes(x=date,y=value,group=variable),alpha=0,color=part_palette[1])+
  geom_line(aes(x=date.y,y=betaa_true),color=true_palette[1])+
  labs(x='Year',y='# adult mosquitoes emerging\nfrom pupae per person per day')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
eir_plot <- ggplot(four_years)+
  geom_line(data=sim_seasonal_hist4plots$EIR_sample,aes(x=date,y=value,group=variable),alpha=alpha_var,color=part_palette[2])+
  geom_line(aes(x=date,y=EIR_true),color=true_palette[2])+
  geom_line(data=sim_seasonal_hist4plots$EIR_history,aes(x=date,y=EIR.median),linetype='dashed',color=true_palette[2])+
  labs(x='Year',y='Infectious mosquito bites\nper person per year')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
eir_truth_plot <- ggplot(four_years)+
  geom_line(data=sim_seasonal_hist4plots$EIR_sample,aes(x=date,y=value,group=variable),alpha=0,color=part_palette[2])+
  geom_line(aes(x=date,y=EIR_true),color=true_palette[2])+
  labs(x='Year',y='Infectious mosquito bites\nper person per year')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

incidence_plot <- ggplot(four_years)+
  geom_line(data=sim_seasonal_hist4plots$inc_sample,aes(x=date,y=value*1000,group=variable),alpha=alpha_var,color=part_palette[3])+
  geom_line(aes(x=date,y=inc_all_true*1000),color=true_palette[3])+
  geom_line(data=sim_seasonal_hist4plots$inc_history,aes(x=date,y=incall.median*1000),linetype='dashed',color=true_palette[3])+
  labs(x='Year',y='Clinical cases\nper 1000 persons')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
incidence_truth_plot <- ggplot(four_years)+
  geom_line(data=sim_seasonal_hist4plots$inc_sample,aes(x=date,y=value*1000,group=variable),alpha=0,color=part_palette[3])+
  geom_line(aes(x=date,y=inc_all_true*1000),color=true_palette[3])+
  labs(x='Year',y='Clinical cases\nper 1000 persons')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
prev05_plot <- ggplot(four_years)+
  geom_line(data=sim_seasonal_hist4plots$prev_sample,aes(x=date,y=value,group=variable),alpha=alpha_var,color=part_palette[4])+
  geom_line(aes(x=date,y=prev05),color=true_palette[4])+
  geom_line(data=sim_seasonal_hist4plots$prev_history,aes(x=date,y=prev.median),linetype='dashed',color=true_palette[4])+
  labs(x='Year',y='Prevalence (<5 years old)')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
prev05_truth_plot <- ggplot(four_years_cis)+
  geom_line(data=sim_seasonal_hist4plots$prev_sample,aes(x=date,y=value,group=variable),alpha=0,color=part_palette[4])+
  geom_line(aes(x=date,y=mean),color=true_palette[4])+
  geom_errorbar(aes(x=date,ymin=lower,ymax=upper),width=0,alpha=0,color=true_palette[4])+
  labs(x='Year',y='Prevalence (<5 years old)')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
prev05_data_only_plot <- ggplot(four_years_cis)+
  geom_line(data=sim_seasonal_hist4plots$prev_sample,aes(x=date,y=value,group=variable),alpha=0,color=part_palette[4])+
  geom_point(aes(x=date,y=mean),color=true_palette[4])+
  geom_errorbar(aes(x=date,ymin=lower,ymax=upper),width=0,color=true_palette[4])+
  labs(x='Year',y='Prevalence (<5 years old)')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
prev05_data_plot <- ggplot(four_years_cis)+
  geom_line(data=sim_seasonal_hist4plots$prev_sample,aes(x=date,y=value,group=variable),alpha=alpha_var,color=part_palette[4])+
  geom_line(data=sim_seasonal_hist4plots$prev_history,aes(x=date,y=prev.median),linetype='dashed',color=true_palette[4])+
  geom_point(aes(x=date,y=mean),color=true_palette[4])+
  geom_errorbar(aes(x=date,ymin=lower,ymax=upper),width=0,color=true_palette[4])+
  labs(x='Year',y='Prevalence (<5 years old)')+
  coord_cartesian(ylim = c(0,NA),xlim = as.Date(c('2023-01-01','2024-06-01')))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

moz_plot <- ggplot(four_years)+
  geom_line(data=sim_seasonal_hist4plots$mv_sample,aes(x=date,y=value,group=variable),alpha=0.5,color='lightgrey')+
  geom_line(aes(x=date,y=mv))+
  geom_line(data=sim_seasonal_hist4plots$mv_history,aes(x=date,y=mv.median),linetype='dashed')+
  coord_cartesian(ylim = c(0,NA))+
  scale_x_date(date_labels = "%b",date_breaks = '4 months',expand = c(0,NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

windows(15,6.5)
truth <- betaa_truth_plot + eir_truth_plot + incidence_truth_plot + prev05_truth_plot + plot_layout(nrow=1)
observed <- betaa_truth_plot + eir_truth_plot + incidence_truth_plot + prev05_data_only_plot + plot_layout(nrow=1,ncol=4)
estimates <- betaa_plot + eir_plot + incidence_plot + prev05_data_plot + plot_layout(nrow=1)
ggsave('true_sim_demo_7.tiff',plot = truth,units = 'cm',width = 6,height=12,dpi=300)
ggsave('true_sim_demo_wide_7.tiff',plot = truth,units = 'cm',width = 17,height=4,dpi=300)
ggsave('estimate_sim_demo_wide_7.tiff',plot = estimates,units = 'cm',width = 17,height=4,dpi=300)
ggsave('observed_sim_demo_wide_7.tiff',plot = observed,units = 'cm',width = 17,height=4,dpi=300)
regr_test_df <- cars
str(cars)
str(datasets::iris)
lm(Species~.,data=iris)

names(sim_seasonal_run_result$history[,1,1])

results <- sim_seasonal_run_result
sim_data <- four_years

