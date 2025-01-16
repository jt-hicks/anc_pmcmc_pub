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
library(hipercow)
source('utils.R')

theme_set(theme_minimal(base_size = 7)+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))

###Set up cluster
obj_sim_trial <- cluster_setup(context_name = 'sim_trial',template='20Core', cores=10)

##Get saved processed data
WA_pg_data_list <- readRDS('Q:/anc_pmcmc/trial/Data/WA_pg_data_list.rds')
WA_sg_data_list <- readRDS('Q:/anc_pmcmc/trial/Data/WA_sg_data_list.rds')
WA_all_data_list <- readRDS('Q:/anc_pmcmc/trial/Data/WA_all_data_list.rds')

EA_pg_data_list <- readRDS('Q:/anc_pmcmc/trial/Data/EA_pg_data_list.rds')
EA_mg_data_list <- readRDS('Q:/anc_pmcmc/trial/Data/EA_mg_data_list.rds')

#MAP estimates for 2010

##Need to convert first year of ANC prev to under 5 prevalence
prevs_wa <- c(`Burkina Faso`=0.583962944,Gambia=0.127724696,Ghana=0.681794701,Mali=0.297692604)
country_ea <- names(EA_pg_data_list)
prevs_ea <- c(0.42,0.391)

trials_pg_data_all_list <- append(WA_pg_data_list,EA_pg_data_list)
trials_mg_data_all_list <- append(WA_sg_data_list,EA_mg_data_list)
trials_pg_data_all <- bind_rows(trials_pg_data_all_list)%>%
  mutate(measure='prev_pg')%>%
  select(site, month, positive, tested, measure)
trials_mg_data_all <- bind_rows(trials_mg_data_all_list)%>%
  mutate(measure=case_when(
    gravidity == 2 ~ "prev_sg",
    grav_cat=='mg' ~ "prev_mg",
    TRUE ~ NA))%>%
  select(site, month, positive, tested, measure)
trials_both_data_all <- bind_rows(trials_pg_data_all,trials_mg_data_all)
trials_both_data_all <- addCIs(trials_both_data_all,Ys=trials_both_data_all$positive,Ns=trials_both_data_all$tested)
trials_both_data_all <- trials_both_data_all%>%
  rename(country=site)
admins_trial <- c(`Burkina Faso` = 'Plateau-Central',
                  Gambia = 'Upper River',
                  Ghana = 'Upper East',
                  Mali = 'Koulikoro',
                  Kenya = 'Siaya',
                  Malawi = 'Blantyre') #Two sites are in Blantyre, one is in neighboring Chikwawa
country_trial <- names(trials_pg_data_all_list)

first_annual_prev_trials <- trials_pg_data_all%>%
  group_by(site)%>%
  slice_head(n=12)%>%
  summarise(annual_prev = sum(positive)/sum(tested))%>%
  select(site,annual_prev)
targetprevs_trials <- sapply(first_annual_prev_trials$annual_prev, function(x) get_u5prev_fromanc(avg_prev=x,comparison='pg'))
plot(targetprevs_trials,first_annual_prev_trials$annual_prev,xlim=c(0,1),ylim=c(0,1))
abline(0,1)
names(targetprevs_trials) <- names(trials_pg_data_all_list)
comparisons_trials <- c(`Burkina Faso` = 'pgsg',
                        Gambia = 'pgsg',
                        Ghana = 'pgsg',
                        Mali = 'pgsg',
                        Kenya = 'pgmg',
                        Malawi = 'pgmg')

obj_sim_trial$login()
trials_short_informed_run_1 <- obj_sim_trial$enqueue_bulk(1:6,function(x,data_raw_pg,data_raw_mg,annual_prev,country,comparisons){
  sifter::run_pmcmc(data_raw_pg=data_raw_pg[[x]],
                    data_raw_mg=data_raw_mg[[x]],
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
                    comparison = comparisons[[country[[x]]]],
                    initial = 'informed',
                    check_flexibility=TRUE
  )
},data_raw_pg=trials_pg_data_all_list,data_raw_mg=trials_mg_data_all_list,annual_prev=targetprevs_trials,country=country_trial,comparisons = comparisons_trials)
trials_short_informed_run_1$status() #gentlewomanly_portuguesemanofwar
trials_short_informed_run_1 <- obj_sim_trial$task_bundle_get('gentlewomanly_portuguesemanofwar')
trials_short_informed_run_1$tasks[[4]]$log()
sapply(trials_pg_data_all_list,function(x) nrow(x))

trials_short_informed_results <- lapply(1:6, function(id){
  trials_short_informed_run_1$tasks[[id]]$result()
})
trials_short_informed_plots <- lapply(1:6,function(i) create_diag_figs(trials_short_informed_results[[i]],
                                                                                country = country_trial[[i]],
                                                                                district = 1,
                                                                                folderpath = './trial_figs/diag',
                                                                                name = 'trials_short_informed'))
trials_short_informed_props <- lapply(1:6, function(id){
  var(trials_short_informed_results[[id]]$pars[500:1000,])
})

trials_short_informed_prepped <- lapply(1:6, function(x) {
  prep_results(results=trials_short_informed_results[[x]],
               sim_data=trials_pg_data_all_list[[x]],
               burnin=0.5,
               country = country_trial[[x]],
               district = NA,
               anc=TRUE)
})
trials_short_informed_summary <- bind_rows(lapply(1:6, function(x){
  trials_short_informed_prepped[[x]]$summary
}))
trials_short_informed_sample <- bind_rows(lapply(1:6, function(x){
  trials_short_informed_prepped[[x]]$sample
}))
volatility_sum_trials <- bind_rows(lapply(1:6, function(x){
  get_param_sum(results=trials_short_informed_results[[x]],
                param_name='volatility',
                posterior_name='log_posterior',
                n_chains=1,
                burnin=0.5,
                country=country_trial[[x]])
}))%>%
  select(measure,country,median,lower_hpd,upper_hpd)
measure_levels <- c('prev_pg','prev_sg','prev_mg','prev_05','incall','EIR','betaa')
ggplot()+
  geom_line(data=trials_short_informed_sample,aes(x=date,y=value,group=variable),color='lightgrey',alpha=0.5)+
  geom_line(data = trials_short_informed_summary, aes(x=date,y=median),color='darkgrey',linewidth=1)+
  geom_point(data=trials_both_data_all,aes(x=as.Date(month,frac=0.5),y=mean))+
  geom_errorbar(data=trials_both_data_all,aes(x=as.Date(month,frac=0.5),ymin=lower,ymax=upper))+
  facet_wrap(factor(country,country_trial)~factor(measure,levels=measure_levels),scales = 'free_y',nrow=6,ncol=6)+
  scale_y_continuous(limits = c(0,1))

table(trials_short_informed_summary$measure)
source('./create_trial_plots_function.R')
trial_rainfall4plot <- bind_rows(trial_rainfall %>%
                                   mutate(month_car=as.character(month))%>%
                                   right_join(trial_dates,by=c('country','month_car'))%>%
                                   group_by(country) %>%
                                   group_split() %>%
                                   lapply(function(x,window_size){
                                     df <- calc_annual_proportion(x)
                                     df$country <- unique(x$country)

                                     return(right_join(df,x,by=join_by(month==month_num,country==country)))}))
trials_prev_pg_plot <- create_dashboard_plots_trial(results=trials_short_informed_prepped,
                                                 observed = trials_pg_data_all,
                                                 rainfall = trial_rainfall,
                                                 var= 'prev_pg',
                                                 title = 'ANC Prevalence\nPrimigravidae')
trials_anc_data_all <- bind_rows(trials_pg_data_all,trials_mg_data_all)%>%
  mutate(measure=ifelse(measure=='prev_pg','prev_pg','prev_mg'))
trials_prev_anc_plot <- create_dashboard_plots_trial(results=trials_short_informed_prepped,
                                                    observed = trials_anc_data_all,
                                                    rainfall = trial_rainfall,
                                                    var= 'prev_anc'
                                                    )
trials_prev_anc_plot_data <- create_dashboard_plots_trial(results=trials_short_informed_prepped,
                                                     observed = trials_anc_data_all,
                                                     rainfall = trial_rainfall,
                                                     var= 'prev_anc',
                                                     show_fits = FALSE
)
trials_prev_anc_plot_2levels <- create_dashboard_plots_trial_2(results=trials_short_informed_prepped,
                                                     observed = trials_anc_data_all,
                                                     rainfall = trial_rainfall,
                                                     var= 'prev_anc'
)
ggsave('./trial_figs/output/trials_prev_anc_plot_2levels.tiff',plot=trials_prev_anc_plot_2levels,width=2,height=5,units='in')

trials_prev_mg_plot <- create_dashboard_plots_trial(results=trials_short_informed_prepped,
                                                    observed = trials_mg_data_all,
                                                    rainfall = trial_rainfall,
                                                    var= 'prev_mg',
                                                    title = 'ANC Prevalence\nSecundi- or Multigravidae')
saveRDS(trials_short_informed_prepped,'trials_short_informed_prepped.rds')
africa_shape <- terra::vect('C:/Users/jthicks/Documents/africa_shape_files/world-administrative-boundaries.shp')
#see markham_index.R for downloading and formatting admin shapes
admin_list <- list(kassena_nankana_shape,oubritenga_shape,kati_san_shape,fulladu_east_shape,
                                 siaya_shape,blant_chik_shape)
library(sf)
area=kassena_nankana_shape
# Calculate centroids for each administrative area
centroids <- lapply(admin_list, function(area) {
  centroid <- st_centroid(st_as_sf(area))
  country <- area$COUNTRY  # Adjust this based on the actual attribute name for country in your shapefiles
  centroid$country <- country
  return(centroid)
})
# Combine centroids into a single spatial object
centroids_sf <- bind_rows(centroids)
africa_sf <- st_as_sf(africa_shape)
country_iso <- c('GHA','BFA','MLI','GMB','KEN','MWI')

africa_sf$include <- ifelse(africa_sf$iso3 %in% country_iso,'yes','no')

library(tidyterra)
trial_color_palette <- c(viridis::viridis(6,begin=0,end=0.95))
names(trial_color_palette) <-  c('Ghana','Burkina Faso', 'Mali', 'Gambia','Malawi','Kenya')
labels_trial_colors <- c('Ghana','Burkina Faso', 'Mali', 'The Gambia','Malawi','Kenya')
names(labels_trial_colors) <- names(trial_color_palette)
included_palette <- c(yes='white',no='lightgrey')
included_countries <- africa_sf[africa_sf$include=='yes',]
trial_sites_map <- ggplot(data = africa_sf) +
  geom_sf(color = "darkgrey",fill='white') +
  geom_sf(data=included_countries,fill='darkgrey')+
  geom_sf(data = centroids_sf, aes(fill = country), size = 2,shape=21,color='black') +
  scale_fill_manual(values=trial_color_palette,labels=labels_trial_colors,guide='none')+
  coord_sf(xlim=c(NA,55))+
  # scale_fill_manual(values=included_palette, guide='none')+
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        panel.spacing.y = unit(3, "mm"))

fit_panel <- trial_sites_map + trials_prev_anc_plot + plot_layout(ncol=2)
ggsave('./trial_figs/output/trials_short_fit_with_map_3.tiff',plot=fit_panel,width=7,height=5,units='in')

ggplot(trials_pg_data_all)+
  geom_point(aes(x=month,y=tested),fill='orange')
observed_incidence <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/navrongo_cases.txt')
observed_incidence$date_aligned <- as.Date(date_decimal(observed_incidence$date))
# max(observed_incidence$cases)
observed_incidence$mean <- observed_incidence$cases/142800
observed_incidence$country <- 'Ghana'
observed_incidence$lower <- NA
observed_incidence$upper <- NA
observed_incidence <- observed_incidence %>%
  filter(date_aligned >= min(as.Date(trials_pg_data_all[trials_pg_data_all$site=='Ghana',]$month)) & date_aligned <= max(as.Date(trials_pg_data_all[trials_pg_data_all$site=='Ghana',]$month)))
trials_incidence_plot <- create_dashboard_plots_trial(results=trials_short_informed_prepped,
                                                    observed = NULL,
                                                    rainfall = trial_rainfall,
                                                    var= 'inc05',
                                                    title = NULL,
                                                    max_value = 25,
                                                    multiplier = 1000,
                                                    rainfall_multiplier = 20)
ggsave('./trial_figs/output/trials_incidence_plot.tiff',plot=trials_incidence_plot,width=2,height=5,units='in')

trials_incidence_plot_pres <- create_dashboard_plots_trial_pres(results=trials_short_informed_prepped,
                                                      observed = NULL,
                                                      rainfall = trial_rainfall,
                                                      var= 'inc05',
                                                      title = NULL,
                                                      max_value = 25,
                                                      multiplier = 1000,
                                                      rainfall_multiplier = 20)
# ggsave('./trial_figs/output/trials_incidence_plot_pres.tiff',plot=trials_incidence_plot_pres,width=2,height=5,units='in')

trials_betaa_plot <- create_dashboard_plots_trial(results=trials_short_informed_prepped,
                                                      observed = NULL,
                                                      var= 'betaa',
                                                     rainfall = trial_rainfall4plot,
                                                     title = 'Mosquito Emergence\nNumber mosquitos per person per day',
                                                      max_value = NA,
                                                      multiplier = 1,
                                                      facet_scales = 'free_y')
trials_plots <- trials_prev_pg_plot + trials_prev_mg_plot + trials_incidence_plot + plot_layout(ncol=3)
ggsave('./trial_figs/output/trials_short_fit_plots_2.tiff',plot=trials_plots,width=7,height=5,units='in')

fit_panel <- trial_sites_map + trials_prev_anc_plot_data + trials_incidence_plot + plot_layout(ncol=3,widths=c(0.75,1,0.5))
ggsave('./trials_prev_anc_plot_data.tiff',plot=trials_prev_anc_plot_data,width=10,height=10,units='cm')
ggsave('./trials_prev_anc_plot_fits.tiff',plot=trials_prev_anc_plot,width=10,height=10,units='cm')
ggsave('./trial_figs/output/trials_short_fit_with_map_data.tiff',plot=fit_panel,width=26,height=13,units='cm')
fit_panel <- trial_sites_map + trials_prev_anc_plot + trials_incidence_plot + plot_layout(ncol=3,widths=c(0.75,1,0.5))
ggsave('./trial_figs/output/trials_short_fit_with_map_withinc.tiff',plot=fit_panel,width=26,height=13,units='cm')

###Longer runs####
windows_authenticate()
hipercow::hipercow_init()
hipercow::hipercow_configure('windows')
hipercow::hipercow_provision()
hipercow_environment_create(packages=c('dplyr','ggplot2'))
orderly2::orderly_init()
trials_short_informed_props

names(trials_pg_data_all_list)

df_trial_submit <- data.frame(name=names(trials_pg_data_all_list),
                              prop=unlist(trials_short_informed_props),
                              comparison_trial = c('pgsg','pgsg','pgsg','pgsg','pgmg','pgmg'))
saveRDS(trials_pg_data_all_list,'./src/run_trial_pmcmc/trials_pg_data_all_list.RDS')
saveRDS(trials_mg_data_all_list,'./src/run_trial_pmcmc/trials_mg_data_all_list.RDS')
resources <- hipercow_resources(cores=32)
long_trial_id <- task_create_bulk_expr(orderly2::orderly_run('run_trial_pmcmc',parameters=list(name=name,
                                                                                    proposal_matrix = prop,
                                                                                    comparison_trial = comparison_trial,
                                                                                    length=10000)),
                                  data=df_trial_submit,
                                  resources=resources)
task_status(long_trial_id$ids) #greenish_ilsamochadegu
task_log_show(long_trial_id$ids[6])
df_trial_submit_2 <- df_trial_submit[c(1,3),]
long_trial_id <- task_create_bulk_expr(orderly2::orderly_run('run_trial_pmcmc',parameters=list(name=name,
                                                                                               proposal_matrix = prop,
                                                                                               comparison_trial = comparison_trial,
                                                                                               length=10000)),
                                       data=df_trial_submit,
                                       resources=resources)
long_trial_id_2 <- task_create_bulk_expr(orderly2::orderly_run('run_trial_pmcmc',parameters=list(name=name,
                                                                                               proposal_matrix = prop,
                                                                                               comparison_trial = comparison_trial,
                                                                                               length=10000,
                                                                                               seed=209)),
                                       data=df_trial_submit_2,
                                       resources=resources)
task_status(long_trial_id_2$ids)
task_log_show(long_trial_id_2$ids[2])

df_trial_submit_3 <- df_trial_submit[c(1),]
long_trial_id_3 <- task_create_bulk_expr(orderly2::orderly_run('run_trial_pmcmc',parameters=list(name=name,
                                                                                                 proposal_matrix = prop,
                                                                                                 comparison_trial = comparison_trial,
                                                                                                 length=10000,
                                                                                                 seed=309)),
                                         data=df_trial_submit_3,
                                         resources=resources)
task_status(long_trial_id_3$ids) #'quasihonourable_abyssiniangroundhornbill'
# task_cancel(long_trial_id_3$ids)

gambia <- trials_short_informed_prepped[[2]]$summary %>%
  filter(country=='Gambia'&measure=='inc05')
ggplot(gambia)+
  geom_point(aes(x=date,y=median))
##Explore fitted betaa
trials_short_informed_betaa <- trials_short_informed_sample%>%
  filter(measure=='betaa')%>%
  group_by(country,variable)%>%
  mutate(last_value=lag(value),
         step=log(value/last_value))%>%
  filter(value!=125)
ggplot(trials_short_informed_betaa)+
  geom_histogram(aes(x=step))+
  facet_wrap(.~country)
plot(sim_seas_long_informed_run_bangon_ghana$tasks[[1]]$result()$mcmc[2001:20000,'volatility'],sim_seas_long_informed_run_bangon_ghana$tasks[[1]]$result()$mcmc[2001:20000,'posterior'])

