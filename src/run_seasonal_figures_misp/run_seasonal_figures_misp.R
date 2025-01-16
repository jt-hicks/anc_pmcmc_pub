orderly2::orderly_shared_resource('prep_results_sim.R')
orderly2::orderly_shared_resource('summarise_smc_results.R')
orderly2::orderly_shared_resource('theme_base.R')

orderly2::orderly_dependency("run_diagnostics_misp", quote(latest(parameter:n_datasets==32)),
                             c('results_list.RDS'))
orderly2::orderly_dependency("run_diagnostics_misp", quote(latest(parameter:n_datasets==16)),
                             c('results_list_control.RDS'='results_list.RDS'))

orderly2::orderly_dependency("create_sim_data", quote(latest()),
                             c('sim_seasonal_dataraw_list.RDS'))

orderly2::orderly_artefact(files=c('misspecification_1mo.tiff','misspecification_1y.tiff'))

source('prep_results_sim.R')
source('summarise_smc_results.R')
source('theme_base.R')

sim_data_raw <- readRDS('sim_seasonal_dataraw_list.RDS')
sim_data_raw <- sim_data_raw[5:8]
results_list <- readRDS('results_list.RDS')
results_list_control <- readRDS('results_list_control.RDS')


misspecification <- c(-0.2,0.2)
admin_list <- c('Tanga','Upper East','Fatick','Equateur')
country_list <- c('Tanzania','Ghana','Senegal','Democratic Republic of Congo')
init_EIR_list <- c(50)
start_pf_time_list <- c(30,90,180,360)
names_list <- c(5:8)

names_key <- expand.grid(misspecification=misspecification,x=c(1:4),start_pf_time=start_pf_time_list)
names_key$name <- names_list[names_key$x]
names_key$admin <- admin_list[names_key$x]
names_key$country <- country_list[names_key$x]
names_key$init_EIR <- 50

sim_data_raw_df <- bind_rows(lapply(1:length(sim_data_raw),function(x){
  df <- sim_data_raw[[x]]
  df$init_EIR <- 50
  return(df)
}))

for(x in 1:length(country_list)){
  admin <- admin_list[x]
  country <- country_list[x]
  init_EIR <- 50
  orderly2::orderly_dependency(name="create_sim_data", query=quote(latest()),
                               c("data/${country}_${init_EIR}.rds" = paste0('sim_data_EIR',init_EIR,'_',admin,'_',country,'.RDS')))
}


prepped_results <- lapply(1:length(results_list),function(x){
  admin <- names_key[x,'admin']
  country <- names_key[x,'country']
  init_EIR <- names_key[x,'init_EIR']
  name <- names_key[x,'name']
  misspecification <- names_key[x,'misspecification']
  start_pf_time <- names_key[x,'start_pf_time']
  sim_data_true <- readRDS(paste0('data/',country,'_',init_EIR,'.rds'))
  sim_data_true$init_EIR <- init_EIR
  prepped <- prep_results_sim(results=results_list[[x]],
               sim_data_raw=sim_data_raw[[names_key[x,'x']]],
               sim_data_true=sim_data_true,
               burnin=0.1,
               site = paste(admin,country,misspecification,start_pf_time,sep='_'),
               anc=FALSE)
  return(prepped)
})

names_key_control <- expand.grid(misspecification=c(0),x=c(1:4),start_pf_time=start_pf_time_list)
names_key_control$name <- names_list[names_key_control$x]
names_key_control$admin <- admin_list[names_key_control$x]
names_key_control$country <- country_list[names_key_control$x]
names_key_control$init_EIR <- 50

prepped_results_control <- lapply(1:length(results_list_control),function(x){
  admin <- names_key_control[x,'admin']
  country <- names_key_control[x,'country']
  init_EIR <- names_key_control[x,'init_EIR']
  name <- names_key_control[x,'name']
  misspecification <- names_key_control[x,'misspecification']
  start_pf_time <- names_key_control[x,'start_pf_time']
  sim_data_true <- readRDS(paste0('data/',country,'_',init_EIR,'.rds'))
  sim_data_true$init_EIR <- init_EIR
  prepped <- prep_results_sim(results=results_list_control[[x]],
                              sim_data_raw=sim_data_raw[[names_key_control[x,'x']]],
                              sim_data_true=sim_data_true,
                              burnin=0.1,
                              site = paste(admin,country,misspecification,start_pf_time,sep='_'),
                              anc=FALSE)
  return(prepped)
})

prepped_results_all <- append(prepped_results,prepped_results_control)

seasonal_short_informed_summary <- bind_rows(lapply(1:length(prepped_results_all), function(x){
  df <- prepped_results_all[[x]]$summary
  names <- unlist(strsplit(unique(df$site),split='_'))
  df$admin <- names[1]
  df$country <- names[2]
  df$misspecification <- as.numeric(names[3])
  df$start_pf_time <- as.integer(names[4])
  df$init_EIR <- 50
  return(df)
}))

measure_palette <- RColorBrewer::brewer.pal(4,name='Set2')

misp_palette <- RColorBrewer::brewer.pal(3,name='Set2')
names(misp_palette) <- c('Annual Prevalence','Underestimated','Overestimated')

results4plot_misp <- seasonal_short_informed_summary %>%
  mutate(month=as.yearmon(date))%>%
  left_join(sim_data_raw_df%>%mutate(date=as.Date(date),month=as.yearmon(date)),by=c('date','admin','country','month','init_EIR'))%>%
  mutate(true_value = ifelse(measure=='prev_05',prev_05,true_value),
         misp_factor = factor(misspecification,levels=c(0,-0.2,0.2),labels=c('Annual Prevalence','Underestimated','Overestimated')),
         country = factor(country,levels=c('Democratic Republic of Congo','Tanzania','Ghana','Senegal'),labels=c('DRC','Tanzania','Ghana','Senegal')))
#&results4plot_misp$date>=as.Date('2021-01-01')
prev_fit <- ggplot(data=results4plot_misp[results4plot_misp$measure=='prev_05'&results4plot_misp$start_pf_time==30&results4plot_misp$date<=as.Date('2021-01-01'),])+
  geom_point(aes(x=date,y=true_value))+
  geom_line(aes(x=date,y=median,color=misp_factor),size=1)+
  # geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=misp_factor),alpha=0.5)+
  # scale_fill_manual(values=misp_palette)+
  scale_color_manual(values=misp_palette)+
  scale_y_continuous(limits = c(0,0.7),expand = c(0,0))+
  facet_grid(.~country)+
  labs(y='RDT Prevalence,\nunder 5 years',
       color='Prevalence used to\ninform baseline burden',
       fill='Prevalence used to\ninform baseline burden')+
  theme(panel.spacing.y = unit(0,'lines'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
inc_fit <- ggplot(data=results4plot_misp[results4plot_misp$measure=='clininc_05'&results4plot_misp$start_pf_time==30&results4plot_misp$date<=as.Date('2021-01-01'),])+
  geom_point(aes(x=date,y=true_value*1000))+
  geom_line(aes(x=date,y=median*1000,color=misp_factor),size=1)+
  # geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=misp_factor),alpha=0.5)+
  scale_y_continuous(limits = c(0,70),expand = c(0,0))+
  # scale_fill_manual(values=misp_palette)+
  scale_color_manual(values=misp_palette)+
  facet_grid(.~country)+
  labs(y='Clinical Cases per\n1000 children under 5 years',
       color='Prevalence used to\ninform baseline burden',
       fill='Prevalence used to\ninform baseline burden')+
  theme(panel.spacing.y = unit(0,'lines'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank())
eir_fit <- ggplot(data=results4plot_misp[results4plot_misp$measure=='EIR'&results4plot_misp$start_pf_time==30&results4plot_misp$date<=as.Date('2021-01-01'),])+
  geom_point(aes(x=date,y=true_value))+
  geom_line(aes(x=date,y=median,color=misp_factor),size=1)+
  # geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=misp_factor),alpha=0.5)+
  # scale_fill_manual(values=misp_palette)+
  scale_color_manual(values=misp_palette)+
  scale_y_log10(breaks=c(0.1,1,10,100),labels=c(0.1,1,10,100))+
  facet_grid(.~country)+
  labs(y='EIR',
       color='Prevalence used to\ninform baseline burden',
       fill='Prevalence used to\ninform baseline burden')+
  theme(panel.spacing.y = unit(0,'lines'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank())
betaa_fit <- ggplot(data=results4plot_misp[results4plot_misp$measure=='betaa'&results4plot_misp$start_pf_time==30&results4plot_misp$date<=as.Date('2021-01-01'),])+
  geom_point(aes(x=date,y=true_value))+
  geom_line(aes(x=date,y=median,color=misp_factor),size=1)+
  # geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=misp_factor),alpha=0.5)+
  # scale_y_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_y_log10(breaks=c(0.01,0.1,1,10),labels=c(0.01,0.1,1,10))+
  # scale_fill_manual(values=misp_palette)+
  scale_color_manual(values=misp_palette)+
  facet_grid(.~country)+
  labs(y='Mosquito Emergence Rate',
       color='Prevalence used to\ninform baseline burden',
       fill='Prevalence used to\ninform baseline burden')+
  scale_x_date(breaks = as.Date(c('2017-01-01','2018-01-01','2019-01-01','2020-01-01','2021-01-01')),labels=c('2017','2018','2019','2020','2021'))+
  theme(panel.spacing.y = unit(0,'lines'),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1,vjust=1.2))
library(patchwork)
composite_misp_1mo <- prev_fit+inc_fit+eir_fit+betaa_fit+plot_layout(ncol=1,nrow=4,guides = 'collect')+
  plot_annotation(tag_levels = c('A','B','C','D'))&
  theme(text = element_text(size=10),
        axis.title.y = element_text(size=12),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0),
        panel.spacing.y = unit(3,'pt'),
        plot.tag = element_text(size=14))
composite_misp_1mo
ggsave(plot=composite_misp_1mo,filename = 'misspecification_1mo.tiff',width=7,height=9,unit='in')

prev_fit_1y <- ggplot(data=results4plot_misp[results4plot_misp$measure=='prev_05'&results4plot_misp$start_pf_time==360&results4plot_misp$date<=as.Date('2021-01-01'),])+
  geom_point(aes(x=date,y=true_value))+
  geom_line(aes(x=date,y=median,color=misp_factor),size=1)+
  # geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=misp_factor),alpha=0.5)+
  # scale_fill_manual(values=misp_palette)+
  scale_color_manual(values=misp_palette)+
  scale_y_continuous(limits = c(0,0.7),expand = c(0,0))+
  facet_grid(.~country)+
  labs(y='RDT Prevalence,\nunder 5 years',
       color='Prevalence used to\ninform baseline burden',
       fill='Prevalence used to\ninform baseline burden')+
  theme(panel.spacing.y = unit(0,'lines'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
inc_fit_1y <- ggplot(data=results4plot_misp[results4plot_misp$measure=='clininc_05'&results4plot_misp$start_pf_time==360&results4plot_misp$date<=as.Date('2021-01-01'),])+
  geom_point(aes(x=date,y=true_value*1000))+
  geom_line(aes(x=date,y=median*1000,color=misp_factor),size=1)+
  # geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=misp_factor),alpha=0.5)+
  scale_y_continuous(limits = c(0,70),expand = c(0,0))+
  # scale_fill_manual(values=misp_palette)+
  scale_color_manual(values=misp_palette)+
  facet_grid(.~country)+
  labs(y='Clinical Cases per\n1000 children under 5 years',
       color='Prevalence used to\ninform baseline burden',
       fill='Prevalence used to\ninform baseline burden')+
  theme(panel.spacing.y = unit(0,'lines'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank())
eir_fit_1y <- ggplot(data=results4plot_misp[results4plot_misp$measure=='EIR'&results4plot_misp$start_pf_time==360&results4plot_misp$date<=as.Date('2021-01-01'),])+
  geom_point(aes(x=date,y=true_value))+
  geom_line(aes(x=date,y=median,color=misp_factor),size=1)+
  # geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=misp_factor),alpha=0.5)+
  # scale_fill_manual(values=misp_palette)+
  scale_color_manual(values=misp_palette)+
  scale_y_log10(breaks=c(0.1,1,10,100),labels=c(0.1,1,10,100))+
  facet_grid(.~country)+
  labs(y='EIR',
       color='Prevalence used to\ninform baseline burden',
       fill='Prevalence used to\ninform baseline burden')+
  theme(panel.spacing.y = unit(0,'lines'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank())
betaa_fit_1y <- ggplot(data=results4plot_misp[results4plot_misp$measure=='betaa'&results4plot_misp$start_pf_time==360&results4plot_misp$date<=as.Date('2021-01-01'),])+
  geom_point(aes(x=date,y=true_value))+
  geom_line(aes(x=date,y=median,color=misp_factor),size=1)+
  # geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=misp_factor),alpha=0.5)+
  # scale_y_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_y_log10(breaks=c(0.01,0.1,1,10),labels=c(0.01,0.1,1,10))+
  # scale_fill_manual(values=misp_palette)+
  scale_color_manual(values=misp_palette)+
  facet_grid(.~country)+
  labs(y='Mosquito Emergence Rate',
       color='Prevalence used to\ninform baseline burden',
       fill='Prevalence used to\ninform baseline burden')+
  scale_x_date(breaks = as.Date(c('2017-01-01','2018-01-01','2019-01-01','2020-01-01','2021-01-01')),labels=c('2017','2018','2019','2020','2021'))+
  theme(panel.spacing.y = unit(0,'lines'),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1,vjust=1.2))
composite_misp_1y <- prev_fit_1y+inc_fit_1y+eir_fit_1y+betaa_fit_1y+plot_layout(ncol=1,nrow=4,guides = 'collect')+
  plot_annotation(tag_levels = c('A','B','C','D'))&
  theme(text = element_text(size=10),
        axis.title.y = element_text(size=12),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0),
        panel.spacing.y = unit(3,'pt'),
        plot.tag = element_text(size=14))
composite_misp_1y
ggsave(plot=composite_misp_1y,filename = 'misspecification_1y.tiff',width=7,height=9,unit='in')
