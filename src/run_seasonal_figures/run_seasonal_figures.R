orderly2::orderly_shared_resource('prep_results_sim.R')
orderly2::orderly_shared_resource('summarise_smc_results.R')
orderly2::orderly_shared_resource('theme_base.R')

orderly2::orderly_dependency("run_diagnostics", quote(latest()),
                             c('results_list.RDS'))
orderly2::orderly_dependency("create_sim_data", quote(latest()),
                             c('sim_seasonal_dataraw_list.RDS'))

orderly2::orderly_artefact(files=c('fit2truth.tiff','prop_in_ci.tiff','correlationwithtrue.tiff'),description = 'Figures for supplement')

source('prep_results_sim.R')
source('summarise_smc_results.R')
source('theme_base.R')

sim_data_raw <- readRDS('sim_seasonal_dataraw_list.RDS')
results_list <- readRDS('results_list.RDS')

admin_list <- c('Tanga','Upper East','Fatick','Equateur')
country_list <- c('Tanzania','Ghana','Senegal','Democratic Republic of Congo')
init_EIR_list <- c(10,50,100,500)

names_key <- expand.grid(x=1:4,y=1:4)
names_key$name <- c(1:16)
names_key$admin <- admin_list[names_key$x]
names_key$country <- country_list[names_key$x]
names_key$init_EIR <- init_EIR_list[names_key$y]

sim_data_raw_df <- bind_rows(lapply(1:length(sim_data_raw),function(x){
  df <- sim_data_raw[[x]]
  df$init_EIR <- names_key[x,'init_EIR']
  return(df)
}))
for(x in 1:length(results_list)){
  admin <- names_key[x,'admin']
  country <- names_key[x,'country']
  init_EIR <- names_key[x,'init_EIR']
  orderly2::orderly_dependency(name="create_sim_data", query=quote(latest()),
                               c("data/${x}.rds" = paste0('sim_data_EIR',init_EIR,'_',admin,'_',country,'.RDS')))
}

prepped_results <- lapply(1:length(results_list),function(x){
  admin <- names_key[x,'admin']
  country <- names_key[x,'country']
  init_EIR <- names_key[x,'init_EIR']
  sim_data_true <- readRDS(paste0('data/',x,'.rds'))
  sim_data_true$init_EIR <- init_EIR
  prepped <- prep_results_sim(results=results_list[[x]],
               sim_data_raw=sim_data_raw[[x]],
               sim_data_true=sim_data_true,
               burnin=0.1,
               site = paste(admin,country,init_EIR,sep='-'),
               anc=FALSE)
  return(prepped)
})


seasonal_short_informed_summary <- bind_rows(lapply(1:length(results_list), function(x){
  df <- prepped_results[[x]]$summary
  names <- unlist(strsplit(unique(df$site),split='-'))
  df$admin <- names[1]
  df$country <- names[2]
  df$init_EIR <- as.integer(names[3])
  return(df)
}))
seasonal_short_r_squared <- seasonal_short_informed_summary%>%
  filter(measure!='prev_05')%>%
  group_by(measure,site,admin,country,init_EIR)%>%
  summarise(r_squared=cor(true_value,median)^2)
seasonal_short_informed_pointest <- bind_rows(lapply(1:length(results_list), function(x){
  df <- prepped_results[[x]]$point_est
  names <- unlist(strsplit(unique(df$site),split='-'))
  df$admin <- names[1]
  df$country <- names[2]
  df$init_EIR <- as.integer(names[3])
  return(df)
}))%>%
  mutate(country=factor(country,levels=c('Democratic Republic of Congo','Tanzania','Ghana','Senegal'),labels=c('DRC','Tanzania','Ghana','Senegal')))
seasonal_short_informed_sample <- bind_rows(lapply(1:length(results_list), function(x){
  df <- prepped_results[[x]]$sample
  names <- unlist(strsplit(unique(df$site),split='-'))
  df$admin <- names[1]
  df$country <- names[2]
  df$init_EIR <- as.integer(names[3])
  return(df)
}))

measure_palette <- RColorBrewer::brewer.pal(4,name='Set2')

results4plot <- seasonal_short_informed_summary %>%
  filter(init_EIR != 500)%>%
  mutate(month=as.yearmon(date))%>%
  left_join(sim_data_raw_df%>%mutate(date=as.Date(date),month=as.yearmon(date)),by=c('date','admin','country','month','init_EIR'))%>%
  mutate(true_value = ifelse(measure=='prev_05',prev_05,true_value),
         country = factor(country,levels=c('Democratic Republic of Congo','Tanzania','Ghana','Senegal'),labels=c('DRC','Tanzania','Ghana','Senegal')))

prev_fit <- ggplot(data=results4plot[results4plot$measure=='prev_05'&results4plot$date>=as.Date('2021-01-01'),])+
  geom_point(aes(x=date,y=true_value))+
  geom_line(aes(x=date,y=median),color=measure_palette[1],size=1)+
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper),fill=measure_palette[1],alpha=0.5)+
  scale_y_continuous(limits = c(0,0.7),expand = c(0,0))+
  facet_grid(init_EIR~country)+
  labs(y='RDT Prevalence,\nunder 5 years')+
  theme(panel.spacing.y = unit(0,'lines'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
inc_fit <- ggplot(data=results4plot[results4plot$measure=='clininc_05'&results4plot$date>=as.Date('2021-01-01'),])+
  geom_point(aes(x=date,y=true_value*1000))+
  geom_line(aes(x=date,y=median*1000),color=measure_palette[2],size=1)+
  geom_ribbon(aes(x=date,ymin=lower*1000,ymax=upper*1000),fill=measure_palette[2],alpha=0.5)+
  scale_y_continuous(limits = c(0,50),expand = c(0,0),breaks=c(0,20,40))+
  facet_grid(init_EIR~country)+
  labs(y='Clinical Cases per\n1000 children under 5 years')+
  theme(panel.spacing.y = unit(0,'lines'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank())
eir_fit <- ggplot(data=results4plot[results4plot$measure=='EIR'&results4plot$date>=as.Date('2021-01-01'),])+
  geom_point(aes(x=date,y=true_value))+
  geom_line(aes(x=date,y=median),color=measure_palette[3],size=1)+
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper),fill=measure_palette[3],alpha=0.5)+
  scale_y_log10(breaks=c(0.01,1,100),labels=c(0.01,1,100))+
  facet_grid(init_EIR~country)+
  labs(y='EIR')+
  theme(panel.spacing.y = unit(0,'lines'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank())
betaa_fit <- ggplot(data=results4plot[results4plot$measure=='betaa'&results4plot$date>=as.Date('2021-01-01'),])+
  geom_point(aes(x=date,y=true_value))+
  geom_line(aes(x=date,y=median),color=measure_palette[4],size=1)+
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper),fill=measure_palette[4],alpha=0.5)+
  # scale_y_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_y_log10(breaks=c(0.001,0.1,10),labels=c(0.001,0.1,10))+
  facet_grid(init_EIR~country)+
  labs(y='Mosquito Emergence Rate')+
  scale_x_date(breaks = as.Date(c('2021-01-01','2022-01-01','2023-01-01','2024-01-01')),labels=c('2021','2022','2023','2024'))+
  theme(panel.spacing.y = unit(0,'lines'),
        axis.title.x = element_blank(),
        strip.text.x = element_blank())
composite <- prev_fit+inc_fit+eir_fit+betaa_fit+plot_layout(ncol=1,nrow=4)+
  plot_annotation(tag_levels = c('A','B','C','D'))&
  theme(text = element_text(size=10),
        axis.title.y = element_text(size=12),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0),
        panel.spacing.y = unit(3,'pt'),
        plot.tag = element_text(size=14))
plot_anno <- grid::textGrob('Initial EIR',rot=-90,hjust=0,vjust=2)
composite_plus_text <- (composite | plot_anno) + plot_layout(widths=c(9,1))+
  plot_annotation(tag_levels = list(c('A','B','C','D'),NA))&
  theme(text = element_text(size=10),
        axis.title.y = element_text(size=12),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0),
        panel.spacing.y = unit(3,'pt'),
        plot.tag = element_text(size=12))
ggsave(plot=composite_plus_text,filename = 'fit2truth.tiff',width=7,height=9,unit='in')

measure_levels <- c('prev_05','clininc_05','EIR','betaa')

results4plot_errors <- results4plot%>%
  mutate(rel_error = (median-true_value)/true_value,
         abs_error = median-true_value)

df_split_errors <- split(results4plot_errors,~measure+country)

df_split <- split(results4plot,~measure+country)
x <- c('DRC','Tanzania','Ghana','Senegal')
y <- c('prev_05','clininc_05','EIR','betaa')
order <- expand.grid(x,y)%>%
  mutate(name=paste(Var2,Var1,sep='.'))
df_split <- df_split[order$name]

plot_fun <- function(x, y) {
  facet_layer <- if (!grepl("Senegal$", y) && grepl("^prev_05", y))
    facet_grid(.~country)
  else if (grepl("Senegal$", y) && grepl("^prev_05", y))
    facet_grid(measure~country)
  else if (grepl("Senegal$", y) && !grepl("^prev_05", y))
    facet_grid(measure~.)

  scale_y_layer <- if (grepl("^EIR", y) || grepl("^betaa",y))
    scale_y_log10(expand = c(0,0))
  else
    scale_y_continuous(expand = c(0,0))

  scale_x_layer <- if (grepl("^EIR", y) || grepl("^betaa",y))
    scale_x_log10()

  measure_labels <- c(prev_05='Prevalence, <5yo',clininc_05='Clinical Incidence, <5yo',EIR='EIR',betaa='Mosquito Emergence')
  x$measure <- factor(x$measure, labels=measure_labels[unique(x$measure)])

  df_error <- if (grepl("^EIR", y) || grepl("^betaa", y)){
    x %>%
      mutate(true_value_log=log(true_value),
             median_log = log(median))%>%
      ungroup()%>%
      group_by(measure,country)%>%
      summarise(r_squared=round(cor(true_value_log,median_log)^2,digits=2),
                sq_error = round(sqrt(sum(((median_log-true_value_log)/true_value_log)^2)/n()),digits=2),
                max_y = max(upper),
                min_x = min(true_value))
  }else{
    x %>%
      ungroup()%>%
      group_by(measure,country)%>%
      summarise(r_squared=round(cor(true_value,median)^2,digits=2),
                sq_error = round(sqrt(sum(((median-true_value)/true_value)^2)/n()),digits=2),
                max_y = max(upper),
                min_x = min(true_value))
  }

  ggplot(data=x)+
    geom_point(aes(x=true_value,y=median,color=factor(init_EIR)))+
    geom_errorbar(aes(x=true_value,ymin=lower,ymax=upper,color=factor(init_EIR)),width=0)+
    geom_abline(linetype='dashed',color='black')+
    geom_text(data=df_error,
              aes(label = paste0('R2 = ',r_squared,'\nRMSE = ',sq_error),x=min_x,y=max_y),
                  size=2,hjust=0,vjust=1.1)+
    scale_y_layer+
    scale_x_layer+
    scale_color_manual(values=init_eir_palette)+
    labs(color='Initial EIR')+
    # scale_y_log10()+
    # facet_grid(factor(measure,levels=measure_levels)~country,scales='free',labeller = label_wrap_gen(multi_line=FALSE))+
    # stat_cor(aes(x=true_value,y=median,label = after_stat(rr.label)), color = "black", size=2)+
    labs(x='True Value',y='Estimated Value') +
    facet_layer
}
init_eir_palette <- brewer.pal(4,'Reds')[2:4]

grid_plots <- purrr::imap(df_split, plot_fun) %>%
  wrap_plots()
grid_plots_final <- grid_plots+plot_layout(guides = 'collect')&
  theme(text = element_text(size=8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
ggsave(plot=grid_plots_final,filename='correlationwithtrue.tiff',width=7,height=6,units = 'in')

# corr_all <- ggplot(data=results4plot)+
#   geom_point(aes(x=true_value,y=median,color=factor(init_EIR)))+
#   geom_errorbar(aes(x=true_value,ymin=lower,ymax=upper,color=factor(init_EIR)),width=0)+
#   geom_abline(linetype='dashed',color='black')+
#   scale_y_continuous(expand = c(0,0))+
#   scale_color_manual(values=init_eir_palette)+
#   # scale_y_log10()+
#   facet_grid(factor(measure,levels=measure_levels)~country,scales='free',labeller = label_wrap_gen(multi_line=FALSE))+
#   stat_cor(aes(x=true_value,y=median,label = after_stat(rr.label)), color = "black", size=2)+
#   labs(x='True Value',y='Estimated Value')

prop_in_ci <- ggplot(data=seasonal_short_informed_pointest[seasonal_short_informed_pointest$measure!='prev_05'&seasonal_short_informed_pointest$init_EIR!=500,])+
  geom_point(aes(x=country,y=truth_captured,color=factor(init_EIR),shape=factor(measure,levels = c('clininc_05','clininc_all','EIR','betaa'),labels = c('Incidence, under 5 years','Incidence, all ages','EIR','Mosquito Emergence'))))+
  scale_color_manual(values=init_eir_palette)+
  labs(y='Proportion of Months with True Value in 95%CI',
       color='Initial EIR',
       shape = 'Transmission\nindicator')+
  scale_y_continuous(limits=c(0,1.05),expand=c(0,0))+
  theme(axis.title.y = element_blank(),
        legend.position="bottom",
        legend.box = "vertical",
        legend.justification = "right",
        axis.title.x = element_text(size=6,margin = margin(t = 1)),
        legend.text = element_text(size=6),
        legend.spacing.y = unit(-6,'pt'))+
  guides(shape=guide_legend(nrow=2))+
  coord_flip()
ggsave(plot=prop_in_ci,filename='prop_in_ci.tiff',width=3,height=3,units='in')

