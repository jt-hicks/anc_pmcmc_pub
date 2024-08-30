install.packages('ggbreak')
library(ggspatial)
library(ggbreak)

trials_short_informed_prepped <- readRDS('trials_short_informed_prepped.rds')
amoah_u5_short_results_prepped <- readRDS('amoah_u5_short_results_prepped_2.rds')
rogerson_prev_short_results_prepped_2 <- readRDS('rogerson_prev_short_results_prepped_2.RDS')

results_prepped_list <- list(trials_short_informed_prepped[[6]],amoah_u5_short_results_prepped,rogerson_prev_short_results_prepped_2)
malawi_sources <- c('Madanitsa, et al. 2016','Amoah, et al. 2021','Rogerson, et al. 2000')

names(results_prepped_list) <- malawi_sources

results_all_sample <- bind_rows(lapply(malawi_sources,function(x){
  df <- results_prepped_list[[x]]$sample
  df$source <- x
  return(df)
}))
results_all_summary <- bind_rows(lapply(malawi_sources,function(x){
  df <- results_prepped_list[[x]]$summary
  df$source <- x
  return(df)
}))

results_prev_summary <- results_all_summary%>%
  mutate(include = case_when(
    source == 'Madanitsa, et al. 2016' & measure %in% c('prev_pg','prev_mg') ~ "yes",
    source == 'Amoah, et al. 2021' & measure == c('prev_05') ~ "yes",
    source == 'Rogerson, et al. 2000' & measure == c('prev_anc_all') ~ "yes",
    TRUE ~ NA))%>%
  filter(include=='yes')
results_prev_sample <- results_all_sample%>%
  mutate(include = case_when(
    source == 'Madanitsa, et al. 2016' & measure %in% c('prev_pg','prev_mg') ~ "yes",
    source == 'Amoah, et al. 2021' & measure == c('prev_05') ~ "yes",
    source == 'Rogerson, et al. 2000' & measure == c('prev_anc_all') ~ "yes",
    TRUE ~ NA))%>%
  filter(include=='yes')
results_prev_sample$groups <- paste0(results_prev_sample$variable,results_prev_sample$measure)
results_eir_summary <- results_all_summary%>%
  filter(measure=='EIR')
results_eir_sample <- results_all_sample%>%
  filter(measure=='EIR')
results_eir_sample$groups <- paste0(results_eir_sample$variable,results_eir_sample$source)



amoah_prev_u5$source <- 'Amoah, et al. 2021'
amoah_prev_u5$measure <- 'prev_05'
rogerson_prev$source <- 'Rogerson, et al. 2000'
rogerson_prev$measure <- 'prev_anc_all'
rogerson_prev_cis <- addCIs(rogerson_prev,rogerson_prev$positive,rogerson_prev$tested)
trials_anc_data_all_cis <- addCIs(trials_anc_data_all,trials_anc_data_all$positive,trials_anc_data_all$tested)
trials_anc_data_all_cis$country <- trials_anc_data_all_cis$site
trials_anc_data_all_cis$date <- as.Date(trials_anc_data_all_cis$month,frac=0.5)
trials_anc_malawi <- trials_anc_data_all_cis%>%
  filter(country=='Malawi')%>%
  mutate(source='Madanitsa, et al. 2016')

malawi_observed_prev <- bind_rows(trials_anc_malawi,amoah_prev_u5,rogerson_prev_cis)
unique(malawi_observed_prev$measure)

tri_vir <- viridis::viridis(3,begin = 0.25,end=0.75)
colors_measures <- c(prev_pg=tri_vir[2],prev_mg=lighten(tri_vir[2],0.5),prev_05=tri_vir[3],prev_anc_all=tri_vir[1])
malawi_rainfall_prev <- trial_rainfall[trial_rainfall$country=='Malawi',]
malawi_rainfall_prev$rainfall_rel <- malawi_rainfall_prev$rainfall*0.8/max(malawi_rainfall_prev$rainfall)

multiplier <- 1
plot_malawi_prev <- ggplot()+
  geom_col(data=malawi_rainfall_prev,aes(x=date,y=rainfall_rel),alpha = 1,fill = 'darkgrey',just=0)+
  scale_y_continuous(expand=c(0,0),limits=c(0,1),sec.axis = sec_axis(~ . /(0.8/max(rainfall$rainfall)), name = "Monthly rainfall (cm)"))+
  geom_line(data=results_prev_sample,aes(x=as.Date(date),y=value*multiplier,color=measure,group=groups),alpha=0.1,linewidth=0.2)+
  geom_line(data=results_prev_summary,aes(x=as.Date(date),y=median*multiplier,color=measure,group=measure),linewidth=0.8)+
  geom_point(data=malawi_observed_prev,aes(x=as.Date(date),y=mean*multiplier,color=measure),pch = 19,position=position_dodge(width=10),size=0.5)+
  geom_errorbar(data=malawi_observed_prev,aes(x=as.Date(date),ymin=lower*multiplier,ymax=upper*multiplier,color=measure),width = 0,position=position_dodge(width=10),linewidth=0.5)+
  scale_color_manual(values=colors_measures)+
  # coord_cartesian(ylim = c(0,max_value),
  #                 xlim = range(results_summary$date))+
  # facet_grid(.~source,scales = 'free_x')+
  scale_x_break(c(as.Date('1998-07-15'),as.Date('2011-08-15')), space = 0.2) +
  scale_x_break(c(as.Date('2012-11-15'),as.Date('2016-04-15')), space = 0.2) +
  scale_x_date(date_labels = "%b '%y", date_breaks = '3 months',limits =c(min(results_prev_summary$date),as.Date('2017-11-15')) )+
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
ggsave('./trial_figs/output/plot_malawi_prev.tiff',plot=plot_malawi_prev,width=7,height=2,units='in')
malawi_rainfall <- trial_rainfall[trial_rainfall$country=='Malawi',]
malawi_rainfall$rainfall_rel <- malawi_rainfall$rainfall*500/max(malawi_rainfall$rainfall)
malawi_sources
colors_sources <- c(`Madanitsa, et al. 2016`=tri_vir[2],`Amoah, et al. 2021`=tri_vir[3],`Rogerson, et al. 2000`=tri_vir[1])

plot_malawi_eir <- ggplot()+
  geom_col(data=malawi_rainfall,aes(x=date,y=rainfall_rel),alpha = 1,fill = 'darkgrey',just=0)+
  scale_y_continuous(expand=c(0,0),limits=c(0,500),sec.axis = sec_axis(~ . /(rainfall_multiplier*500/max(rainfall$rainfall)), name = "Monthly rainfall (cm)"))+
  theme(strip.placement = "outside")+
  geom_line(data=results_eir_sample,aes(x=as.Date(date),y=value*multiplier,color=source,group=groups),alpha=0.1,linewidth=0.2)+
  geom_line(data=results_eir_summary,aes(x=as.Date(date),y=median*multiplier,color=source,group=source),linewidth=0.8)+
  scale_color_manual(values=colors_sources)+
  # facet_grid(.~source,scales = 'free_x')+
  coord_cartesian(ylim=c(0,500))+
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
ggsave('./trial_figs/output/plot_malawi_eir.tiff',plot=plot_malawi_eir,width=7,height=2,units='in')

malawi_timeline <- plot_malawi_prev / plot_malawi_eir
ggsave('./trial_figs/output/malawi_comparison_fits.tiff',width=6,height=5,units='in')

malawi_rainfall_rolling <- malawi_rainfall %>%
  mutate(season_year = ifelse(month_num<=7,year-1,year))%>%
  group_by(season_year) %>%
  filter(!season_year %in% c(1996,2019))%>%
  mutate(total_rainfall = sum(rainfall),
         rainfall_perc = rainfall / total_rainfall,
         rainfall_perc_2mo = rollapply(rainfall_perc,width=2,FUN=sum,align='left',fill=NA,partial=TRUE),
         rainfall_perc_3mo = rollapply(rainfall_perc,width=3,FUN=sum,align='left',fill=NA,partial=TRUE),
         rainfall_perc_4mo = rollapply(rainfall_perc,width=4,FUN=sum,align='left',fill=NA,partial=TRUE))

ggplot(malawi_rainfall_rolling)+
  geom_col(aes(x=date,y=rainfall_perc))+
  facet_wrap(season_year~.,scales='free_x')+
  scale_x_date(date_labels='%m',date_breaks = 'month')

one_month_malawi <- malawi_rainfall_rolling %>%
  filter(rainfall_perc>0.6)
two_month_malawi <- malawi_rainfall_rolling %>%
  filter(rainfall_perc_2mo>0.6)
three_month_malawi <- malawi_rainfall_rolling %>%
  filter(rainfall_perc_3mo>0.6)%>%
  filter(!season_year %in% c(1996,2019))%>%
  group_by(season_year)%>%
  filter(rainfall_perc_3mo==max(rainfall_perc_3mo))
table(three_month_malawi$season_year)
four_month_malawi <- malawi_rainfall_rolling %>%
  filter(rainfall_perc_4mo>0.6)%>%
  filter(!season_year %in% c(1996,2019))%>%
  group_by(season_year)%>%
  filter(rainfall_perc_4mo==max(rainfall_perc_4mo))
table(four_month_malawi$season_year)
# Applying the function to the dataset
result <- find_consecutive_months(malawi_rainfall)

##Malawi Maps
malawi_shape<-terra::vect('C:/Users/jthicks/Documents/africa_shape_files/gadm41_MWI_shp/gadm41_MWI_2.shp')

malawi_border <- terra::aggregate(malawi_shape,by='COUNTRY')
malawi_admin1 <- terra::aggregate(malawi_shape,by='NAME_1')
malawi_admin1$interest <- ifelse(malawi_admin1$NAME_1 %in% c('Blantyre,Chikwawa'),'Yes','No')
blant_chik_shape <- subset(malawi_shape,malawi_shape$NAME_1 %in% c('Blantyre','Chikwawa'))
blant_chik_shape$interest <- 'yes'
blant_chik_shape <- terra::aggregate(blant_chik_shape,by='interest')
blant_shape <- subset(malawi_shape,malawi_shape$NAME_1 %in% c('Blantyre'))
blant_shape <- terra::aggregate(blant_shape,by='NAME_1')

chik_shape <- subset(malawi_shape,malawi_shape$NAME_1 %in% c('Chikwawa'))
chik_shape <- terra::aggregate(chik_shape,by='NAME_1')


malawi_sites <- data.frame(source=c(rep(malawi_sources[1],3),rep(malawi_sources[2],3),malawi_sources[3]),
                           latitude = c(-15.874831,-15.967982,-16.023803671884767,-15.830859,-15.935394,-15.987641,-15.803091),
                           longitude = c(34.956554,34.908439,34.791939287730976,34.483789,34.759846,34.481616,35.021586))
coordinates_sf <- st_as_sf(malawi_sites, coords = c("longitude", "latitude"), crs = 4326)
joined_data <- st_join(coordinates_sf, malawi_shape, join = st_within)

malawi_map <- ggplot() +
  geom_sf(data=malawi_border,color='darkgrey')+
  annotate("rect", xmin = 34, xmax = 35.4, ymin = -17, ymax = -15,alpha = .6,fill = '#C3A995FF')+
  # geom_sf(data=included_countries,fill='darkgrey')+
  # geom_sf(data = centroids_sf, aes(fill = country), size = 2,shape=21,color='black') +
  # scale_fill_manual(values=trial_color_palette,labels=labels_trial_colors,guide='none')+
  # scale_fill_manual(values=included_palette, guide='none')+
  theme(legend.title = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        panel.spacing.y = unit(3, "mm"))
malawi_sites_map <- ggplot() +
  geom_sf(data=malawi_border,color='darkgrey')+
  geom_sf(data=blant_shape,fill='darkgrey')+
  geom_sf(data=chik_shape,fill='darkgrey')+
  geom_sf(data = coordinates_sf, aes(fill = source), size = 1.5,shape=21,color='black')+
  # geom_sf(data=included_countries,fill='darkgrey')+
  # geom_sf(data = centroids_sf, aes(fill = country), size = 2,shape=21,color='black') +
  # scale_fill_manual(values=trial_color_palette,labels=labels_trial_colors,guide='none')+
  annotation_scale(location = "br", width_hint = 0.3,bar_cols = c("darkgrey", "white")) +
  coord_sf(xlim=c(34,35.4),ylim=c(-15,-17))+
  scale_fill_manual(values=colors_sources,guide='none')+
  theme(legend.title = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        panel.spacing.y = unit(3, "mm"))

mal_map_winsert <- malawi_sites_map + inset_element(malawi_map, left = -.05, bottom = 0, right = 0.26,top = 0.41)
ggsave('./trial_figs/output/mal_map_winsert.tiff',plot=mal_map_winsert,width=3,height=3,units='in')

prev_plus_smc <- plot_malawi_prev / plot_malawi_smc
timeseries_smc <- ggarrange(plot_facet_timeseries,prev_plus_smc,heights=c(1,1.5),ncol = 1)
map_plus_smc <- ggarrange(mal_map_winsert,timeseries_smc,widths=c(1,2))
ggsave('./trial_figs/output/malawimap_plus_smc.tiff',plot=map_plus_smc,width=7,height=5,units='in')
