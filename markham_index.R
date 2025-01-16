library(parallel)
library(scales)
library(forcats)
source('./utils.R')
theme_set(theme_minimal(base_size = 8)+
            theme(text = element_text(size=6),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  axis.title.x = element_text(size=6,margin = margin(t = -1)),
                  axis.title.y = element_text(size=6),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom',
                  legend.box.spacing = unit(0,'pt'),
                  legend.key.spacing.y = unit(-8, 'pt'),
                  legend.margin = margin(t=0,r=0,b=0,l=0),
                  legend.box.margin = margin(t=0,r=0,b=0,l=0),
                  legend.text = element_text(size=7,margin=margin(0,0,0,0)),
                  plot.margin = margin(t=4,l=4,r=0,b=0))
          )
theme_set(theme_minimal(base_size = 8)+
            theme(text = element_text(size=6),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  axis.title.x = element_text(size=6,margin = margin(t = -2)),
                  axis.title.y = element_text(size=6),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom',
                  legend.box.spacing = unit(0,'pt'),
                  legend.key.spacing.y = unit(-8, 'pt'),
                  legend.margin = margin(t=0,r=0,b=0,l=0),
                  legend.box.margin = margin(t=0,r=0,b=0,l=0),
                  legend.text = element_text(size=7,margin=margin(0,0,0,0)),
                  plot.margin = margin(t=4,l=4,r=0,b=0))
)

#Markhan Seasonality Index
month_arc <- data.frame(date = seq.Date(from=as.Date('2010-01-01'),to=as.Date('2010-12-31'),by='days'),
                        degree = 0:364)%>%
  mutate(radians = degree*pi/180,
         mid_month = as.character(as.Date(as.yearmon(date),frac=0.5)),
         month=month(date))%>%
  filter(date==as.Date(mid_month))%>%
  select(month,degree,radians)

##Import and format comparison data
cairns_angle_observed <- readxl::read_excel('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/msi_angles_cairns.xlsx',
                                            sheet='Sheet1')

siaya_total <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/white_paper_siaya.txt')
siaya_total$date <- as.Date(date_decimal(siaya_total$date_decimal))
siaya_total$year <- year(siaya_total$date)
siaya_total$month <- month(siaya_total$date)
siaya_total$value <- siaya_total$cases_total

#Roca-Feltrer, et al. 2009 https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-8-276
blantyre_total <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/blantyre_hosp.txt')
blantyre_total$date <- as.Date(date_decimal(blantyre_total$date_decimal))
blantyre_total$year <- year(blantyre_total$date)
blantyre_total$month <- month(blantyre_total$date)
blantyre_total$value <- blantyre_total$cases_total
blantyre_total$location <- 'Blantyre'

#Roca-Feltrer, et al. 2009 https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-8-276
cereb_total <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/blantyre_cereb.txt')
cereb_total$date <- as.Date(date_decimal(cereb_total$date_decimal))
cereb_total$year <- year(cereb_total$date)
cereb_total$month <- month(cereb_total$date)
cereb_total$location <- 'Blantyre'

#Tizifia, et al 2021 https://malariajournal.biomedcentral.com/articles/10.1186/s12936-021-04013-5
chikwawa_total <- readxl::read_excel('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/chikwawa incidence.xlsx',
                                     sheet='Sheet1')%>%
  filter(month!=as.Date('2019-05-01'))

chikwawa_total$date <- as.Date(chikwawa_total$month)
chikwawa_total$year <- year(chikwawa_total$date)
chikwawa_total$month <- month(chikwawa_total$date)
chikwawa_total$value <- chikwawa_total$incidence
chikwawa_total$source <- 'Tizifia - Cohort Clinical Incidence'
chikwawa_total$location <- 'Chikwawa'
chikwawa_total_ave <- chikwawa_total%>%
  group_by(month)%>%
  summarise(value=mean(value))%>%
  mutate(date=as.Date(as.yearmon(paste0('2010-',month))))
cereb_total$source <- 'Roca-Feltrer - Cerebral malaria'
blantyre_total$source <- 'Roca-Feltrer - All severe malaria'
chikwawa_total_ave$source <- 'Tizifia - Cohort Clinical Incidence'


basse_total <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/gambia_incidence.txt')
basse_total$date <- as.Date(date_decimal(basse_total$date_decimal))
basse_total$year <- year(basse_total$date)
basse_total$month <- month(basse_total$date)
basse_total <- bind_rows(data.frame(cases=0,
                                    year=2011,
                                    month=c(1:6),
                                    date = seq.Date(as.Date('2011-01-01'),as.Date('2011-06-01'),by='month')),
                         basse_total)
basse_total$value <- basse_total$cases

##Calculate MSI and vectors for observed data
siaya_msi <- return_msi(dataframe=siaya_total) %>%
  mutate(x_coord = 0 + rk * cos(theta_k),
         y_coord = 0 + rk * sin(theta_k),
         index=1)%>%
  bind_rows(data.frame(year=unique(siaya_msi$year),
                       index=0,
                       x_coord=0,
                       y_coord=0))%>%
  arrange(year,index)
siaya_vectors <- return_vectors(dataframe = siaya_total)%>%
  mutate(xstart = lag(x_coord),ystart = lag(y_coord))%>%
  na.omit()

ggplot()+
  geom_segment(data=siaya_vectors,aes(x=xstart,y=ystart,xend=x_coord,yend=y_coord),size=1,arrow = arrow(type='closed',length = unit(0.1,'inches')))+
  geom_path(data=siaya_msi,aes(x=x_coord,y=y_coord),size=1.5,color='darkred',arrow = arrow(type='closed',length = unit(0.1,'inches')))+
  facet_wrap(.~year)
ggplot()+
  geom_segment(data=siaya_vectors,aes(x=ystart,y=xstart,xend=y_coord,yend=x_coord),size=1,arrow = arrow(type='closed',length = unit(0.1,'inches')))+
  geom_path(data=siaya_msi,aes(x=y_coord,y=x_coord),size=1.5,color='darkred',arrow = arrow(type='closed',length = unit(0.1,'inches')))+
  facet_wrap(.~year)

siaya_total %>%
  group_by(year, month) %>%
  summarise(max_value = max(value), .groups = 'drop') %>%
  group_by(year) %>%
  slice(which.max(max_value)) %>%
  ungroup()
mean(siaya_msi$msi)
sample_test <- trials_short_informed_prepped[[1]]$sample

blantyre_msi <- return_msi(dataframe=blantyre_total,single_year = TRUE) %>%
  mutate(x_coord = 0 + rk * cos(theta_k),
         y_coord = 0 + rk * sin(theta_k))
blantyre_msi$source <- 'Roca-Feltrer - All severe malaria'
blantyre_vectors <- return_vectors(dataframe = blantyre_total,single_year = TRUE)%>%
  mutate(xstart = lag(x_coord),ystart = lag(y_coord))%>%
  na.omit()
ggplot()+
  geom_segment(data=blantyre_vectors,aes(x=ystart,y=xstart,xend=y_coord,yend=x_coord),size=1,arrow = arrow(type='closed',length = unit(0.1,'inches')))+
  geom_segment(data=blantyre_msi,aes(x=0,y=0,xend=y_coord,yend=x_coord),size=1.5,color='darkred',arrow = arrow(type='closed',length = unit(0.1,'inches')))
basse_msi <- return_msi(dataframe=basse_total,single_year = TRUE) %>%
  mutate(x_coord = 0 + rk * cos(theta_k),
         y_coord = 0 + rk * sin(theta_k))
basse_vectors <- return_vectors(dataframe = basse_total,single_year = TRUE)%>%
  mutate(xstart = lag(x_coord),ystart = lag(y_coord))%>%
  na.omit()
ggplot()+
  geom_segment(data=basse_vectors,aes(x=ystart,y=xstart,xend=y_coord,yend=x_coord),size=1,arrow = arrow(type='closed',length = unit(0.1,'inches')))+
  geom_segment(data=basse_msi,aes(x=0,y=0,xend=y_coord,yend=x_coord),size=1.5,color='darkred',arrow = arrow(type='closed',length = unit(0.1,'inches')))

cereb_msi <- return_msi(dataframe=cereb_total,single_year = TRUE) %>%
  mutate(x_coord = 0 + rk * cos(theta_k),
         y_coord = 0 + rk * sin(theta_k))
cereb_msi$source <- 'Roca-Feltrer - Cerebral malaria'
cereb_vectors <- return_vectors(dataframe = cereb_total,single_year = TRUE)%>%
  mutate(xstart = lag(x_coord),ystart = lag(y_coord))%>%
  na.omit()
ggplot()+
  geom_segment(data=cereb_vectors,aes(x=ystart,y=xstart,xend=y_coord,yend=x_coord),size=1,arrow = arrow(type='closed',length = unit(0.1,'inches')))+
  geom_segment(data=cereb_msi,aes(x=0,y=0,xend=y_coord,yend=x_coord),size=1.5,color='darkred',arrow = arrow(type='closed',length = unit(0.1,'inches')))

nrow(chikwawa_total)
chikwawa_msi <- return_msi(dataframe=chikwawa_total_ave,single_year = TRUE) %>%
  mutate(x_coord = 0 + rk * cos(theta_k),
         y_coord = 0 + rk * sin(theta_k))
chikwawa_msi$source <- 'Tizifia - Cohort Clinical Incidence'
chikwawa_vectors <- return_vectors(dataframe = chikwawa_total_ave,single_year = TRUE)%>%
  mutate(xstart = lag(x_coord),ystart = lag(y_coord))%>%
  na.omit()
ggplot()+
  geom_segment(data=chikwawa_vectors,aes(x=ystart,y=xstart,xend=y_coord,yend=x_coord),size=1,arrow = arrow(type='closed',length = unit(0.1,'inches')))+
  geom_segment(data=chikwawa_msi,aes(x=0,y=0,xend=y_coord,yend=x_coord),size=1.5,color='darkred',arrow = arrow(type='closed',length = unit(0.1,'inches')))

##Get ANC-estimated incidence
incidence <- trials_short_informed_sample[trials_short_informed_sample$measure=='inc05',]%>%
  filter(!(country=='Gambia'&time==13))
##Gambia incidence to plot Gambia time series
gambia_incidence <- incidence[incidence$country=='Gambia',]

trial_color_palette <- c(viridis::viridis(6,begin=0,end=0.85))


show_col(trial_color_palette)
,'#8D5B4CFF'
names(trial_color_palette) <-  c('Gambia','Mali','Burkina Faso','Ghana','Malawi','Kenya')
labels_trial_colors <- c('The Gambia','Mali','Burkina Faso','Ghana','Malawi','Kenya')
names(labels_trial_colors) <- names(trial_color_palette)

sample_msi_list <- incidence %>%
  group_by(country,variable) %>%
  group_split() %>%
  lapply(process_group)
sample_msi <- bind_rows(sample_msi_list)

sample_msi_summary <- sample_msi%>%
  group_by(country)%>%
  summarise(msi_median=median(msi),
            msi_lower=quantile(msi,0.025),
            msi_upper=quantile(msi,0.975),
            angle_median=median(theta_k_deg),
            angle_lower=quantile(theta_k_deg,0.025),
            angle_upper=quantile(theta_k_deg,0.975))

msi_observed <- data.frame(country=c('Ghana','Burkina Faso', 'Mali', 'Gambia','Kenya','Malawi'),
                           msi=c(0.545,0.742,0.823,basse_msi$msi,0.111,cereb_msi$msi),
                           angle=c(cairns_angle_observed[cairns_angle_observed$country=='Ghana',]$angle,
                                   cairns_angle_observed[cairns_angle_observed$country=='Burkina Faso',]$angle,
                                   cairns_angle_observed[cairns_angle_observed$country=='Mali','angle',]$angle,
                                   basse_msi$theta_k_deg,
                                   NA,
                                   cereb_msi$theta_k_deg))
msi_both <- merge(sample_msi_summary,msi_observed,by='country') %>%
  mutate(angle_date_median=as.Date(ifelse(msi_median>=0.25,as.Date(date_decimal(2010+angle_median/360)),as.Date('2009-12-15'))),
         angle_date_lower=as.Date(ifelse(msi_median>=0.25,as.Date(date_decimal(2010+angle_lower/360)),as.Date('2009-12-15'))),
         angle_date_upper=as.Date(ifelse(msi_median>=0.25,as.Date(date_decimal(2010+angle_upper/360)),as.Date('2009-12-15'))),
         angle_date=as.Date(ifelse(msi>=0.25,as.Date(date_decimal(2010+angle/360)),as.Date('2009-12-15'))))


msi_plot <- ggplot(msi_both)+
  geom_point(aes(x=msi,y=msi_median,color=country))+
  geom_errorbar(aes(x=msi,ymin=msi_lower,ymax=msi_upper,color=country),width=0)+
  geom_abline(linetype='dashed')+
  # geom_text(aes(x=msi_obs,y=median,label=country),hjust = 0, nudge_x = 0.02)+
  scale_color_manual(name='Country',values=trial_color_palette,labels=labels_trial_colors,breaks=names(trial_color_palette))+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  labs(x='Observed Seasonality',y='Model Seasonality')+
  # labs(x='Observed Incidence* MSI',y='Estimated Incidence MSI')+
  theme(legend.title = element_blank(),
        legend.position = 'none')
ggsave('./trial_figs/output/msi_final_plot.tiff',plot=msi_plot,width=3.5,height=3.5,units='in')
date_breaks <- seq.Date(from=as.Date('2009-12-15'),as.Date('2010-12-15'),by='month')
date_labels <- c('No Peak','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
msi_angle_plot <- ggplot(msi_both)+
  annotate("rect", xmin = as.Date('2009-12-01'), xmax = as.Date('2009-12-31'), ymin = as.Date('2009-12-01'), ymax = as.Date('2010-12-15'),alpha = .4,fill = '#C3A995FF')+
  annotate("rect", ymin = as.Date('2009-12-01'), ymax = as.Date('2009-12-31'), xmin = as.Date('2009-12-01'), xmax = as.Date('2010-12-15'),alpha = .4,fill = '#C3A995FF')+
  geom_point(aes(x=angle_date,y=angle_date_median,color=country))+
  geom_errorbar(aes(x=angle_date,ymin=angle_date_lower,ymax=angle_date_upper,color=country),width=0)+
  geom_abline(linetype='dashed')+
  # geom_text(aes(x=msi_obs,y=median,label=country),hjust = 0, nudge_x = 7)+
  scale_color_manual(name='Country',values=trial_color_palette,labels=labels_trial_colors,breaks=names(trial_color_palette))+
  scale_x_date(limits = as.Date(c('2009-12-01','2010-12-15')),labels = date_labels,breaks = date_breaks)+
  scale_y_date(limits = as.Date(c('2009-12-01','2010-12-15')),labels = date_labels,breaks = date_breaks)+
  labs(x='Observed Incidence Peak',y='Estimated Incidence Peak')+
  theme(legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45,vjust=1,hjust=1))
ggsave('./trial_figs/output/peak_plot_final.tiff',plot=msi_angle_plot,width=3.5,height=3.5,units='in')
ggsave('./trial_figs/output/peak_plot_labeled.tiff',plot=msi_angle_plot,width=7,height=7,units='in')

##Rain-based MSI
##Get rainfall for sites
##Download rainfall data
source('./trial_figs/rainfall_data.R')
devtools::install_github("mrc-ide/umbrella")
rainfall_folder<-"C:/Users/jthicks/Documents/rainfall_data_chirps_2000-2019"
unique_months <- data.frame(date=seq.Date(as.Date('1997-01-01'),as.Date('2023-12-01'),by='month'))%>%
  mutate(month=month(date),
         year=year(date))
sapply(seq_along(unique_months$month), function(x){
  get_chirps_monthly(year=unique_months$year[x],month=unique_months$month[x],folder=rainfall_folder)
})
earlier_months <- data.frame(date=seq.Date(as.Date('1997-01-01'),as.Date('1999-12-01'),by='month'))%>%
  mutate(month=month(date),
         year=year(date))
sapply(seq_along(earlier_months$month), function(x){
  get_chirps_monthly(year=earlier_months$year[x],month=earlier_months$month[x],folder=rainfall_folder)
})
my_raster_files<-list.files(rainfall_folder)
my_raster <- terra::rast(paste(rainfall_folder,my_raster_files,sep="/"))

##Get site shape files
ghana_shape <- terra::vect('C:/Users/jthicks/Documents/africa_shape_files/gadm41_GHA_shp/gadm41_GHA_2.shp')
burkinafaso_shape <-terra::vect('C:/Users/jthicks/Documents/africa_shape_files/gadm41_BFA_shp/gadm41_BFA_2.shp')
mali_shape<-terra::vect('C:/Users/jthicks/Documents/africa_shape_files/gadm41_MLI_shp/gadm41_MLI_2.shp')
gambia_shape<-terra::vect('C:/Users/jthicks/Documents/africa_shape_files/gadm41_GMB_shp/gadm41_GMB_2.shp')
kenya_shape<-terra::vect('C:/Users/jthicks/Documents/africa_shape_files/gadm41_KEN_shp/gadm41_KEN_2.shp')
malawi_shape<-terra::vect('C:/Users/jthicks/Documents/africa_shape_files/gadm41_MWI_shp/gadm41_MWI_2.shp')

##Subet shape files to only districts of interest
kassena_nankana_shape <- subset(ghana_shape,ghana_shape$NAME_2=='Kasena Nankana East')
oubritenga_shape <- subset(burkinafaso_shape,burkinafaso_shape$NAME_2=='Oubritenga')
kati_san_shape <- subset(mali_shape,mali_shape$NAME_2 %in% c('Kati','San'))
kati_san_shape$interest <- 'yes'
kati_san_shape <- terra::aggregate(kati_san_shape,by='interest')
fulladu_east_shape <- subset(gambia_shape,gambia_shape$NAME_2=='Fulladu East')
siaya_shape <- subset(kenya_shape,kenya_shape$NAME_1=='Siaya')
siaya_shape <- terra::aggregate(siaya_shape,by='NAME_1')
blant_chik_shape <- subset(malawi_shape,malawi_shape$NAME_1 %in% c('Blantyre','Chikwawa'))
blant_chik_shape$interest <- 'yes'
blant_chik_shape <- terra::aggregate(blant_chik_shape,by='interest')
shape_list <- list(kassena_nankana_shape,oubritenga_shape,kati_san_shape,fulladu_east_shape,
                   siaya_shape,blant_chik_shape)
country_list_shapes <- c('Ghana','Burkina Faso','Mali','Gambia','Kenya','Malawi')
trial_rainfall <- bind_rows(lapply(1:6,function(x){
  return_rainfall_df(shape=shape_list[[x]],
                     country=country_list_shapes[x])
}))
trial_rainfall$year <- year(trial_rainfall$month)
trial_rainfall$month_num <- month(trial_rainfall$month)

trial_rainfall$date <- as.Date(trial_rainfall$month)
trial_rainfall$value <- trial_rainfall$rainfall
trial_rainfall$variable <- 1

ggplot(trial_rainfall)+
  geom_line(aes(x=month,y=rainfall))+
  facet_grid(country~.)
ggplot(trial_rainfall[trial_rainfall$country=='Gambia',])+
  geom_line(aes(x=month_num,y=rainfall))+
  facet_grid(year~.)
malawi_rainfall <- trial_rainfall[trial_rainfall$country=='Malawi',]%>%
  mutate(month_car=as.character(month),
         season_year = ifelse(month_num<=7,year,year-1),
         seas_month_num = ifelse(month_num<=6,month_num+6,month_num-6),
         rainfall_rel = rainfall*1/max(rainfall))
ggplot(trial_rainfall[trial_rainfall$country=='Malawi'&trial_rainfall$year>=2008,]%>%
         mutate(month_car=as.character(month),
                season_year = ifelse(month_num<=7,year,year-1),
                seas_month_num = ifelse(month_num<=6,month_num+6,month_num-6))%>%
         filter(season_year>=2009,season_year<2019))+
  geom_line(aes(x=seas_month_num,y=rainfall),linewidth=1)+
  facet_grid(season_year~.)+
  scale_x_continuous(breaks=1:12,labels=c('Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'))

ggsave('./trial_figs/output/malawi_rain_seasons_plot.tiff',width=3.5,height=7,units='in')
trial_dates <- incidence %>%
  group_by(country,date)%>%
  summarise(value=1)%>%
  select(!value)%>%
  mutate(month_car=as.character(as.yearmon(date)))%>%
  select(!date)
trial_rainfall_subset <-trial_rainfall %>%
  mutate(month_car=as.character(month))%>%
  right_join(trial_dates,by=c('country','month_car'))
rain_msi <- bind_rows(trial_rainfall %>%
                        mutate(month_car=as.character(month))%>%
                        right_join(trial_dates,by=c('country','month_car'))%>%
                        group_by(country) %>%
                        group_split() %>%
                        lapply(process_group))
malawi_rain_msi <- bind_rows(trial_rainfall[trial_rainfall$country=='Malawi'&trial_rainfall$year>=2008,] %>%
                               mutate(month_car=as.character(month),
                                      season_year = ifelse(month_num<=6,year,year-1))%>%
                               filter(season_year>=2009 & season_year<2019)%>%
                               group_by(season_year) %>%
                               group_split() %>%
                               lapply(process_group))%>%
  mutate(angle_date=as.Date(ifelse(msi>=0.25,as.Date(date_decimal(2010+theta_k_deg/360)),as.Date('2009-12-15'))),
         year=2009:2018)

ggplot(malawi_rain_msi)

msi_both_rain <- merge(sample_msi_summary,rain_msi,by='country') %>%
  mutate(angle_date_median=as.Date(ifelse(msi_median>=0.25,as.Date(date_decimal(2010+angle_median/360)),as.Date('2009-12-15'))),
         angle_date_lower=as.Date(ifelse(msi_median>=0.25,as.Date(date_decimal(2010+angle_lower/360)),as.Date('2009-12-15'))),
         angle_date_upper=as.Date(ifelse(msi_median>=0.25,as.Date(date_decimal(2010+angle_upper/360)),as.Date('2009-12-15'))),
         angle_date=as.Date(ifelse(msi>=0.25,as.Date(date_decimal(2010+theta_k_deg/360)),as.Date('2009-12-15'))))

msi_rain_plot <- ggplot(msi_both_rain)+
  geom_point(aes(x=msi,y=msi_median,color=country))+
  geom_errorbar(aes(x=msi,ymin=msi_lower,ymax=msi_upper,color=country),width=0)+
  geom_abline(linetype='dashed')+
  # geom_text(aes(x=msi_obs,y=median,label=country),hjust = 0, nudge_x = 0.02)+
  scale_color_manual(name='Country',values=trial_color_palette,labels=labels_trial_colors,breaks=names(trial_color_palette))+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  # labs(x='Rainfall MSI',y='Estimated Incidence MSI')+
  labs(x='Rainfall Seasonality',y='Model Seasonality')+
  theme(legend.title = element_blank(),
        legend.position = 'none')
ggsave('./trial_figs/output/msi_rain_plot.tiff',plot=msi_rain_plot,width=3.5,height=3.5,units='in')

msi_angle_rain_plot <- ggplot(msi_both_rain)+
  annotate("rect", xmin = as.Date('2009-12-01'), xmax = as.Date('2009-12-31'), ymin = as.Date('2009-12-01'), ymax = as.Date('2010-12-15'),alpha = .4,fill = '#C3A995FF')+
  annotate("rect", ymin = as.Date('2009-12-01'), ymax = as.Date('2009-12-31'), xmin = as.Date('2009-12-01'), xmax = as.Date('2010-12-15'),alpha = .4,fill = '#C3A995FF')+
  geom_point(aes(x=angle_date,y=angle_date_median,color=country))+
  geom_errorbar(aes(x=angle_date,ymin=angle_date_lower,ymax=angle_date_upper,color=country),width=0)+
  geom_abline(linetype='dashed')+
  # geom_text(aes(x=msi_obs,y=median,label=country),hjust = 0, nudge_x = 7)+
  scale_color_manual(name='Country',values=trial_color_palette,labels=labels_trial_colors,breaks=names(trial_color_palette))+
  scale_x_date(limits = as.Date(c('2009-12-01','2010-12-15')),labels = date_labels,breaks = date_breaks)+
  scale_y_date(limits = as.Date(c('2009-12-01','2010-12-15')),labels = date_labels,breaks = date_breaks)+
  labs(x='Rainfall Peak',y='Estimated Incidence Peak')+
  theme(legend.title=element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45,vjust=1,hjust=1))
ggsave('./trial_figs/output/peak_plot_rainfall.tiff',plot=msi_angle_rain_plot,width=3.5,height=3.5,units='in')
msi_plots_combo <- msi_plot + msi_angle_plot + msi_rain_plot + msi_angle_rain_plot + plot_layout(ncol = 2, guides = 'collect')
msi_plots_norain <- msi_plot + msi_angle_plot  + plot_layout(nrow = 2, guides = 'collect')

ggsave('./trial_figs/output/msi_rain_combos.tiff',plot=msi_plots_combo,width=4,height=4,units='in')
ggsave('./msi_rain_combos_pres.tiff',plot=msi_plots_combo,width=3.2,height=3.2,units='in')

msi_comp_plusinc <- ggarrange(trials_incidence_plot,msi_plots_combo,ncol=2,widths = c(1, 2))
ggsave('./trial_figs/output/msi_rain_combos_withinc.tiff',plot=msi_comp_plusinc,width=6,height=3.5,units='in')

msi_comp_plusprevinc <- ggarrange(trials_prev_anc_plot_2levels,trials_incidence_plot,msi_plots_norain,ncol=3,widths = c(1.1,1.1,1))
ggsave('./trial_figs/output/msi_withprev&inc.tiff',plot=msi_comp_plusprevinc,width=5.5,height=2.75,units='in')

trials_incidence_plot | ((msi_plot + msi_angle_plot)/(msi_rain_plot + msi_angle_rain_plot))
trials_incidence_plot | ((msi_plot + msi_angle_plot)/(msi_rain_plot + msi_angle_rain_plot))
##Peak Proportion
obs_list <- list(Gambia=basse_total,Kenya=siaya_total,Malawi=chikwawa_total)
obs_list_malaw <- list(Gambia=basse_total,Kenya=siaya_total,Malawi=cereb_total,Malawi=blantyre_total,Malawi=chikwawa_total)

peak_obs <- bind_rows(lapply(c('Gambia','Kenya','Malawi'),function(country){
  prop_max_3mo <- calc_case_conc(obs_list[[country]],3)
  prop_max_3mo$period <- '3 months'
  prop_max_4mo <- calc_case_conc(obs_list[[country]],4)
  prop_max_4mo$period <- '4 months'
  prop_max_both <- bind_rows(prop_max_3mo,prop_max_4mo)
  prop_max_both$country <- country
  return(prop_max_both)
}))
color_list <- c('black','black','darkred','darkorange','gold')
country_list <- c('Gambia','Kenya','Malawi','Malawi','Malawi')
peak_obs_malaw <- bind_rows(lapply(1:5,function(x){
  prop_max_3mo <- calc_case_conc(obs_list_malaw[[x]],3)
  prop_max_3mo$period <- '3 months'
  prop_max_3mo$country <- country_list[x]
  prop_max_3mo$color <- color_list[x]
  return(prop_max_3mo)
}))
sample_peak_list <- incidence %>%
  group_by(country,variable) %>%
  group_split() %>%
  lapply(process_group_peak)
sample_peak <- bind_rows(sample_peak_list)
sample_peak_summary <- sample_peak%>%
  group_by(country,period)%>%
  summarise(median=median(prop),
            lower=quantile(prop,0.025),
            upper=quantile(prop,0.975))
peak_both <- merge(sample_peak_summary,peak_obs,by=c('country','period'))
peak_both$color <- ifelse(peak_both$country=='Malawi','darkred','black')
peak_plot <- ggplot(peak_both[peak_both$period=='3 months'&peak_both$country!='Gambia',])+
  geom_point(aes(x=prop,y=median,color=color))+
  geom_errorbar(aes(x=prop,ymin=lower,ymax=upper,color=color),width=0)+
  geom_abline(linetype='dashed')+
  scale_color_manual(values=named_color_list,guide='none')+
  # geom_text(aes(x=prop,y=median,label=country),hjust = 0, nudge_x = 0.02)+
  # facet_wrap(.~period)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  labs(x='Observed proportion',y='Estimated proportion')
ggsave('./trial_figs/output/peak_prop_final_plot.tiff',plot=peak_plot,width=14,height=7,units='in')
ggsave('./trial_figs/output/peak_prop_plot_labeled.tiff',plot=peak_plot,width=14,height=7,units='in')
peak_both_malaw <- merge(sample_peak_summary,peak_obs_malaw,by=c('country','period'))
peak_malaw_plot <- ggplot(peak_both_malaw)+
  geom_point(aes(x=prop,y=median,color=color))+
  geom_errorbar(aes(x=prop,ymin=lower,ymax=upper,color=color))+
  geom_abline(linetype='dashed')+
  # geom_text(aes(x=prop,y=median,label=country),hjust = 0, nudge_x = 0.02)+
  scale_color_manual(values=named_color_list,guide='none')+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  labs(x='Observed proportion',y='Estimated proportion')
ggsave('./trial_figs/output/peak_prop_plot-malaw.tiff',plot=peak_malaw_plot,width=3.5,height=3.5,units='in')
malaw_combo_plot <- msi_plot + msi_angle_plot + peak_malaw_plot + plot_layout(ncol = 3)
ggsave('./trial_figs/output/malaw-combo.tiff',plot=malaw_combo_plot,width=11,height=3.5,units='in')
rain_peak <- bind_rows(trial_rainfall %>%
                         mutate(month_car=as.character(month))%>%
                         right_join(trial_dates,by=c('country','month_car'))%>%
                         group_by(country) %>%
                         group_split() %>%
                         lapply(process_group_peak))

peak_both_rain <- merge(sample_peak_summary,rain_peak,by=c('country','period'))
peak_rain_plot <- ggplot(peak_both_rain[peak_both_rain$period=='3 months',])+
  geom_point(aes(x=prop,y=median))+
  geom_errorbar(aes(x=prop,ymin=lower,ymax=upper),width=0)+
  geom_abline(linetype='dashed')+
  geom_text(aes(x=prop,y=median,label=country),hjust = 0, nudge_x = 0.02)+
  # facet_wrap(.~period)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  labs(x='Rainfall proportion',y='Estimated proportion')
ggsave('./trial_figs/output/rainfall_prop_compare.tiff',plot=peak_rain_plot,width=3.5,height=3.5,units='in')
all_combo_plot <- msi_plot + msi_angle_plot + peak_rain_plot + plot_layout(ncol = 3)
ggsave('./trial_figs/output/all_seasonality_comparisons.tiff',plot=all_combo_plot,width=11,height=3.5,units='in')

##Malawi Comparison
amoah_eir <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/Amoah_elife_2021.txt')
amoah_eir$date <- as.Date(date_decimal(amoah_eir$date_decimal))
amoah_eir$year <- year(amoah_eir$date)
amoah_eir$month <- month(amoah_eir$date)
amoah_eir$eir <- ifelse(amoah_eir$eir<0,0,amoah_eir$eir)
amoah_eir$value <- amoah_eir$eir
amoah_eir$source <- 'Amoah - EIR'
amoah_eir$site <- 'Malawi'
amoah_eir_ave <- amoah_eir%>%
  group_by(month)%>%
  summarise(value=mean(value))%>%
  mutate(date=as.Date(as.yearmon(paste0('2010-',month))))
amoah_eir_ave$source <- 'Amoah - EIR'
amoah_eir$location <- 'Chikwawa'

amoah_ento <- readxl::read_excel('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/Amoah_ento_data.xlsx',
                                            sheet='ento_data')
table(amoah_ento$time.mod)
month_match <- data.frame(t=1:38,
                          date=seq.Date(from=as.Date('2015-04-15'),length.out=38,by='month'))
amoah_ento_sum <- amoah_ento %>%
  mutate(moz_count = arabss_female + funss_female)%>%
  left_join(month_match,by=join_by(time.mod==t))%>%
  group_by(date)%>%
  summarise(moz_count=sum(moz_count),
            arabss_female = sum(arabss_female),
            funss_female = sum(funss_female))

amoah_prev <- readxl::read_excel('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/Amoah_parasitaemia_data.xlsx',
                                 sheet='parasitaemia_data')
amoah_prev$age.years <- amoah_prev$age.months/12
hist(amoah_prev[amoah_prev$age.years<10,]$age.years)
max(amoah_prev[amoah_prev$age.years<10,]$age.years)

amoah_prev_u5 <- amoah_prev %>%
  filter(age.years<=5)%>%
  group_by(time.mod)%>%
  summarise(tested = n(),
            positive = sum(RDT.result))%>%
  left_join(month_match,by=join_by(time.mod==t))
amoah_prev_u5 <- addCIs(amoah_prev_u5,Ys=amoah_prev_u5$positive,Ns=amoah_prev_u5$tested)
amoah_prev_u5$month <- as.yearmon(amoah_prev_u5$date)
str(amoah_prev_u5)
windows_authenticate()


amoah_u5_short <- task_create_expr({
  mamasante::run_pmcmc(data_raw=amoah_prev_u5,
                       init_EIR = 100,
                       n_particles=200,
                       proposal_matrix = matrix(1),
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
},
parallel = hipercow_parallel('parallel',cores_per_process = 32),
resources=resources)
hipercow::task_log_show(amoah_u5_short)
hipercow::task_log_watch(amoah_u5_short)
task_status(amoah_u5_short)
task_info(amoah_u5_short)
amoah_u5_short_results <- task_result(amoah_u5_short)
combine_elements <- function(e1, e2) {
  if (is.vector(e1) && is.vector(e2)) {
    return(c(e1, e2))
  } else if (is.data.frame(e1) && is.data.frame(e2)) {
    return(rbind(e1, e2))
  } else {
    stop("Elements must be of the same type (both vectors or both dataframes).")
  }
}

create_diag_figs(amoah_u5_short_results,
                 country = 'Amoah',
                 district = 1,
                 folderpath = './trial_figs/diag',
                 name = 'amoah_u5_short_results')
amoah_u5_short_results_prepped <- prep_results(results=amoah_u5_short_2_results,
             sim_data=amoah_prev_u5,
             burnin=0.1,
             country = 'Malawi',
             district = NA,
             anc=FALSE)
saveRDS(amoah_u5_short_results_prepped,'amoah_u5_short_results_prepped.rds')

amoah_u5_short_results$history[,1,1]
amoah_prev_u5$site <- 'Malawi'
amoah_prev_u5_plot <- create_dashboard_plots_trial(results=amoah_u5_short_results_prepped,
                                                    observed = amoah_prev_u5,
                                                    rainfall = trial_rainfall[trial_rainfall$country=='Malawi',],
                                                    single_site = TRUE,
                                                    var= 'prev_05',
                                                    title = 'Under 5 Prevalence\nAmoah, et al.')
amoah_eir_normalized <- amoah_eir %>%
  mutate(value = 100*((value - min(value)) / (max(value) - min(value))))

amoah_eir_plot <- create_dashboard_plots_trial(results=amoah_u5_short_results_prepped,
                                               observed = amoah_eir_normalized,
                                               rainfall = trial_rainfall[trial_rainfall$country=='Malawi',],
                                               max_value = 150,
                                               multiplier = 1,
                                               rainfall_multiplier = 150,
                                               single_site = TRUE,
                                               var= 'EIR',
                                               title = 'EIR\nAmoah, et al.')

amoah_plots_inital <- amoah_prev_u5_plot + amoah_eir_plot + plot_layout(ncol=2)
ggsave('./trial_figs/output/amoah_initial_fit_2.tiff',plot=amoah_plots_inital,width=7,height=3.5,units='in')
resources_8 <- hipercow_resources(cores=8)
amoah_prev_year1 <- amoah_prev_u5%>%
  slice_head(n=12)%>%
  summarise(annual_prev = sum(positive)/sum(tested))

amoah_u5_short_2 <- task_create_expr({
  mamasante::run_pmcmc(data_raw=amoah_prev_u5,
                       init_EIR = 100,
                       n_particles=200,
                       proposal_matrix = matrix(1),
                       target_prev = 0.321,
                       target_prev_group = 'u5',
                       max_param=125,
                       prop_treated = 0.4,
                       n_steps = 2000,
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
},
parallel = hipercow_parallel('parallel',cores_per_process = 32),
resources=resources)
task_status(amoah_u5_short_2)
task_info(amoah_u5_short_2)
amoah_u5_short_2_results <-task_result(amoah_u5_short_2)
create_diag_figs(amoah_u5_short_2_results,
                 country = 'Amoah',
                 district = 1,
                 folderpath = './trial_figs/diag',
                 name = 'amoah_u5_short_results_2')
amoah_u5_short_results_prepped_2 <- prep_results(results=amoah_u5_short_2_results,
                                               sim_data=amoah_prev_u5,
                                               burnin=0.1,
                                               country = 'Malawi',
                                               district = NA,
                                               anc=FALSE)
saveRDS(amoah_u5_short_results_prepped_2,'amoah_u5_short_results_prepped_2.rds')
amoah_prev_u5_plot <- create_dashboard_plots_trial(results=amoah_u5_short_results_prepped_2,
                                                   observed = amoah_prev_u5,
                                                   rainfall = trial_rainfall[trial_rainfall$country=='Malawi',],
                                                   single_site = TRUE,
                                                   var= 'prev_05',
                                                   title = 'Under 5 Prevalence\nAmoah, et al.')

amoah_eir_plot <- create_dashboard_plots_trial(results=amoah_u5_short_results_prepped_2,
                                               observed = amoah_eir_normalized,
                                               rainfall = trial_rainfall[trial_rainfall$country=='Malawi',],
                                               max_value = 150,
                                               multiplier = 1,
                                               rainfall_multiplier = 150,
                                               single_site = TRUE,
                                               var= 'EIR',
                                               title = 'EIR\nAmoah, et al.')

# task_cancel(amoah_u5_short_2)

task_info(amoah_u5_short_2)
task_log_show(amoah_u5_short_2)

# task_cancel(amoah_u5_short)
ggplot(amoah_prev_u5)+
  geom_line(aes(x=date,y=mean))+
  scale_y_continuous(limits = c(0,1))
ggplot(amoah_ento_sum)+
  geom_col(data=malawi_rainfall,aes(x=date,y=rainfall_rel*60),fill='lightgrey')+
  geom_line(aes(x=date,y=moz_count),linewidth=1)+
  scale_x_date(date_breaks='month',limits = range(amoah_ento_sum$date))
ggplot(amoah_ento_sum)+
  geom_col(data=malawi_rainfall,aes(x=date,y=rainfall_rel*60),fill='lightgrey')+
  geom_line(aes(x=date,y=arabss_female),linewidth=1)+
  scale_x_date(date_breaks='month',limits = range(amoah_ento_sum$date))
str(amoah_eir)

rogerson_prev <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/rogerson_2000.txt')%>%
  mutate(month_num = round(month),
         prev_anc = prev/100,
         tested = 397,
         positive = round(397*prev/100))%>%
  mutate(year = ifelse(month_num >= 8,1997,1998))%>%
  mutate(month = as.yearmon(paste0(year,'-',month_num),'%Y-%m'))%>%
  mutate(date = as.Date(month,frac=0.5))%>%
  arrange(date)
rogerson_prev_year1 <- rogerson_prev%>%
  summarise(annual_prev = sum(positive)/sum(tested))
targetprevs_trials <- get_u5prev_fromanc(avg_prev=rogerson_prev_year1$annual_prev,comparison='ancall')

rogerson_prev_short <- task_create_expr({
  mamasante::run_pmcmc(data_raw=rogerson_prev,
                       init_EIR = 100,
                       n_particles=200,
                       proposal_matrix = matrix(1),
                       target_prev = 0.5805684,
                       target_prev_group = 'u5',
                       max_param=125,
                       prop_treated = 0.4,
                       n_steps = 2000,
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
},
parallel = hipercow_parallel('parallel',cores_per_process = 32),
resources=resources)
task_status(rogerson_prev_short)
task_log_show(rogerson_prev_short)

rogerson_prev_short_results <- task_result(rogerson_prev_short)
rogerson_prev_short_results$run_time/3600

create_diag_figs(rogerson_prev_short_results,
                 country = 'Rogerson',
                 district = 1,
                 folderpath = './trial_figs/diag',
                 name = 'rogerson_short_results')
rogerson_prev_short_results_2$history[,1,1]
rogerson_prev_short_results_prepped <- prep_results(results=rogerson_prev_short_results,
                                               sim_data=rogerson_prev,
                                               burnin=0.1,
                                               country = 'Malawi',
                                               district = NA,
                                               anc=TRUE)
saveRDS(rogerson_prev_short_results_prepped,'rogerson_prev_short_results_prepped.rds')
rogerson_prev$site <- 'Malawi'
rogerson_prev_plot <- create_dashboard_plots_trial(results=rogerson_prev_short_results_prepped,
                                                   observed = rogerson_prev,
                                                   rainfall = trial_rainfall[trial_rainfall$country=='Malawi',],
                                                   single_site = TRUE,
                                                   var= 'prev_05',
                                                   title = 'Under 5 Prevalence\nRogerson, et al.')
rogerson_eir_plot <- create_dashboard_plots_trial(results=rogerson_prev_short_results_prepped,
                                               observed = NULL,
                                               rainfall = trial_rainfall[trial_rainfall$country=='Malawi',],
                                               max_value = 150,
                                               multiplier = 1,
                                               rainfall_multiplier = 150,
                                               single_site = TRUE,
                                               var= 'EIR',
                                               title = 'EIR\nRogerson, et al.')

windows_authenticate()
rogerson_prev_short_2 <- task_create_expr({
  mamasante::run_pmcmc(data_raw=rogerson_prev,
                       init_EIR = 100,
                       n_particles=200,
                       proposal_matrix = matrix(1),
                       target_prev = 0.5805684,
                       target_prev_group = 'u5',
                       max_param=125,
                       prop_treated = 0.4,
                       n_steps = 2000,
                       n_threads = 10,
                       n_chains = 1,
                       n_workers = 1,
                       state_check = 0,## Run equilibrium checks
                       seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                       seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                       seed = 1L,
                       start_pf_time = 30*12,
                       particle_tune = FALSE,
                       comparison = 'ancall',
                       initial = 'informed')
},
parallel = hipercow_parallel('parallel',cores_per_process = 32),
resources=resources)
task_info(rogerson_prev_short_2)
task_log_show(rogerson_prev_short_2)
rogerson_prev_short_2_results <- task_result(rogerson_prev_short_2)
rogerson_prev_short_2_results$history[,1,1]
create_diag_figs(rogerson_prev_short_2_results,
                 country = 'Rogerson',
                 district = 1,
                 folderpath = './trial_figs/diag',
                 name = 'rogerson_short_results_2')

rogerson_prev_short_results_prepped_2 <- prep_results(results=rogerson_prev_short_2_results,
                                                    sim_data=rogerson_prev,
                                                    burnin=0.1,
                                                    country = 'Malawi',
                                                    district = NA,
                                                    anc=TRUE)
saveRDS(rogerson_prev_short_results_prepped_2,'rogerson_prev_short_results_prepped_2.rds')
View(rogerson_prev_short_results_prepped_2$summary)
rogerson_prev_plot_2 <- create_dashboard_plots_trial(results=rogerson_prev_short_results_prepped_2,
                                                   observed = rogerson_prev,
                                                   rainfall = trial_rainfall[trial_rainfall$country=='Malawi',],
                                                   single_site = TRUE,
                                                   var= 'prev_anc_all',
                                                   title = 'ANC Prevalence\nRogerson, et al.')
rogerson_eir_plot_2 <- create_dashboard_plots_trial(results=rogerson_prev_short_results_prepped_2,
                                                  observed = NULL,
                                                  rainfall = trial_rainfall[trial_rainfall$country=='Malawi',],
                                                  max_value = 500,
                                                  multiplier = 1,
                                                  rainfall_multiplier = 150,
                                                  single_site = TRUE,
                                                  var= 'EIR',
                                                  title = 'EIR\nRogerson, et al.')

feasey_inc <- readxl::read_xlsx('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/Feasey_S1_Table.xlsx')
str(feasey_inc)
feasey_inc$date <- as.Date(as.yearmon(paste0(feasey_inc$Year,'-',feasey_inc$Month),'%Y-%B'),frac=0.5)
feasey_inc$year <- year(feasey_inc$date)
feasey_inc$month <- month(feasey_inc$date)
feasey_inc$value <- feasey_inc$malaria_cases
feasey_inc$source <- 'Feasey - Routine Data'
feasey_inc$location <- 'Blantyre'
feasey_inc_ave <- feasey_inc%>%
  group_by(month)%>%
  summarise(value=mean(value))%>%
  mutate(date=as.Date(as.yearmon(paste0('2010-',month))))
feasey_inc_ave$source <- 'Feasey - Routine Data'


incidence_anc_malawi <- incidence %>%
  filter(country=='Malawi')%>%
  mutate(month=month(date))%>%
  group_by(date)%>%
  summarise(value=median(value),
            lower=quantile(value,c(0.025,0.975)[[1]]),
            upper=quantile(value,c(0.025,0.975)[[2]]))%>%
  mutate(source='ANC-estimated Incidence')

incidence_malawi_all <- bind_rows(blantyre_total,
                                  chikwawa_total,
                                  cereb_total,
                                  amoah_eir,
                                  feasey_inc)%>%
  mutate(source=fct_rev(factor(source,levels=c('Tizifia - Cohort Clinical Incidence',
                                               'Amoah - EIR',
                                               'Roca-Feltrer - Cerebral malaria',
                                               'Roca-Feltrer - All severe malaria',
                                               'Feasey - Routine Data'))))%>%
  group_by(source)%>%
  mutate(norm_value = (value - min(value)) / (max(value) - min(value)))%>%
  ungroup()

source_palette <- trial_color_palette
names(source_palette) <- unique(incidence_malawi_all$source)
date_breaks_ts <- seq.Date(from=as.Date('2010-01-01'),as.Date('2010-12-01'),by='month')
date_labels_ts <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

plot_timeseries <- ggplot(incidence_malawi_all)+
  geom_line(aes(x=date,y=norm_value,color=source),size=1)+
  scale_color_manual(name='Source',values=source_palette,breaks=names(source_palette))+
  scale_x_date(limits = as.Date(c('2010-01-01','2010-12-15')),labels = date_labels_ts,breaks = date_breaks_ts)+
  labs(y='Normalized value')+
  theme(axis.title.x = element_blank())
colors_sources_blantchik <- colors_sources[2:3]
names(colors_sources_blantchik) <- c('Chikwawa','Blantyre')
plot_facet_timeseries <- ggplot(incidence_malawi_all)+
  geom_line(data=malawi_rainfall ,aes(x=date,y=rainfall_rel/500),alpha = 1,color = 'darkgrey')+
  geom_line(aes(x=date,y=norm_value,color=location))+
  facet_grid(source~.)+
  coord_cartesian(xlim=c(as.Date('2001-01-01'),NA))+
  scale_color_manual(name='Source',values=colors_sources_blantchik,breaks=names(colors_sources_blantchik),guide='none')+
  labs(y='Normalized value')+
  theme(axis.title.x = element_blank())
ggsave('./trial_figs/output/malawi_timeseries_3.tiff',plot=plot_facet_timeseries,width=7,height=6,units='in')

amoah_msi <- return_msi(dataframe=amoah_eir_ave,single_year = TRUE) %>%
  mutate(x_coord = 0 + rk * cos(theta_k),
         y_coord = 0 + rk * sin(theta_k))
amoah_msi$source <- 'Amoah - EIR'
amoah_vectors <- return_vectors(dataframe = amoah_eir_ave,single_year = TRUE)%>%
  mutate(xstart = lag(x_coord),ystart = lag(y_coord))%>%
  na.omit()

feasey_msi <- return_msi(dataframe=feasey_inc_ave,single_year = TRUE) %>%
  mutate(x_coord = 0 + rk * cos(theta_k),
         y_coord = 0 + rk * sin(theta_k))
feasey_msi$source <- 'Feasey - Routine Data'
amoah_vectors <- return_vectors(dataframe = feasey_inc_ave,single_year = TRUE)%>%
  mutate(xstart = lag(x_coord),ystart = lag(y_coord))%>%
  na.omit()
anc_malawi_msi <- sample_msi_summary%>%
  filter(country=='Malawi')%>%
  mutate(source = 'ANC-estimated Incidence')%>%
  rename(theta_k_deg=angle_median,
         msi=msi_median)
malawi_msi_compare <- bind_rows(anc_malawi_msi,
                                blantyre_msi,
                                chikwawa_msi,
                                cereb_msi,
                                amoah_msi,
                                feasey_msi)%>%
  mutate(source=fct_rev(factor(source,levels=c('ANC-estimated Incidence',
                                       'Amoah - EIR',
                                       'Tizifia - Cohort Clinical Incidence',
                                       'Roca-Feltrer - Cerebral malaria',
                                       'Feasey - Routine Data',
                                       'Roca-Feltrer - All severe malaria'))))%>%
  mutate(angle_date_median=as.Date(ifelse(theta_k_deg>=0.25,as.Date(date_decimal(2010+theta_k_deg/360)),as.Date('2009-12-15'))),
         angle_date_lower=as.Date(ifelse(theta_k_deg>=0.25,as.Date(date_decimal(2010+angle_lower/360)),as.Date('2009-12-15'))),
         angle_date_upper=as.Date(ifelse(theta_k_deg>=0.25,as.Date(date_decimal(2010+angle_upper/360)),as.Date('2009-12-15'))))

malawi_msi_plot <- ggplot(malawi_msi_compare)+
  geom_point(aes(x=msi,y=source,color=source))+
  geom_errorbarh(aes(xmin=msi_lower,xmax=msi_upper,y=source,color=source),height=0)+
  coord_cartesian(xlim = c(0,1))+
  scale_color_manual(name='Source',values=source_palette,breaks=names(source_palette))+
  labs(x='MSI')+
  theme(axis.title.y = element_blank())

malawi_rain_msi_plot <- ggplot(malawi_rain_msi)+
  geom_point(aes(y=msi,x=year))+
  coord_cartesian(ylim = c(0,1))+
  labs(y='MSI')+
  scale_x_continuous(labels=c(2009:2019),breaks=c(2009:2019))+
  theme(axis.title.x = element_blank())
malawi_rain_msi_plot_2 <- ggplot(malawi_rain_msi)+
  geom_point(aes(y=msi,x=year))+
  coord_cartesian(ylim = c(0,1))+
  labs(y='MSI')+
  scale_x_continuous(labels=c(2009:2019),breaks=c(2009:2019))+
  theme(axis.title.x = element_blank())

malawi_rain_peak_plot <- ggplot(malawi_rain_msi)+
  # annotate("rect", xmin = as.Date('2009-12-01'), xmax = as.Date('2009-12-31'), ymin = as.Date('2009-12-01'), ymax = as.Date('2010-12-15'),alpha = .4,fill = '#C3A995FF')+
  geom_point(aes(y=angle_date,x=year))+
  scale_x_continuous(labels=c(2009:2019),breaks=c(2009:2019))+
  scale_y_date(limits = as.Date(c('2009-12-01','2010-12-15')),labels = date_labels,breaks = date_breaks)+
  labs(y='Estimated Peak')+
  theme(axis.title.x = element_blank())
malawi_rain_peak_plot_2 <- ggplot(malawi_rain_msi)+
  # annotate("rect", xmin = as.Date('2009-12-01'), xmax = as.Date('2009-12-31'), ymin = as.Date('2009-12-01'), ymax = as.Date('2010-12-15'),alpha = .4,fill = '#C3A995FF')+
  geom_point(aes(y=angle_date,x=year))+
  scale_x_continuous(labels=c(2009:2019),breaks=c(2009:2019))+
  scale_y_date(limits = as.Date(c('2009-12-01','2010-12-15')),labels = date_labels,breaks = date_breaks)+
  labs(y='Estimated Peak')+
  theme(axis.title.x = element_blank())

malawi_rain_plots_2 <- malawi_rain_msi_plot_2 + malawi_rain_peak_plot_2 + plot_layout(ncol=2)
str(malawi_msi_compare)
ggsave('./trial_figs/output/malawi_rain_plots_2.tiff',plot=malawi_rain_plots_2,width=7,height=3.5,units='in')

malawi_peak_plot <- ggplot(malawi_msi_compare)+
  # annotate("rect", xmin = as.Date('2009-12-01'), xmax = as.Date('2009-12-31'), ymin = as.Date('2009-12-01'), ymax = as.Date('2010-12-15'),alpha = .4,fill = '#C3A995FF')+
  geom_point(aes(x=angle_date_median,y=source,color=source))+
  geom_errorbarh(aes(y=source,xmin=angle_date_lower,xmax=angle_date_upper,color=source),height=0)+
  scale_x_date(limits = as.Date(c('2009-12-01','2010-12-15')),labels = date_labels,breaks = date_breaks)+
  scale_color_manual(name='Source',values=source_palette,breaks=names(source_palette))+
  labs(x='Estimated Peak')+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())

malawi_plots <- plot_facet_timeseries / (malawi_msi_plot + malawi_peak_plot) & theme(legend.position = 'none')
ggsave('./trial_figs/output/malawi_plots.tiff',plot=malawi_plots,width=6,height=5,units='in')


##Chikwawa prospective cohort
chikwawa_fu_dates <- readxl::read_xlsx('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/chikwawa_follow-up_time.xlsx')
chikwawa_event_dates <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/chikwawa_event_times.csv')
chikwawa_all <- left_join(chikwawa_fu_dates,chikwawa_event_dates,by=join_by(participant_id==openhdsindividualId))
chikwawa_all$malaria_final <- sapply(1:nrow(chikwawa_all), function(x){ifelse(any(chikwawa_all[x,]$Malaria,chikwawa_all[x,]$Malaria1,chikwawa_all[x,]$Malaria2),TRUE,FALSE)})
chikwawa_all$visit_date_1 <- as.Date(chikwawa_all$visit_date,format="%d/%m/%Y")
chikwawa_all <- chikwawa_all %>%
  mutate(malaria_final = ifelse(is.na(malaria_final),FALSE,malaria_final))
num_cases_perchild <- data.frame(table(chikwawa_all$participant_id))
hist(num_cases_perchild$Freq)
chikwawa_cases_per_month <- chikwawa_all %>%
  mutate(month=month(visit_date_1),
         year=year(visit_date_1),
         yearmon = as.yearmon(visit_date_1))%>%
  group_by(yearmon)%>%
  summarise(case_count=sum(malaria_final))
str(chikwawa_all)

##Zomba hospitalizations
zomba_1 <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/zomba_hosp_marks_1.txt')
zomba_2 <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/zomba_hosp_marks_2.txt')
zomba_3 <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/zomba_hosp_marks_3.txt')
zomba_4 <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/zomba_hosp_marks_4.txt')
zomba_5 <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/zomba_hosp_marks_5.txt')
zomba_6 <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/zomba_hosp_marks_6.txt')
zomba_7 <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/zomba_hosp_marks_7.txt')
zomba_all <- bind_rows(zomba_1,
                       zomba_2,
                       zomba_3,
                       zomba_4,
                       zomba_5,
                       zomba_6,
                       zomba_7)
zomba_all$date <- as.Date(date_decimal(zomba_all$date_decimal))

zomba_total <- zomba_all %>%
  group_by(date)%>%
  summarise(value=mean(value))%>%
  ungroup()%>%
  mutate(year = year(date),
         month = month(date))%>%
  group_by(year,month)%>%
  summarise(value=mean(value))%>%
  ungroup()%>%
  group_by(month)%>%
  summarise(value=mean(value))
zomba_total_year <- zomba_all %>%
  group_by(date)%>%
  summarise(value=mean(value))%>%
  ungroup()%>%
  mutate(year = year(date),
         month = month(date))%>%
  group_by(year,month)%>%
  summarise(value=mean(value))%>%
  ungroup()
zomba_total$date <- seq.Date(as.Date('2000-01-01'),as.Date('2000-12-01'),by='month')
plot(zomba_total$value)
lines(zomba_total$value)
zomba_msi <- return_msi(dataframe=zomba_total,single_year = TRUE) %>%
  mutate(x_coord = 0 + rk * cos(theta_k),
         y_coord = 0 + rk * sin(theta_k))
zomba_msi_years <- return_msi(dataframe=zomba_total_year,single_year = FALSE) %>%
  mutate(x_coord = 0 + rk * cos(theta_k),
         y_coord = 0 + rk * sin(theta_k))
zomba_vectors <- return_vectors(dataframe = zomba_total,single_year = TRUE)%>%
  mutate(xstart = lag(x_coord),ystart = lag(y_coord))%>%
  na.omit()

##Gambia inlay
basse_proportion <- basse_total%>%
  mutate(prop=value/sum(value))

incidence_filtered <- incidence%>%
  filter(country=='Gambia'&time!=13)
sample_gambia_prop <- incidence_filtered %>%
  mutate(month=month(as.Date(date)))%>%
  group_by(variable) %>%
  mutate(prop=value/sum(value))

gambia_prop_summary <- sample_gambia_prop%>%
  group_by(month)%>%
  summarise(median=median(prop),
            lower=quantile(prop,0.025),
            upper=quantile(prop,0.975))
gambia_prop_both <- merge(basse_proportion,gambia_prop_summary,by='month')
gambia_prop_both$date <- seq.Date(as.Date('2000-01-01'),as.Date('2000-12-01'),by='month')
gambia_prop_both$prop_2 <- ifelse(gambia_prop_both$prop==0,NA,gambia_prop_both$prop)

gambia_plot <- ggplot(gambia_prop_both)+
  geom_point(aes(x=date,y=median),color='darkorange',size=1.5)+
  geom_errorbar(aes(x=date,ymin=lower,ymax=upper),color='darkorange',width=0)+
  geom_point(aes(x=date,y=prop_2),color='darkblue',size=1.5)+
  scale_x_date(limits = as.Date(c('2000-01-01','2000-12-01')),labels = date_format("%b"),breaks = "1 month")+
  scale_y_continuous(limits=c(0,0.75),expand=c(0,0))+
  labs(x='Month',
         y='Monthly proportion of annual cases')
ggsave('./trial_figs/output/gambia_comparison_timeseries.tiff',plot=gambia_plot,width=3.5,height=3.5,units='in')

