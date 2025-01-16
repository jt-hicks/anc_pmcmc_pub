create_dashboard_plots_trial_2 <- function(results,
                                         observed,
                                         rainfall=NULL,
                                         prev_all=NULL,
                                         var,
                                         title=NULL,
                                         max_value=1,
                                         multiplier=1,
                                         rainfall_multiplier=0.75,
                                         facet_scales = 'fixed',
                                         single_site = FALSE,
                                         show_fits = TRUE){
  ## fn to return prevalence from log_odds
  country <- c("Burkina Faso","Gambia","Ghana","Mali","Kenya","Malawi")
  country_labels <- c("Burkina Faso","The Gambia","Ghana","Mali","Kenya","Malawi")

  names(country_labels) <- country
  smc_times <- data.frame(country=factor(c(rep('Ghana',8),rep('Burkina Faso',6),rep('Mali',6),rep('Gambia',3)),
                                         levels=c('Gambia','Mali','Burkina Faso','Ghana','Malawi','Kenya')),
                          smc_times=as.Date(c('2010-07-01','2010-08-01','2010-09-01','2010-10-01',
                                              '2011-07-01','2011-08-01','2011-09-01','2011-10-01',
                                              '2010-08-01','2010-09-01','2010-10-01',
                                              '2011-08-01','2011-09-01','2011-10-01',
                                              '2010-08-01','2010-09-01','2010-10-01',
                                              '2011-08-01','2011-09-01','2011-10-01',
                                              '2010-09-01','2010-10-01','2010-11-01')))

  shift_trial_dates <- function(df,single_site) {
    if(single_site){
      df$date_aligned <- df$date
      return(df)
    }
    df %>%
      mutate(
        date_aligned = case_when(
          country == "Kenya" ~ as.Date(date) %m-% years(2),
          country == "Malawi" ~ as.Date(date) %m-% years(1),
          TRUE ~ date                   # Default case: keep date
        )
      )
  }
  colors <-c(viridis::viridis(6,begin=0,end=0.85))
  names(colors) <- c('Gambia','Mali','Burkina Faso','Ghana','Malawi','Kenya')
  twolevels <- c(colors, lighten(colors,0.5))
  names(twolevels) <- c(paste0('prev_pg-',names(colors)),paste0('prev_mg-',names(colors)))

  if(single_site & length(results)>1){
    results <- list(results)
  }
  if(var=='prev_anc'){
    measure_labels <- c('Primigravidae','Secundi- or Multigravidae')
    names(measure_labels) <- c('prev_pg','prev_mg')

    results_summary <- bind_rows(lapply(1:length(results), function(x){
      results[[x]]$summary[results[[x]]$summary$measure %in% c('prev_pg','prev_mg'),]
    }))
    results_sample <- bind_rows(lapply(1:length(results), function(x){
      results[[x]]$sample[results[[x]]$sample$measure%in% c('prev_pg','prev_mg'),]
    }))
    results_summary$measure <- factor(results_summary$measure, levels=c('prev_pg','prev_mg'))
    results_sample$measure <- factor(results_sample$measure, levels=c('prev_pg','prev_mg'))
  }else{
    measure_labels <- NULL
    results_summary <- bind_rows(lapply(1:length(results), function(x){
      results[[x]]$summary[results[[x]]$summary$measure==var,]
    }))
    results_sample <- bind_rows(lapply(1:length(results), function(x){
      results[[x]]$sample[results[[x]]$sample$measure==var,]
    }))
  }

  if(!is.null(observed) & grepl('prev', var, fixed = TRUE)){
    observed_cis <- addCIs(observed,observed$positive,observed$tested)
    observed_cis$country <- observed_cis$site
    observed_cis$date <- as.Date(observed_cis$month,frac=0.5)
    observed_cis <- shift_trial_dates(df=observed_cis,single_site)
  } else if(!is.null(observed)){
    observed_cis <- observed
    if('site' %in% names(observed)){
      observed_cis$country <- observed_cis$site
    }
    if(!('mean' %in% names(observed))){
      observed_cis$mean <- observed_cis$value
    }
    if(!('lower' %in% names(observed))){
      observed_cis$lower <- NA
      observed_cis$upper <- NA
    }
    observed_cis <- shift_trial_dates(observed_cis,single_site)
  } else {
    observed_cis <- data.frame(country=character(),
                                     date_aligned = Date(),
                                     mean = numeric(),
                                     lower = numeric(),
                                     upper = numeric())}
  results_summary$date_aligned <- results_summary$date
  results_summary <- shift_trial_dates(results_summary,single_site)
  results_sample$date_aligned <- results_sample$date
  results_sample <- shift_trial_dates(results_sample,single_site)

  results_summary$country <- factor(results_summary$country,levels=c('Gambia','Mali','Burkina Faso','Ghana','Malawi','Kenya'))
  results_summary$twolevels <- paste0(results_summary$measure,'-',results_summary$country)
  results_sample$country <- factor(results_sample$country,levels=c('Gambia','Mali','Burkina Faso','Ghana','Malawi','Kenya'))
  results_sample$twolevels <- paste0(results_sample$measure,'-',results_sample$country)
  observed_cis$country <- factor(observed_cis$country,levels=c('Gambia','Mali','Burkina Faso','Ghana','Malawi','Kenya'))
  observed_cis$twolevels <- paste0(observed_cis$measure,'-',observed_cis$country)
  y_axis_label <- NULL
  if(var=='prev_anc'){
    observed_cis$measure <- factor(observed_cis$measure, levels=c('prev_pg','prev_mg'))
    y_axis_label <- 'ANC prevalence'
  }else if(var=='inc05'){
    y_axis_label <- 'Daily clinical cases per 1000 children under 5 years old'
  }

  plot_base <- ggplot()+
    scale_y_continuous(expand=c(0,0))

  if(!is.null(rainfall)){
    rainfall$date_aligned <- as.Date(rainfall$month,frac=0.5)
    rainfall <- shift_trial_dates(rainfall,single_site)
    rainfall$rainfall_rel <- rainfall$rainfall*rainfall_multiplier/max(rainfall$rainfall)
    rainfall$country <- factor(rainfall$country,levels=c('Gambia','Mali','Burkina Faso','Ghana','Malawi','Kenya'))

    plot_base <- plot_base +
      geom_col(data=rainfall,aes(x=date_aligned,y=rainfall_rel),alpha = 1,fill = 'darkgrey',just=0)+
      # scale_y_continuous(expand=c(0,0),sec.axis = sec_axis(~ . /(rainfall_multiplier*10/max(rainfall$rainfall)), name = "Monthly rainfall (cm)"))+
      theme(strip.placement = "outside")

  }
  if(show_fits){
    plot_base <- plot_base+
      # geom_line(data=results_sample,aes(x=as.Date(date_aligned),y=value*multiplier,color=country,group=variable),alpha=0.1,linewidth=0.2)+
      geom_line(data=results_summary,aes(x=as.Date(date_aligned),y=median*multiplier,color=twolevels,group=twolevels),linewidth=0.8)
  }

  plot <- plot_base+
    geom_point(data=observed_cis,aes(x=as.Date(date_aligned),y=mean*multiplier,color=twolevels,group=twolevels),pch = 19,position=position_dodge(width=10),size=0.5)+
    geom_errorbar(data=observed_cis,aes(x=as.Date(date_aligned),ymin=lower*multiplier,ymax=upper*multiplier,color=twolevels,group=twolevels),width = 0,position=position_dodge(width=10),linewidth=0.5)+
    scale_color_manual(values=twolevels)+
    scale_x_date(date_labels = "%b", date_breaks = '3 months')+
    scale_y_continuous(expand=c(0,0),breaks=c(0,0.25,0.5,0.75))+
    coord_cartesian(ylim = c(0,max_value),
                    xlim = range(results_summary$date_aligned))+
    labs(title = title,
         y = y_axis_label)+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          # axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none',
          panel.spacing.y = unit(0, "mm")
    )
  if(var!='prev_anc'){
    plot <- plot +
      theme(strip.text.x = element_blank())
  }
  if(!single_site){
    plot <- plot +
      geom_vline(data=smc_times,aes(xintercept=smc_times),linetype='dashed',color='black',linewidth=0.5)+
      facet_grid(country~.,scales=facet_scales,labeller = labeller(country = country_labels, measure=measure_labels))
  }
  return(plot)
}
