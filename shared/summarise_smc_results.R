create_dashboard_plots_smc <- function(results,
                                       observed,
                                       analysis_type,
                                       rainfall=NULL,
                                       prev_all=NULL,
                                       smc_times=NULL,
                                       var,
                                       title=NULL,
                                       max_value=1,
                                       multiplier=1,
                                       rainfall_multiplier=0.75,
                                       facet_scales = 'fixed',
                                       show_fits = TRUE,
                                       zoom_in = NULL){


  # colors <-c(viridis::viridis(length(results),begin=0,end=0.95))
  colors <-c("#377EB8","#FF7F00","#984EA3")
  # colors <-c("#2F6B8EFF",'#EB5E28',"#E41A1C")
  # colors <-c("green",'#EB5E28',"#BEAED4")

    # names(colors) <- names(results)


  if(var=='prev_anc'){
    measure_labels <- c('Primigravidae','Secundi- or Multigravidae')
    names(measure_labels) <- c('prev_pg','prev_mg')

    results_summary <- dplyr::bind_rows(lapply(1:length(results), function(x){
      results[[x]]$summary[results[[x]]$summary$measure %in% c('prev_pg','prev_mg'),]
    }))
    results_sample <- dplyr::bind_rows(lapply(1:length(results), function(x){
      results[[x]]$sample[results[[x]]$sample$measure%in% c('prev_pg','prev_mg'),]
    }))
    results_summary$measure <- factor(results_summary$measure, levels=c('prev_pg','prev_mg'))
    results_sample$measure <- factor(results_sample$measure, levels=c('prev_pg','prev_mg'))
  } else {
    measure_labels <- NULL
    results_summary <- dplyr::bind_rows(lapply(1:length(results), function(x){
      results[[x]]$summary[results[[x]]$summary$measure==var,]
    }))
    results_sample <- dplyr::bind_rows(lapply(1:length(results), function(x){
      results[[x]]$sample[results[[x]]$sample$measure==var,]
    }))
  }

  plot_corr <- NA

  if(!is.null(observed) & grepl('prev', var, fixed = TRUE)){
    observed_cis <- addCIs(observed,observed$positive,observed$tested)
    if(analysis_type == 'SMC'){
      observed_cis$site <- observed_cis$Council
      if(!is.null(zoom_in)){
        observed_cis$site <- factor(observed_cis$site,levels=zoom_in)
      }

    }else {
      observed_cis$site <- observed_cis$Region
    }
    observed_cis$date <- zoo::as.Date(observed_cis$yearmon,frac=0.5)
  } else if(!is.null(observed)){
    observed_cis <- observed
    if(analysis_type == 'SMC'){
      ##Convert to daily rate
      observed_cis$mean <- observed_cis$mean/30
      observed_cis$lower <- observed_cis$lower/30
      observed_cis$upper <- observed_cis$upper/30

      observed_cis$site <- observed_cis$council
      if(!is.null(zoom_in)){
        observed_cis$site <- factor(observed_cis$site,levels=zoom_in)
      }
      observed_cis$date <- zoo::as.Date(zoo::as.yearmon(zoo::as.Date(paste0(observed_cis$month,'-01-',observed_cis$year),format='%B-%d-%Y')),frac=0.5)

      inc_comp <- left_join(results_summary,observed_cis,by=join_by(site==site,date==date),suffix=c('.est','.obs'))
    #   plot_corr <- ggplot(data=inc_comp)+
    #     geom_point(aes(x=mean.obs*multiplier,y=median*multiplier,color=site))+
    #     geom_errorbar(aes(x=mean.obs*multiplier,ymin=lower.est*multiplier,ymax=upper.est*multiplier,color=site))+
    #     geom_errorbarh(aes(xmin=lower.obs*multiplier,xmax=upper.obs*multiplier,y=median*multiplier,color=site))+
    #     geom_abline(linetype='dashed')+
    #     facet_wrap(site~.,scales='fixed')+
    #     scale_y_continuous(expand=c(0,0),limits=c(0,25))+
    #     scale_x_continuous(expand=c(0,0),limits=c(0,25))+
    #     scale_color_manual(values=colors)+
    #     labs(x='Daily Clinical Case Reports per 1000 children under 5',
    #          y='Daily Estimated Clinical Incidence per 1000 children under 5')+
    #     theme(legend.title = element_blank(),
    #           axis.text.x=element_text(angle=45, hjust=1, vjust=1),
    #           axis.ticks.x = element_line(linewidth = 0.5),
    #           axis.ticks.length = unit(3, "pt"),
    #           legend.position = 'none',
    #           panel.spacing.y = unit(3, "mm"))


    } else {
      observed_cis$site <- observed_cis$Region
    }
    if(!('mean' %in% names(observed))){
      observed_cis$mean <- observed_cis$value
    }
    if(!('lower' %in% names(observed))){
      observed_cis$lower <- NA
      observed_cis$upper <- NA
    }
  } else {
    observed_cis <- data.frame(site=character(),
                               date = as.Date(integer()),
                               mean = numeric(),
                               lower = numeric(),
                               upper = numeric())}
  if(!is.null(zoom_in)){
    results_sample$site <- factor(results_sample$site,levels=zoom_in)
    results_summary$site <- factor(results_summary$site,levels=zoom_in)
  }

  # print(str(observed_cis))
  # print(str(results_sample))
  # print(str(results_summary))

  y_axis_label <- NULL
  if(var=='prev_anc_all'){
    y_axis_label <- 'ANC prevalence'
  }else if(var=='inc05'){
    y_axis_label <- 'Daily clinical cases per 1000 children under 5 years old'
  }

  plot_base <- ggplot()+
    scale_y_continuous(expand=c(0,0))

  if(!is.null(rainfall)){
    rainfall$date <- zoo::as.Date(zoo::as.yearmon(rainfall$yearmon),frac=0.5)
    rainfall$rainfall_rel <- rainfall$Rainfall*rainfall_multiplier/max(rainfall$Rainfall)
    if(analysis_type == 'SMC'){
      rainfall$site <- rainfall$Council
      if(!is.null(zoom_in)){
        rainfall$site <- factor(rainfall$site,levels=zoom_in)
      }

    }else {
      rainfall$site <- rainfall$Region
    }
    # print(str(rainfall))

    plot_base <- plot_base +
      geom_col(data=rainfall,aes(x=date,y=rainfall_rel),alpha = 1,fill = 'darkgrey',just=0)+
      scale_y_continuous(expand=c(0,0),sec.axis = sec_axis(~ . /(rainfall_multiplier*10/max(rainfall$Rainfall)), name = "Monthly rainfall (cm)"))+
      theme(strip.placement = "outside")

  }
  if(show_fits& grepl('inc', var, fixed = TRUE)){
    plot <- plot_base+
      geom_line(data=results_sample,aes(x=as.Date(date),y=value*multiplier,group=variable),alpha=0.1,linewidth=0.2,color=colors[1])+
      geom_line(data=results_summary,aes(x=as.Date(date),y=median*multiplier,group=site),linewidth=0.8,color=colors[1])+
      geom_point(data=observed_cis,aes(x=as.Date(date),y=mean*multiplier,group=site),color=colors[2],pch = 19,position=position_dodge(width=10),size=0.8,alpha=1)+
      geom_errorbar(data=observed_cis,aes(x=as.Date(date),ymin=lower*multiplier,ymax=upper*multiplier,group=site),color=colors[2],width = 0,position=position_dodge(width=10),linewidth=0.5,alpha=1)

  }
  if(!show_fits& grepl('inc', var, fixed = TRUE)){
    plot <- plot_base+
      geom_point(data=observed_cis,aes(x=as.Date(date),y=mean*multiplier,group=site),color=colors[2],pch = 19,position=position_dodge(width=10),size=0.8,alpha=1)+
      geom_errorbar(data=observed_cis,aes(x=as.Date(date),ymin=lower*multiplier,ymax=upper*multiplier,group=site),color=colors[2],width = 0,position=position_dodge(width=10),linewidth=0.5,alpha=1)

  }


  if(show_fits& grepl('prev', var, fixed = TRUE)){
    plot <- plot_base+
      geom_point(data=observed_cis,aes(x=as.Date(date),y=mean*multiplier,group=site),color=colors[3],pch = 19,position=position_dodge(width=10),size=0.8,alpha=1)+
      geom_errorbar(data=observed_cis,aes(x=as.Date(date),ymin=lower*multiplier,ymax=upper*multiplier,group=site),color=colors[3],width = 0,position=position_dodge(width=10),linewidth=0.5,alpha=1)+
      geom_line(data=results_sample,aes(x=as.Date(date),y=value*multiplier,group=variable),alpha=0.1,linewidth=0.2,color=colors[1])+
      geom_line(data=results_summary,aes(x=as.Date(date),y=median*multiplier,group=site),linewidth=0.8,color=colors[1])
  }
  if(!show_fits& grepl('prev', var, fixed = TRUE)){
    plot <- plot_base+
      geom_point(data=observed_cis,aes(x=as.Date(date),y=mean*multiplier,group=site),color=colors[3],pch = 19,position=position_dodge(width=10),size=0.8,alpha=1)+
      geom_errorbar(data=observed_cis,aes(x=as.Date(date),ymin=lower*multiplier,ymax=upper*multiplier,group=site),color=colors[3],width = 0,position=position_dodge(width=10),linewidth=0.5,alpha=1)
  }

  plot <- plot +
    # scale_color_manual(values=colors)+
    scale_x_date(date_labels = "%b %y", date_breaks = '1 year')+
    coord_cartesian(ylim = c(0,max_value),
                    xlim = range(results_summary$date))+
    labs(title = title,
         y = y_axis_label)+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          # axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          # axis.ticks.x = element_line(linewidth = 0.5),
          # axis.ticks.length = unit(3, "pt"),
          legend.position = 'none',
          panel.spacing.y = unit(3, "mm")
    )
  if(analysis_type == 'SMC'){
    plot <- plot +
      geom_vline(data=smc_times,aes(xintercept=smc_times),linetype='dashed',color='black',linewidth=0.5)+
      facet_wrap(site~.,ncol=1,scales=facet_scales)
  } else {
    province_grid <- read.csv('Province_grid.csv')

    province_grid$name=province_grid$code=province_grid$NAME_1
    province_grid[province_grid$NAME_1=='Dar Es Salaam',]$NAME_1 <- 'DAR'
    province_grid[province_grid$name=='Dar Es Salaam',]$name <- 'DAR'
    province_grid[province_grid$code=='Dar Es Salaam',]$code <- 'DAR'

    plot <- plot +
      geom_vline(data=smc_times,aes(xintercept=smc_times),linetype='dashed',color='black',linewidth=0.5)+
      geofacet::facet_geo(~ sites, grid = province_grid%>%
                  select(row,col,code,name))
  }

  return(list(plot,plot_corr))
}
