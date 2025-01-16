create_diag_figs <- function(result,name,country,district,burnin=0){
  length <- length(result$mcmc[,1])
  start_chain <- round(burnin*length)+1
  mcmc <- result$mcmc[start_chain:length,]
  print('acceptance rate')
  ar <- 1 - coda::rejectionRate(coda::as.mcmc(mcmc))
  print(ar)
  print('effective size')
  ess <- coda::effectiveSize(coda::as.mcmc(mcmc))
  print(ess)

  title <- paste0('Diagnostic plots for seasonal model - ',district,', ',country)
  pars_list <- names(mcmc)
  trace_plots <- lapply(pars_list, function(x){
    bayesplot::mcmc_trace(mcmc,pars = x) +
      ggtitle(paste0(x,' / AR: ',round(ar[x],3),' / ESS: ',round(ess[x],1)))+
      theme(title = element_text(size=6),
            axis.title.y = element_blank())
  })
  dense_plots <- lapply(pars_list, function(x){
    bayesplot::mcmc_dens(mcmc,pars = x) +
      ggtitle(paste0(x,' / AR: ',round(ar[x],2),' / ESS: ',round(ess[x],1)))+
      theme(title = element_text(size=6),
            axis.title.x = element_blank())
  })
  diag <- (trace_plots[[1]]+dense_plots[[1]])/
    (trace_plots[[2]]+dense_plots[[2]])/
    (trace_plots[[3]]+dense_plots[[3]])/
    (trace_plots[[4]]+dense_plots[[4]]) +
    plot_layout(guides = "collect") + plot_annotation(title = title)

  ggsave(filename=paste0(name,'-',district,'-',country,'.png'),plot=diag,dpi=300,height = 11,width=8,units = 'in')
  return(diag)
}
