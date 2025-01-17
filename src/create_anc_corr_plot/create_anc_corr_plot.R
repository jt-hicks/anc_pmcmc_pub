orderly2::orderly_dependency(name='run_anc_corr',query = 'latest()',
                             files = c('pgsgmg_mcmc_run.rds'))
orderly2::orderly_artefact(files='mcmc_sim_summary.rds',description='MCMC summary')
orderly2::orderly_artefact(files='correlation_pgmgdiff_plot.tiff',description='Correlation plot')


run_model_pgsgmg <- readRDS('pgsgmg_mcmc_run.rds')
##attach the output
attach.bugs(run_model_pgsgmg)

prev_child <- seq(0.01,max(all_data_pgsgmg4model$upper.cs),by=0.01)
logodds_child <-logit(prev_child)

mcmc_sim_summary <- dplyr::bind_rows(lapply(1:length(logodds_child),function(i){
  prev_preg_pg <- get_prev_from_log_odds(logodds_child[i]+intercept_pg+gradient_pg*(logodds_child[i]-av_lo_child))
  prev_preg_sg <- get_prev_from_log_odds(logit(prev_preg_pg)+intercept_sg+gradient_sg*(logodds_child[i]-av_lo_child))
  prev_preg_mg <- get_prev_from_log_odds(logit(prev_preg_pg)+intercept_mg+gradient_mg*(logodds_child[i]-av_lo_child))
  #Primigrav-specific gradient
  log_OR_pp_v_c <-intercept_pg+gradient_pg*(logodds_child[i]-av_lo_child)
  #Secundigrav-specific gradient
  log_OR_ps_v_pp<-intercept_sg+gradient_sg*(logodds_child[i]-av_lo_child)
  #Multigrav-specific gradient
  log_OR_pm_v_pp <-intercept_mg+gradient_mg*(logodds_child[i]-av_lo_child)

  prev_preg_pg_quant <- quantile(prev_preg_pg,c(0.025,0.5,0.975))
  prev_preg_sg_quant <- quantile(prev_preg_sg,c(0.025,0.5,0.975))
  prev_preg_mg_quant <- quantile(prev_preg_mg,c(0.025,0.5,0.975))
  log_OR_pp_v_c_quant <- quantile(log_OR_pp_v_c,c(0.025,0.5,0.975))
  log_OR_ps_v_pp_quant <- quantile(log_OR_ps_v_pp,c(0.025,0.5,0.975))
  log_OR_pm_v_pp_quant <- quantile(log_OR_pm_v_pp,c(0.025,0.5,0.975))

  data.frame(prev_child = get_prev_from_log_odds(logodds_child[i]),
             prev_preg_pg_median = prev_preg_pg_quant[[2]],
             prev_preg_pg_lower = prev_preg_pg_quant[[1]],
             prev_preg_pg_upper = prev_preg_pg_quant[[3]],
             prev_preg_sg_median = prev_preg_sg_quant[[2]],
             prev_preg_sg_lower = prev_preg_sg_quant[[1]],
             prev_preg_sg_upper = prev_preg_sg_quant[[3]],
             prev_preg_mg_median = prev_preg_mg_quant[[2]],
             prev_preg_mg_lower = prev_preg_mg_quant[[1]],
             prev_preg_mg_upper = prev_preg_mg_quant[[3]],
             log_OR_pp_v_c_median = log_OR_pp_v_c_quant[[2]],
             log_OR_pp_v_c_lower = log_OR_pp_v_c_quant[[1]],
             log_OR_pp_v_c_upper = log_OR_pp_v_c_quant[[3]],
             log_OR_ps_v_pp_median = log_OR_ps_v_pp_quant[[2]],
             log_OR_ps_v_pp_lower = log_OR_ps_v_pp_quant[[1]],
             log_OR_ps_v_pp_upper = log_OR_ps_v_pp_quant[[3]],
             log_OR_pm_v_pp_median = log_OR_pm_v_pp_quant[[2]],
             log_OR_pm_v_pp_lower = log_OR_pm_v_pp_quant[[1]],
             log_OR_pm_v_pp_upper = log_OR_pm_v_pp_quant[[3]])
}))

saveRDS(mcmc_sim_summary,'mcmc_sim_summary.rds')

colors_nosmc <- c(viridis(4,begin=0.2,end=0.9))

gravpgsgmg_pg <- ggplot(all_data_pgsgmg4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc.pg*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc.pg*100,ymax=upper.anc.pg*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc.pg*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=mcmc_sim_summary,aes(x=prev_child*100,ymin=prev_preg_pg_lower*100,ymax=prev_preg_pg_upper*100),alpha=0.2)+
  geom_line(data=mcmc_sim_summary,aes(x=prev_child*100,y=prev_preg_pg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='First pregnancy',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')


gravpgsgmg_mg <- ggplot(all_data_pgsgmg4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc.mg*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc.mg*100,ymax=upper.anc.mg*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc.mg*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=mcmc_sim_summary,aes(x=prev_child*100,ymin=prev_preg_mg_lower*100,ymax=prev_preg_mg_upper*100),alpha=0.2)+
  geom_line(data=mcmc_sim_summary,aes(x=prev_child*100,y=prev_preg_mg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='Second or later pregnancy',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')




diff_pgsgmg_pgmg <- gravpgsgmg_pg + gravpgsgmg_mg+ plot_layout(ncol=2,guides = 'collect')*theme(legend.position = 'bottom')
ggsave('correlation_pgmgdiff_plot.tiff',plot=diff_pgsgmg_pgmg,units='cm',height=12,width=20)


