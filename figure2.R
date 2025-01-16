library(R2OpenBUGS)
library(coda)
library(car)
library(MASS)
library(plotrix)
library(binom)
library(scales)
library(zoo)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(readxl)
source('utils.R')
library(patchwork)
library(bayesplot)
library(tidyverse)
library(excel.link)
library(zoo)
library(gridExtra)

##Copy over data needed for correlation
dir.create('src/process_nnp_data/data_raw')
dir.create('src/process_nnp_data/data_out')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/ANC-based surveillance/data_anc_mother_nigeria.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_child_nigeria_2020.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_hh_nigeria_2020.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_child_nigeria_2021.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_hh_nigeria_2021.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_child_nigeria_2022.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_hh_nigeria_2022.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Mozambique/ANC-based surveillance/ANC_Mozambique_2022.09_FINAL.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Mozambique/Cross-sectional survey/Moz Baseline CSS Datasets 2020.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Mozambique/Cross-sectional survey/Moz Midline CSS Datasets 2021.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Mozambique/Cross-sectional survey/Moz Endline CSS Datasets 2022.xlsx',
          to='src/process_nnp_data/data_raw/')

dir.create('src/process_mipmon_data/data_raw')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Mozambique/isglobal_cism_data/mipmon_merged.csv',
          to='src/process_mipmon_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Mozambique/isglobal_cism_data/cross_merged.csv',
          to='src/process_mipmon_data/data_raw/')

dir.create('src/process_vaneijk_data/data_raw')
file.copy(from='Q:/anc_pmcmc/nnp/data/paper_data_strat_G.xlsx',
          to='src/process_vaneijk_data/data_raw/')


##Process NNP data
nnp_data_id <- orderly2::orderly_run('process_nnp_data')

##Process mipmon data
mipmon_data_id <- orderly2::orderly_run('process_mipmon_data')

##Process Van Eijk data
vaneijk_data_id <- orderly2::orderly_run('process_vaneijk_data')

##Combine and process all data
combine_data_id <- orderly2::orderly_run('process_data_for_correlation')



N_sites=nrow(all_data_total4model)

####All pregnancies####
windows(height=5,width=10)
windows(height=5,width=3.3)
par(mar=c(5,5,1,1))

### now define the model in bugs terms to fit to the data
preg_child_model<-function(){
  for (i in 1:N) {
    pos_child[i] ~ dbin(p_child[i],total_child[i])
    pos_preg[i] ~ dbin(p_preg[i],total_preg[i])
    logit(p_child[i]) <- log_odds_child[i]
    logit(p_preg[i]) <- log_odds_child[i]+log_OR_p_v_c[i]
    log_odds_child[i]~dnorm(av_lo_child,sigma_c_inv)
    log_OR_p_v_c[i]<-RE_intercept[i]+gradient*(log_odds_child[i]-av_lo_child)
    RE_intercept[i]~dnorm(intercept,sigma_int_inv)
  }
  intercept~dnorm(0,0.001)
  gradient ~ dnorm(0,0.001)
  av_lo_child~dnorm(0,0.001)
  sigma_c_inv~dgamma(0.001,0.001)
  sigma_int_inv~dgamma(0.001,0.001)
  sigma_child<-1/sqrt(sigma_c_inv)
  sigma_intercept<-1/sqrt(sigma_int_inv)
}


## write to directory
model.file <- file.path(tempdir(),"model.txt")
write.model(preg_child_model, model.file)

## put in bugs data format
data_list<-list(pos_child=all_data_total4model$positive.cs,
                total_child=all_data_total4model$total.cs,
                pos_preg=all_data_total4model$positive.anc,
                total_preg=all_data_total4model$total.anc,
                N=N_sites)

## fn for random initial conditions
inits<-function(){
  list(av_lo_child=rnorm(1),intercept=rnorm(1),gradient=rnorm(1),sigma_c_inv=runif(1),sigma_int_inv=runif(1))
}
## params of interest
params<-c("av_lo_child","intercept","gradient","sigma_child","sigma_intercept")

## run the model
run_model <- bugs(data_list, inits, params,
                  model.file, n.iter=50000,n.burnin=10000)

##attach the output
attach.bugs(run_model)
### should recapture params - something slightly funny with gradient- need to check SD formats and priors
quantile(gradient,c(0.025,0.5,0.975))
quantile(intercept,c(0.025,0.5,0.975))
quantile(av_lo_child,c(0.025,0.5,0.975))
get_prev_from_log_odds(quantile(av_lo_child,c(0.025,0.5,0.975)))
## lets plot the relationship (apols for horrible base R)
prev_child=seq(0.01,0.99,by=0.01)
logodds_child=log(get_odds_from_prev(prev_child))
prev_preg_lower=array(dim=length(logodds_child))
prev_preg_upper=array(dim=length(logodds_child))
for(i in 1:length(logodds_child)){
  prev_preg_lower[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.025))
  prev_preg_upper[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.975))
}
prev_preg_median=get_prev_from_log_odds(logodds_child+median(gradient)*(logodds_child-median(av_lo_child))+median(intercept))
pred_all <- data.frame(prev_child,prev_preg_lower,prev_preg_upper,prev_preg_median)

##Try to make similar graph in ggplot
colors <- c(viridis(4,begin=0.2,end=0.9),"#D95F02")
"#414487FF" "#27808EFF" "#35B779FF" "#BBDF27FF"
colors_nosmc <- c(viridis(4,begin=0.2,end=0.9))

grav_all <- ggplot(all_data_total4fig)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_all,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_all,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='All pregnancies',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
windows(7,7)
grav_all
##plot without smc
grav_all_nosmc <- ggplot(all_data_total4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_all,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_all,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='All pregnancies',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
grav_all_nosmc

########################################
##Primigrav, Secundigrav and Multigrav simultaneously
########################################
N_sites=nrow(all_data_pgsgmg4model)

### now define the model in bugs terms to fit to the data
preg_pgsgmg_child_model<-function(){
  ##For each site (N)

  ##For no secundigrav
  for (i in 1:N) {
    #Probability distribution of prevalence measures
    pos_child[i] ~ dbin(p_child[i],total_child[i])
    pos_preg_pg[i] ~ dbin(p_preg_pg[i],total_preg_pg[i])
    pos_preg_mg[i] ~ dbin(p_preg_mg[i],total_preg_mg[i])

    #Log odds conversions from probability
    logit(p_child[i]) <- log_odds_child[i]
    logit(p_preg_pg[i]) <- log_odds_child[i]+log_OR_pp_v_c[i]
    logit(p_preg_mg[i]) <- logit(p_preg_pg[i])+log_OR_pm_v_pp[i]


    log_odds_child[i]~dnorm(av_lo_child,sigma_c_inv)

    #Primigrav-specific gradient
    log_OR_pp_v_c[i]<-RE_intercept_pg[i]+gradient_pg*(log_odds_child[i]-av_lo_child)
    #Multigrav-specific gradient
    log_OR_pm_v_pp[i]<-RE_intercept_mg[i]+gradient_mg*(log_odds_child[i]-av_lo_child)

    #Random-effects intercept for each gravidity type
    RE_intercept_pg[i]~dnorm(intercept_pg,sigma_int_inv)
    RE_intercept_mg[i]~dnorm(intercept_mg,sigma_int_inv)

    ##Add secundi
    #Prior on total number of secundigrav - based on unknown probability of a multigrav woman being a secundigrav
    # total_preg_sg[i] ~ dbin(prob_sg[i],total_preg_mg[i])
    # logit(prob_sg[i]) <- log_odds_sg[i]
    # log_odds_sg[i]~dnorm(av_lo_sg,sigma_sg_inv)

    pos_preg_sg[i] ~ dbin(p_preg_sg[i],total_preg_sg[i])##Make NA total_preg_sg to 1
    logit(p_preg_sg[i]) <- logit(p_preg_pg[i])+log_OR_ps_v_pp[i]
    log_OR_ps_v_pp[i]<-RE_intercept_sg[i]+gradient_sg*(log_odds_child[i]-av_lo_child)
    RE_intercept_sg[i]~dnorm(intercept_sg,sigma_int_inv)
  }
  intercept_pg~dnorm(0,0.001)
  intercept_sg~dnorm(0,0.001)
  intercept_mg~dnorm(0,0.001)
  gradient_pg ~ dnorm(0,0.001)
  gradient_sg ~ dnorm(0,0.001)
  gradient_mg ~ dnorm(0,0.001)
  # av_lo_child~dnorm(0,0.001)
  # av_lo_sg~dnorm(0,0.001)
  av_lo_child~dunif(-3.663562,3.663562)
  # av_lo_sg~dunif(-3.663562,3.663562)
  sigma_c_inv~dgamma(0.001,0.001)
  sigma_int_inv~dgamma(0.001,0.001)
  # sigma_sg_inv~dgamma(0.001,0.001)
  sigma_child<-1/sqrt(sigma_c_inv)
  sigma_intercept<-1/sqrt(sigma_int_inv)
  # sigma_sg<-1/sqrt(sigma_sg_inv)
}

## write to directory
model.file <- file.path(tempdir(),"model.txt")
write.model(preg_pgsgmg_child_model, model.file)

## put in bugs data format
data_list<-list(pos_child=all_data_pgsgmg4model$positive.cs,
                total_child=all_data_pgsgmg4model$total.cs,
                pos_preg_pg=all_data_pgsgmg4model$positive.anc.pg,
                total_preg_pg=all_data_pgsgmg4model$total.anc.pg,
                pos_preg_sg=all_data_pgsgmg4model$positive.anc.sg,
                total_preg_sg=ifelse(is.na(all_data_pgsgmg4model$total.anc.sg),1,all_data_pgsgmg4model$total.anc.sg),
                pos_preg_mg=all_data_pgsgmg4model$positive.anc.mg,
                total_preg_mg=all_data_pgsgmg4model$total.anc.mg,
                N=N_sites)

## fn for random initial conditions
inits<-function(){
  inits <- list(av_lo_child=rnorm(1),
                intercept_pg=rnorm(1),
                intercept_sg=rnorm(1),
                intercept_mg=rnorm(1),
                gradient_pg=rnorm(1),
                gradient_sg=rnorm(1),
                gradient_mg=rnorm(1),
                sigma_c_inv=runif(1),
                sigma_int_inv=runif(1))
  print(inits)
  return(inits)
}

## params of interest
params<-c("av_lo_child",
          "intercept_pg",
          "intercept_sg",
          "intercept_mg",
          "gradient_pg",
          "gradient_sg",
          "gradient_mg",
          "sigma_child",
          "sigma_intercept")

# av_lo_sg_sample <- rnorm(1000)
# sigma_sg_inv_sample <- runif(1000)
# log_odds_sg_sample <- sapply(1:1000, function(i){rnorm(1,av_lo_sg_sample[i],sigma_sg_inv_sample[i])})
# prob_sg_sample <- get_prev_from_log_odds(log_odds_sg_sample)
# total_preg_mg=all_data_pgsgmg4model$total.anc.mg
# i=1
# x=1
# total_preg_sg_sample <- bind_rows(lapply(1:1000,function(i){
#   c(sample = sapply(1:length(total_preg_mg),function(x,i){
#     sample_binom = rbinom(1,size=total_preg_mg[x],prob=prob_sg_sample[i])
#   },i=i))}))
#
# pos_preg_sg_sample <- bind_rows(lapply(1:1000,function(i){
#   c(sample = sapply(1:length(total_preg_mg),function(x,i){
#     sample_binom = rbinom(1,prob=prob_sg_sample[i],size=total_preg_sg_sample[[i,x]])
#   },i=i))}))
#
#
#
#
# log_odds_sg[i]~dnorm(av_lo_sg,sigma_sg_inv)
## run the model
run_model_pgsgmg <- bugs(data_list, inits, params,
                         model.file, n.iter=50000,n.burnin=10000)
run_model_pgsgmg_test <- bugs(data_list, inits, params,
                              model.file, n.iter=50000,n.burnin=10000)
pgsgmg_corr_sample <- sample_n(as.data.frame(run_model_pgsgmg$sims.matrix),1000)
saveRDS(pgsgmg_corr_sample,'./nnp/Corr/pgsgmg_mcmc_sample_290524.rds')
saveRDS(run_model_pgsgmg,'./nnp/Corr/pgsgmg_mcmc_run290524.rds')
run_model_pgsgmg_2 <- readRDS('./nnp/Corr/pgsgmg_mcmc_run070523.rds')
run_model_pgsgmg$summary
print(run_model_pgsgmg)
plot(run_model_pgsgmg)
traceplot(as.mcmc(run_model_pgsgmg$sims.matrix))
pgsgmg_mcmc <- as.mcmc(run_model_pgsgmg$sims.matrix)
traceplot(pgsgmg_mcmc)
get_prev_from_log_odds(run_model_pgsgmg$summary['av_lo_child',])
exp(run_model_pgsgmg$summary[c('intercept_pg','intercept_sg','intercept_mg'),c('2.5%','50%','97.5%')])
pgsgmg_mcmc_df <- as.data.frame(pgsgmg_mcmc)
ggplot(pgsgmg_mcmc_df)+
  geom_density(aes(x=get_prev_from_log_odds(av_lo_child)))
ggplot(pgsgmg_mcmc_df)+
  geom_density(aes(x=get_prev_from_log_odds(av_lo_sg)))
ggplot(pgsgmg_mcmc_df)+
  geom_density(aes(x=gradient_pg))
ggplot(pgsgmg_mcmc_df)+
  geom_density(aes(x=gradient_sg))
ggplot(pgsgmg_mcmc_df)+
  geom_density(aes(x=gradient_mg))
ggplot(pgsgmg_mcmc_df)+
  geom_density(aes(x=intercept_pg))
ggplot(pgsgmg_mcmc_df)+
  geom_density(aes(x=intercept_sg))
ggplot(pgsgmg_mcmc_df)+
  geom_density(aes(x=intercept_mg))
ggplot(pgsgmg_mcmc_df)+
  geom_density(aes(x=sigma_child))
ggplot(pgsgmg_mcmc_df)+
  geom_density(aes(x=sigma_intercept))
ggplot(pgsgmg_mcmc_df)+
  geom_density(aes(x=sigma_sg))

##attach the output
attach.bugs(run_model_pgsgmg)

### should recapture params - something slightly funny with gradient- need to check SD formats and priors
gradient_pg_quant <- quantile(gradient_pg,c(0.025,0.5,0.975))
gradient_sg_quant <- quantile(gradient_sg,c(0.025,0.5,0.975))
gradient_mg_quant <- quantile(gradient_mg,c(0.025,0.5,0.975))
intercept_pg_quant <- quantile(intercept_pg,c(0.025,0.5,0.975))
intercept_sg_quant <- quantile(intercept_sg,c(0.025,0.5,0.975))
intercept_mg_quant <- quantile(intercept_mg,c(0.025,0.5,0.975))
av_lo_child_quant <- quantile(av_lo_child,c(0.025,0.5,0.975))
get_prev_from_log_odds(quantile(av_lo_child,c(0.025,0.5,0.975)))
get_prev_from_log_odds(quantile(av_lo_sg,c(0.025,0.5,0.975)))

prev_child <- seq(0.01,0.99,by=0.01)
logodds_child <-logit(prev_child)

prev_preg_lower <- array(dim=length(logodds_child))
prev_preg_upper <- array(dim=length(logodds_child))
logodds_child[i]+intercept_pg+gradient_pg*(logodds_child[1]-av_lo_child)
ggplot()+
  geom_density(aes(x=prev_preg_pg),color='darkblue')+
  geom_density(aes(x=prev_preg_sg),color='darkred')+
  geom_density(aes(x=prev_preg_mg),color='darkgreen')

ggplot()+
  geom_density(aes(x=prev_preg_sg))
i <- 1
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

mcmc_sim_summary_4lucy <- mcmc_sim_summary[,1:10]
saveRDS(mcmc_sim_summary_4lucy,'./nnp/Corr/mcmc_sim_summary_4lucy.rds')
colors_pgsgmg <- c(viridis(3,begin=0.1,end=0.9))
names(colors_pgsgmg) <- c('Primigrav','Secundigrav','Multigrav')
compare_prevs <- ggplot(mcmc_sim_summary)+
  geom_abline(color='darkgrey',linetype='dashed',size=1)+
  geom_line(aes(x=prev_child,y=prev_preg_pg_median,color='Primigrav'),size=1)+
  geom_line(aes(x=prev_child,y=prev_preg_sg_median,color='Secundigrav'),size=1)+
  geom_line(aes(x=prev_child,y=prev_preg_mg_median,color='Multigrav'),size=1)+
  scale_color_manual(name='Gravidity',values=colors_pgsgmg,breaks=names(colors_pgsgmg))+
  labs(x='Cross-section Prevalence (<5 yo)',y='Median ANC Prevalence')
windows(7,7)
ggsave(plot=compare_prevs,filename = './nnp/Corr/pgsgmg_prevdiff_070524.tiff',width=6,height=6,units='in')
colors_pgsgmg_or <- c(viridis(3,begin=0.1,end=0.9))
names(colors_pgsgmg_or) <- c('PG vs child','SG vs PG','MG vs PG')
compare_ors <- ggplot(mcmc_sim_summary)+
  geom_hline(yintercept = 1.0,linetype='dashed',size=1,color='darkgrey')+
  geom_line(aes(x=prev_child,y=exp(log_OR_pp_v_c_median),color='PG vs child'),size=1)+
  geom_line(aes(x=prev_child,y=exp(log_OR_ps_v_pp_median),color='SG vs PG'),size=1)+
  geom_line(aes(x=prev_child,y=exp(log_OR_pm_v_pp_median),color='MG vs PG'),size=1)+
  scale_color_manual(name='Gravidity',values=colors_pgsgmg_or,breaks=names(colors_pgsgmg_or))+
  # scale_y_log10(limits=c(0.3,NA))+
  scale_y_continuous(limits=c(0,NA))+
  labs(x='Cross-section Prevalence (<5 yo)',y='Median Odds Ratio')
windows(7,7)
compare_ors
ggsave(plot=colors_pgsgmg_or,filename = './nnp/Corr/pgsgmg_ordiff_070524.tiff',width=6,height=6,units='in')
prev_preg_pg <- bind_cols(lapply(1:length(logodds_child),function(i){
  logodds_child[i]+intercept_pg+gradient_pg*(logodds_child[i]-av_lo_child)
}))
for(i in 1:length(logodds_child)){
  prev_preg_lower[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept_pg+gradient_pg*(logodds_child[i]-av_lo_child),0.025))
  prev_preg_upper[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept_pg+gradient_pg*(logodds_child[i]-av_lo_child),0.975))
}
prev_preg_median=get_prev_from_log_odds(logodds_child+median(gradient_pg)*(logodds_child-median(av_lo_child))+median(intercept_pg))
pred_pgmg_pg <- data.frame(prev_child,prev_preg_lower,prev_preg_upper,prev_preg_median)
prev_preg_lower=array(dim=length(logodds_child))
prev_preg_upper=array(dim=length(logodds_child))
p_preg_mg[i] <- logit(p_preg_pg[i])+RE_intercept_mg[i]+gradient_mg*(log_odds_child[i]-av_lo_child)


for(i in 1:length(logodds_child)){
  prev_preg_lower[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept_mg+gradient_mg*(logodds_child[i]-av_lo_child),0.025))
  prev_preg_upper[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept_mg+gradient_mg*(logodds_child[i]-av_lo_child),0.975))
}
prev_preg_median=get_prev_from_log_odds(logodds_child+median(gradient_mg)*(logodds_child-median(av_lo_child))+median(intercept_mg))
pred_pgmg_mg <- data.frame(prev_child,prev_preg_lower,prev_preg_upper,prev_preg_median)

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
windows(7,7)
gravpgsgmg_pg

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
windows(7,7)
gravpgsgmg_mg

gravpgsgmg_sg <- ggplot(all_data_pgsgmg4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc.sg*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc.sg*100,ymax=upper.anc.sg*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc.sg*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=mcmc_sim_summary,aes(x=prev_child*100,ymin=prev_preg_sg_lower*100,ymax=prev_preg_sg_upper*100),alpha=0.2)+
  geom_line(data=mcmc_sim_summary,aes(x=prev_child*100,y=prev_preg_sg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='Second pregnancy',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
windows(7,7)
gravpgsgmg_sg

gravpgsgmg_pgsg <- ggplot(all_data_pgsgmg4model)+
  geom_point(aes(x=mean.anc.pg*100,y=mean.anc.sg*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.anc.pg*100,ymin=lower.anc.sg*100,ymax=upper.anc.sg*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc.sg*100,xmin=lower.anc.pg*100,xmax=upper.anc.pg*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=mcmc_sim_summary,aes(x=prev_preg_pg_median*100,ymin=prev_preg_sg_lower*100,ymax=prev_preg_sg_upper*100),alpha=0.2)+
  geom_line(data=mcmc_sim_summary,aes(x=prev_preg_pg_median*100,y=prev_preg_sg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='PG v SG pregnancy',x='PG ANC prevalence',y='SG ANC Prevalence')

gravpgsgmg_pgmg <- ggplot(all_data_pgsgmg4model)+
  geom_point(aes(x=mean.anc.pg*100,y=mean.anc.mg*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.anc.pg*100,ymin=lower.anc.mg*100,ymax=upper.anc.mg*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc.mg*100,xmin=lower.anc.pg*100,xmax=upper.anc.pg*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=mcmc_sim_summary,aes(x=prev_preg_pg_median*100,ymin=prev_preg_mg_lower*100,ymax=prev_preg_mg_upper*100),alpha=0.2)+
  geom_line(data=mcmc_sim_summary,aes(x=prev_preg_pg_median*100,y=prev_preg_mg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='PG v MG pregnancy',x='PG ANC prevalence',y='MG ANC Prevalence')

diff_gravs <- gravpgsgmg_pgsg + gravpgsgmg_pgmg + plot_layout(ncol=2,guides='collect') & theme(legend.position = 'bottom')
ggsave('Corr/correlation_gravsdiff_070524.tiff',plot=diff_gravs,units='cm',height=12,width=20)

diff_pgsgmg_model <- gravpgsgmg_pg + gravpgsgmg_sg+ gravpgsgmg_mg+ plot_layout(ncol=3,guides = 'collect')
ggsave('Corr/correlation_pgsgmgdiff_070524.tiff',plot=diff_pgsgmg_model,units='cm',height=12,width=30)
ggsave('Corr/correlation_pgsgmgdiff_data_070524.tiff',plot=diff_pgsgmg_model,units='cm',height=12,width=30)

diff_pgsgmg_pgmg <- gravpgsgmg_pg + gravpgsgmg_mg+ plot_layout(ncol=2,guides = 'collect')
ggsave('Corr/correlation_pgmgdiff_data_070524.tiff',plot=diff_pgsgmg_pgmg,units='cm',height=12,width=20)
ggsave('Corr/correlation_pgmgdiff_070524.tiff',plot=diff_pgsgmg_pgmg,units='cm',height=12,width=20)


gravpgmg_mg <- ggplot(all_data_pgsgmg4model)+
  geom_point(aes(x=mean.anc.pg*100,y=mean.anc.mg*100,col=country),size=3)+
  geom_text(aes(x=mean.anc.pg*100,y=mean.anc.mg*100,label=site),hjust=0.1, vjust=0.1)+
  geom_errorbar(aes(x=mean.anc.pg*100,ymin=lower.anc.mg*100,ymax=upper.anc.mg*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc.mg*100,xmin=lower.anc.pg*100,xmax=upper.anc.pg*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(x='Primigrav prevalence',y='Multigrav prevalence')
windows(7,7)
gravpgmg_mg


gravpgmg_sg <- ggplot(all_data_pgsgmg4model)+
  geom_point(aes(x=mean.anc.pg*100,y=mean.anc.sg*100,col=country),size=3)+
  geom_text(aes(x=mean.anc.pg*100,y=mean.anc.sg*100,label=site),hjust=0.1, vjust=0.1)+
  geom_errorbar(aes(x=mean.anc.pg*100,ymin=lower.anc.sg*100,ymax=upper.anc.sg*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc.sg*100,xmin=lower.anc.pg*100,xmax=upper.anc.pg*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(x='Primigrav prevalence',y='Secundigrav prevalence')
windows(7,7)
gravpgmg_sg

gravpgmg_sgmg <- ggplot(all_data_pgsgmg4model)+
  geom_point(aes(x=mean.anc.sg*100,y=mean.anc.mg*100,col=country),size=3)+
  geom_text(aes(x=mean.anc.sg*100,y=mean.anc.mg*100,label=site),hjust=0.1, vjust=0.1)+
  geom_errorbar(aes(x=mean.anc.sg*100,ymin=lower.anc.mg*100,ymax=upper.anc.mg*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc.mg*100,xmin=lower.anc.sg*100,xmax=upper.anc.sg*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(x='Secundigrav prevalence',y='Multigrav prevalence')
windows(7,7)
gravpgmg_sgmg
diff_pgsgmg <- gravpgmg_mg + gravpgmg_sg+ gravpgmg_sgmg+ plot_layout(ncol=3,guides = 'collect')
ggsave('Corr/correlation_pgsgmgdiff_190224.pdf',plot=diff_pgsgmg,height=5.5,width=13)
diff_pgsgmg <- gravpgmg_mg + gravpgmg_sg+ gravpgmg_sgmg+ plot_layout(ncol=3,guides = 'collect')
ggsave('Corr/correlation_pgsgmgdiff_190224_wlabs.pdf',plot=diff_pgsgmg,height=5.5,width=13)
diff_pgsgmg <- gravpgmg_mg + gravpgmg_sg+ gravpgmg_sgmg+ plot_layout(ncol=3,guides = 'collect')
ggsave('Corr/correlation_pgsgmgdiff_190224_wlabs2.pdf',plot=diff_pgsgmg,height=5.5,width=13)

or_diff_pgmg <- gravpgmg_pg + gravpgmg_mg+
  plot_annotation(title = 'OR difference between PG and MG') + plot_layout(ncol=2,guides = 'collect')
ggsave('Corr/correlation_ordiff_130224.pdf',plot=or_diff_pgmg,height=5.5,width=30)
