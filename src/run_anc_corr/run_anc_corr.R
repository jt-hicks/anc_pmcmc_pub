orderly2::orderly_dependency(name='process_data_for_correlation',query = 'latest()',
                             files = c('all_data_pgsgmg4model.rds'))
orderly2::orderly_artefact(files='model.txt',description = 'Regression model')
orderly2::orderly_artefact(files='pgsgmg_mcmc_run.rds',description = 'Regression model output')
orderly2::orderly_artefact(files=c('trace_plot_1.png','trace_plot_2.png','trace_plot_3.png','trace_plot_4.png'),description = 'Regression model trace plots')

all_data_pgsgmg4model <- readRDS('all_data_pgsgmg4model.rds')

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
model.file <- ("model.txt")
R2OpenBUGS::write.model(preg_pgsgmg_child_model, model.file)

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


run_model_pgsgmg <- bugs(data_list, inits, params,
                         model.file, n.iter=50000,n.burnin=10000)
png('trace_plot_%d.png')
trace_plot <- plot(coda::as.mcmc(run_model_pgsgmg$sims.matrix))
dev.off()

saveRDS(run_model_pgsgmg,'pgsgmg_mcmc_run.rds')



