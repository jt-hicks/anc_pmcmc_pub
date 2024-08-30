library(SBC)

source('utils.R')
theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom',
                  title = element_text(size=12)))

cache_dir <- "./_implementing_backends_SBC_cache"
if(!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}

set.seed(54882236)
n_sims <- 100  # Number of SBC iterations to run

sim_generator <- SBC_generator_function(datasim4pmcmc, N = 12)
sim_dataset <- generate_datasets(
  sim_generator,
  n_sims)
# sim_dataset_100 <- sim_dataset
sim_dataset_100$variables
sim_dataset$generated
sim_dataset$generated

sim_dataset_100_combined <- bind_rows(lapply(1:100,function(i){
  df <- data.frame(sim_dataset_100$generated[i])
  df$draw <- i
  df$volatility <- sim_dataset_100$variables[[i,'volatility']]
  df$init_EIR <- sim_dataset_100$variables[[i,'init_EIR']]
  return(df)
}))

sim_dataset_100_2_combined <- bind_rows(lapply(1:100,function(i){
  df <- data.frame(sim_dataset$generated[i])
  df$draw <- i
  df$volatility <- sim_dataset$variables[[i,'volatility']]
  df$init_EIR <- sim_dataset$variables[[i,'init_EIR']]
  return(df)
}))

seasonal_peak_example <- data.frame(sim_dataset$generated[11])
seasonal_peak_example$betaa_true

ggplot(sim_dataset_100_2_combined)+
  geom_point(aes(x=month,y=mv))+
  facet_wrap(.~draw,scales = 'free_y')

SBC_backend_pmcmc <- function(...) {
  message('SBC_backend_pmcmc')
  args = list(...)
  if(any(names(args) == "data_raw")) {
    stop(paste0("Argument 'data' cannot be provided when defining a backend",
                " as it needs to be set by the SBC package"))
  }

  structure(list(args = args), class = "SBC_backend_pmcmc")
}

SBC_fit.SBC_backend_pmcmc <- function(backend, generated, cores) {
  message('SBC_fit.SBC_backend_pmcmc')
  args_with_data <- backend$args
  args_with_data$data_raw <- generated

  fit <- do.call(run_pmcmc, args_with_data)
  structure(fit,class = 'SBC_backend_pmcmc')
}

SBC_fit_to_draws_matrix.SBC_backend_pmcmc <- function(fit) {
  message('SBC_fit_to_draws_matrix.pmcmc')
  posterior::as_draws_matrix(fit$mcmc)
}

backend_pmcmc <- SBC_backend_pmcmc(init_EIR = 100,
                                   n_particles=200,
                                   proposal_matrix = diag(0.5,2),
                                   max_param=125,
                                   prop_treated = 0.4,
                                   n_steps = 5000,
                                   n_threads = 8,
                                   state_check = 0,## Run equilibrium checks
                                   seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                                   seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                                   seed = 1L,
                                   start_pf_time = 30,
                                   particle_tune = FALSE,
                                   comparison = 'u5',
                                   initial = 'fitted')
SBC_backend_iid_draws.SBC_backend_pmcmc <- function(backend) {
  FALSE
}
res_pmcmc <- compute_SBC(sim_dataset,
                           backend_pmcmc,
                           thin_ranks = 1,
                           cache_mode = "results",
                           cache_location = file.path(cache_dir,"test")
                         # cache_mode = "none"
)
res_pmcmc$default_diagnostics
res_pmcmc$fits[[1]]$run_time
res_pmcmc$messages
mcmc_pairs(res_pmcmc$fits[[1]]$mcmc)
mcmc_trace(res_pmcmc$fits[[1]]$mcmc)
recompute_SBC_statistics(res_pmcmc,datasets = sim_dataset)
sim_1 <- sim_dataset$generated[[1]] %>%
  mutate(date=as.Date(date))
history.df.prev <- as.data.frame(t(res_pmcmc$fits[[1]]$history[1, 1001:5000, -1]))
prev_history <- history.df.prev%>%
  dplyr::mutate(date=sim_dataset$generated[[1]]$date)%>%
  melt(id='date')%>%
  group_by(date)%>%
  dplyr::summarise(prev.median=median(value),
                   prev.mean=mean(value),
                   prev.upper=quantile(value,probs=0.975),
                   prev.lower=quantile(value,probs=0.025))%>%
  mutate(date=as.Date(date))
ggplot()+
  geom_point(data=sim_1,aes(x=date,y=prev_all_true))+
  geom_line(data = prev_history, aes(x=date,y=prev.median),color='darkgrey',linewidth=1)+
  geom_ribbon(data=prev_history,aes(x=date,ymin=prev.lower,ymax=prev.upper),fill='lightgrey',alpha=0.5)+
  coord_cartesian(ylim = c(0,NA))

history.df.incall <- as.data.frame(t(res_pmcmc$fits[[1]]$history[4, 1001:5000, -1]))
inc_history <- history.df.incall%>%
  dplyr::mutate(date=sim_1$date)%>%
  melt(id='date')%>%
  group_by(date)%>%
  dplyr::summarise(incall.median=median(value),
                   incall.mean=mean(value),
                   incall.upper=quantile(value,probs=0.975),
                   incall.lower=quantile(value,probs=0.025))
ggplot()+
  geom_point(data=sim_1,aes(x=date,y=inc_all_true))+
  geom_line(data = inc_history, aes(x=date,y=incall.median),color='darkgrey')+
  geom_ribbon(data=inc_history,aes(x=date,ymin=incall.lower,ymax=incall.upper),fill='lightgrey',alpha=0.5)+
  coord_cartesian(ylim = c(0,NA))

history.df.mv <- as.data.frame(t(res_pmcmc$fits[[1]]$history[16, 1001:5000, -1]))
mv_history <- history.df.mv%>%
  dplyr::mutate(date=sim_1$date)%>%
  melt(id='date')%>%
  group_by(date)%>%
  dplyr::summarise(mv.median=median(value),
                   mv.mean=mean(value),
                   mv.upper=quantile(value,probs=0.975),
                   mv.lower=quantile(value,probs=0.025))
ggplot()+
  geom_point(data=sim_1,aes(x=date,y=mv))+
  geom_line(data = mv_history, aes(x=date,y=mv.median),color='darkgrey')+
  geom_ribbon(data=mv_history,aes(x=date,ymin=mv.lower,ymax=mv.upper),fill='lightgrey',alpha=0.5)+
  coord_cartesian(ylim = c(0,NA))

history.df.spzrate <- as.data.frame(t(res_pmcmc$fits[[1]]$history[15, 1001:5000, -1]))
spz_history <- history.df.spzrate%>%
  dplyr::mutate(date=sim_1$date)%>%
  melt(id='date')%>%
  group_by(date)%>%
  dplyr::summarise(spzrate.median=median(value),
                   spzrate.mean=mean(value),
                   spzrate.upper=quantile(value,probs=0.975),
                   spzrate.lower=quantile(value,probs=0.025))
ggplot()+
  geom_point(data=sim_1,aes(x=date,y=spz_rate))+
  geom_line(data = spz_history, aes(x=date,y=spzrate.median),color='darkgrey')+
  geom_ribbon(data=spz_history,aes(x=date,ymin=spzrate.lower,ymax=spzrate.upper),fill='lightgrey',alpha=0.5)+
  coord_cartesian(ylim = c(0,NA))
t(res_pmcmc$fits[[1]]$history[, 5000, -1])
