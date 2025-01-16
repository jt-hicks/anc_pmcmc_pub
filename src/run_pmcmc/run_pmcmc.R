orderly2::orderly_parameters(name=NULL,
                             proposal_matrix=1,
                             length=NULL,
                             workers=1,
                             chain=1,
                             seed=1L,
                             start_pf_time=360,
                             misspecification = 0)
orderly2::orderly_dependency("create_sim_data", quote(latest()),
                             c('sim_seasonal_dataraw_list.RDS'))
data_list <- readRDS('sim_seasonal_dataraw_list.RDS')

if(proposal_matrix!=1){
  orderly2::orderly_dependency(name="run_diagnostics", query='latest(parameter:length==1000)',
                               c("prop.rds" = "props.RDS"))
  prop_file <- readRDS('prop.rds')
  proposal_matrix=prop_file[name,'prop']
}
##Calculate first year prevalence and convert to children under 5
first_annual_prev <- data_list[[name]]%>%
  slice_head(n=12)

annual_prev = (sum(first_annual_prev$positive)/sum(first_annual_prev$tested))+misspecification

# hipercow::hipercow_init()
# hipercow::hipercow_configure('windows')
# hipercow::hipercow_provision()
# hipercow::hipercow_environment_create(packages=c('dplyr','ggplot2'))
# resources <- hipercow::hipercow_resources(cores=1)

result <- mamasante::run_pmcmc(data_raw=data_list[[name]],
                               init_EIR = 100,
                               n_particles=200,
                               proposal_matrix = matrix(proposal_matrix),
                               target_prev = annual_prev,
                               target_prev_group = 'u5',
                               max_param=125,
                               prop_treated = 0.4,
                               n_steps = length,
                               n_threads = hipercow::hipercow_parallel_get_cores(),
                               n_chains = chain,
                               n_workers = workers,
                               state_check = 0,## Run equilibrium checks
                               seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                               seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                               seed = seed,
                               start_pf_time = start_pf_time,
                               particle_tune = FALSE,
                               comparison = 'u5',
                               initial = 'informed')

saveRDS(result,'result.rds')
