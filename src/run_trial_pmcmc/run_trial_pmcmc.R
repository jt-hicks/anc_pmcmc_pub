orderly2::orderly_shared_resource('get_u5prev_fromanc.R')
orderly2::orderly_parameters(name=NULL,
                             proposal_matrix=1,
                             comparison_trial=NULL,
                             length=NULL,
                             workers=1,
                             chain=1,
                             seed=1L)
orderly2::orderly_resource('trials_pg_data_all_list.RDS')
orderly2::orderly_resource('trials_mg_data_all_list.RDS')
orderly2::orderly_artefact(description='pmcmc results',file='result.rds')

source('get_u5prev_fromanc.R')
trials_pg_data_all_list <- readRDS('trials_pg_data_all_list.RDS')
trials_mg_data_all_list <- readRDS('trials_mg_data_all_list.RDS')

print(name)

##Calculate first year prevalence and convert to children under 5
first_annual_prev <- trials_pg_data_all_list[[name]]%>%
  slice_head(n=12)
annual_prev = sum(first_annual_prev$positive)/sum(first_annual_prev$tested)

result <- mamasante::run_pmcmc(data_raw_pg=trials_pg_data_all_list[[name]],
                  data_raw_mg=trials_mg_data_all_list[[name]],
                  init_EIR = 100,
                  n_particles=200,
                  proposal_matrix = matrix(proposal_matrix),
                  target_prev = annual_prev,
                  target_prev_group='u5',
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
                  start_pf_time = 30*12,
                  particle_tune = FALSE,
                  comparison = comparison_trial,
                  initial = 'informed',
                  check_flexibility=TRUE)

saveRDS(result,'result.rds')
