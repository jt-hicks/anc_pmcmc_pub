source('utils.R')
source('submit_sbc.R')

obj_sbc <- cluster_setup(context_name = 'sbc',template='20Core', cores=10)

+install.packages("didehpc")
test <- sbc_submit(sim_dataset)
test1 <- obj_sbc$enqueue(sbc_submit(sim_dataset))
obj_sbc$task_list()

test1$status()
test1$log()
# obj_sbc$unsubmit("a360c99d0a175c81d419ed302cfccb67")
test1_result <- test1$result()
test1_result$messages
test1_result$errors
test1_result$outputs
obj_sbc$login()
test_data <- obj_sbc$enqueue(print(sim_data))
test_data$status()
test_data$result()
test_data_source <- sim_dataset$generated[[1]]
test_workers <- obj_sbc$enqueue(sifter::run_pmcmc(data_raw=test_data_source,
                                          init_EIR = 100,
                                          n_particles=10,
                                          proposal_matrix = diag(0.5,2),
                                          max_param=125,
                                          prop_treated = 0.4,
                                          n_steps = 10,
                                          n_threads = 8,
                                          n_chains = 2,
                                          n_workers = 2,
                                          state_check = 0,## Run equilibrium checks
                                          seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                                          seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                                          seed = 1L,
                                          start_pf_time = 30,
                                          particle_tune = FALSE,
                                          comparison = 'u5',
                                          initial = 'fitted'))
test_workers$status()
test_workers$log()
test_workers_result <- test_workers$result()
test_workers_result$mcmc

sim10_runs <- obj_sbc$enqueue_bulk(1:10,function(x,data){
  generated_data <-  data$generated[[x]]
  sifter::run_pmcmc(data_raw=generated_data,
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = diag(0.5,2),
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 8,
                    n_chains = 2,
                    n_workers = 2,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'fitted')
},data=sim_dataset_100)
sim10_runs$status() #cynophobic_drafthorse
result_1 <- sim10_runs$tasks[[1]]$result()
mcmc_pairs(result_1$mcmc)
mcmc_trace(result_1$mcmc)

result_2 <- sim10_runs$tasks[[2]]$result()
mcmc_pairs(result_2$mcmc)
mcmc_trace(result_2$mcmc)
hist2 <- result_2$history

obj_sbc$config
sim20_runs <- obj_sbc$enqueue_bulk(11:30,function(x,data){
  generated_data <-  data$generated[[x]]
  sifter::run_pmcmc(data_raw=generated_data,
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = diag(0.5,2),
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 8,
                    n_chains = 2,
                    n_workers = 2,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'fitted')
},data=sim_dataset_100)
sim20_runs$status() #concretionary_earwig

sim50_runs <- obj_sbc$enqueue_bulk(31:50,function(x,data){
  generated_data <-  data$generated[[x]]
  sifter::run_pmcmc(data_raw=generated_data,
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = diag(0.5,2),
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 8,
                    n_chains = 2,
                    n_workers = 2,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'fitted')
},data=sim_dataset_100)
sim50_runs$status() #sticky_agouti

sim_seasonal_run <- obj_sbc$enqueue_bulk(1,function(x,data){
  sifter::run_pmcmc(data_raw=data,
                    init_EIR = 100,
                    n_particles=200,
                    proposal_matrix = diag(0.5,2),
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 8,
                    n_chains = 1,
                    n_workers = 1,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'fitted')
},data=four_years)
sim_seasonal_run$status() #foolish_arrowworm
sim_seasonal_run_result <- sim_seasonal_run$tasks[[1]]$result()
sim_seasonal_run_result <- readRDS('T:/jth/sbc/db/data/790c9c8250091785442fa3c18f0fb266.rds')
sim_seasonal_result()
sim_seasonal_run_result
ar <- as.data.frame(t(1 - coda::rejectionRate(as.mcmc(sim_seasonal_run_result$mcmc))))
ess <- as.data.frame(t(coda::effectiveSize(as.mcmc(sim_seasonal_run_result$mcmc[101:1000,]))))
mcmc_trace(sim_seasonal_run_result$mcmc)

data_t <- diff(sim_data$t)
results_t <- diff(results$times)

sim10_informed_runs <- obj_sbc$enqueue_bulk(1:10,function(x,data){
  generated_data <-  data$generated[[x]]
  sifter::run_pmcmc(data_raw=generated_data,
                    init_EIR = data$variables[[x,'init_EIR']],
                    n_particles=200,
                    proposal_matrix = diag(0.5,2),
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 8,
                    n_chains = 2,
                    n_workers = 2,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed')
},data=sim_dataset_100)
sim10_informed_runs$status() #ununpentium_eider

sim10_informed_longerflex_runs <- obj_sbc$enqueue_bulk(1:10,function(x,data){
  generated_data <-  data$generated[[x]]
  sifter::run_pmcmc(data_raw=generated_data,
                    init_EIR = data$variables[[x,'init_EIR']],
                    n_particles=200,
                    proposal_matrix = diag(0.5,2),
                    max_param=125,
                    prop_treated = 0.4,
                    n_steps = 1000,
                    n_threads = 8,
                    n_chains = 2,
                    n_workers = 2,
                    state_check = 0,## Run equilibrium checks
                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30*4,
                    particle_tune = FALSE,
                    comparison = 'u5',
                    initial = 'informed')
},data=sim_dataset_100)
sim10_informed_longerflex_runs$status() #loveydovey_africanwilddog

##Sim for illustration fit, using hipercow##
library(hipercow)
hipercow_init()
hipercow_configure(driver='windows')
id <- task_create_expr(sessionInfo())
id
windows_authenticate()
task_status(id)
task_result(id)
hipercow_provision()
resources <- hipercow_resources(cores=32)
sim_run_id <- task_create_expr(mamasante::run_pmcmc(data_raw=four_years,
                                                    init_EIR = 100,
                                                    n_particles=200,
                                                    proposal_matrix = diag(0.5,2),
                                                    max_param=125,
                                                    prop_treated = 0.4,
                                                    n_steps = 100,
                                                    n_threads = 8,
                                                    n_chains = 1,
                                                    n_workers = 1,
                                                    state_check = 0,## Run equilibrium checks
                                                    seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                                                    seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                                                    seed = 2507,
                                                    start_pf_time = 30,
                                                    particle_tune = FALSE,
                                                    comparison = 'u5',
                                                    initial = 'fitted'),
                               parallel = hipercow_parallel('parallel',cores_per_process = 32),
                               resources=resources)
# task_cancel(sim_run_id)
task_info(sim_run_id)
task_result(sim_run_id)
task_log_show(sim_run_id)
task_log_watch(sim_run_id)
sim_seasonal_run_result <- task_result(sim_run_id)

ar <- as.data.frame(t(1 - coda::rejectionRate(as.mcmc(sim_seasonal_run_result$mcmc))))
ess <- as.data.frame(t(coda::effectiveSize(as.mcmc(sim_seasonal_run_result$mcmc[101:1000,]))))
mcmc_trace(sim_seasonal_run_result$mcmc)
plot(four_years$t,four_years$prev05)
