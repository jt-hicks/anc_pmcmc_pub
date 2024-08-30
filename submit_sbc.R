sbc_submit <- function(sim_data){
  # options(mc.cores = parallel::detectCores())
  # plan(multisession)
  print(getwd())
  cache_dir <- "./_implementing_backends_SBC_cache"
  if(!dir.exists(cache_dir)) {
    dir.create(cache_dir)
  }

  SBC_backend_pmcmc <- function(...) {
    args = list(...)
    if(any(names(args) == "data_raw")) {
      stop(paste0("Argument 'data' cannot be provided when defining a backend",
                  " as it needs to be set by the SBC package"))
    }

    structure(list(args = args), class = "SBC_backend_pmcmc")
  }

  SBC_fit.SBC_backend_pmcmc <- function(backend, generated, cores) {
    args_with_data <- backend$args
    args_with_data$data_raw <- generated

    fit <- do.call(sifter::run_pmcmc, args_with_data)
    print('what fit looks like:')
    print(fit)
    print(str(fit))
    structure(fit,class = 'SBC_backend_pmcmc')
    # do.call(sifter::run_pmcmc, args_with_data)
  }

  SBC_fit_to_draws_matrix.SBC_backend_pmcmc <- function(fit) {
    posterior::as_draws_matrix(fit$mcmc)
  }

  backend_pmcmc <- SBC_backend_pmcmc(init_EIR = 100,
                                     n_particles=10,
                                     proposal_matrix = diag(0.5,2),
                                     max_param=125,
                                     prop_treated = 0.4,
                                     n_steps = 10,
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
  SBC_backend_iid_draws.SBC_backend_pmcmc <- function(backend) {
    FALSE
  }
  res_pmcmc <- SBC::compute_SBC(sim_data,
                           backend_pmcmc,
                           thin_ranks = 1,
                           # cache_mode = "results",
                           # cache_location = file.path(cache_dir,"test"),
                           cache_mode = "none",
                           globals = c('SBC_fit_to_draws_matrix.SBC_backend_pmcmc',
                                       'SBC_backend_iid_draws.SBC_backend_pmcmc',
                                       'SBC_fit.SBC_backend_pmcmc')
  )
  return(res_pmcmc)
}
