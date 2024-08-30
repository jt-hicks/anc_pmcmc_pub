data_gen_piecewise <- function(init_EIR=100,
                               betaa_sample,
                               prop_treated = 0.4,
                               init_age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80),
                               het_brackets = 5,
                               plot_instance = TRUE){

  # model_file <- system.file("odin", "MiP_odin_model_nodelay.R", package = "mamasante")
  model_file <- 'Q:/mamasante/inst/odin/MiP_odin_model_nodelay.R'




  ##Convert init_EIR into an init_betaa
  mpl_initial <- mamasante::model_param_list_create(init_EIR = init_EIR,
                                                    init_ft = prop_treated,
                                                    comparison='u5'
  )

  pars_initial <- mamasante::equilibrium_init_create_stripped(age_vector = init_age,
                                                              init_EIR = init_EIR,
                                                              ft = prop_treated,
                                                              model_param_list = mpl_initial,
                                                              het_brackets = het_brackets)

  ################## generate the data ######################
  ##set up the simulation for the simualted data
  time<- max(betaa_sample$t) + 30
  out_step=1

  mpl <- mamasante::model_param_list_create(init_EIR = pars$init_EIR,
                                            init_ft = prop_treated,
                                            betaa_times=c(0,betaa_sample$t),
                                            betaa_vals=c(pars$betaa_eq, unlist(betaa_sample[,2])),
                                            lag_rates = 10,
                                            comparison = 'u5'
  )

  pars <- mamasante::equilibrium_init_create_stripped(age_vector = init_age,
                                                      init_EIR = init_EIR,
                                                      ft = prop_treated,
                                                      model_param_list = mpl,
                                                      het_brackets = het_brackets)
  ##The malaria model
  generator <- odin(model_file)
  state_use <- pars[names(pars) %in% coef(generator)$name]

  # create model with initial values
  mod <- generator(user = state_use, use_dde = TRUE)
  tt <- seq(0, time, out_step)

  # run the simulation to base the data
  mod_run <- mod$run(tt, step_max_n = 1e7,
                     atol = 1e-5,
                     rtol = 1e-5)

  # shape output
  out <- mod$transform_variables(mod_run)

  # plot data and generate data
  # plot(out$t,out$betaa_out,type='l',col="red",ylim=c(0,125))
  # lines(out$t,out$prev*100,col="blue",lwd=4)
  out_df <- data.frame(t=out$t,
                       prev_05=out$prev,
                       EIR=out$EIR_out,
                       betaa=out$betaa_out,
                       inc05=out$inc05,
                       inc_all = out$inc,
                       prev_all = out$prev_all,
                       mv = out$Sv+out$Ev+out$Iv)
  return(out_df)
}
