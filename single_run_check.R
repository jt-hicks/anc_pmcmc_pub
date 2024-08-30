devtools::install_github('hyunjimoon/SBC')
devtools::install_github('jt-hicks/sifter@cleaning_1223')

library(sifter)
library(SBC)
library('tidyverse')
library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(didehpc)
library(pkgdepends)
library("coda")
library(binom)
library(ggplot2)
library(bayesplot)
library(reshape2)
library(ggpubr)
library(zoo)
library(patchwork)
library(RColorBrewer)
source('utils.R')

##
vol_test <-data.frame(volatility=rgamma(1000, shape = 3.4, rate = 3.1))
ggplot(data=vol_test,aes(x=volatility))+geom_histogram()

eir_test <-data.frame(log_init_EIR=rnorm(1000, mean = 4, sd = 3))
eir_test$init_EIR <- exp(eir_test$log_init_EIR)
ggplot(data=eir_test,aes(x=init_EIR))+geom_histogram()+scale_x_log10()
ggplot(data=eir_test,aes(x=log_init_EIR))+geom_histogram()
median(eir_test$init_EIR)
quantile(eir_test$init_EIR)

sim_data <- datasim4pmcmc(N=12, max_param = 125)

test <- run_pmcmc(data_raw=sim_data$generated,
                      init_EIR = 100,
                      n_particles=10,
                      proposal_matrix = diag(0.5,2),
                      max_param=125,
                      prop_treated = 0.4,
                      n_steps = 10,
                      n_threads = 1,
                      state_check = 0,## Run equilibrium checks
                      seasonality_on = FALSE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                      seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                      seed = 1L,
                      start_pf_time = 30,
                      particle_tune = FALSE,
                      comparison = 'u5',
                      initial = 'fitted')
posterior::as_draws_matrix(test$probs)
test$pars
