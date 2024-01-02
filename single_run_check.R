library(sifter)
library(odin)
library(zoo)
library(tidyr)
library(dplyr)
source('utils.R')

sim_data <- datasim4pmcmc(volatility=1,init_EIR=100, max_param = 125)
