orderly2::orderly_shared_resource('get_seasonal_profile.R',
                                  'get_seasonal_sim.R',
                                  'odin_model_stripped_seasonal.R')

source('get_seasonal_profile.R')
source('get_seasonal_sim.R')

admin_units_seasonal <- mamasante::load_file("admin_units_seasonal.rds")

admin_list <- list('Tanga','Upper East','Fatick','Equateur')
country_list <- list('Tanzania','Ghana','Senegal','Democratic Republic of Congo')
init_EIR_list <- c(10,50,100,500)
model_file <- "odin_model_stripped_seasonal.R"

seasonal_profiles <- bind_rows(lapply(1:4,function(x){
  get_seasonal_profile(admin=admin_list[[x]],country=country_list[[x]])
}))

combi <- expand.grid(x=1:4,y=1:4)

sim_seasonal_data_list <- apply(combi,1,function(i){
  print(i['x'])
    sim <- gen_seasonal_sim(init_EIR=init_EIR_list[i['y']],
                                  max_param=125,
                                  model_file= model_file,
                                  country = country_list[[i['x']]],
                                  admin_unit = admin_list[[i['x']]],
                                  sim_length = 8)
    sim$true_data$country <- country_list[[i['x']]]
    sim$true_data$admin <- admin_list[[i['x']]]
    sim$true_data$site <- paste0(admin_list[[i['x']]],', ',country_list[[i['x']]])
    sim$data_raw$country <- country_list[[i['x']]]
    sim$data_raw$admin <- admin_list[[i['x']]]
    sim$data_raw$site <- paste0(admin_list[[i['x']]],', ',country_list[[i['x']]])
    saveRDS(sim$true_data,paste0('sim_data_EIR',init_EIR_list[i['y']],'_',admin_list[[i['x']]],'_',country_list[[i['x']]],'.RDS'))
    return(sim$data_raw)
})
saveRDS(sim_seasonal_data_list,'sim_seasonal_dataraw_list.RDS')
