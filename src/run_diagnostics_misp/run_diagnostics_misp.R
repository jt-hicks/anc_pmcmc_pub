orderly2::orderly_shared_resource('create_diag_figs.R')
orderly2::orderly_parameters(proposal_matrix=1,
                             n_datasets = 32,
                             length=NULL,
                             workers=1,
                             chain=1,
                             seed=1L,
                             start_pf_time=360)

source('create_diag_figs.R')

results_list <- list()
prop_df <- data.frame(name=integer(),
                      prop=numeric())

names_list <- c(1:n_datasets)

if(n_datasets==32){
  misspecification <- c(-0.2,0.2)
  admin_list <- c('Tanga','Upper East','Fatick','Equateur')
  country_list <- c('Tanzania','Ghana','Senegal','Democratic Republic of Congo')
  init_EIR_list <- c(50)
  start_pf_time_list <- c(30,90,180,360)
  names_list <- c(5:8)

  names_key <- expand.grid(misspecification=misspecification,x=c(1:4),start_pf_time=start_pf_time_list)
  names_key$name <- names_list[names_key$x]
  names_key$admin <- admin_list[names_key$x]
  names_key$country <- country_list[names_key$x]
  names_key$init_EIR <- 50
}else if(n_datasets==16){
  misspecification <- c(0)
  admin_list <- c('Tanga','Upper East','Fatick','Equateur')
  country_list <- c('Tanzania','Ghana','Senegal','Democratic Republic of Congo')
  init_EIR_list <- c(50)
  start_pf_time_list <- c(30,90,180,360)
  names_list <- c(5:8)

  names_key <- expand.grid(misspecification=misspecification,x=c(1:4),start_pf_time=start_pf_time_list)
  names_key$name <- names_list[names_key$x]
  names_key$admin <- admin_list[names_key$x]
  names_key$country <- country_list[names_key$x]
  names_key$init_EIR <- 50

}


# for(name in names_list){
for(x in c(1:n_datasets)){
  name <- names_key[x,'name']
  misspecification <- names_key[x,'misspecification']
  start_pf_time <- names_key[x,'start_pf_time']
  orderly2::orderly_dependency(name="run_pmcmc", query=paste0('latest(parameter:name == environment:name && parameter:misspecification == environment:misspecification && parameter:start_pf_time == environment:start_pf_time)'),
                               c("data/${x}.rds" = "result.rds"))
  result <- readRDS(paste0('data/',x,".rds"))
  result <- append(result,list(sim_pars=names_key[x,c('admin','country','init_EIR','misspecification','start_pf_time')]))

  temp_df <- data.frame(name=x,
                        prop=as.numeric(var(result$pars)))
  prop_df <- dplyr::bind_rows(prop_df,temp_df)

  temp_list <- list(result)
  names(temp_list) <- name
  results_list <- append(results_list,temp_list)

  plot <- create_diag_figs(result,
                   country = paste(names_key[x,'country'],names_key[x,'misspecification']),
                   district = names_key[x,'start_pf_time'],
                   name = x)
  ggsave(filename=paste0('plots/',x,'.png'),plot,width=5,height=7,units='in')
}

saveRDS(prop_df,'props.RDS')
saveRDS(results_list,'results_list.RDS')
