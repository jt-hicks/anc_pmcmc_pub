orderly2::orderly_shared_resource('create_diag_figs.R')
orderly2::orderly_parameters(proposal_matrix=1,
                             n_datasets = 16,
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

admin_list <- c('Tanga','Upper East','Fatick','Equateur')
country_list <- c('Tanzania','Ghana','Senegal','Democratic Republic of Congo')
init_EIR_list <- c(10,50,100,500)

names_key <- expand.grid(x=1:4,y=1:4)
names_key$name <- c(1:16)
names_key$admin <- admin_list[names_key$x]
names_key$country <- country_list[names_key$x]
names_key$init_EIR <- init_EIR_list[names_key$y]

# for(name in names_list){
for(x in c(1:16)){
  # print(x)
  # orderly2::orderly_search(name="run_pmcmc", expr='latest(parameter:name == environment:x)')
  # que <- orderly_query(name="run_pmcmc", expr='latest(parameter:name == environment:x)')
  orderly2::orderly_dependency(name="run_pmcmc", query=paste0('latest(parameter:name == ',x,')'),
                               c("data/${x}.rds" = "result.rds"))

  result <- readRDS(paste0('data/',x,".rds"))
  result <- append(result,list(sim_pars=names_key[names_key$name==x,c('admin','country','init_EIR')]))

  temp_df <- data.frame(name=x,
                        prop=as.numeric(var(result$pars)))
  prop_df <- dplyr::bind_rows(prop_df,temp_df)

  temp_list <- list(result)
  names(temp_list) <- x
  results_list <- append(results_list,temp_list)

  create_diag_figs(result,
                   country = paste(names_key[x,'country'],names_key[x,'init_EIR']),
                   district = names_key[x,'admin'],
                   name = x)
}

saveRDS(prop_df,'props.RDS')
saveRDS(results_list,'results_list.RDS')
