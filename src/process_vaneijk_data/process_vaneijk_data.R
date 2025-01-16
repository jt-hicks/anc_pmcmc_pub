orderly2::orderly_resource("data_raw/")
orderly2::orderly_artefact(c('vaneijk_raw.RDS'),
                           description='Van Eijk processed data')

vaneijk_data <- read_excel('data_raw/paper_data_strat_G.xlsx')

saveRDS(vaneijk_data,'vaneijk_raw.RDS')
