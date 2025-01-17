library(R2OpenBUGS)
library(coda)
library(car)
library(MASS)
library(plotrix)
library(binom)
library(scales)
library(zoo)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(readxl)
source('utils.R')
library(patchwork)
library(bayesplot)
library(tidyverse)
library(excel.link)
library(zoo)
library(gridExtra)
library(hipercow)

##Copy over data needed for correlation
dir.create('src/process_nnp_data/data_raw')
dir.create('src/process_nnp_data/data_out')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/ANC-based surveillance/data_anc_mother_nigeria.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_child_nigeria_2020.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_hh_nigeria_2020.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_child_nigeria_2021.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_hh_nigeria_2021.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_child_nigeria_2022.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_hh_nigeria_2022.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Mozambique/ANC-based surveillance/ANC_Mozambique_2022.09_FINAL.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Mozambique/Cross-sectional survey/Moz Baseline CSS Datasets 2020.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Mozambique/Cross-sectional survey/Moz Midline CSS Datasets 2021.xlsx',
          to='src/process_nnp_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Mozambique/Cross-sectional survey/Moz Endline CSS Datasets 2022.xlsx',
          to='src/process_nnp_data/data_raw/')

dir.create('src/process_mipmon_data/data_raw')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Mozambique/isglobal_cism_data/mipmon_merged.csv',
          to='src/process_mipmon_data/data_raw/')
file.copy(from='C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Mozambique/isglobal_cism_data/cross_merged.csv',
          to='src/process_mipmon_data/data_raw/')

dir.create('src/process_vaneijk_data/data_raw')
file.copy(from='Q:/anc_pmcmc/nnp/data/paper_data_strat_G.xlsx',
          to='src/process_vaneijk_data/data_raw/')


##Process NNP data
nnp_data_id <- orderly2::orderly_run('process_nnp_data')

##Process mipmon data
mipmon_data_id <- orderly2::orderly_run('process_mipmon_data')

##Process Van Eijk data
vaneijk_data_id <- orderly2::orderly_run('process_vaneijk_data')

##Combine and process all data
combine_data_id <- orderly2::orderly_run('process_data_for_correlation')

##Run correlation model
run_corr_id <- orderly2::orderly_run('run_anc_corr')

##Summarise model output and produce plot
create_plot_id <- orderly2::orderly_run('create_anc_corr_plot')
