orderly2::orderly_resource("data_raw/")
orderly2::orderly_artefact(c('NG_CS_all_grouped_site.rds','NG_ANC_mother_grouped_sitegrav.rds','NG_ANC_mother_grouped_sitegrav_sg.rds',
                             'MZ_CS_all_grouped_site_180123.rds','MZ_ANC_mother_grouped_sitegrav_0822.rds','MZ_ANC_mother_grouped_sitegrav_0822_sg.rds'),
                           description='NNP processed data')

#######################
########NIGERIA########
#######################
#Read in Nigeria data and rename most used variables
NG_ANC_mother <- xl.read.file('data_raw/data_anc_mother_nigeria.xlsx')%>%
  dplyr::rename(primigrav = q1_preg,
                prev_pregs = q2_preg,
                mal_symp = q1_mal,
                rdt = q2_mal) %>%
  dplyr::mutate(month = as.yearmon(month),
                grav = prev_pregs+1,
                ward = toupper(ward))

#Read in Nigeria data and rename most used variables
NG_CS_child_2020 <- xl.read.file('data_raw/data_nnp_survey_child_nigeria_2020.xlsx')%>%
  dplyr::rename(rdt = q85c_result,
                mal_symp = q86a_symptoms) %>%
  mutate(age = ifelse(age==0,age_months/12,age))

NG_CS_hh_2020 <- xl.read.file('data_raw/data_nnp_survey_hh_nigeria_2020.xlsx')%>%
  dplyr::rename(individual_id = hh_id,
                date = start) %>%
  mutate(date = as.Date(date))

#Merge child to hh to get date of survey
NG_CS_2020 <- merge(NG_CS_child_2020,NG_CS_hh_2020,by='submission_id',all.x=TRUE,suffixes=c('.child','.hh'))

##Read in 2021 data
NG_CS_child_2021 <- xl.read.file('data_raw/data_nnp_survey_child_nigeria_2021.xlsx')%>%
  dplyr::rename(rdt = q90c_result_6to59months,
                mal_symp = q91_a_symptoms) %>%
  mutate(age = ifelse(age==0,age_months/12,age))

NG_CS_hh_2021 <- xl.read.file('data_raw/data_nnp_survey_hh_nigeria_2021.xlsx')%>%
  dplyr::rename(individual_id = hh_id,
                date = start) %>%
  mutate(date = as.Date(date))

#Merge child to hh to get date of survey
NG_CS_2021 <- merge(NG_CS_child_2021,NG_CS_hh_2021,by='submission_id',all.x=TRUE,suffixes=c('.child','.hh'))

##Read in 2022 data
NG_CS_child_2022 <- xl.read.file('data_raw/data_nnp_survey_child_nigeria_2022.xlsx')%>%
  dplyr::rename(rdt = q90c_result_6to59months,
                mal_symp = q91a_symptoms) %>%
  mutate(age = ifelse(age==0,age_months/12,age))

NG_CS_hh_2022 <- xl.read.file('data_raw/data_nnp_survey_hh_nigeria_2022.xlsx')%>%
  dplyr::rename(individual_id = hh_id,
                date = start) %>%
  mutate(date = as.Date(date))

#Merge child to hh to get date of survey
NG_CS_2022 <- merge(NG_CS_child_2022,NG_CS_hh_2022,by='submission_id',all.x=TRUE,suffixes=c('.child','.hh'))

#Combine 3 years of data
NG_CS_all <- plyr::rbind.fill(NG_CS_2020,NG_CS_2021,NG_CS_2022)
NG_CS_all$district_code <- as.integer(substr(NG_CS_all$cluster_n,1,1))
NG_CS_all$district_n <- unlist(lapply(NG_CS_all$district_code, FUN = function(x) switch(x,'Ejigbo','Ife North','Asa','Moro')))

#By site and gravidity#
NG_CS_all_grouped_site <- NG_CS_all %>%
  dplyr::rename(site=district_n)%>%
  mutate(month=as.yearmon(date),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA)))
  )%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(site))%>%
  group_by(site,month,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
NG_CS_all_grouped_site <- addCIs(NG_CS_all_grouped_site,NG_CS_all_grouped_site$positive,NG_CS_all_grouped_site$total)
saveRDS(NG_CS_all_grouped_site,'NG_CS_all_grouped_site.rds')

#Group by grav
NG_ANC_mother_grouped_sitegrav <- NG_ANC_mother %>%
  dplyr::rename(site = lga) %>%
  dplyr::mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")))%>%
  filter(!is.na(site))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
NG_ANC_mother_grouped_sitegrav <- addCIs(NG_ANC_mother_grouped_sitegrav,NG_ANC_mother_grouped_sitegrav$positive,NG_ANC_mother_grouped_sitegrav$total)
saveRDS(NG_ANC_mother_grouped_sitegrav,'NG_ANC_mother_grouped_sitegrav.rds')

##secundigrav
NG_ANC_mother_grouped_sitegrav_sg <- NG_ANC_mother %>%
  dplyr::rename(site = lga) %>%
  dplyr::mutate(grav_cat=cut(grav,breaks=c(0,1,2,Inf),labels=c("Gravidities 1","Gravidities 2","Gravidities 3+")))%>%
  filter(!is.na(site))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
NG_ANC_mother_grouped_sitegrav_sg <- addCIs(NG_ANC_mother_grouped_sitegrav_sg,NG_ANC_mother_grouped_sitegrav_sg$positive,NG_ANC_mother_grouped_sitegrav_sg$total)
saveRDS(NG_ANC_mother_grouped_sitegrav_sg,'NG_ANC_mother_grouped_sitegrav_sg.rds')


#######################
######MOZAMBIQUE#######
#######################
#Read in Mozambique data and rename most used variables
MZ_ANC_mother <- xl.read.file('data_raw/ANC_Mozambique_2022.09_FINAL.xlsx',
                              xl.sheet = 'ANC')%>%
  dplyr::rename(primigrav = q4_primagravidae,
                prev_pregs = q5_num_pregnancy,
                mal_symp = q6_mal_symptoms,
                rdt = q7_rdt_result,
                age = q2_age_n) %>%
  dplyr::mutate(month = as.yearmon(as.Date(date_interview_n)),
                grav = ifelse(primigrav=='Yes',1,ifelse(prev_pregs==0,NA,prev_pregs+1)))

#Read in Mozambique data and rename most used variables
MZ_CS_rdt_base <- xl.read.file('data_raw/Moz Baseline CSS Datasets 2020.xlsx',
                               xl.sheet = 'RDT')%>%
  mutate(date = as.Date(date))
MZ_CS_rdt_mid <- xl.read.file('data_raw/Moz Midline CSS Datasets 2021.xlsx',
                              xl.sheet = 'RDT')%>%
  mutate(date = as.Date(date,"%m/%d/%Y"))
MZ_CS_rdt_end <- xl.read.file('data_raw/Moz Endline CSS Datasets 2022.xlsx',
                              xl.sheet = 'RDT')%>%
  mutate(date_interview = as.Date(date_interview))%>%
  dplyr::rename(date = date_interview)
MZ_CS_rdt <- rbind(MZ_CS_rdt_base,MZ_CS_rdt_mid,MZ_CS_rdt_end)

#Group CS data by site and calculate prevalence and CI by month
MZ_CS_all_grouped_site <- MZ_CS_rdt %>%
  dplyr::rename(site = district) %>%
  mutate(month=as.yearmon(date),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA)))
  )%>%
  filter(!is.na(rdt)) %>%
  filter(site %in% c('Changara','Chemba','Guro'))%>%
  mutate(month=as.yearmon(ifelse(month=='Sep 2020',as.yearmon('Oct 2020'),as.yearmon(month))))%>% ##Group Sep and Oct 2020 together for the pre-intervention CX survey
  group_by(site,month,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
MZ_CS_all_grouped_site <- addCIs(MZ_CS_all_grouped_site,MZ_CS_all_grouped_site$positive,MZ_CS_all_grouped_site$total)
saveRDS(MZ_CS_all_grouped_site,'MZ_CS_all_grouped_site_180123.rds')

#Group ANC data by grav and site
MZ_ANC_mother_grouped_sitegrav <- MZ_ANC_mother %>%
  rename(site = district_n) %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA))))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age>=12&age<50) %>%
  filter(!is.na(grav))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
MZ_ANC_mother_grouped_sitegrav <- addCIs(MZ_ANC_mother_grouped_sitegrav,MZ_ANC_mother_grouped_sitegrav$positive,MZ_ANC_mother_grouped_sitegrav$total)
saveRDS(MZ_ANC_mother_grouped_sitegrav,'MZ_ANC_mother_grouped_sitegrav_0822.rds')

### secundigrav
MZ_ANC_mother_grouped_sitegrav_sg <- MZ_ANC_mother %>%
  dplyr::rename(site = district_n) %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,2,Inf),labels=c("Gravidities 1","Gravidities 2","Gravidities 3+")),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA))))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age>=12&age<50) %>%
  filter(!is.na(grav))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
MZ_ANC_mother_grouped_sitegrav_sg <- addCIs(MZ_ANC_mother_grouped_sitegrav_sg,MZ_ANC_mother_grouped_sitegrav_sg$positive,MZ_ANC_mother_grouped_sitegrav_sg$total)
saveRDS(MZ_ANC_mother_grouped_sitegrav,'MZ_ANC_mother_grouped_sitegrav_0822_sg.rds')



