orderly2::orderly_dependency("process_nnp_data", quote(latest()),
                             c('NG_CS_all_grouped_site.rds','NG_ANC_mother_grouped_sitegrav.rds','NG_ANC_mother_grouped_sitegrav_sg.rds',
                               'MZ_CS_all_grouped_site_180123.rds','MZ_ANC_mother_grouped_sitegrav_0822.rds','MZ_ANC_mother_grouped_sitegrav_0822_sg.rds'))
orderly2::orderly_dependency("process_mipmon_data", quote(latest()),
                             c('mipmon4corr.RDS','mipmon4corr_sg.RDS'))
orderly2::orderly_dependency("process_vaneijk_data", quote(latest()),
                             c('vaneijk_raw.RDS'))

orderly2::orderly_artefact('all_data_pgsgmg4model.rds',description='Processed data for correlation')

###NNP###
NG_cs <- readRDS('NG_CS_all_grouped_site.rds')
NG_anc <- readRDS('NG_ANC_mother_grouped_sitegrav.rds')
MZ_cs <- readRDS('MZ_CS_all_grouped_site_180123.rds')
MZ_anc <- readRDS('MZ_ANC_mother_grouped_sitegrav_0822.rds')

MZ_anc_sg <- readRDS('MZ_ANC_mother_grouped_sitegrav_0822_sg.rds')
NG_anc_sg <- readRDS('NG_ANC_mother_grouped_sitegrav_sg.rds')

NG_cs$country <- 'Nigeria'
NG_anc$country <- 'Nigeria'
MZ_cs$country <- 'Mozambique'
MZ_anc$country <- 'Mozambique'

NG_anc_sg$country <- 'Nigeria'
MZ_anc_sg$country <- 'Mozambique'

all_cs <- rbind(NG_cs,MZ_cs)
all_anc <- rbind(NG_anc,MZ_anc)
all_cs$month_adj <- all_cs$month
all_anc_sg <- rbind(NG_anc_sg,MZ_anc_sg)
##Create an adjusted month variable so that the first month of CX data matches the first month of ANC data
all_cs[all_cs$country=='Mozambique'&(all_cs$month=='Oct 2020'),]$month_adj <- as.yearmon('Dec 2020')
all_cs[all_cs$country=='Nigeria'&(all_cs$month=='Oct 2020'),]$month_adj <- as.yearmon('Nov 2020')
all_cs <- all_cs[,-2] %>%
  dplyr::rename(month=month_adj)

all_both <- merge(all_cs,all_anc, by = c('site','month','country'), suffixes = c('.cs','.anc'))
all_both_total <- all_both %>%
  group_by(country,site,month,.drop=FALSE)%>%
  dplyr::summarise(positive.anc=sum(positive.anc),total.anc=sum(total.anc),positive.cs=mean(positive.cs),total.cs=mean(total.cs),
                   mean.cs=mean(mean.cs),upper.cs=mean(upper.cs),lower.cs=mean(lower.cs),
                   mean.anc=mean(mean.anc),upper.anc=mean(upper.anc),lower.anc=mean(lower.anc))%>%
  mutate(grav_cat = 'All pregnancies')
all_both <- rbind(all_both,all_both_total)

all_both_sg <- merge(all_cs,all_anc_sg, by = c('site','month','country'), suffixes = c('.cs','.anc'))
all_both_total_sg <- all_both_sg %>%
  group_by(country,site,month,.drop=FALSE)%>%
  dplyr::summarise(positive.anc=sum(positive.anc),total.anc=sum(total.anc),positive.cs=mean(positive.cs),total.cs=mean(total.cs),
                   mean.cs=mean(mean.cs),upper.cs=mean(upper.cs),lower.cs=mean(lower.cs),
                   mean.anc=mean(mean.anc),upper.anc=mean(upper.anc),lower.anc=mean(lower.anc))%>%
  mutate(grav_cat = 'All pregnancies')
all_both_sg <- rbind(all_both_sg,all_both_total_sg)
all_nnp_total <- all_both_total %>%
  dplyr::mutate(country=dplyr::recode(country,`Nigeria`='NNP - Nigeria',
                                      `Mozambique`='NNP - Mozambique'))
all_nnp_pg <- all_both[all_both$grav_cat=='Gravidities 1',]%>%
  dplyr::mutate(country=dplyr::recode(country,`Nigeria`='NNP - Nigeria',
                                      `Mozambique`='NNP - Mozambique'))
all_nnp_mg <- all_both[all_both$grav_cat %in% c('Gravidities 2-3','Gravidities 4+'),]%>%
  group_by(country,site,month)%>%
  dplyr::summarise(positive.cs=mean(positive.cs),
                   total.cs=mean(total.cs),
                   positive.anc=sum(positive.anc),
                   total.anc=sum(total.anc))%>%
  dplyr::mutate(country=dplyr::recode(country,`Nigeria`='NNP - Nigeria',
                                      `Mozambique`='NNP - Mozambique'))
all_nnp_mg <- addCIs(all_nnp_mg,all_nnp_mg$positive.cs,all_nnp_mg$total.cs)%>%
  dplyr::rename(mean.cs=mean,
                upper.cs=upper,
                lower.cs=lower)
all_nnp_mg <- addCIs(all_nnp_mg,all_nnp_mg$positive.anc,all_nnp_mg$total.anc)%>%
  dplyr::rename(mean.anc=mean,
                upper.anc=upper,
                lower.anc=lower)

all_nnp_sg <- all_both_sg[all_both_sg$grav_cat=='Gravidities 2',]%>%
  dplyr::mutate(country=dplyr::recode(country,`Nigeria`='NNP - Nigeria',
                                      `Mozambique`='NNP - Mozambique'))


###Van Eijk###
vaneijk <- readRDS('vaneijk_raw.RDS')
primi_ve <- vaneijk[vaneijk$subG == 'P',3:6]%>%
  dplyr::rename(positive.cs=child_Y,
                total.cs=child_N,
                positive.anc=preg_Y,
                total.anc=preg_N)
multi_ve <- vaneijk[vaneijk$subG == 'M',3:6]%>%
  dplyr::rename(positive.cs=child_Y,
                total.cs=child_N,
                positive.anc=preg_Y,
                total.anc=preg_N)
all_ve <- vaneijk %>%
  subset(subG %in% c('P','M')) %>%
  group_by(Site)%>%
  dplyr::summarise(inf_prev_n=mean(child_Y),
                   inf_N=mean(child_N),
                   ANC_prev_n=sum(preg_Y),
                   ANC_N=sum(preg_N))%>%
  dplyr::select(-Site) %>%
  dplyr::rename(positive.cs=inf_prev_n,
                total.cs=inf_N,
                positive.anc=ANC_prev_n,
                total.anc=ANC_N)
primi_ve <- addCIs_anc(primi_ve,primi_ve$positive.cs,primi_ve$total.cs,primi_ve$positive.anc,primi_ve$total.anc)%>%
  mutate(country='Van Eijk')
multi_ve <- addCIs_anc(multi_ve,multi_ve$positive.cs,multi_ve$total.cs,multi_ve$positive.anc,multi_ve$total.anc)%>%
  mutate(country='Van Eijk')
all_ve <- addCIs_anc(all_ve,all_ve$positive.cs,all_ve$total.cs,all_ve$positive.anc,all_ve$total.anc)%>%
  mutate(country='Van Eijk')


###MiPMon###
mipmon <- readRDS('mipmon4corr.RDS')
all_mipmon_total <- mipmon[mipmon$grav=='All',]%>%
  group_by(site,month)%>%
  dplyr::summarise(positive.cs=mean(positive.cs),
                   total.cs=mean(total.cs),
                   positive.anc=sum(positive.anc),
                   total.anc=sum(total.anc))
all_mipmon_pg <- mipmon[mipmon$grav=='primi',]
all_mipmon_mg <- mipmon[mipmon$grav=='multi',]
mipmon_sg <- readRDS('mipmon4corr_sg.RDS')
all_mipmon_sg <- mipmon_sg[mipmon_sg$grav_sg=='secundi',]

all_mipmon_total <- addCIs_anc(all_mipmon_total,all_mipmon_total$positive.cs,all_mipmon_total$total.cs,all_mipmon_total$positive.anc,all_mipmon_total$total.anc)%>%
  mutate(country='MiPMon')
all_mipmon_pg <- addCIs_anc(all_mipmon_pg,all_mipmon_pg$positive.cs,all_mipmon_pg$total.cs,all_mipmon_pg$positive.anc,all_mipmon_pg$total.anc)%>%
  mutate(country='MiPMon')
all_mipmon_mg <- addCIs_anc(all_mipmon_mg,all_mipmon_mg$positive.cs,all_mipmon_mg$total.cs,all_mipmon_mg$positive.anc,all_mipmon_mg$total.anc)%>%
  mutate(country='MiPMon')
all_mipmon_sg <- addCIs_anc(all_mipmon_sg,all_mipmon_sg$positive.cs,all_mipmon_sg$total.cs,all_mipmon_sg$positive.anc,all_mipmon_sg$total.anc)%>%
  mutate(country='MiPMon')

###Combine all data
all_data_total <- plyr::rbind.fill(all_nnp_total,all_ve,all_mipmon_total)
all_data_pg <- plyr::rbind.fill(all_nnp_pg,primi_ve,all_mipmon_pg)
all_data_mg <- plyr::rbind.fill(all_nnp_mg,multi_ve,all_mipmon_mg)
all_data_sg <- plyr::rbind.fill(all_nnp_sg,all_mipmon_sg)

all_data_total4model <- all_data_total %>%
  filter(!(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month >= as.yearmon('Nov 2021')))
all_data_pg4model <- all_data_pg %>%
  filter(!(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month >= as.yearmon('Nov 2021')))
all_data_mg4model <- all_data_mg %>%
  filter(!(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month >= as.yearmon('Nov 2021')))
all_data_pgmg4model <- merge(all_data_pg4model,all_data_mg4model, by = c('site','month','country','positive.cs','total.cs','mean.cs','upper.cs','lower.cs'), suffixes = c('.pg','.mg'))

all_data_sg4model <- all_data_sg %>%
  filter(!(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month >= as.yearmon('Nov 2021')))

all_data_pgsgmg4model <- dplyr::left_join(all_data_pgmg4model,all_data_sg4model, by = c('site','month','country','positive.cs','total.cs','mean.cs','upper.cs','lower.cs'), suffix = c('.pgmg','.sg'))%>%
  dplyr::rename(positive.anc.sg = positive.anc,
                total.anc.sg = total.anc,
                mean.anc.sg = mean.anc,
                upper.anc.sg = upper.anc,
                lower.anc.sg = lower.anc)

saveRDS(all_data_pgsgmg4model,'all_data_pgsgmg4model.rds')

