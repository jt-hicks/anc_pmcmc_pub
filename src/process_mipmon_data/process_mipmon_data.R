orderly2::orderly_resource("data_raw/")
orderly2::orderly_artefact(c('mipmon4corr.RDS','mipmon4corr_sg.RDS'),
                           description='MiPMon processed data')

mipmon <- read.csv("data_raw/mipmon_merged.csv")

mipmon.all <- mipmon %>%
  mutate(date = as.Date(visdate,tryFormats='%d/%m/%Y'),
         year = format(date, format='%Y'),
         week = week(date),
         month = as.yearmon(date),
         N=1,
         rdt=ifelse(pcrpos=='PCR-',
                    0,
                    ifelse(pcrpos=='PCR+',
                           ifelse(density<100&!is.na(density),
                                  0,
                                  ifelse(!is.na(density),
                                         1,
                                         NA)
                           ),
                           NA)
         ),
         site=recode(posto_code,
                     `Ilha-Josina`='Ilha Josina',
                     `Magude-sede`='Magude Sede',
                     `MOtzae`='Magude Sede',
                     `Manhica-Sede`='Manhica',
                     .default = NA_character_),
         grav=ifelse(gestnum==1,'primi',ifelse(!is.na(gestnum),'multi',NA)),
         grav_sg = ifelse(grav=='primi','primi',ifelse(gestnum==2,'secundi',ifelse(gestnum > 2,'multi',NA))))%>%
  filter(visit=='PN'&!is.na(site)&!is.na(gestnum)&!is.na(month)&!is.na(rdt))

mipmon.grouped <- mipmon.all%>%
  group_by(site, month, grav) %>%
  summarise(positive.anc=sum(rdt),
            total.anc=sum(N))
mipmon.grouped_sg <- mipmon.all%>%
  group_by(site, month, grav_sg) %>%
  summarise(positive.anc=sum(rdt),
            total.anc=sum(N))

mipmon.grouped.all <- mipmon.all%>%
  group_by(site, month) %>%
  summarise(positive.anc=sum(rdt),
            total.anc=sum(N))

mipmon.grouped.all2 <- rbind(mipmon.grouped,mipmon.grouped.all)%>%
  mutate(grav=ifelse(is.na(grav),'All',grav))

mipmon.grouped.all2_sg <- rbind(mipmon.grouped_sg,mipmon.grouped.all)%>%
  mutate(grav_sg=ifelse(is.na(grav_sg),'All',grav_sg))

#read in cross-sectional data
cross <- read.csv('data_raw/cross_merged.csv')%>%
  mutate(date = as.Date(visdate),
         year = format(date, format='%Y'),
         week = week(date),
         month = as.yearmon(date),
         N=1,
         rdt=ifelse(rdt=='Negative',0,ifelse(rdt=='Positive',1,NA)),
         site=recode(area,`Ilha Josina`='Ilha Josina',
                     `Magude Sede`='Magude Sede',
                     `Motaze`='Magude Sede',
                     `ManhiÃ§a`='Manhica',
                     .default = NA_character_))%>%
  filter(!is.na(site)&!is.na(month)&!is.na(rdt))

cross.grouped <- cross%>%
  group_by(site, month) %>%
  summarise(positive.cs=sum(rdt),
            total.cs=sum(N))

mipmon.combined <- merge(cross.grouped,mipmon.grouped.all2, by = c('site','month'), suffixes = c('.cs','.anc'))
saveRDS(mipmon.combined,'mipmon4corr.RDS')
mipmon.combined_sg <- merge(cross.grouped,mipmon.grouped.all2_sg, by = c('site','month'), suffixes = c('.cs','.anc'))
saveRDS(mipmon.combined_sg,'mipmon4corr_sg.RDS')

