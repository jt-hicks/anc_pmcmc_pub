library(zoo)
library(lubridate)
library(dplyr)

#Markhan Seasonality Index
month_arc <- data.frame(date = seq.Date(from=as.Date('2010-01-01'),to=as.Date('2010-12-31'),by='days'),
                        degree = 0:364)%>%
  mutate(radians = degree*pi/180,
         mid_month = as.character(as.Date(as.yearmon(date),frac=0.5)),
         month=month(date))%>%
  filter(date==as.Date(mid_month))%>%
  select(month,degree,radians)


siaya_total <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/Writing/Manuscripts/pMCMC/Real_comparisons/white_paper_siaya.txt')
siaya_total$date <- as.Date(date_decimal(siaya_total$date_decimal))
siaya_total$year <- year(siaya_total$date)
siaya_total$month <- month(siaya_total$date)

return_msi <- function(dataframe,single_year=FALSE){
  dataframe$year <- year(dataframe$date)
  dataframe$month <- month(dataframe$date)
  if(single_year){
    if(nrow(dataframe)!=12){
      stop('There need to be 12 months of data to calculate MSI for a single year.')
    }
    dataframe$year = 1
  }
  output <- dataframe %>%
    left_join(month_arc, by=join_by(month))%>%
    mutate(sin = value*sin(radians),
           cos = value*cos(radians))%>%
    group_by(year)%>%
    summarise(rk = sqrt(sum(sin)^2 + sum(cos)^2),
              theta_k = atan2(sum(sin),sum(cos)),
              theta_k_deg = (theta_k*180/pi + 360) %% 360,
              msi = rk/sum(value),
              month_num=n()) %>%
    filter(month_num==12)

  return(output)
}
siaya_total$value <- siaya_total$cases_total
siaya_msi <- return_msi(siaya_total)
mean(siaya_msi$msi)
