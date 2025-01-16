get_u5prev_fromanc <- function(avg_prev,comparison=NULL){
  ## Gravidity prevalence conversion coefficients:
  coefs_pgsgmg_df <- apply(mamasante::load_file('pgsgmg_corr_sample.RDS'),2,median)
  coefs_all_df <- apply(mamasante::load_file('all_corr_sample.RDS'),2,median)

  av_lo_child <- coefs_pgsgmg_df[['av_lo_child']]
  gradient_pg <- coefs_pgsgmg_df[['gradient_pg']]
  intercept_pg <- coefs_pgsgmg_df[['intercept_pg']]
  gradient_sg <- coefs_pgsgmg_df[['gradient_sg']]
  intercept_sg <- coefs_pgsgmg_df[['intercept_sg']]
  gradient_mg <- coefs_pgsgmg_df[['gradient_mg']]
  intercept_mg <- coefs_pgsgmg_df[['intercept_mg']]
  av_lo_child_all <- coefs_all_df[['av_lo_child']]
  gradient_all <- coefs_all_df[['gradient']]
  intercept_all <- coefs_all_df[['intercept']]

  #Determine average log-odds of childhood prevalence depending on first year of available data
  if(comparison=='ancall'){
    log_odds_pall <- log(mamasante::get_odds_from_prev(avg_prev))
    log_odds_child <- ((log_odds_pall - intercept_all) + av_lo_child_all*gradient_all)/(gradient_all + 1)
    prev_u5 <- mamasante::get_prev_from_log_odds(log_odds_child)
  }else {
    log_odds_pg <- log(mamasante::get_odds_from_prev(avg_prev[1]))
    log_odds_child <- ((log_odds_pg - intercept_pg) + av_lo_child*gradient_pg)/(gradient_pg + 1)
    prev_u5 <- mamasante::get_prev_from_log_odds(log_odds_child)
  }

  return(prev_u5)
}
