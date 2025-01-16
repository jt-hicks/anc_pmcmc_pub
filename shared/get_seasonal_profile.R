get_seasonal_profile <- function(admin,country,years=1){
  admin_matches <- mamasante::admin_match(admin_unit = admin, country = country,
                                          admin_units_seasonal = admin_units_seasonal)

  ssa0 <- admin_units_seasonal$a0[admin_matches]
  ssa1 <- admin_units_seasonal$a1[admin_matches]
  ssa2 <- admin_units_seasonal$a2[admin_matches]
  ssa3 <- admin_units_seasonal$a3[admin_matches]
  ssb1 <- admin_units_seasonal$b1[admin_matches]
  ssb2 <- admin_units_seasonal$b2[admin_matches]
  ssb3 <- admin_units_seasonal$b3[admin_matches]
  theta_c <- admin_units_seasonal$theta_c[admin_matches]

  t <- c(1:365*years)
  return(data.frame(rainfall=pmax((ssa0+ssa1*cos(2*pi*-t/365)+ssa2*cos(2*2*pi*-t/365)+ssa3*cos(3*2*pi*-t/365)+ssb1*sin(2*pi*-t/365)+ssb2*sin(2*2*pi*-t/365)+ ssb3*sin(3*2*pi*-t/365))/theta_c,0.001),
                    t=t,
                    admin=admin,
                    country=country))
}

