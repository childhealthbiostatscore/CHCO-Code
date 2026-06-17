# StO2 cortex
sto2_cor <- function(vb_cor, hct, r2star_cor) {
  gamma = 267500000
  deltaxi0 = 2.64E-07
  pi_4_3 = (4/3)*pi
  b0 = 3
  k = gamma * deltaxi0 * pi_4_3 * b0
  r2_cor = 7.35
  hct_90 = hct*(0.009)

  sto2 <- 1 - ((r2star_cor - r2_cor) / (k * vb_cor * hct_90))
  
  return(sto2)
}

# StO2 medulla
sto2_med <- function(vb_med, hct, r2star_med) {
  gamma = 267500000
  deltaxi0 = 2.64E-07
  pi_4_3 = (4/3)*pi
  b0 = 3
  k = gamma * deltaxi0 * pi_4_3 * b0
  r2_med = 6.31
  hct_90 = hct*(0.009)
  
  sto2 <- 1 - ((r2star_med - r2_med) / (k * vb_med * hct_90))
  
  return(sto2)
}
