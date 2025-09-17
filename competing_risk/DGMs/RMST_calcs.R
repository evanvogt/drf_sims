################################################################################
# naiive calculation not accounting for competing risks

# work through an example first

# discharge parameters
b0 <- 4 # baseline discharge time
bW <- -1 # treatment shortens ICU stay
b1 <- -0.05 # prognostic - female
b2 <- 0.7 # prognostic - APACHE ish

# mortality parameters
c0 <- 8 # baseline mortality time
cW <- 1 # treatment prolongs mortality
c1 <- 0.05 # prognostic - female
c2 <- -1 # prognostic

# use a weibull distribution to simulate the event times


# found the SSRMST package - simulation based

install.packages("SSRMST")
library(SSRMST)

ssrmst(ac_rate = NULL, ac_period = 30, ac_number = 100, tot_time = 90, tau = 10, shape0 = 1.5, scale0 = 7.7, shape1 = 1.5, scale1 = 9.7, seed = 1998)
