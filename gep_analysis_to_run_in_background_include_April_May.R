# Run the GEP analysis 
# Include all months
# Initial time is June
 
library(tidyverse)
library(lubridate)
library(gnnlab)
library(readxl)
library(ggpubr)

library(brms)
library(tidybayes)
library(bayesplot)
library(multidplyr)

getwd()
dset = readRDS("prepared_brms_data.rds")
# dset = dset %>% mutate(PPFD = PPFD * sd_ppfd + mean_ppfd,TEMP = TEMP * sd_temp + mean_temp)

WARMUP  = 5000
SAMPLES = 2000
SEED    = 2021
CHAINS = CORES = 4

CTRL = list(adapt_delta = 0.9999999999, max_treedepth = 30)

dsetZ = dset %>% filter(str_detect(location, "Z"))
dsetS = dset %>% filter(str_detect(location, "S"))

ggplot(dsetZ) +
  geom_point(aes(x= month, y = PPFD, group = month)) +
  facet_grid(rows = vars(state))

################################################################################
# Model testig -----------------------------------------------------------------
## Check the knots.
## For time-series, set the knots outside the bounds to 
## enforce a good fit with the cyclic cubic spline.
## Also I set the knots to 4 so that the curvature in
## May and April match. 

# outT = mgcv::gam(TEMP~s(month, bs = "cc", k = 6, by = state), knots  = list(month = c(0.5, 12.5)), data = dsetZ)
# outP = mgcv::gam(PPFD~s(month, bs = "cc", k = 4, by = state, id = 0), knots  = list(month = c(0.5, 12.5)), data = dsetZ, 
#                  family = gaussian())
# 
# dsetZ %>% expand(state,  month = seq(0.5, 12.5, by = 0.5)) %>% 
#   mutate(pred = predict(outP, newdata = ., type = "response")) %>% 
#   ggplot() + 
#   geom_line(aes(x = month, y = pred, color = state)) +
#   geom_point(aes(x = month, y = PPFD, color = state), data =dset, alpha = 0.2)
# 
# 
# 
# dsetZ %>% expand(state,  month = seq(0.5, 12.5, by = 0.5)) %>% 
#   mutate(pred = predict(outT, newdata = .)) %>% 
#   ggplot() + 
#   geom_line(aes(x = month, y = pred, color = state)) +
#   geom_point(aes(x = month, y = TEMP, color = state), data =dset, alpha = 0.2)
# 

ppfdmodel = bf(PPFD  ~ s(month, k = 4, by = state, bs = "cc", id = 0) + state,
               sigma ~ s(month, k = 6, by = state, id = 0) + state ) + gaussian()

# testZ = brm(ppfdmodel, data = dsetZ, seed = SEED, chains = CHAINS, cores = CORES, knots = list(month = c(0.5, 12.5)), control = list(adapt_delta = 0.999, max_treedepth = 20))
# testS = update(testZ, newdata = dsetS, seed = SEED, chains = CHAINS, cores = CORES, knots = list(month = c(0.5, 12.5)), control = list(adapt_delta = 0.999, max_treedepth = 20))
# brms::pp_check(testZ)
# brms::pp_check(testS)

################################################################################
tempmodel = bf(TEMP  ~ s(month, k = 6, by = state, bs = "cc", id = 0) + state,
               sigma ~ s(month, k = 6, by = state, id = 0) + state) + gaussian()
gepmodel  = bf(GEP ~ PPFD * TEMP * state) + Gamma("log")
# gepmodel  = bf(GEP ~ PPFD * TEMP  + (0 + PPFD + TEMP || state) + state) + Gamma("log") 

completemodel = gepmodel + ppfdmodel + tempmodel + set_rescor(FALSE)

get_prior(completemodel, data = dset)

PRIORS = 
  set_prior("student_t(3, 3, 2)", class = "Intercept", resp = "GEP") +
  set_prior("student_t(3, 0, 1)", class = "b", coef = "PPFD", resp = "GEP") +
  set_prior("student_t(3, 0, 1)", class = "b", coef = "TEMP", resp = "GEP") +
  set_prior("student_t(3, 0, 1)", class = "b", coef = "PPFD:TEMP", resp = "GEP") +
  set_prior("student_t(3, 0, 1)", class = "b", coef = "stateVegetated", resp = "GEP") +
  set_prior("normal(0, 1)", class = "shape", resp = "GEP") +
  
  set_prior("student_t(3, 0, 2)", class = "b",         resp = "PPFD") +
  set_prior("student_t(3, 0, 1)",  class = "sds",       resp = "PPFD") +
  set_prior("normal(0, 1)",      class = "Intercept", resp = "PPFD") +
  set_prior("student_t(3, 0, 1)",  class = "b",         resp = "PPFD", dpar = "sigma") +
  set_prior("student_t(3, 0, 1)",  class = "sds",       resp = "PPFD", dpar = "sigma") +
  set_prior("normal(0, 2)",        class = "Intercept", resp = "PPFD", dpar = "sigma") +  
  
  set_prior("student_t(3, 0, 1)",  class = "b",         resp = "TEMP") +
  set_prior("student_t(3, 0, 1)",  class = "sds",       resp = "TEMP") +
  set_prior("normal(0, 2)",        class = "Intercept", resp = "TEMP") +
  set_prior("student_t(3, 0, 1)",  class = "b",         resp = "TEMP", dpar = "sigma") +
  set_prior("student_t(3, 0, 1)",  class = "sds",       resp = "TEMP", dpar = "sigma") +
  set_prior("normal(0, 2)",       class = "Intercept", resp = "TEMP", dpar = "sigma") 

z1 = lubridate::now()
boutZ = brm(completemodel, 
            data = dsetZ, 
            seed = SEED,
            chains = CHAINS, cores = CORES, iter = SAMPLES + WARMUP,
            warmup = WARMUP,
            control = CTRL,
            prior = PRIORS,
            backend = "cmdstanr",
            threads = threading(4),
            refresh = 3000,
            knots = list(month = c(0.5, 12.5)))
z2 = lubridate::now()

summary(boutZ)
boutZ %>% tidy_draws() %>% pull(treedepth__) %>% range() # Range of treedepth
boutZ %>% tidy_draws() %>% pull(divergent__) %>% sum()


s1 = lubridate::now()
boutS = update(boutZ, 
            newdata = dsetS, 
            seed = SEED,
            chains = CHAINS, cores = CORES, iter = SAMPLES + WARMUP,
            warmup = WARMUP,
            control = CTRL,
            backend = "cmdstanr",
            threads = threading(4),
            refresh = 3000,
            knots = list(month = c(0.5, 12.5)))
s2 = lubridate::now()


saveRDS(boutZ, "boutZ_include_April_May.rds")
saveRDS(boutS, "boutS_include_April_May.rds")

################################################################################

# boutZ = readRDS( "boutZ_include_April_May.rds")
# boutS = readRDS( "boutS_include_April_May.rds")

aseries = function(n, floor = TRUE) {
  # Function to calculate ISO216 A Series
  n = as.integer(n)
  wodd = function(n) {
    1 / (2 ^ ((n + 1)/2)) * 1000 * sqrt(sqrt(2))
  }
  hodd = function(n) {
    1 / (2 ^ ((n - 1)/2)) * 1000 / sqrt(sqrt(2))
  }
  weven = function(n) {
    1 / (2 ^ (n / 2)) * 1000 / sqrt(sqrt(2))
  }
  heven = function(n) {
    1 / (2 ^ (n / 2)) * 1000 * sqrt(sqrt(2))
  }
  if(n %% 2)  {
    w = wodd(n)
    h = hodd(n)
  } else {
    w = weven(n)
    h = heven(n)
  }
  if(floor) {return(floor(c(w,h)))}
  return(c(w,h))
}

NSAMPLES = 50

# Zostera
boutZ %>% tidy_draws() %>% pull(treedepth__) %>% range() # Range of treedepth
boutS %>% tidy_draws() %>% pull(treedepth__) %>% range()
boutS %>% tidy_draws() %>% pull(divergent__) %>% sum()
boutZ %>% tidy_draws() %>% pull(divergent__) %>% sum()   # Number of divergent transitions

grepz = posterior_predict(boutZ, resp = "GEP", nsamples  = NSAMPLES)
gz = boutZ$data$GEP
trepz = posterior_predict(boutZ, resp = "TEMP", nsamples  = NSAMPLES)
tz = boutZ$data$TEMP
prepz = posterior_predict(boutZ, resp = "PPFD", nsamples  = NSAMPLES)
pz = boutZ$data$PPFD

# Sargassum

greps = posterior_predict(boutS, resp = "GEP", nsamples  = NSAMPLES)
gs = boutS$data$GEP

treps = posterior_predict(boutS, resp = "TEMP", nsamples  = NSAMPLES)
ts = boutS$data$TEMP

preps = posterior_predict(boutS, resp = "PPFD", nsamples  = NSAMPLES)
ps = boutS$data$PPFD

p1s = ppc_dens_overlay(gs, greps)
p2s = ppc_dens_overlay(ps, preps)
p3s = ppc_dens_overlay(ts, treps)

p1z = ppc_dens_overlay(gz, grepz)
p2z = ppc_dens_overlay(pz, prepz)
p3z = ppc_dens_overlay(tz, trepz)

diags = cowplot::plot_grid(p1s,p2s,p3s,
                   p1z, p2z, p3z, ncol = 3,
                   labels = c("GS", "PS", "TS", 
                              "GZ", "PZ", "TZ"))


wh = aseries(4);wh
DPI = 300
ggsave("PPCplot.png", plot = diags, width = wh[2], height = wh[1], dpi = DPI, units = "mm")




