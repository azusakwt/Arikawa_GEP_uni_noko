# Run the GEP analysis 
# Include all months
# Initial time is June
# Production version.
# 2021 March 09
################################################################################ 
library(tidyverse)
library(lubridate)
library(gnnlab)
library(readxl)
library(ggpubr)

library(brms)
library(tidybayes)
library(bayesplot)
library(multidplyr)

dset = readRDS("prepared_brms_data.rds")
dset = dset %>% mutate(PPFD = PPFD * sd_ppfd + mean_ppfd,TEMP = TEMP * sd_temp + mean_temp)

dset %>% 
  summarise(across(c(TEMP, PPFD), list(min, max)))

dset = dset %>% mutate(TEMP = TEMP / 30,
                       PPFD = PPFD / 60) %>% 
  mutate(amonth = sqrt(abs(month - 6))) %>% 
  mutate(location = factor(location),
         sl = interaction(state, location))

WARMUP  = 1500
SAMPLES = 2000
SEED    = 2021
CHAINS = CORES = 4

CTRL = list(adapt_delta = 0.9999999, max_treedepth = 20)

dsetZ = dset %>% filter(str_detect(location, "Z"))
dsetS = dset %>% filter(str_detect(location, "S"))

# dset %>% group_by(location, state, month) %>% 
#   mutate(TEMP = TEMP - mean(TEMP)) %>% 
#   ggplot(aes(x = TEMP)) +
#     geom_histogram(bins = 30) +
#   facet_grid(cols = vars(location),
#              rows = vars(state))
# 
# ggplot(dset) +
#   geom_histogram(aes(x = (TEMP))) +
#   facet_grid(cols = vars(location, state))
# 
# ggplot(dset) +
#   geom_histogram(aes(x = PPFD)) +
#   facet_grid(cols = vars(location))
# 
# 
# ggplot(dset) +
#   geom_histogram(aes(x = GEP)) +
#   facet_grid(cols = vars(location))



################################################################################
ppfdmodel = bf(PPFD  ~ s(month, k = 4, bs = "cc") + s(month, by = state) + state + (0+1|state) +  (0+1|month/state)) + Gamma("log")
tempmodel = bf(TEMP  ~ s(month, k = 6, bs = "cc") + s(month, by = state) + s(month, by = location) + state + location + (0+1|state/location) + (0 + 1|month/state/location)) + gaussian()
gepmodel  = bf(GEP ~ PPFD + TEMP + state + location + (0+PPFD+TEMP|state/location) + (0 + 1|month/state/location)) + Gamma("log")
# gepmodel  = bf(GEP ~ PPFD * TEMP  + (0 + PPFD + TEMP || state) + state) + Gamma("log") 

completemodel = gepmodel + ppfdmodel + tempmodel + set_rescor(FALSE)

get_prior(completemodel, data = dset)

PRIORS = get_prior(completemodel, data = dset) %>% 
  mutate(prior = ifelse(str_detect(class, "b|Intercept|b"),
                        "normal(0, 1)", prior)) %>% 
  mutate(prior = ifelse(str_detect(class, "sd"),
                        "normal(0, 1)", prior)) %>% 
  mutate(prior = ifelse(str_detect(class, "shape|sigma"),
                        "normal(0, 1)", prior)) %>% 
  filter(str_detect(coef, "^s\\(", negate = T))

y = readRDS("inverse_mass_matrix.rds")

z1 = lubridate::now()
bout = brm(completemodel, 
           data = dset, 
           seed = SEED,
           chains = CHAINS, 
           cores = CORES, 
           iter = SAMPLES + WARMUP,
           warmup = WARMUP,
           control = CTRL,
           prior = PRIORS,
           backend = "cmdstanr",
           threads = 4,
           refresh = 100,
           knots = list(month = c(0.5, 12.5)),
           output_dir = "sampleroutput",
           save_pars = save_pars(all = TRUE),
           
           inv_metric = y[[1]])

z2 = lubridate::now()
saveRDS(bout, "bout_hierarchical.rds")

sprintf("Running time was %f hours", (as.numeric(z2) -as.numeric(z1))/3600)

pp_check(bout, resp = "GEP")
pp_check(bout, resp = "PPFD")
pp_check(bout, resp = "TEMP")

# 
# 
# x = readRDS("bout_hierarchical.rds")
# 
# library("cmdstanr")
# fname = dir("sampleroutput/", pattern = "file36", full =T)
# x = lapply(fname, read_cmdstan_csv)
# z = lapply(x, function(x) {x$inv_metric})
# y = list(z[[1]][[1]],
#          z[[2]][[2]],
#          z[[3]][[3]],
#          z[[4]][[4]])
# 
# saveRDS(y, "inverse_mass_matrix.rds")
# 
# 
# 
# 
# 
# 