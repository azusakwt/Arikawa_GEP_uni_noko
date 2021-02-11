# Run the GEP analysis 
# Exclude April and May
# Initial time is June
 
library(tidyverse)
library(lubridate)
library(gnnlab)
library(readxl)
library(ggpubr)

library(brms)
library(tidybayes)
library(multidplyr)

dset = readRDS("prepared_brms_data.rds")

dset = dset %>% filter(str_detect(month.abb, "Apr|May", negate = T))


WARMUP  = 2000
SAMPLES = 2000
SEED    = 2020
CHAINS = CORES = 4

CTRL = list(adapt_delta = 0.99999, max_treedepth = 20)

dsetZ = dset %>% filter(str_detect(location, "Z"))
dsetS = dset %>% filter(str_detect(location, "S"))

ppfdmodel = bf(PPFD ~ s(month, k = 6, by = state) + state) + gaussian()
tempmodel = bf(TEMP ~ s(month, k = 6, by = state) + state) + gaussian()
gepmodel  = bf(GEP ~ PPFD * TEMP * state) + Gamma("log")
gepmodel  = bf(GEP ~ PPFD * TEMP  + (0 + PPFD + TEMP || state) + state) + Gamma("log") 

completemodel = gepmodel + ppfdmodel + tempmodel + set_rescor(FALSE)

get_prior(completemodel, data = dset)

PRIORS = 
  set_prior("normal(0, 1)", class = "b", resp = "TEMP") +
  set_prior("normal(0, 1)", class = "b", resp = "PPFD") +
  set_prior("normal(0, 1)", class = "b", coef = "PPFD", resp = "GEP") +
  set_prior("normal(0, 1)", class = "b", coef = "TEMP", resp = "GEP") +
  set_prior("normal(0, 1)", class = "b", coef = "PPFD:TEMP", resp = "GEP") +
  set_prior("normal(0, 1)", class = "sd", resp = "GEP") +
  set_prior("normal(0, 1)", class = "b", coef = "stateVegetated", resp = "GEP") +
  
  set_prior("student_t(3, 3, 3)", class = "Intercept", resp = "GEP") +
  set_prior("student_t(3, 0, 3)", class = "Intercept", resp = "PPFD") +
  set_prior("student_t(3, 0, 3)", class = "Intercept", resp = "TEMP") +
  
  set_prior("student_t(3, 0, 2)", class = "sds", resp = "PPFD") +
  set_prior("student_t(3, 0, 2)", class = "sds", resp = "TEMP") +
  set_prior("normal(0, 1)", class = "shape", resp = "GEP") +
  set_prior("normal(0, 1)", class = "sigma", resp = "PPFD") +
  set_prior("normal(0, 1)", class = "sigma", resp = "TEMP") 

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
            knots = list(month = c(0.5, 12.5)))
z2 = lubridate::now()

s1 = lubridate::now()
boutS = update(boutZ, 
            newdata = dsetS, 
            seed = SEED,
            chains = CHAINS, cores = CORES, iter = SAMPLES + WARMUP,
            warmup = WARMUP,
            control = CTRL,
            backend = "cmdstanr",
            threads = threading(4),
            knots = list(month = c(0.5, 12.5)))
s2 = lubridate::now()

saveRDS(boutZ, "boutZ_exclude_April_May.rds")
saveRDS(boutS, "boutS_exclude_April_May.rds")

