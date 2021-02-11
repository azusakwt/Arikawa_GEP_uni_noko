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
library(multidplyr)

dset = readRDS("prepared_brms_data.rds")

WARMUP  = 2000
SAMPLES = 2000
SEED    = 2020
CHAINS = CORES = 4

CTRL = list(adapt_delta = 0.99999, max_treedepth = 20)

dsetZ = dset %>% filter(str_detect(location, "Z"))
dsetS = dset %>% filter(str_detect(location, "S"))

dsetZ %>% t.test(GEP ~ state, data = .)
dsetS %>% t.test(GEP ~ state, data = .)

# dsetZ = dsetZ %>% group_by(state, year, month) %>% summarise(across(c(GEP, PPFD, TEMP), mean, na.rm=T))
# dsetS = dsetS %>% group_by(state, year, month) %>% summarise(across(c(GEP, PPFD, TEMP), mean, na.rm=T))

ppfdmodel = bf(PPFD ~ s(month, k = 6, by = state) + state) + gaussian()
tempmodel = bf(TEMP ~ s(month, k = 6, by = state) + state) + gaussian()
gepmodel  = bf(GEP ~ PPFD * TEMP * state) + Gamma("log")
gepmodel  = bf(GEP ~ PPFD + TEMP  +  state) + gaussian()

completemodel = gepmodel + ppfdmodel + tempmodel + set_rescor(TRUE)

get_prior(completemodel, data = dset)

PRIORS = 
  set_prior("lkj(1)", class = "rescor") + 
  set_prior("normal(0, 1)", class = "b", resp = "TEMP") +
  set_prior("normal(0, 1)", class = "b", resp = "PPFD") +
  set_prior("normal(0, 1)", class = "b", resp = "GEP") +
  
  set_prior("student_t(3, 3, 3)", class = "Intercept", resp = "GEP") +
  set_prior("student_t(3, 0, 3)", class = "Intercept", resp = "PPFD") +
  set_prior("student_t(3, 0, 3)", class = "Intercept", resp = "TEMP") +
  
  set_prior("student_t(3, 0, 2)", class = "sds", resp = "PPFD") +
  set_prior("student_t(3, 0, 2)", class = "sds", resp = "TEMP") +
  set_prior("normal(0, 1)", class = "sigma", resp = "GEP") +
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
            refresh = 500,
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
            refresh = 500,
            knots = list(month = c(0.5, 12.5)))
s2 = lubridate::now()


summary(boutZ)
summary(boutS)

boutS$data %>% 
  ggplot() + 
  geom_boxplot(aes(x = state, y = GEP, fill = state))



boutZ$data %>% 
  t.test(GEP ~ state, data = .)


















