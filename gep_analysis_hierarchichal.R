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
library(furrr)
library(patchwork)

# Functions
arrhenius = function(x, x0) {
  eV = 0.65
  k = 8.617e-5
  iK = (1 / (273.15 + x0)) - (1/ (273.15 + x))
  exp((eV / k) * iK)
}

################################################################################

dset = read_rds("prepared_brms_data.rds")
# dset = dset %>% mutate(PPFD = PPFD * sd_ppfd + mean_ppfd,TEMP = TEMP * sd_temp + mean_temp)

# dset %>% summarise(across(c(TEMP, PPFD), list(min=min, mean=mean, median=median, max=max)))

## GEPの水温補正  
## Don't normalize GEP with temperature, since having temperature in the response and the
## observation does not make any sense.
dset = dset %>% 
  # mutate(GEPoriginal = GEP, GEP = GEP / arrhenius(TEMP, 20)) %>% 
  mutate(TEMP = TEMP / 30, PPFD = PPFD / 60) %>% 
  mutate(location = factor(location))


################################################################################
ppfdmodel = bf(PPFD  ~ s(month, k = 6, bs = "cc") + s(month, by = state) + s(month, by = location) + state + location  + (0 + 1|month/state/location)) + Gamma("log")
tempmodel = bf(TEMP  ~ s(month, k = 6, bs = "cc") + s(month, by = state) + s(month, by = location) + state + location  + (0 + 1|month/state/location)) + gaussian()
gepmodel  = bf(GEP ~ PPFD * TEMP * state * location + (0 + PPFD + TEMP||state/location) + (0 + 1|month/state/location)) + Gamma("log")
# gepmodel  = bf(GEP ~ PPFD * TEMP  + (0 + PPFD + TEMP || state) + state) + Gamma("log") 

completemodel = gepmodel + ppfdmodel + tempmodel + set_rescor(FALSE)

# get_prior(completemodel, data = dset)

PRIORS =
  get_prior(completemodel, data = dset) %>% 
  mutate(prior = ifelse(str_detect(class, "b") & str_detect(coef, "[A-z]"), "normal(0, 1)", prior)) %>% 
  mutate(prior = ifelse(str_detect(class, "^sd$")  & (str_detect(coef, "[A-z]") | str_detect(group, "[A-z]")),"normal(0, 1)", prior)) %>%
  mutate(prior = ifelse(str_detect(class, "^sds$"),"normal(0, 1)", prior)) %>% 
  mutate(prior = ifelse(str_detect(class, "shape|sigma"), "normal(0, 2)", prior)) 

get_prior(completemodel, data = dset) 
  

PRIORS =
  prior(normal(0,1), resp = GEP, class = b) + 
  prior(normal(0,1), resp = PPFD, class = b) + 
  prior(normal(0,1), resp = TEMP, class = b) + 
  prior(normal(0,1), resp = GEP, class = sd) +
  prior(normal(0,1), resp = PPFD, class = sd) + 
  prior(normal(0,1), resp = TEMP, class = sd) + 
  prior(normal(0,1), resp = PPFD, class = sds) +
  prior(normal(0,1), resp = TEMP, class = sds) + 
  prior(normal(0,2), resp = GEP, class = shape) +
  prior(normal(0,2), resp = PPFD, class = shape) + 
  prior(normal(0,2), resp = TEMP, class = sigma) 

make_stancode(completemodel, data = dset, prior = PRIORS, backend = "cmdstanr", threads = 4)
# y = readRDS("inverse_mass_matrix.rds")

WARMUP  = 1500
SAMPLES = 2000
SEED    = 2021
CHAINS = CORES = 4

CTRL = list(adapt_delta = 0.999999999, max_treedepth = 20)
# CTRL = list(adapt_delta = 0.95, max_treedepth = 15)
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
           refresh = 1000,
           knots = list(month = c(0.5, 12.5)),
           output_dir = "sampleroutput",
           save_pars = save_pars(all = TRUE))

z2 = lubridate::now()
sprintf("Running time was %f hours", (as.numeric(z2) -as.numeric(z1))/3600)

FNAME = "bout_hierarchical_final_20210312.rds"
write_rds(bout, FNAME)

aseries = function(n, floor = TRUE) {
  # Function to calculate ISO216 A Series
  n = as.integer(n)
  if(n %% 2) {
    wh = c(w = 1 / (2 ^ ((n + 1)/2)) * 1000 * sqrt(sqrt(2)),
           h = 1 / (2 ^ ((n - 1)/2)) * 1000 / sqrt(sqrt(2)))
  }  else {
    wh = c(w = 1 / (2 ^ (n / 2)) * 1000 / sqrt(sqrt(2)),
           h = 1 / (2 ^ (n / 2)) * 1000 * sqrt(sqrt(2)))
  }
  if(floor) {return(floor(wh))} else {return(wh)}
}
bout = read_rds(FNAME)

wh = aseries(4)

p1 = pp_check(bout, resp = "GEP" , nsamples = 50)
p2 = pp_check(bout, resp = "PPFD", nsamples = 50)
p3 = pp_check(bout, resp = "TEMP", nsamples = 50)
pout = p1 / p2 / p3

ggsave("boutHierarchical_ppcheck.png", plot = pout, width = wh[2], height = wh[1], dpi = 300, units = "mm")

p1 = pp_check(bout, type = "stat_grouped", resp = "GEP", stat = "mean", group = "state")
p2 = pp_check(bout, type = "stat_grouped", resp = "GEP", stat = "mean", group = "location")
pout = p1/p2

ggsave("boutHierarchical_ppcheck_gep.png", plot = pout, width = wh[2], height = wh[1], dpi = 300, units = "mm")

p1 = pp_check(bout, type = "stat_grouped", resp = "PPFD", stat = "mean", group = "state")
p2 = pp_check(bout, type = "stat_grouped", resp = "PPFD", stat = "mean", group = "location")
pout = p1/p2

ggsave("boutHierarchical_ppcheck_ppfd.png", plot = pout, width = wh[2], height = wh[1], dpi = 300, units = "mm")

p1 = pp_check(bout, type = "stat_grouped", resp = "TEMP", stat = "mean", group = "state")
p2 = pp_check(bout, type = "stat_grouped", resp = "TEMP", stat = "mean", group = "location")
pout = p1/p2

ggsave("boutHierarchical_ppcheck_temp.png", plot = pout, width = wh[2], height = wh[1], dpi = 300, units = "mm")






