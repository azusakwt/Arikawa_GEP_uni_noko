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

dset = readRDS("prepared_brms_data.rds")
dset = dset %>% mutate(PPFD = PPFD * sd_ppfd + mean_ppfd,TEMP = TEMP * sd_temp + mean_temp)

dset %>% 
  summarise(across(c(TEMP, PPFD), list(min, max)))

dset = dset %>% mutate(TEMP = TEMP / 30,
                       PPFD = PPFD / 60) %>% 
  mutate(amonth = sqrt(abs(month - 6))) %>% 
  mutate(location = factor(location),
         sl = interaction(state, location))

library(gamm4)
library(lme4)
outT = gamm4(TEMP ~ s(month, bs = "cc", k = 6) + 
               s(month,  by = sl) + 
               sl, 
             random = ~(1|sl),
             knots  = list(month = seq(0.5, 12.5, length = 6)), 
             data = dset)

outP = gamm4(PPFD ~ s(month, bs = "cc", k = 6) + 
               s(month,  by = sl) + 
               sl, 
             random = ~(1|sl),
             knots  = list(month = seq(0.5, 12.5, length = 6)), 
             data = dset,
             family = Gamma("log"))

dset = dset %>% 
  mutate(TEMP2 = predict(outT$gam, newdata = ., type = "response")) %>% 
  mutate(PPFD2 = predict(outP$gam, newdata = ., type = "response")) 

ggplot(dset) + 
  geom_line(aes(x = month, y = TEMP2)) +
  facet_grid(cols = vars(state),
             rows = vars(location))


ggplot(dset) + 
  geom_line(aes(x = month, y = PPFD2)) +
  facet_grid(cols = vars(state),
             rows = vars(location))

ggplot(dset) + 
  geom_point(aes(x = PPFD2, y = log(GEP))) +
  facet_grid(cols = vars(state),
             rows = vars(location))


ggplot(dset) + 
  geom_point(aes(x = TEMP, y = log(GEP))) +
  facet_grid(cols = vars(state),
             rows = vars(location))

gout = brm(GEP ~ PPFD + TEMP + state + location + (0 + PPFD + TEMP || state / location) + (0 + 1|month/state/location),
                     data = dset,
                     family = Gamma("log"),
                     chains = 4, cores = 4,
                     control = list(adapt_delta = 0.999999,
                                    max_treedepth = 15))
gout

summary(gout)
pp_check(gout)
fixef(gout)
ranef(gout)

X = dset %>% 
  select(PPFD = PPFD2, TEMP = TEMP2, state, location, month) %>% 
  distinct() %>% 
  tidybayes::add_predicted_draws(gout, type = "response")

X %>% 
  group_by(month, state, location) %>% 
  mean_hdci(.prediction) %>% 
  ggplot() + 
  geom_ribbon(aes(x = month, ymin = .lower, ymax = .upper), alpha = 0.2) + 
  geom_line(aes(x = month, y = .prediction)) + 
  geom_point(aes(x = month, y = GEP), data = dset) +
  facet_grid(cols = vars(state),
             rows = vars(location))







dset %>% 
  select(PPFD = PPFD2, TEMP = TEMP2, state, location, month) %>% 
  distinct() %>% 
  tidybayes::add_linpred_draws(gout, type = "response") %>% 
  group_by(month, state, location) %>% 
  select(state, location, month, .draw, .value) %>% 
  pivot_wider(names_from = state,
              values_from = .value) %>% 
  mutate(VD = Vegetated - Desertified) %>% 
  group_by(location, month) %>% 
  mean_hdci(VD) %>% 
  ggplot() +
  geom_ribbon(aes(x = month, ymin = .lower, ymax = .upper), alpha = 0.2) +
  geom_line(aes(x = month, y = VD)) +
  geom_hline(yintercept = 0) +
  facet_grid(rows = vars(location))







