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
library(patchwork)

# Functions
arrhenius = function(x, x0) {
  eV = 0.65
  k = 8.617e-5
  iK = (1 / (273.15 + x0)) - (1/ (273.15 + x))
  exp((eV / k) * iK)
}

# Page size functions -----
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

bseries = function(n) {
  # Function to calculate ISO216 BS Series
  n = as.integer(n)
  if(n == 0) {
    wh = c(1000, 1414)
  } else {
    wh  = aseries(n, floor = F)
    wh2 = aseries(n-1, floor = F)
    w = sqrt(wh[1] * wh2[1])
    h = sqrt(wh[2] * wh2[2])
  }
  return(floor(c(w,h)))
}

aseries(4)
################################################################################
# (1) create_dataset.R 
# (2) prepare_brms_data.R 
dset = readRDS("prepared_brms_data.rds")
# dset = dset %>% mutate(PPFD = PPFD * sd_ppfd + mean_ppfd,TEMP = TEMP * sd_temp + mean_temp)

dset %>% summarise(across(c(TEMP, PPFD), list(min=min, mean=mean, median=median, max=max)))

## GEPの水温補正  
dset = dset %>% 
  # mutate(GEPoriginal = GEP, GEP = GEP / arrhenius(TEMP, 20)) %>% 
  mutate(TEMP = TEMP / 30,
         PPFD = PPFD / 60) %>% 
  mutate(location = factor(location),
         state = factor(state))


################################################################################
WARMUP  = 1500
SAMPLES = 2000
SEED    = 2021
CHAINS = CORES = 4

gepmodel  = bf(GEP ~ PPFD * TEMP * state *　location) + Gamma("log")
CTRL = list(adapt_delta = 0.9995, max_treedepth = 15)

get_prior(gepmodel, data = dset)

PRIOR = c(prior(normal(0, 1), class = b),
          prior(normal(0, 1), class = Intercept),
          prior(normal(0, 1), class = shape))

z1 = lubridate::now()
bout = brm(gepmodel, 
           data = dset, 
           seed = SEED,
           chains = CHAINS, 
           cores = CORES, 
           prior = PRIOR, 
           iter = SAMPLES + WARMUP,
           warmup = WARMUP,
           control = CTRL,
           backend = "cmdstanr",
           threads = 4,
           refresh = 500,
           save_pars = save_pars(all = TRUE))
z2 = lubridate::now()
summary(bout)
fixef(bout)

write_rds(bout, "brms_simplemodel_20210311.rds")

# pp_check(bout, nsamples = 50)

# conditional_effects(bout)

###############################################################
# Read Urchin and Sargassum data ----
sargassum = read_csv("~/Lab_Data/kamigotoR2/sargassum/sargassum.csv")
urchins =   read_csv("~/Lab_Data/kamigotoR2/urchins/urchins.csv")

urchins %>% 
  drop_na() %>% 
  mutate(urchin = recode(urchin, 
                         "hcrassispina" = "Heliocidaris crassispina",
                         "dsetosum" = "Diadema setosum")) %>% 
  group_by(year) %>% 
  summarise_at(vars(individual),
               list(sd = sd, 
                    mean = mean,
                    n = length)) %>% 
  mutate(se = sd/sqrt(n - 1))

sargassum %>% 
  mutate(count = count / 3) %>% 
  group_by(year) %>% 
  summarise(mean = mean (count, na.rm =T))


urchins = urchins %>% 
  drop_na() %>% 
  mutate(urchin = recode(urchin, 
                         "hcrassispina" = "Heliocidaris crassispina",
                         "dsetosum" = "Diadema setosum")) %>% 
  group_by(year, urchin) %>% 
  summarise_at(vars(individual),
               list(sd = sd, 
                    mean = mean,
                    n = length)) %>% 
  mutate(se = sd/sqrt(n - 1))




sargassum = sargassum %>% 
  drop_na() %>% 
  group_by(year) %>% 
  summarise_at(vars(count),
               list(sd = sd, 
                    mean = mean,
                    n = length)) %>% 
  mutate(se = sd/sqrt(n - 1))

# 最低活動水温
# イスズミ Kyphosus lembus
# ノトイスズミ Kyphosys bigibus
kyphosys = 15
# アイゴ Signanus fuscescens
siganus = 19

# ガンガゼ Diadema setosum
diadema = 15
# ムラサキウニ Heliocidaris crassipina
# heliocidaris =
urchins = urchins %>% rename(type = urchin)
sargassum = sargassum %>% mutate(type = "Sargassum macrocarpum")

benthos = full_join(urchins, sargassum)

# dset = readRDS("prepared_brms_data.rds")

# Read BRMS fit ----------------------------------------------------------------
bout = read_rds("brms_simplemodel_20210311.rds")
summary(bout)
# Make predictions of GEP, PPFD, Temperature -----------------------------------

se = function(x, na.rm = T) {
  N = sum(!is.na(x))
  sd(x, na.rm = TRUE) / sqrt(N-1)
}

dset_mean = dset %>% 
  select(location, state, PPFD, TEMP, GEP, month) %>% 
  group_by(location, state, month) %>% 
  summarise(across(c(PPFD, TEMP, GEP), list(mean = mean, se = se)))
  

# 予測区間
predictions = dset_mean %>% 
  mutate(PPFD = PPFD_mean, TEMP = TEMP_mean) %>% 
  ungroup() %>% 
  add_predicted_draws(bout) %>% 
  mean_hdci() %>% 
  filter(str_detect(location, "S"))

# 期待値の区間
fitted = dset_mean %>% ungroup() %>% 
  mutate(PPFD = PPFD_mean, TEMP = TEMP_mean) %>% 
  add_fitted_draws(bout) %>% 
  mean_hdci()%>% 
  filter(str_detect(location, "S"))

dset_mean = dset_mean %>%
  filter(str_detect(location, "S"))

################################################################################

# GEP vs month
CLRS = as.vector(palette.colors(palette = "Okabe-Ito"))
FONTSIZE = 15
DODGE = 0.2
xlabel = "Month"
scale_month = 
  scale_x_continuous("Month", limits = c(0.5, 12.5), 
                     breaks = 1:12,
                     labels = str_sub(month.abb[c(6:12, 1:5)], 1, 1))

ylabel = "GEP~(g~O[2]~m^{-2}~d^{-1})"
ggplot() + 
  geom_ribbon(aes(x = month, ymin = .lower, ymax = .upper, fill = state), data = predictions, alpha = 0.2,
              position = position_dodge(DODGE)) +
  geom_ribbon(aes(x = month, ymin = .lower, ymax = .upper, fill = state), data = fitted, alpha = 0.4,
              position = position_dodge(DODGE)) +
  geom_line(aes(x = month, y = .value, color = state), data = fitted,
            position = position_dodge(DODGE)) +
  geom_point(aes(x = month, y = GEP_mean, color = state),
             position = position_dodge(DODGE),
             data = dset_mean) +
  geom_errorbar(aes(x = month, 
                    ymin = GEP_mean - GEP_se,
                    ymax = GEP_mean + GEP_se, 
                    color = state),
                width = 0,
                position = position_dodge(DODGE),
             data = dset_mean) +
  scale_month +
  scale_y_continuous(name = parse(text = ylabel), limits = c(0, 50), breaks = seq(0, 50, by = 10)) + 
  scale_color_manual(values = CLRS[c(7,3)], labels = c("Desertified [D]", "Vegetated [V]")) +
  scale_fill_manual(values =  CLRS[c(7,3)], labels = c("Desertified [D]", "Vegetated [V]")) +
  ggpubr::theme_pubr(base_size = FONTSIZE) +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"))

wh = aseries(5);wh
DPI = 300
ggsave("GEPvsMonth.png", width = wh[2], height = wh[1], dpi = DPI, units = "mm")

# GEPv - GEPd vs month

dset_vd = dset_mean %>% ungroup() %>% 
  mutate(PPFD = PPFD_mean, TEMP = TEMP_mean) %>% 
  add_fitted_draws(bout) %>% 
  filter(str_detect(location, "S")) %>% 
  ungroup() %>% 
  select(state, month, .value) %>% 
  pivot_wider(names_from = state, 
              values_from = .value,
              values_fn = list) %>% 
  drop_na() %>% 
  unnest(everything()) %>% 
  mutate(VD = Vegetated - Desertified) %>% 
  group_by(month) %>% 
  mean_hdci(VD)
  
ylabel = "Delta*GEP~(g~O[2]~m^{-2}~d^{-1})"

deltaGEPplot = ggplot(dset_vd) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_ribbon(aes(x = month, ymin = .lower, ymax = .upper), 
              alpha = 0.2) +
  geom_line(aes(x = month, y = VD)) +
  scale_month + 
  scale_y_continuous(name = parse(text = ylabel),
                     limits = c(-10, 21)) + 
  ggpubr::theme_pubr(base_size = FONTSIZE) +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white")); deltaGEPplot

wh = aseries(5);wh
DPI = 300
ggsave("GEPVDvsMonth.png", width = wh[2], height = wh[1], plot = deltaGEPplot, dpi = DPI, units = "mm")

# GEP vs temperature

xlabel = "Temperature (°C)"
ylabel = "GEP~(g~O[2]~m^{-2}~d^{-1})"

fitted_temperature = 
  dset_mean %>% ungroup() %>% 
  mutate(PPFD = PPFD_mean, TEMP = TEMP_mean) %>% 
  expand(TEMP = seq(min(TEMP), max(TEMP), length = 9),
         PPFD = seq(min(PPFD), max(PPFD), length = 3),
         month, location, state) %>% 
  add_fitted_draws(bout, allow_new_levels = TRUE) %>% ungroup() %>% 
  group_by(TEMP, PPFD, location, state) %>% 
  mean_hdci(.value)%>% 
  filter(str_detect(location, "S"))


plabel = fitted_temperature%>% ungroup() %>%  select(PPFD) %>% distinct() %>% 
  mutate(PPFD = PPFD * 60) %>%
  mutate(PPFD = round(PPFD, 0)) %>% 
  pull(PPFD)

plabel = str_glue("{plabel}~mol~photons~m^{{-2}}~d^{{-1}}")

xlabel = "Temperature (°C)"
ylabel = "GEP~(g~O[2]~m^{-2}~d^{-1})"
ggplot() + 
  geom_line(aes(x = TEMP, y = .value, color = state, linetype = factor(PPFD)), size = 1, data = fitted_temperature) +
  geom_point(aes(x = TEMP_mean, y = GEP_mean, color = state), data = dset_mean) +
  geom_errorbar(aes(x = TEMP_mean, ymin = GEP_mean - GEP_se,ymax = GEP_mean + GEP_se,  color = state), width = 0, data = dset_mean) +
  scale_x_continuous(xlabel, limits =   c(0.4, 1), breaks = seq(0.4, 1, by  = 0.1), labels = seq(0.4, 1, by  = 0.1) * 30) +
  scale_y_continuous(name = parse(text = ylabel), limits = c(0, 80)) + 
  scale_color_manual(values = CLRS[c(7,3)], labels = c("Desertified [D]", "Vegetated [V]")) +
  scale_linetype_manual(values = 1:3, labels = parse(text = plabel)) +
  ggpubr::theme_pubr(base_size = FONTSIZE) +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5), color = NA),
        strip.text = element_text(size = 14, color = "white"))

wh = aseries(5);wh
DPI = 300
ggsave("GEPvsTemperature.png", width = wh[2], height = wh[1], dpi = DPI, units = "mm")

# GEP vs PPFD

xlabel = "PPFD~(mol~photons~m^{-2}~d^{-1})"
ylabel = "GEP~(g~O[2]~m^{-2}~d^{-1})"

fitted_ppfd = dset_mean %>% ungroup() %>% 
  mutate(PPFD = PPFD_mean, TEMP = TEMP_mean) %>% 
  expand(PPFD = seq(min(PPFD), max(PPFD), length = 9),
         TEMP = seq(min(TEMP), max(TEMP), length = 3),
         month, location, state) %>% 
  add_fitted_draws(bout, allow_new_levels = TRUE) %>% ungroup() %>% 
  group_by(TEMP, PPFD, location, state) %>% 
  mean_hdci(.value)%>% 
  filter(str_detect(location, "S"))

tlabel = fitted_ppfd %>% ungroup() %>%  select(TEMP) %>% distinct() %>% 
  mutate(TEMP = TEMP * 30) %>%
  mutate(TEMP = round(TEMP, 1)) %>% 
  pull(TEMP)

tlabel = str_glue("{tlabel}°C")

ggplot() + 
  geom_line(aes(x =   PPFD, y = .value, color = state, linetype = factor(TEMP)),  size = 1, data = fitted_ppfd) +
  geom_point(aes(x =  PPFD_mean, y = GEP_mean, color = state), data = dset_mean) +
  geom_errorbar(aes(x = PPFD_mean,  ymin = GEP_mean - GEP_se, ymax = GEP_mean + GEP_se, color = state), width = 0, data = dset_mean) +
  scale_x_continuous(name = parse(text = xlabel), limits =   c(0.0, 0.3),
                     breaks = seq(0.0, 0.3, by  = 0.1),
                     labels = seq(0.0, 0.3, by  = 0.1) * 60) +
  scale_y_continuous(name = parse(text = ylabel), limits = c(0, 40)) + 
  scale_color_manual(values = CLRS[c(7,3)], labels = c("Desertified [D]", "Vegetated [V]")) +
  scale_linetype_manual(values = 1:3,  labels = tlabel) +
  ggpubr::theme_pubr(base_size = FONTSIZE) +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5), color = NA),
        strip.text = element_text(size = 14, color = "white"))

wh = aseries(5);wh
DPI = 300
ggsave("GEPvsPPFD.png", width = wh[2], height = wh[1], dpi = DPI, units = "mm")


# GEP vs state

# This uses the expectation
# conditional_effects(bout, "state:location")

fitted_state = bout$data %>% 
  expand(PPFD = mean(PPFD),
         TEMP = mean(TEMP),
         GEP = mean(GEP),
         state = unique(state), 
         location = unique(location)) %>% 
  add_fitted_draws(bout) %>% 
  filter(location == "S") %>% 
  group_by(state, location) %>%  
  summarise(fit = mean(.value),
            sd = sd(.value)) %>% 
  mutate(lower = fit - sd,
         upper = fit + sd)
fitted_state
  
xlabel = "State"
ylabel = "GEP[20]~(g~O[2]~m^{-2}~d^{-1})"

fitted_state %>% 
ggplot() + 
  geom_pointrange(aes(x =  state, y = fit, 
                      ymin = lower, ymax = upper,
                 color = state)) +
  scale_x_discrete(name = "State") +
  scale_y_continuous(name = parse(text = ylabel),
                     limits = c(20, 23)) + 
  scale_color_manual("", 
                     values = CLRS[(c(7,3))]) +
  ggpubr::theme_pubr(base_size = FONTSIZE) +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"))

wh = aseries(5);wh
DPI = 300
ggsave("GEPvsState.png", width = wh[2], height = wh[1], dpi = DPI, units = "mm")

# GEP vs state

# This uses the expectation
# conditional_effects(bout, "state:location")


fitted_state = dset %>% select(state, location) %>% distinct()
tmp = dset %>% ungroup() %>% 
  group_by(month) %>% 
  summarise(PPFD = mean(PPFD),
            TEMP = mean(TEMP))

fitted_state = fitted_state %>% mutate(data = list(tmp))  %>% 
  unnest(data) %>% 
  add_fitted_draws(bout) %>% 
  filter(location == "S") %>% 
  ungroup() %>% 
  select(-.row) %>% 
  pivot_wider(names_from = state,
              values_from = .value) %>% 
 mutate(DV = 100*( 1 - Desertified / Vegetated) )


xlabel = "Percent decrease (%)"
fitted_state %>% mutate(month = factor(month,
                                       levels = c(6:12, 1:5),
                                       labels = month.abb[c(6:12, 1:5)])) %>% 
ggplot() + 
  ggdist::stat_halfeye(aes(x = DV, y = month),
                       fill = CLRS[3]) +
  geom_vline(xintercept = 0 , linetype = "dashed",color = "grey")+
  scale_x_continuous(name = xlabel) +
  scale_y_discrete("Month") +
  scale_color_manual("", 
                     values = CLRS[(c(7,3))]) +
  ggpubr::theme_pubr(base_size = FONTSIZE) +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"))

wh = aseries(5);wh
DPI = 300
ggsave("Difference.png", width = wh[2], height = wh[1], dpi = DPI, units = "mm")



## GEP Percent decrease (only Sargassum)

fitted_state = dset %>% select(state, location) %>% distinct()
tmp = dset %>% ungroup() %>% 
  group_by(month) %>% 
  summarise(PPFD = mean(PPFD), TEMP = mean(TEMP))

fitted_state = fitted_state %>% mutate(data = list(tmp))  %>% 
  unnest(data) %>% 
  add_fitted_draws(bout) %>% 
  filter(location == "S") %>% 
  ungroup() %>% 
  select(-.row) %>% 
  pivot_wider(names_from = state,
              values_from = .value) %>% 
  mutate(DV = 100*( 1 - Desertified / Vegetated)) 

xlabel = "Percent decrease (%)"
ggplot(fitted_state) + 
  ggdist::stat_eye(aes(y = DV, x = month), fill = CLRS[3]) +
  geom_hline(yintercept = 0 , linetype = "dashed",color = "grey")+
  scale_month + 
  scale_y_continuous(name = xlabel) +
  scale_color_manual("",  values = CLRS[(c(7,3))]) +
  ggpubr::theme_pubr(base_size = FONTSIZE) +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5), color = NA),
        strip.text = element_text(size = 14, color = "white"))

ggsave("Percent_decrease.png", width = wh[2], height = wh[1], dpi = DPI, units = "mm")


















