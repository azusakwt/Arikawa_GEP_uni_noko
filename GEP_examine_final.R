# GEP
# Azusa Kawate, Greg
# 2021/01/09, 01/10

#####################################
#package-------------
library(tidyverse)
library(lubridate)
library(gnnlab)
library(readxl)
library(ggpubr)
library(tidybayes)
library(brms)

###############################################################
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


###############################################################
# Read Urchin and Sargassum data ----
sargassum = read_csv("sargassum.csv")
urchins = read_csv("urchins.csv")

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

dset = readRDS("prepared_brms_data.rds")

# Read BRMS fit ----------------------------------------------------------------
dsetZ = dset %>% filter(str_detect(location, "Z")) # Observations
dsetS = dset %>% filter(str_detect(location, "S")) # Observations

bout = tibble(fname = dir(pattern = "bout[SZ]_include.*rds", full = T)) %>% mutate(data = map(fname, readRDS))

# Make predictions of GEP, PPFD, Temperature -----------------------------------
make_pred = function(bout) {
  preddata = bout$data %>% 
    expand(month= seq(1, 12, by = 0.5), 
           state,
           TEMP = 0,
           PPFD = 0)
  
  predtemp = preddata  %>%
    add_linpred_draws(bout, resp = c("TEMP"), scale = "response") %>% 
    group_by(state, month) %>% 
    mean_hdci(.value) %>% 
    select(state, month, temp = .value, temp.lower = .lower, temp.upper = .upper) 
  
  predppfd = preddata %>%
    add_linpred_draws(bout, resp = c("PPFD"), scale = "response") %>% 
    group_by(state, month) %>% 
    mean_hdci(.value) %>% 
    select(state, month, ppfd = .value, ppfd.lower = .lower, ppfd.upper = .upper)
  
  predgep = full_join(predtemp %>% select(state, month, TEMP = temp),
                      predppfd %>% select(state, month, PPFD = ppfd)) %>% 
    add_predicted_draws(bout, resp = "GEP") %>% 
    group_by(state, month, TEMP, PPFD) %>% 
    mean_hdci(.prediction) %>% ungroup() %>% 
    select(state, month, gep = .prediction, gep.lower = .lower, gep.upper = .upper) %>% ungroup()
  
  meangep = full_join(predtemp %>% select(state, month, TEMP = temp),
                      predppfd %>% select(state, month, PPFD = ppfd)) %>% 
    add_linpred_draws(bout, resp = "GEP", scale = "response") %>% 
    group_by(state, month, TEMP, PPFD) %>% 
    mean_hdci(.value) %>% ungroup() %>% 
    select(state, month, e = .value, elower = .lower, eupper = .upper)  %>% ungroup() 
  
  predtemp = preddata  %>%
    add_predicted_draws(bout, resp = c("TEMP"), scale = "response") %>% 
    group_by(state, month) %>% 
    mean_hdci(.prediction) %>% 
    select(state, month, temp = .prediction, temp.lower = .lower, temp.upper = .upper) 
  
  predppfd = preddata %>%
    add_predicted_draws(bout, resp = c("PPFD"), scale = "response") %>% 
    group_by(state, month) %>% 
    mean_hdci(.prediction) %>% 
    select(state, month, ppfd = .prediction, ppfd.lower = .lower, ppfd.upper = .upper)
  
  full_join(predtemp, predppfd, by = c("state", "month")) %>% full_join(predgep, by = c("state", "month")) %>% 
    full_join(meangep, by = c("state", "month"))
}


bout = bout %>% mutate(predictions = map(data, make_pred))
################################################################################

# GEP Plot ---------------------------------------------------------------------
predictions = bout %>% 
  mutate(location = ifelse(str_detect(fname, "boutS"), "Sargassum", "Zostera")) %>% 
  select(location, predictions) %>% unnest(predictions)

dset = dset %>% mutate(location = recode(location, "Z" = "Zostera",
                                         "S" = "Sargassum"))

plotlabel = dset %>% group_by(location) %>% 
  summarise(gep = max(GEP)) %>% 
  mutate(month = 12, 
         gep = max(gep),
         l = str_glue("({LETTERS[1:2]})"),
         l2 = c("Sargassum", "Zostera"))

ylabel = "GEP~(g~O[2]~m^{-2}~d^{-1})"

p1 = ggplot() + 
  geom_point(aes(x = month, y = GEP, color = state),
             data = dset, alpha = 0.2,
             position = position_jitterdodge(0.1, 0, 0.3)) + 
  geom_ribbon(aes(x = month, ymin = gep.lower, ymax = gep.upper, fill = state),
              data = predictions,
              alpha  = 0.1,
              position = position_jitterdodge(0.1, 0, 0.3)) +
  geom_ribbon(aes(x = month, ymin = elower, ymax = eupper, fill = state), 
              data = predictions,
              alpha  = 0.4,
              position = position_jitterdodge(0.1, 0, 0.3)) +
  geom_line(aes(x = month, y = e, color = state),
            position = position_jitterdodge(0.1, 0, 0.3),
            data = predictions) + 
  geom_text(aes(x = month, y = gep, label = l2),
            vjust = 1, hjust = 1,
            data = plotlabel) +
  scale_y_continuous(parse(text = ylabel), labels = label_parsed) + 
  scale_x_continuous("Month",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = str_sub(month.abb[c(6:12, 1:5)],1,1)) +
  scale_color_manual(values = viridis::viridis(3), labels = c("Desertified [D]", "Vegetated [V]")) +
  scale_fill_manual(values = viridis::viridis(3), labels = c("Desertified [D]", "Vegetated [V]")) +
  lemon::facet_rep_grid(cols = vars(location)) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(0.6,1),
        legend.justification = c(0.0,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

wh = aseries(5);wh
DPI = 300
ggsave("GEP_plot.png", plot = p1, width = wh[2], height = wh[1], dpi = DPI, units = "mm")
################################################################################


make_data2 = function(z) {
  
  preddata = z$data %>% 
    expand(month= seq(1, 12, by = 0.5), 
           state,
           TEMP = 0,
           PPFD = 0)
  
  predtemp = preddata  %>%
    add_linpred_draws(z, resp = c("TEMP")) %>% 
    group_by(state, month) %>% 
    mean_hdci(.value) %>% 
    select(state, month, temp = .value, temp.lower = .lower, temp.upper = .upper) 
  
  predppfd = preddata %>%
    add_linpred_draws(z, resp = c("PPFD")) %>% 
    group_by(state, month) %>% 
    mean_hdci(.value) %>% 
    select(state, month, ppfd = .value, ppfd.lower = .lower, ppfd.upper = .upper)
  
  predgep = full_join(predtemp %>% select(state, month, TEMP = temp),
                      predppfd %>% select(state, month, PPFD = ppfd)) %>% 
    add_linpred_draws(z, resp = "GEP", scale = "response")
  
  predgep = predgep %>% 
    ungroup() %>% 
    select(state, month, .value, .draw) %>% 
    pivot_wider(values_from = .value,
                names_from = state) %>% 
    mutate(dGEP = Vegetated - Desertified) %>% 
    group_by(month) %>% 
    mean_hdci(dGEP) %>% ungroup() %>% 
    select(month, g=dGEP, g.lower=.lower, g.upper=.upper) %>% ungroup()
  
  predtemp = 
    preddata  %>%
    add_linpred_draws(z, resp = c("TEMP")) %>% ungroup() %>% 
    select(state, month, .value, .draw) %>% 
    pivot_wider(values_from = .value,
                names_from = state) %>% 
    mutate(dTEMP = Vegetated - Desertified) %>% 
    group_by(month) %>% 
    mean_hdci(dTEMP) %>%  
    select(month, t=dTEMP, t.lower=.lower, t.upper=.upper) 
  
  predppfd = 
    preddata  %>%
    add_linpred_draws(z, resp = c("PPFD")) %>% ungroup() %>% 
    select(state, month, .value, .draw) %>% 
    pivot_wider(values_from = .value,
                names_from = state) %>% 
    mutate(dPPFD = Vegetated - Desertified) %>% 
    group_by(month) %>% 
    mean_hdci(dPPFD) %>%  
    select(month, p=dPPFD, p.lower = .lower, p.upper = .upper) 
  
  full_join(predtemp, predppfd) %>% full_join(predgep) 
}


bout = bout %>% mutate(diff = map(data, make_data2))

# Make GEP and delta GEP plots -------------------------------------------------

differences = bout %>% 
  mutate(location = ifelse(str_detect(fname, "boutS"), "Sargassum", "Zostera")) %>% 
  select(location, diff) %>% unnest(diff)
dset = dset %>% mutate(location = recode(location, "Z" = "Zostera", "S" = "Sargassum"))

plotlabel = dset %>% group_by(location) %>% 
  summarise(gep = max(GEP)) %>% 
  mutate(month = 12, 
         gep = 15,
         l = str_glue("({LETTERS[1:2]})"),
         l2 = c("Sargassum", "Zostera"))

ylabel = "ΔGEP == GEP[V] - GEP[D]"

p2 = ggplot() + 
  geom_ribbon(aes(x = month, ymin = g.lower, ymax = g.upper), 
              data = differences,
              alpha  = 0.4) +
  geom_line(aes(x = month, y = g),
            data = differences) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  # geom_text(aes(x = month, y = 10, label = l2),
  #           vjust = 1, hjust = 1,
  #           data = plotlabel) +
  scale_y_continuous(parse(text = ylabel), labels = label_parsed,
                     limits = c(-5, 10)) + 
  
  scale_x_continuous("Month",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = str_sub(month.abb[c(6:12, 1:5)],1,1)) +
  scale_color_manual(values = viridis::viridis(3)) +
  scale_fill_manual(values = viridis::viridis(3)) +
  lemon::facet_rep_grid(cols = vars(location)) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(0.6,1),
        legend.justification = c(0.0,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

#

p21 = cowplot::plot_grid(p1,p2,align = "hv", ncol = 1)

wh = aseries(5);wh
DPI = 300
ggsave("GEP_dGEP_plot.png", plot = p21, width = wh[2], height = wh[1], dpi = DPI, units = "mm")

################################################################################
# Delta PPFD and Delta Temperature ---------------------------------------------
ylabel  = "Delta*PPFD == Delta*P[(V-D)]"
p3 = ggplot() + 
  geom_ribbon(aes(x = month, ymin = p.lower, ymax = p.upper), 
              data = differences, alpha  = 0.4) +
  geom_line(aes(x = month, y = p),
            data = differences) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  # geom_text(aes(x = month, y = gep, label = l),
  #           vjust = 1, hjust = 1,
  #           data = plotlabel) +
  # geom_text(aes(x = 1, y = gep, label = l2),
  #           vjust = 1, hjust = 0,
  #           data = plotlabel) +
  scale_y_continuous(parse(text = ylabel), labels = label_parsed,
                     limits = c(-1, 1.1),
                     breaks = c(-1,0,1))+
  scale_x_continuous("Month",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = str_sub(month.abb[c(6:12, 1:5)],1,1)) +
  scale_color_manual(values = viridis::viridis(3)) +
  scale_fill_manual(values = viridis::viridis(3)) +
  lemon::facet_rep_grid(cols = vars(location)) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(0.6,1),
        legend.justification = c(0.0,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

ylabel  = "Temperature == Delta*T[(V - D)]"
p4 = ggplot() + 
  geom_ribbon(aes(x = month, ymin = t.lower, ymax = t.upper), 
              data = differences, alpha  = 0.4) +
  geom_line(aes(x = month, y = t),
            data = differences) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  # geom_text(aes(x = month, y = gep, label = l),
  #           vjust = 1, hjust = 1,
  #           data = plotlabel) +
  geom_text(aes(x = 12, y = 2, label = l2),
            vjust = 1, hjust = 1,
            data = plotlabel) +
  scale_y_continuous(parse(text = ylabel), labels = label_parsed,
                     limits = c(-1, 2),
                     breaks = c(-1,0,1, 2))+
  scale_x_continuous("Month",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = str_sub(month.abb[c(6:12, 1:5)],1,1)) +
  scale_color_manual(values = viridis::viridis(3)) +
  scale_fill_manual(values = viridis::viridis(3)) +
  lemon::facet_rep_grid(cols = vars(location)) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(0.6,1),
        legend.justification = c(0.0,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

p43 = cowplot::plot_grid(p4,p3,align = "hv", ncol = 1)
wh = aseries(5);wh
DPI = 300
ggsave("dPPFD_dTEMP_plot.png", plot = p43, width = wh[2], height = wh[1], dpi = DPI, units = "mm")
################################################################################
# Temperature and PPFD ---------------------------------------------------------
dscaled = dset %>% 
  select(location, state, starts_with("mean"), starts_with("sd")) %>% distinct() %>% 
  mutate(location = recode(location,
                           "S" = "Sargassum",
                           "Z" = "Zostera")) %>% print()

rescaled = full_join(dscaled, predictions) %>% 
  mutate(across(starts_with("ppfd"), ~. * sd_ppfd + mean_ppfd)) %>% 
  mutate(across(starts_with("temp"), ~. * sd_temp + mean_temp)) 


ylabel = "Temperature (°C)"  
p5 = ggplot(rescaled) + 
  geom_ribbon(aes(x = month, ymin = temp.lower, ymax = temp.upper,
                  fill = state), 
              alpha  = 0.4) +
  geom_line(aes(x = month, y = temp, color = state)) +
  geom_point(aes(x = month, y = TEMP, color = state),
             alpha = 0.2,
             data = mutate(dset, location = recode(location,
                                                   "S" = "Sargassum",
                                                   "Z" = "Zostera")) %>% 
               mutate(across(starts_with("TEMP"), ~. * sd_temp + mean_temp)) ) +
  # geom_text(aes(x = month, y = gep, label = l),
  #           vjust = 1, hjust = 1,
  #           data = plotlabel) +
  geom_text(aes(x = 12, y = 30, label = l2),
            vjust = 1, hjust = 1,
            data = plotlabel) +
  scale_y_continuous( ylabel, limits = c(0, 31)) +
  scale_x_continuous("Month",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = str_sub(month.abb[c(6:12, 1:5)],1,1)) +
  scale_color_manual(values = viridis::viridis(3)) +
  scale_fill_manual(values = viridis::viridis(3)) +
  lemon::facet_rep_grid(cols = vars(location)) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

ylabel = "PPFD~(mol~m^{-2}~d^{-1})"

p6 = ggplot(rescaled) + 
  geom_ribbon(aes(x = month, ymin = ppfd.lower, ymax = ppfd.upper,
                  fill = state), 
              alpha  = 0.4) +
  geom_line(aes(x = month, y = ppfd, color = state)) +
  geom_point(aes(x = month, y = PPFD, color = state),
             alpha = 0.2,
             data = mutate(dset, location = recode(location,
                                                   "S" = "Sargassum",
                                                   "Z" = "Zostera")) %>% 
               mutate(across(starts_with("PPFD"), ~. * sd_ppfd + mean_ppfd)) ) +
  # geom_text(aes(x = month, y = gep, label = l),
  #           vjust = 1, hjust = 1,
  #           data = plotlabel) +
  # geom_text(aes(x = 1, y = gep, label = l2),
  #           vjust = 1, hjust = 0,
  #           data = plotlabel) +
  scale_y_continuous( parse(text = ylabel)) +
  scale_x_continuous("Month",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = str_sub(month.abb[c(6:12, 1:5)],1,1)) +
  scale_color_manual(values = viridis::viridis(3)) +
  scale_fill_manual(values = viridis::viridis(3)) +
  lemon::facet_rep_grid(cols = vars(location)) +
  ggpubr::theme_pubr() +
  theme(legend.position = "none",
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

p65 = cowplot::plot_grid(p5,p6,align = "hv", ncol = 1)
wh = aseries(5);wh
DPI = 300
ggsave("PPFD_TEMP_plot.png", plot = p65, width = wh[2], height = wh[1], dpi = DPI, units = "mm")



################################################################################

bout %>% mutate(summary = map(data, summary)) %>% pull(summary)

x = bout %>% slice(2) %>% pull(data) %>% pluck(1)

# prior_summary(x)
# get_variables(x)


get_cfs = function(x) {
  
  varstoget = get_variables(x) %>% str_subset("b_GEP")
  tidy_draws(x) %>% 
    select(varstoget) %>% 
    gather() %>% 
    group_by(key) %>% 
    # median_qi(.width = c(0.5, 0.8, 0.95, 0.99)) %>% 
    mutate(key = str_remove(key, "b_GEP_")) %>% 
    mutate(key = factor(key, levels = c("Intercept", 
                                        "stateVegetated",
                                        "TEMP", 
                                        "PPFD",
                                        "PPFD:TEMP",
                                        "PPFD:stateVegetated",
                                        "TEMP:stateVegetated",
                                        "PPFD:TEMP:stateVegetated"),
                        labels = c("Intercept [D]",
                                   "Intercept [V - D]",
                                   "Temperature [T]", 
                                   "PPFD [P]", 
                                   "P × T ",
                                   "P × State[V]",
                                   "T × State[V]",
                                   "T × P × State[V]"))) %>% 
    ungroup()
}

bout = bout %>% mutate(cfs = map(data, get_cfs))

cfs = bout %>% 
  mutate(location = ifelse(str_detect(fname, "boutS"), "Sargassum", "Zostera")) %>% 
  select(location, cfs) %>% 
  mutate(location = recode(location,
                           "S" = "Sargassum",
                           "Z" = "Zostera")) %>% 
  ungroup() %>% 
  unnest(cfs) 

cfs = cfs %>% filter(!str_detect(as.character(key), "Intercept \\[D\\]"))

plotlabel = cfs %>% expand(location) %>% 
  mutate(l = str_glue("({LETTERS[1:2]})"),
         l2 = c("Sargassum", "Zostera"))

p7 = ggplot(cfs) + 
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.2) +
   ggdist::stat_dist_halfeye(aes(y = key, 
                                 dist = distributional::dist_student_t(3,0,1),
                                 fill = "Prior"), 
                             data = . %>% tidyr::expand(key)%>% filter(!str_detect(as.character(key), "Intercept \\[D\\]")),
                             show_interval = F,
                             alpha = 0.5) + 
  ggdist::stat_halfeye(aes(x = value, y = key, 
                           fill = "Posterior"),
                       show_interval = F,
                       alpha = 0.5) +
  geom_text(aes(x = 0, y = 0.8, label = l2),
            hjust = 0.5,
            vjust = 0,
            data = plotlabel) +
  scale_x_continuous("Value", breaks = seq(-1, 5, by = 1)) +
  scale_y_discrete("Population coefficient") +
  scale_fill_manual(values = viridis::viridis(3)) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  coord_cartesian(xlim = c(-1,1)) +
  facet_grid(cols = vars(location)) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())
wh = aseries(5);wh
DPI = 300
ggsave("Population_coefficients.png", plot = p7, width = wh[2], height = wh[1], dpi = DPI, units = "mm")

################################################################################


get_pt_epsilon = function(x,z) {
  x$data %>% 
    expand(state, month = 1:12, 
           TEMP = 0, PPFD = 0) %>% 
    add_fitted_draws(x, 
                     resp = z, 
                     dpar = "sigma", n = 1000) %>%  
    ungroup() 
}

bout = bout %>% 
  mutate(st = map(data, get_pt_epsilon, z = "TEMP")) %>% 
  mutate(sp = map(data, get_pt_epsilon, z = "PPFD"))

p1 = bout$st[[1]] %>% 
  ggplot() + 
  ggdist::stat_dist_halfeye(aes(x = month, 
                                dist = distributional::dist_normal(0,1),
                                fill = "Prior"), 
                            show_interval = F,
                            alpha = 0.5,
                            data = . %>% 
                              expand(state, month = 1:12, 
                                     TEMP = 0, PPFD = 0)) + 
  ggdist::stat_halfeye(aes(x = month, y = sigma, 
                           fill = "Posterior"),
                       show_interval = F,
                       alpha = 0.5) +
  scale_x_continuous("Month",
                     breaks = 1:12,
                     labels = str_sub(month.abb[c(6:12, 1:5)],1,1)) +
  scale_y_continuous(parse(text = "sigma[T]")) +
  scale_fill_manual(values = viridis::viridis(3)) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  coord_cartesian(ylim = c(0, 1)) + 
  facet_grid(rows = vars(state)) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.title = element_blank()) +
  facet_grid(rows = vars(state))

p2 = bout$st[[2]] %>% 
  ggplot() + 
  ggdist::stat_dist_halfeye(aes(x = month, 
                                dist = distributional::dist_normal(0,1),
                                fill = "Prior"), 
                            show_interval = F,
                            alpha = 0.5,
                            data = . %>% 
                              expand(state, month = 1:12, 
                                     TEMP = 0, PPFD = 0)) + 
  ggdist::stat_halfeye(aes(x = month, y = sigma, 
                           fill = "Posterior"),
                       show_interval = F,
                       alpha = 0.5) +
  scale_x_continuous("Month",
                     breaks = 1:12,
                     labels = str_sub(month.abb[c(6:12, 1:5)],1,1)) +
  scale_y_continuous(parse(text = "sigma[T]")) +
  scale_fill_manual(values = viridis::viridis(3)) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  coord_cartesian(ylim = c(0, 1)) + 
  facet_grid(rows = vars(state)) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.title = element_blank()) +
  facet_grid(rows = vars(state))

pout1 = cowplot::plot_grid(p2,p1,align = "hv", ncol = 1,
                          label_x = c(0.5,0.5),
                          hjust = c(0.5,0.5),
                          labels = c("Zostera", "Sargassum"))

p1 = bout$sp[[1]] %>% 
  ggplot() + 
  ggdist::stat_dist_halfeye(aes(x = month, 
                                dist = distributional::dist_normal(0,1),
                                fill = "Prior"), 
                            show_interval = F,
                            alpha = 0.5,
                            data = . %>% 
                              expand(state, month = 1:12, 
                                     TEMP = 0, PPFD = 0)) + 
  ggdist::stat_halfeye(aes(x = month, y = sigma, 
                           fill = "Posterior"),
                       show_interval = F,
                       alpha = 0.5) +
  scale_x_continuous("Month",
                     breaks = 1:12,
                     labels = str_sub(month.abb[c(6:12, 1:5)],1,1)) +
  scale_y_continuous(parse(text = "sigma[P]")) +
  scale_fill_manual(values = viridis::viridis(3)) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  coord_cartesian(ylim = c(0, 2)) + 
  facet_grid(rows = vars(state)) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.title = element_blank()) +
  facet_grid(rows = vars(state))

p2 = bout$sp[[2]] %>% 
  ggplot() + 
  ggdist::stat_dist_halfeye(aes(x = month, 
                                dist = distributional::dist_normal(0,1),
                                fill = "Prior"), 
                            show_interval = F,
                            alpha = 0.5,
                            data = . %>% 
                              expand(state, month = 1:12, 
                                     TEMP = 0, PPFD = 0)) + 
  ggdist::stat_halfeye(aes(x = month, y = sigma, 
                           fill = "Posterior"),
                       show_interval = F,
                       alpha = 0.5) +
  scale_x_continuous("Month",
                     breaks = 1:12,
                     labels = str_sub(month.abb[c(6:12, 1:5)],1,1)) +
  scale_y_continuous(parse(text = "sigma[P]")) +
  scale_fill_manual(values = viridis::viridis(3)) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  coord_cartesian(ylim = c(0, 2)) + 
  facet_grid(rows = vars(state)) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.title = element_blank()) +
  facet_grid(rows = vars(state))


pout2 = cowplot::plot_grid(p2,p1,align = "hv", ncol = 1,
                          label_x = c(0.5,0.5),
                          hjust = c(0.5,0.5),
                          labels = c("Zostera", "Sargassum"))

pout12 = cowplot::plot_grid(pout1, pout2, ncol = 2)

wh = aseries(4);wh
DPI = 300
ggsave("Sigma.png", plot = pout12, width = wh[2], height = wh[1], dpi = DPI, units = "mm")


################################################################################

get_epsilon = function(x) {
  varstoget = get_variables(x) %>% str_subset("^shape")
  tidy_draws(x) %>% 
    select(varstoget) %>% 
    gather() %>% 
    group_by(key) %>% 
    # median_qi(.width = c(0.5, 0.8, 0.95, 0.99)) %>% 
    mutate(key = factor(key, levels = c("shape_GEP"),
                        labels = c("shape[G]"))) %>% ungroup()
}

bout = bout %>% mutate(sigma = map(data, get_epsilon))

sigma = bout %>% 
  mutate(location = ifelse(str_detect(fname, "boutS"), "Sargassum", "Zostera")) %>% 
  select(location, sigma) %>% 
  mutate(location = recode(location,
                           "S" = "Sargassum",
                           "Z" = "Zostera")) %>% 
  ungroup() %>% 
  unnest(sigma) 

plotlabel = cfs %>% expand(location) %>% 
  mutate(l = str_glue("({LETTERS[1:2]})"),
         l2 = c("Sargassum", "Zostera"))

p8 = ggplot(sigma) + 
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.2) +
  ggdist::stat_dist_halfeye(aes(y = key, 
                                dist = distributional::dist_normal(0,1),
                                fill = "Prior"), 
                            data = sigma %>% expand(key),
                            show_interval = F,
                            alpha = 0.5) + 
  ggdist::stat_halfeye(aes(x = value, y = key, 
                           fill = "Posterior"),
                       show_interval = F,
                       alpha = 0.5,
                       data = sigma) +
  geom_text(aes(x = 5, y = 0.8, label = l2),
            hjust = 1,
            vjust = 0,
            data = plotlabel) +
  scale_x_continuous("Value", breaks = seq(-1, 7, by = 1)) +
  scale_y_discrete("Epsilon coefficient",
                   labels = label_parsed) +
  scale_fill_manual(values = viridis::viridis(3)) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  coord_cartesian(xlim = c(-1, 8)) + 
  facet_grid(rows = vars(location)) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())
wh = aseries(5);wh
DPI = 300
ggsave("Epsilon_coefficients.png", plot = p8, width = wh[2], height = wh[1], dpi = DPI, units = "mm")











