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
# Read Urchin and ガラモ場 data ----
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
dset = dset %>% mutate(PPFD = PPFD * sd_ppfd + mean_ppfd,TEMP = TEMP * sd_temp + mean_temp)

dset = dset %>%
  mutate(location = factor(location, levels = c("S", "Z"),
                           labels = c("ガラモ場", "アマモ場")))



ylabel = "ウニ~(ind.~m^{-2})~とノコギリモク~(ind.~50~cm^{-2})~の個体群密度"

urchins %>% 
  mutate(type = factor(type, 
                       levels = c(
                         "hcrassispina",
                         "dsetosum"
                       ))) %>% 
  
  ggplot() + 
  geom_col(aes(x = year, y = mean, fill = type), position = position_dodge(1)) +
  geom_line(aes(x = year, y = mean/4,
                color = "sarg"), data = sargassum,
            size = 2) +
  geom_point(aes(x = year, y = mean/4,
                color = "sarg"), data = sargassum,
            size = 3) +
  scale_y_continuous(parse(text = ylabel), 
                     labels = label_parsed,
                     breaks = seq(0, 50, by = 25), limits = c(0, 50)) + 
  scale_x_continuous("年", breaks = c(2018, 2019, 2020)) +
  scale_fill_manual( values = viridis::viridis(4)[c(1,2)],
                     labels = parse(text = c(
                       "ムラサキウニ~(italic('Diadema setosum'))",
                       "ガンガゼ~(italic('Heliocidaris crassipina'))"
                     ))) +
  scale_color_manual(values = viridis::viridis(4)[3],
                     labels = parse(text ="ノコギリモク~(italic('Sargassum macrocarpum'))"))+
  ggpubr::theme_pubr() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "pt"),
        legend.box.spacing = unit(0, "pt"),
        legend.key.size = unit(0.5, "cm"),
  
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"))

wh = aseries(5);wh
DPI = 300
ggsave("kyogikai/urchins.png", width = wh[1], height = wh[1], dpi = DPI, units = "mm")


  # Read BRMS fit ----------------------------------------------------------------
bout = readRDS("bout_hierarchical.rds")

summary(bout)

bout %>% tidy_draws() %>% select(treedepth__, accept_stat__, divergent__) %>% 
  summarise(across(everything(), list(min = min, max = max)))

# Make predictions of GEP, PPFD, Temperature -----------------------------------
make_pred = function(bout) {
  preddata = bout$data %>% 
    select(month, state, location) %>% 
    distinct() %>% 
    mutate(TEMP = 0) %>% 
    mutate(PPFD = 0)
  
  predtemp = preddata  %>%
    add_linpred_draws(bout, resp = c("TEMP"), scale = "response") %>% 
    group_by(state, location, month) %>% 
    mean_hdci(.value) %>% 
    select(state, location, month, temp = .value, temp.lower = .lower, temp.upper = .upper) 
  
  predppfd = preddata %>%
    add_linpred_draws(bout, resp = c("PPFD"), scale = "response") %>% 
    group_by(state, location, month) %>% 
    mean_hdci(.value) %>% 
    select(state, location, month, ppfd = .value, ppfd.lower = .lower, ppfd.upper = .upper)
  
  predgep = full_join(predtemp %>% select(state, location, month, TEMP = temp),
                      predppfd %>% select(state, location, month, PPFD = ppfd)) %>% 
    add_predicted_draws(bout, resp = "GEP") %>% 
    group_by(state, location, month, TEMP, PPFD) %>% 
    mean_hdci(.prediction) %>% ungroup() %>% 
    select(state, location, month, gep = .prediction, gep.lower = .lower, gep.upper = .upper) %>% ungroup()
  
  meangep = full_join(predtemp %>% select(state, location, month, TEMP = temp),
                      predppfd %>% select(state, location, month, PPFD = ppfd)) %>% 
    add_linpred_draws(bout, resp = "GEP", scale = "response") %>% 
    group_by(state, location, month, TEMP, PPFD) %>% 
    mean_hdci(.value) %>% ungroup() %>% 
    select(state, location, month, e = .value, elower = .lower, eupper = .upper)  %>% ungroup() 
  
  predtemp = preddata  %>%
    add_predicted_draws(bout, resp = c("TEMP"), scale = "response") %>% 
    group_by(state, location, month) %>% 
    mean_hdci(.prediction) %>% 
    select(state, location, month, temp = .prediction, temp.lower = .lower, temp.upper = .upper) 
  
  predppfd = preddata %>%
    add_predicted_draws(bout, resp = c("PPFD"), scale = "response") %>% 
    group_by(state, location, month) %>% 
    mean_hdci(.prediction) %>% 
    select(state, location, month, ppfd = .prediction, ppfd.lower = .lower, ppfd.upper = .upper)
  
  full_join(predtemp, predppfd, by = c("state", "location","month")) %>% full_join(predgep, by = c("state", "location", "month")) %>% 
    full_join(meangep, by = c("state", "location", "month"))
}


predictions = make_pred(bout)
################################################################################

# GEP Plot ---------------------------------------------------------------------
predictions = predictions %>% 
  mutate(location = factor(location, levels = c("S", "Z"),
                           labels = c("ガラモ場", "アマモ場"))) 

ylabel = "'生態系総一次生産量, '*GEP~(g~O[2]~m^{-2}~d^{-1})"


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
  scale_y_continuous(parse(text = ylabel), labels = label_parsed) + 
  scale_x_continuous("月",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = c(6:12, 1:5)) +
  scale_color_manual("植生の状態", values = viridis::viridis(3)c(6:12, 1:5)) +
  scale_fill_manual( "植生の状態", values = viridis::viridis(3),  labels = c("繁茂 (2019以前)", "遷移 (2019以後)")) +
  lemon::facet_rep_wrap(facets = vars(location), ncol = 1) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"));p1

wh = aseries(5);wh
DPI = 300
ggsave("kyogikai/GEP_plot.png", plot = p1, width = wh[1], height = wh[1], dpi = DPI, units = "mm")
################################################################################


make_data2 = function(z) {
  
  preddata = z$data %>% 
    select(month, state, location) %>% 
    distinct() %>% 
    mutate(TEMP = 0) %>% 
    mutate(PPFD = 0)
  
  
  predtemp = preddata  %>%
    add_linpred_draws(z, resp = c("TEMP")) %>% 
    group_by(state, location, month) %>% 
    mean_hdci(.value) %>% 
    select(state, location, month, temp = .value, temp.lower = .lower, temp.upper = .upper) 
  
  predppfd = preddata %>%
    add_linpred_draws(z, resp = c("PPFD")) %>% 
    group_by(state, location, month) %>% 
    mean_hdci(.value) %>% 
    select(state, location, month, ppfd = .value, ppfd.lower = .lower, ppfd.upper = .upper)
  
  predgep = full_join(predtemp %>% select(state, location, month, TEMP = temp),
                      predppfd %>% select(state, location, month, PPFD = ppfd)) %>% 
    add_linpred_draws(z, resp = "GEP", scale = "response")
  
  predgep = predgep %>% 
    ungroup() %>% 
    select(state, location, month, .value, .draw) %>% 
    pivot_wider(values_from = .value,
                names_from = state) %>% 
    mutate(dGEP = Vegetated - Desertified) %>%
    mutate(pGEP = (Vegetated - Desertified)/Vegetated) %>% 
    group_by(month, location) %>% 
    mean_hdci(dGEP, pGEP) %>% ungroup() %>% print() %>% 
    select(location, month, 
           g=dGEP, 
           g.lower=dGEP.lower, 
           g.upper=dGEP.upper,
           pg=pGEP,
           pg.lower = pGEP.lower,
           pg.upper = pGEP.upper) %>% ungroup()
  
  predtemp = 
    preddata  %>%
    add_linpred_draws(z, resp = c("TEMP")) %>% ungroup() %>% 
    select(state, location, month, .value, .draw) %>% 
    pivot_wider(values_from = .value,
                names_from = state) %>% 
    mutate(dTEMP = Vegetated - Desertified) %>% 
    group_by(month, location) %>% 
    mean_hdci(dTEMP) %>%  
    select(location, month, t=dTEMP, t.lower=.lower, t.upper=.upper) 
  
  predppfd = 
    preddata  %>%
    add_linpred_draws(z, resp = c("PPFD")) %>% ungroup() %>% 
    select(state, location, month, .value, .draw) %>% 
    pivot_wider(values_from = .value,
                names_from = state) %>% 
    mutate(dPPFD = Vegetated - Desertified) %>% 
    group_by(month, location) %>% 
    mean_hdci(dPPFD) %>%  
    select(location, month, p=dPPFD, p.lower = .lower, p.upper = .upper) 
  
  full_join(predtemp, predppfd) %>% full_join(predgep) 
}


differences = make_data2(bout)

differences = differences %>% 
  mutate(across(c(t, t.lower, t.upper), ~. * 30)) %>% 
  mutate(across(c(p, p.lower, p.upper), ~. * 60))
# Make GEP and delta GEP plots -------------------------------------------------

differences = differences %>% 
  mutate(location = factor(location, levels = c("S", "Z"),
                           labels = c("ガラモ場", "アマモ場"))) 

ylabel = "繁茂と遷移状態の生態系総一次生産量の差"

p2 = ggplot() + 
  geom_ribbon(aes(x = month, ymin = g.lower, ymax = g.upper), 
              data = differences,
              alpha  = 0.4) +
  geom_line(aes(x = month, y = g),
            data = differences) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_y_continuous(ylabel,
                     limits = c(-10, 15)) + 
  
  scale_x_continuous("月",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = c(6:12, 1:5)) +
  scale_color_manual(values = viridis::viridis(3)) +
  scale_fill_manual(values = viridis::viridis(3)) +
  lemon::facet_rep_wrap(facets = vars(location), ncol = 1) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(0.6,1),
        legend.justification = c(0.0,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5), color = NA),
        strip.text = element_text(size = 14, color = "white"))

wh = aseries(5);wh
DPI = 300
ggsave("kyogikai/dGEP_plot.png", plot = p2, width = wh[1], height = wh[1], dpi = DPI, units = "mm")


################################################################################

p2b = ggplot() + 
  geom_ribbon(aes(x = month, ymin = 100*pg.lower, ymax = 100*pg.upper), 
              data = differences,
              alpha  = 0.4) +
  geom_line(aes(x = month, y = 100*pg),
            data = differences) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_y_continuous("生態系総一次生産の減少率 (%)",
                     breaks = seq(-100, 50, by= 25)) + 
  scale_x_continuous("月",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = c(6:12, 1:5)) +
  scale_color_manual(values = viridis::viridis(3)) +
  scale_fill_manual(values = viridis::viridis(3)) +
  lemon::facet_rep_wrap(facets = vars(location), ncol = 1) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(0.6,1),
        legend.justification = c(0.0,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"))
p2b
wh = aseries(5);wh
DPI = 300
ggsave("kyogikai/pGEP_plot.png", plot = p2b, width = wh[1], height = wh[1], dpi = DPI, units = "mm")


################################################################################
# Delta PPFD and Delta Temperature ---------------------------------------------
ylabel  = "繁茂と遷移状態の光量子量の差~(mol~m^{-2}~d^{-1})"
p3 = ggplot() + 
  geom_ribbon(aes(x = month, ymin = p.lower, ymax = p.upper), 
              data = differences, alpha  = 0.4) +
  geom_line(aes(x = month, y = p),
            data = differences) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_y_continuous(parse(text = ylabel), labels = label_parsed,
                     breaks = c(-10,0,10, 20, 30),
                     limits = c(-10, 30))+
  scale_x_continuous("月",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = c(6:12, 1:5)) +
  scale_color_manual(values = viridis::viridis(3)) +
  scale_fill_manual(values = viridis::viridis(3)) +
  lemon::facet_rep_wrap(facets = vars(location), ncol = 1) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(0.6,1),
        legend.justification = c(0.0,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"));p3

ylabel  = "繁茂と遷移状態の平均水温の差~('°C')"
p4 = ggplot() + 
  geom_ribbon(aes(x = month, ymin = t.lower, ymax = t.upper), 
              data = differences, alpha  = 0.4) +
  geom_line(aes(x = month, y = t),
            data = differences) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_y_continuous(parse(text = ylabel), labels = label_parsed,
                     breaks = c(-3,0,3, 6),
                     limits = c(-3, 6))+
  scale_x_continuous("月",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = c(6:12, 1:5)) +
  scale_color_manual(values = viridis::viridis(3)) +
  scale_fill_manual(values = viridis::viridis(3)) +
  lemon::facet_rep_wrap(facets = vars(location), ncol = 1) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(0.6,1),
        legend.justification = c(0.0,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"))

wh = aseries(5);wh
DPI = 300
ggsave("kyogikai/dPPFD_plot.png", plot = p3, width = wh[1], height = wh[1], dpi = DPI, units = "mm")
ggsave("kyogikai/dTEMP_plot.png", plot = p4, width = wh[1], height = wh[1], dpi = DPI, units = "mm")
################################################################################
# Temperature and PPFD ---------------------------------------------------------
dscaled = dset %>% 
  select(location, state, starts_with("mean"), starts_with("sd")) %>% distinct() 

rescaled = predictions %>% 
  mutate(across(starts_with("ppfd"), ~. * 60)) %>% 
  mutate(across(starts_with("temp"), ~. * 30)) 


ylabel = "平均水温 (°C)"  
p5 = ggplot(rescaled) + 
  geom_ribbon(aes(x = month, ymin = temp.lower, ymax = temp.upper,
                  fill = state), 
              alpha  = 0.2) +
  geom_line(aes(x = month, y = temp, color = state)) +
  geom_point(aes(x = month, y = TEMP, color = state),
             alpha = 0.2,
             data = dset) +
  scale_y_continuous( ylabel, limits = c(10, 31)) +
  scale_x_continuous("月",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = c(6:12, 1:5)) +
  scale_color_manual(values = viridis::viridis(3),  labels = c("繁茂 (2019以前)", "遷移 (2019以後)")) +
  scale_fill_manual(values = viridis::viridis(3),  labels = c("繁茂 (2019以前)", "遷移 (2019以後)")) +
  lemon::facet_rep_wrap(facets = vars(location), ncol = 1) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white")); p5

ylabel = "光量子量~(mol~m^{-2}~d^{-1})"

p6 = ggplot(rescaled) + 
  geom_ribbon(aes(x = month, ymin = ppfd.lower, ymax = ppfd.upper,
                  fill = state), 
              alpha  = 0.2) +
  geom_line(aes(x = month, y = ppfd, color = state)) +
  geom_point(aes(x = month, y = PPFD, color = state),
             alpha = 0.2,
             data = dset) +
  scale_y_continuous( parse(text = ylabel)) +
  scale_x_continuous("月",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = c(6:12, 1:5)) +
  scale_color_manual(values = viridis::viridis(3),  labels = c("繁茂 (2019以前)", "遷移 (2019以後)")) +
  scale_fill_manual(values = viridis::viridis(3),  labels = c("繁茂 (2019以前)", "遷移 (2019以後)")) +
  lemon::facet_rep_wrap(facets = vars(location), ncol = 1) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(0.8,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"));p6

wh = aseries(5);wh
DPI = 300
ggsave("kyogikai/PPFD_plot.png", plot = p6, width = wh[1], height = wh[1], dpi = DPI, units = "mm")
ggsave("kyogikai/TEMP_plot.png", plot = p5, width = wh[1], height = wh[1], dpi = DPI, units = "mm")


################################################################################

# prior_summary(x)
# get_variables(x)

bout %>% tidy_draws()
varstoget = get_variables(bout) %>% str_subset("b_GEP")

cfs = tidy_draws(bout) %>% 
  select(all_of(varstoget)) %>% 
  mutate(I_D_S = b_GEP_Intercept,
         I_V_S = b_GEP_Intercept + b_GEP_stateVegetated,
         I_D_Z = b_GEP_Intercept + b_GEP_locationZ,
         I_V_Z = b_GEP_Intercept + b_GEP_locationZ + b_GEP_stateVegetated) %>% 
  select(starts_with("I"),
         PPFD = b_GEP_PPFD,
         TEMP = b_GEP_TEMP) %>% 
  gather() %>% 
  group_by(key) %>% 
  mutate(key = factor(key,
                      levels = unique(.$key),
                      labels = c("ガラモ場・遷移",
                                 "ガラモ場・繁茂",
                                 "アマモ場・遷移",
                                 "アマモ場・繁茂",
                                 "光量子量", 
                                 "平均水温")))

p7 = ggplot(cfs) + 
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.2) +
  ggdist::stat_dist_halfeye(aes(y = key, 
                                dist = distributional::dist_normal(0, 1),
                                fill = "事後分布"), 
                            data = . %>% tidyr::expand(key),
                            show_interval = F,
                            alpha = 0.5) + 
  ggdist::stat_halfeye(aes(x = value, y = key, 
                           fill = "事前分布 [N(0,1)]"),
                       show_interval = F,
                       alpha = 0.5) +
  scale_x_continuous("母系数の値", breaks = seq(-5, 5, by = 1)) +
  scale_y_discrete("モデル母系数") +
  scale_fill_manual(values = viridis::viridis(3)) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  coord_cartesian(xlim = c(-4,4)) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())
p7
wh = aseries(5);wh
DPI = 300
ggsave("kyogikai/Population_coefficients.png", plot = p7, width = wh[1], height = wh[1], dpi = DPI, units = "mm")

################################################################################

cfs2 = tidy_draws(bout) %>% 
  select(varstoget) %>% 
  mutate(I_S = b_GEP_stateVegetated,
         I_Z = b_GEP_locationZ + b_GEP_stateVegetated) %>% 
  select(starts_with("I")) %>% 
  gather() %>% 
  group_by(key) %>% 
  mutate(key = factor(key,
                      levels = unique(.$key),
                      labels = c("ガラモ場",
                                 "アマモ場")))

p8 = 
  cfs2 %>% mutate(value = exp(value)) %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.2) +
  ggdist::stat_halfeye(aes(x = value, y = key),
                       show_interval = F,
                       alpha = 0.8) +
  scale_x_continuous("繁茂と遷移状態の生態系総一次生産量の差", breaks = seq(0, 8, by = 2)) +
  scale_y_discrete("生態系") +
  scale_fill_manual(values = viridis::viridis(3)) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())
p8
wh = aseries(5);wh
DPI = 300
ggsave("kyogikai/Posterior_difference_intercept.png", plot = p8, width = wh[1], height = wh[1], dpi = DPI, units = "mm")




################################################################################

cfs3 = 
  tidy_draws(bout) %>% 
  select(varstoget) %>% 
  mutate(across(everything(), exp)) %>% 
  mutate(I_S = 1 - b_GEP_Intercept / (b_GEP_Intercept + b_GEP_stateVegetated),
         I_Z = 1 - (b_GEP_Intercept + b_GEP_locationZ) / (b_GEP_Intercept + b_GEP_locationZ + b_GEP_stateVegetated)) %>% 
  mutate(across(c(I_S, I_Z), ~. * 100)) %>% 
  select(starts_with("I")) %>% 
  gather() %>% 
  group_by(key) %>% 
  mutate(key = factor(key,
                      levels = unique(.$key),
                      labels = c("ガラモ場",
                                 "アマモ場")))

cfs3l = cfs3 %>% group_by(key) %>% 
  mean_hdci() %>% 
  mutate(l = sprintf("%2.0f %% (%2.0f - %2.0f %%)", 
                     value, .lower, .upper))

p9 = cfs3 %>% 
  ggplot() + 
  ggdist::stat_halfeye(aes(x = value, y = key),
                       show_interval = F,
                       alpha = 0.8) +
  geom_text(aes(x = value, y = key, label = l),
            data = cfs3l,
            size = 6,
            vjust = 1.2) +
  scale_x_continuous("生態系総一次生産の減少率 (%)", 
                     breaks = seq(0, 100, by = 10),
                     limits = c(0, 100)) +
  scale_y_discrete("") +
  scale_fill_manual(values = viridis::viridis(3)) +
  guides(fill = F) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()); p9
p9
wh = aseries(5);wh
DPI = 300
ggsave("kyogikai/Percent_difference_intercept.png", plot = p9, width = wh[1], height = wh[1], dpi = DPI, units = "mm")



rescaled %>% 
  mutate(kyp = temp > kyphosys,
         sig = temp > siganus,
         dia = temp > diadema) %>% 
  ungroup() %>% 
  group_by(state, location) %>% 
  summarise(across(c(kyp, sig, dia), list(s = sum, l = length))) %>% 
  mutate(kyp = kyp_s / kyp_l,
         sig = sig_s / sig_l,
         dia = dia_s / dia_l)




differences %>% 
  select(location, month, t, p, g)

rescaled %>% 
  select(state, location, month, temp) %>% 
  pivot_wider(names_from = state, 
              values_from = temp) %>% ungroup() %>% 
  filter(month %in% c(8 ,9)) %>% 
  mutate(diff = Vegetated - Desertified) %>% 
  summarise(mean(diff),
            mean(Desertified                 ))


rescaled %>% 
  select(state, location, month, ppfd) %>% 
  pivot_wider(names_from = state, 
              values_from = ppfd) %>% ungroup() %>% 
  filter(month %in% c(1 ,2 ,3)) %>% 
  mutate(diff = Vegetated - Desertified) %>% print() %>% 
  summarise(mean(diff),
            mean(Desertified                 ))

