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
dset = dset %>% mutate(PPFD = PPFD * sd_ppfd + mean_ppfd,TEMP = TEMP * sd_temp + mean_temp)

dset = dset %>%
  mutate(location = factor(location, levels = c("S", "Z"),
                           labels = c("Sargassum", "Zostera")))


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
                           labels = c("Sargassum", "Zostera"))) 

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
  # geom_text(aes(x = month, y = gep, label = l2),
  #           vjust = 1, hjust = 1,
  #           data = plotlabel) +
  scale_y_continuous(parse(text = ylabel), labels = label_parsed) + 
  scale_x_continuous("Month",
                     limits = c(0.5, 12.5),
                     breaks = 1:12,
                     labels = str_sub(month.abb[c(6:12, 1:5)],1,1)) +
  scale_color_manual(values = viridis::viridis(3), labels = c("Desertified [D]", "Vegetated [V]")) +
  scale_fill_manual(values = viridis::viridis(3), labels = c("Desertified [D]", "Vegetated [V]")) +
  lemon::facet_rep_grid(cols = vars(location)) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"))

wh = aseries(5);wh
DPI = 300
ggsave("GEP_plot.png", plot = p1, width = wh[2], height = wh[1], dpi = DPI, units = "mm")
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
                           labels = c("Sargassum", "Zostera"))) 

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
                     limits = c(-10, 15)) + 
  
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
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"))

p21 = cowplot::plot_grid(p1,p2,align = "hv", 
                         axis = "tblhr",
                         ncol = 1)

wh = aseries(5);wh
DPI = 300
ggsave("GEP_dGEP_plot.png", plot = p21, width = wh[2], height = wh[1], dpi = DPI, units = "mm")


################################################################################

p2b = ggplot() + 
  geom_ribbon(aes(x = month, ymin = pg.lower, ymax = pg.upper), 
              data = differences,
              alpha  = 0.4) +
  geom_line(aes(x = month, y = pg),
            data = differences) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  # geom_text(aes(x = month, y = 10, label = l2),
  #           vjust = 1, hjust = 1,
  #           data = plotlabel) +
  scale_y_continuous(parse(text = "(GEP[V]- GEP[D])/GEP[V]~('%')"), labels = label_parsed) + 
  
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
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"))

p21b = cowplot::plot_grid(p1,p2b,align = "hv", 
                         axis = "tblhr",
                         ncol = 1)

wh = aseries(5);wh
DPI = 300
ggsave("GEP_pGEP_plot.png", plot = p21b, width = wh[2], height = wh[1], dpi = DPI, units = "mm")


################################################################################
# Delta PPFD and Delta Temperature ---------------------------------------------
ylabel  = "Delta*PPFD[(V-D)]~(mol~m^{-2}~d^{-1})"
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
                     breaks = c(-10,0,10, 20, 30),
                     limits = c(-10, 30))+
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
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"))

ylabel  = "Delta*Temperature[(V - D)]~('°C')"
p4 = ggplot() + 
  geom_ribbon(aes(x = month, ymin = t.lower, ymax = t.upper), 
              data = differences, alpha  = 0.4) +
  geom_line(aes(x = month, y = t),
            data = differences) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  # geom_text(aes(x = month, y = gep, label = l),
  #           vjust = 1, hjust = 1,
  #           data = plotlabel) +
  # geom_text(aes(x = 12, y = 2, label = l2),
  #           vjust = 1, hjust = 1,
  #           data = plotlabel) +
  scale_y_continuous(parse(text = ylabel), labels = label_parsed,
                     breaks = c(-3,0,3, 6),
                     limits = c(-3, 6))+
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
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"))

p43 = cowplot::plot_grid(p4,p3,align = "hv", ncol = 1, axis = "tblhr")
wh = aseries(5);wh
DPI = 300
ggsave("dPPFD_dTEMP_plot.png", plot = p43, width = wh[2], height = wh[1], dpi = DPI, units = "mm")

################################################################################
# Temperature and PPFD ---------------------------------------------------------
dscaled = dset %>% 
  select(location, state, starts_with("mean"), starts_with("sd")) %>% distinct() 

rescaled = predictions %>% 
  mutate(across(starts_with("ppfd"), ~. * 60)) %>% 
  mutate(across(starts_with("temp"), ~. * 30)) 


ylabel = "Temperature (°C)"  
p5 = ggplot(rescaled) + 
  geom_ribbon(aes(x = month, ymin = temp.lower, ymax = temp.upper,
                  fill = state), 
              alpha  = 0.4) +
  geom_line(aes(x = month, y = temp, color = state)) +
  geom_point(aes(x = month, y = TEMP, color = state),
             alpha = 0.2,
             data = dset) +
  # geom_text(aes(x = month, y = gep, label = l),
  #           vjust = 1, hjust = 1,
  #           data = plotlabel) +
  # geom_text(aes(x = 12, y = 30, label = l2),
  #           vjust = 1, hjust = 1,
  #           data = plotlabel) +
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
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"))

ylabel = "PPFD~(mol~m^{-2}~d^{-1})"

p6 = ggplot(rescaled) + 
  geom_ribbon(aes(x = month, ymin = ppfd.lower, ymax = ppfd.upper,
                  fill = state), 
              alpha  = 0.4) +
  geom_line(aes(x = month, y = ppfd, color = state)) +
  geom_point(aes(x = month, y = PPFD, color = state),
             alpha = 0.2,
             data = dset) +
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
        strip.background = element_rect(fill = grey(0, 0.5),
                                        color = NA),
        strip.text = element_text(size = 14, color = "white"))

p65 = cowplot::plot_grid(p5,p6,align = "hv", ncol = 1, axis = "tblhr")
wh = aseries(5);wh
DPI = 300
ggsave("PPFD_TEMP_plot.png", plot = p65, width = wh[2], height = wh[1], dpi = DPI, units = "mm")



################################################################################

# prior_summary(x)
# get_variables(x)

bout %>% tidy_draws()
varstoget = get_variables(bout) %>% str_subset("b_GEP")

cfs = tidy_draws(bout) %>% 
  select(varstoget) %>% 
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
                      labels = c("Sargassum:Desertified",
                                 "Sargassum:Vegetated",
                                 "Zostera:Desertified",
                                 "Zostera:Vegetated",
                                 "PPFD", 
                                 "Temperature")))

p7 = ggplot(cfs) + 
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.2) +
   ggdist::stat_dist_halfeye(aes(y = key, 
                                 dist = distributional::dist_normal(0, 1),
                                 fill = "Prior"), 
                             data = . %>% tidyr::expand(key),
                             show_interval = F,
                             alpha = 0.5) + 
  ggdist::stat_halfeye(aes(x = value, y = key, 
                           fill = "Posterior"),
                       show_interval = F,
                       alpha = 0.5) +
  scale_x_continuous("Value", breaks = seq(-5, 5, by = 1)) +
  scale_y_discrete("Population coefficient") +
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

wh = aseries(5);wh
DPI = 300
ggsave("Population_coefficients.png", plot = p7, width = wh[2], height = wh[1], dpi = DPI, units = "mm")

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
                      labels = c("Sargassum",
                                 "Zostera")))

p8 = 
  cfs2 %>% mutate(value = exp(value)) %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.2) +
  ggdist::stat_halfeye(aes(x = value, y = key, 
                           fill = "Vegetated - Desertified"),
                       show_interval = F,
                       alpha = 0.5) +
  scale_x_continuous(parse(text = "GEP[V] - GEP[D]"), breaks = seq(-5, 5, by = 1)) +
  scale_y_discrete("Ecosystem") +
  scale_fill_manual(values = viridis::viridis(3)) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

wh = aseries(5);wh
DPI = 300
ggsave("Posterior_difference_intercept.png", plot = p8, width = wh[2], height = wh[1], dpi = DPI, units = "mm")




################################################################################

cfs3 = 
  tidy_draws(bout) %>% 
  select(varstoget) %>% 
  mutate(across(everything(), exp)) %>% 
  mutate(I_S = 1 - b_GEP_Intercept / (b_GEP_Intercept + b_GEP_stateVegetated),
         I_Z = 1 - (b_GEP_Intercept + b_GEP_locationZ) / (b_GEP_Intercept + b_GEP_locationZ + b_GEP_stateVegetated)) %>% 
  select(starts_with("I")) %>% 
  gather() %>% 
  group_by(key) %>% 
  mutate(key = factor(key,
                      levels = unique(.$key),
                      labels = c("Sargassum",
                                 "Zostera")))

cfs3l = cfs3 %>% group_by(key) %>% 
  mean_hdci() %>% 
    mutate(l = sprintf("%0.2f (%0.2f - %0.2f) %%", 
                       value, .lower, .upper))

p9 = cfs3 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.2) +
  ggdist::stat_halfeye(aes(x = value, y = key, 
                           fill = "Vegetated - Desertified"),
                       show_interval = F,
                       alpha = 0.5) +
  geom_text(aes(x = value, y = key, label = l),
            data = cfs3l,
            vjust = -0.5) +
  scale_x_continuous(parse(text = "1 - (GEP[D]/GEP[V])~~('%')"), breaks = seq(-1, 1, by = 0.1),
                     limits = c(0, 1)) +
  scale_y_discrete("Ecosystem") +
  scale_fill_manual(values = viridis::viridis(3)) +
  guides(fill = F) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()); p9

wh = aseries(5);wh
DPI = 300
ggsave("Percent_difference_intercept.png", plot = p9, width = wh[2], height = wh[1], dpi = DPI, units = "mm")










