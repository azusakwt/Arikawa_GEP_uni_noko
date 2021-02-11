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
# make figure---------------
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

dset %>% select(month, month.abb)

dset %>% 
  ggplot() +
  geom_point(aes(x = month, y = GEP, color = state),
             position = position_jitterdodge(0.1, 0, 0.3)) +
  geom_smooth(aes(x = month, y = GEP, color = state),
              method = "gam", formula = y~s(x, bs = "cc", k = 4),
              method.args = Gamma("log")) +
  facet_grid(cols = vars(location)) +
  scale_x_continuous(breaks = 1:12,
                     labels = month.abb[c(6:12, 1:5)],
                     limit = c(0.5, 12.5))



dset %>% ungroup() %>% 
  mutate(PPFD= (PPFD - mean(PPFD))/sd(PPFD)) %>% 
  ggplot() +
  geom_point(aes(x = month, y = PPFD, color = state),
             position = position_jitterdodge(0.1, 0, 0.3)) +
  geom_smooth(aes(x = month, y = PPFD, color = state),
              method = "gam", formula = y~s(x, bs = "cs", k = 6))  +
  scale_x_continuous(breaks = 1:12,
                     labels = month.abb[c(6:12, 1:5)],
                     limit = c(0.5, 12.5))

dset %>% 
  ggplot() +
  geom_point(aes(x = month, y = TEMP, color = state),
             position = position_jitterdodge(0.1, 0, 0.3)) +
  geom_smooth(aes(x = month, y = TEMP, color = state),
              method = "gam", formula = y~s(x, bs = "cs", k = 6),
              method.args = gaussian())  +
  scale_x_continuous(breaks = 1:12,
                     labels = month.abb[c(6:12, 1:5)],
                     limit = c(0.5, 12.5))


dsetZ = dset %>% filter(str_detect(location, "Z"))
dsetS = dset %>% filter(str_detect(location, "S"))


boutZ = tibble(fname = dir(pattern = "boutZ.*rds", full = T)) %>%
  mutate(data = map(fname, readRDS))

boutS = tibble(fname = dir(pattern = "boutS.*rds", full = T)) %>%
  mutate(data = map(fname, readRDS))

# Make predictions
make_pred = function(bout) {
  preddata = bout$data %>% 
    expand(month= seq(1, 12, by = 0.5), 
           state,
           TEMP = 0,
           PPFD = 0)
  
  predtemp = preddata  %>%
    add_linpred_draws(bout, resp = c("TEMP")) %>% 
    group_by(state, month) %>% 
    mean_hdci(.value) %>% 
    select(state, month, temp = .value, temp.lower = .lower, temp.upper = .upper) 
  
  predppfd = preddata %>%
    add_linpred_draws(bout, resp = c("PPFD")) %>% 
    group_by(state, month) %>% 
    mean_hdci(.value) %>% 
    select(state, month, ppfd = .value, ppfd.lower = .lower, ppfd.upper = .upper)
  
  predgep = full_join(predtemp %>% select(state, month, TEMP = temp),
                      predppfd %>% select(state, month, PPFD = ppfd)) %>% 
    add_predicted_draws(bout, resp = "GEP") %>% 
    group_by(state, month, TEMP, PPFD) %>% 
    mean_hdci(.prediction) %>% 
    select(state, month, gep = .prediction, gep.lower = .lower, gep.upper = .upper)
  
  
  full_join(predtemp, predppfd) %>% full_join(predgep)
}


boutZ = boutZ %>% mutate(predictions = map(data, make_pred))
boutS = boutS %>% mutate(predictions = map(data, make_pred))

make_temp_plot = function(X) {
  ggplot(X) + 
    geom_ribbon(aes(x = month, ymin = temp.lower, ymax = temp.upper, fill = state), alpha = 0.2) +
    geom_line( aes(x = month, y = temp, color = state), size = 2) +
    scale_fill_manual(values = viridis::viridis(3)) +
    scale_color_manual(values = viridis::viridis(3)) +
    scale_x_continuous(breaks = 1:12,
                       labels = month.abb[c(6:12, 1:5)],
                       limit = c(0.5, 12.5))  +
    scale_y_continuous("Scaled temperature")
}

make_ppfd_plot = function(X) {
  ggplot(X) + 
    geom_ribbon(aes(x = month, ymin = ppfd.lower, ymax = ppfd.upper, fill = state), alpha = 0.2) +
    geom_line( aes(x = month, y = ppfd, color = state), size = 2) +
    scale_fill_manual(values = viridis::viridis(3)) +
    scale_color_manual(values = viridis::viridis(3)) +
    scale_x_continuous(breaks = 1:12,
                       labels = month.abb[c(6:12, 1:5)],
                       limit = c(0.5, 12.5)) +
    scale_y_continuous("Scaled PPFD")
}  


make_gep_plot = function(X, Z){
  ggplot(X) + 
    geom_ribbon(aes(x = month, ymin = gep.lower, ymax = gep.upper, fill = state),
                alpha = 0.2) +
    geom_line( aes(x = month, y = gep, color = state), size = 2) +
    geom_point(aes(x = month, y = GEP, color = state), data = Z$data,
               alpha = 0.5,
               position = position_jitterdodge(0.2, 0, 0.3))  +
    scale_fill_manual(values = viridis::viridis(3)) +
    scale_color_manual(values = viridis::viridis(3)) +
    scale_x_continuous(breaks = 1:12,
                       labels = month.abb[c(6:12, 1:5)],
                       limit = c(0.5, 12.5)) +
    scale_y_continuous("GEP", limits = c(0, 80))
}


boutZ = boutZ %>% 
  mutate(p1 = map(predictions, make_temp_plot),
         p2 = map(predictions, make_ppfd_plot),
         p3 = map2(predictions, data, make_gep_plot))


make_gep_env_plot = function(bout,m = 1){
  preddata = bout$data
    
  preddata = preddata %>%
    filter(near(month, m)) %>% 
    group_by(month, state) %>% 
    summarise(PPFD = mean(PPFD),
             TEMP  = mean(TEMP))
  
  predgep = preddata %>% 
    add_predicted_draws(bout, resp = "GEP") %>% 
    ungroup() %>%
    group_by(state, month) %>% 
    mean_hdci(.prediction) %>% 
    select(state, month, 
           gep = .prediction, gep.lower = .lower, gep.upper = .upper) %>% print()
  
  ggplot(predgep) + 
    geom_errorbar(aes(x = state, ymin = gep.lower, ymax = gep.upper, color = state),
                width = 0, 
                size = 2,
                position = position_dodge(0.2)) +
    geom_point( aes(x = state, y = gep, color = state), 
                size = 4, position = position_dodge(0.2)) +
    
    scale_fill_manual(values = viridis::viridis(3)) +
    scale_color_manual(values = viridis::viridis(3)) +
    scale_x_discrete() +
    scale_y_continuous("GEP", limits = c(0, 80)) 
}

boutS %>%
  mutate(p4 = map(data, make_gep_env_plot, m = 8)) %>% 
  mutate(exp = str_extract(fname, "exclude|include")) %>% 
  select(exp, p4) %>% 
  pivot_wider(names_from = exp, values_from = p4) %>% 
  pull(include) %>% 
  ggpubr::ggarrange(plotlist = ., ncol = 1, align = "hv",
                    common.legend = TRUE,
                    legend = "bottom")



boutZ %>%
  mutate(p4 = map(data, make_gep_env_plot)) %>% 
  mutate(exp = str_extract(fname, "exclude|include")) %>% 
  select(exp, p4) %>% 
  pivot_wider(names_from = exp, values_from = p4) %>% 
  pull(include) %>% 
  ggpubr::ggarrange(plotlist = ., ncol = 1, align = "hv",
                    common.legend = TRUE,
                    legend = "bottom")


################################################################################

exclude = boutZ %>%
  mutate(exp = str_extract(fname, "exclude|include")) %>% 
  select(exp, p1, p2, p3) %>% 
  pivot_longer(cols = c(p1,p2,p3)) %>% 
  pivot_wider(names_from = exp, values_from = value) %>% 
  pull(exclude) %>% 
  ggpubr::ggarrange(plotlist = ., ncol = 1, align = "hv",
                    common.legend = TRUE,
                    legend = "bottom")

include = boutZ %>%
  mutate(exp = str_extract(fname, "exclude|include")) %>% 
  select(exp, p1, p2, p3) %>% 
  pivot_longer(cols = c(p1,p2,p3)) %>% 
  pivot_wider(names_from = exp, values_from = value) %>% 
  pull(include) %>% 
  ggpubr::ggarrange(plotlist = ., ncol = 1, align = "hv",
                    common.legend = TRUE,
                    legend = "bottom")

pout = ggpubr::ggarrange(exclude, include, ncol = 2,
                         labels = c("Exclude Apr and May",
                                    "Include Apr and May"))

pout
fname = str_glue("Output_Zostera_{today()}.png")
ggsave(filename = fname, 
       plot = pout, 
       width = 300, height = 400, 
       units = "mm")


boutS = boutS %>% 
  mutate(p1 = map(predictions, make_temp_plot),
         p2 = map(predictions, make_ppfd_plot),
         p3 = map2(predictions, data, make_gep_plot))


exclude = boutS %>%
  mutate(exp = str_extract(fname, "exclude|include")) %>% 
  select(exp, p1, p2, p3) %>% 
  pivot_longer(cols = c(p1,p2,p3)) %>% 
  pivot_wider(names_from = exp, values_from = value) %>% 
  pull(exclude) %>% 
  ggpubr::ggarrange(plotlist = ., ncol = 1, align = "hv",
                    common.legend = TRUE,
                    legend = "bottom")

include = boutS %>%
  mutate(exp = str_extract(fname, "exclude|include")) %>% 
  select(exp, p1, p2, p3) %>% 
  pivot_longer(cols = c(p1,p2,p3)) %>% 
  pivot_wider(names_from = exp, values_from = value) %>% 
  pull(include) %>% 
  ggpubr::ggarrange(plotlist = ., ncol = 1, align = "hv",
                    common.legend = TRUE,
                    legend = "bottom")

pout = ggpubr::ggarrange(exclude, include, ncol = 2,
                         labels = c("Exclude Apr and May",
                                    "Include Apr and May"))

pout
fname = str_glue("Output_Sargassum_{today()}.png")
ggsave(filename = fname, 
       plot = pout, 
       width = 300, height = 400, 
       units = "mm")
###############################################################################













################################################################################
summary(x)



x = boutZ %>% slice(2) %>% pull(data) %>% pluck(1)
summary(x)

get_variables(x)
varstoget = get_variables(x) %>% 
  str_subset("b_GEP")

tidy_draws(x) %>% 
  select(varstoget) %>% 
  gather() %>% 
  group_by(key) %>% 
  median_qi(.width = c(0.5, 0.8, 0.95, 0.99)) %>% 
  ggplot() + 
  geom_interval(aes(x = value, xmin = .lower, xmax = .upper, y = key)) +
  geom_vline(xintercept = 0, linetype = "dashed")



tidy_draws(x) %>% 
  dplyr::select(starts_with("zs_")) %>% 
  gather() %>% 
  group_by(key) %>% 
  median_qi(.width = c(0.5, 0.8, 0.95, 0.99)) %>% 
  ggplot() + 
  geom_interval(aes(x = value, xmin = .lower, xmax = .upper, y = key)) +
  geom_vline(xintercept = 0, linetype = "dashed")




varstoget = get_variables(x) %>% 
  str_subset("b_PPFD|b_TEMP|_stateVegetated") %>% print()

tidy_draws(x) %>% 
  select(varstoget) %>% 
  gather() %>% 
  group_by(key) %>% 
  median_qi(.width = c(0.5, 0.8, 0.95, 0.99)) %>% 
  ggplot() + 
  geom_interval(aes(x = value, xmin = .lower, xmax = .upper, y = key)) +
  geom_vline(xintercept = 0, linetype = "dashed")

