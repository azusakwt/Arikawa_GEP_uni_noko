library(tidyverse)
library(lubridate)
library(gnnlab)
library(lemon)
library(ggfortify)
library(bayesplot)
library(tidybayes)
library(zoo)

library(FactoMineR)
library(factoextra)
library(brms)

# set system ----
Sys.setlocale("LC_TIME","en_US.UTF-8") # ここでやると、ファイルの時刻データに問題がでる（午前・午後）

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


natural = readRDS("matsuda_natural_experiment.rds") %>% ungroup() %>% 
  rename(mean_temp = temp.ave, diff_depth = dif_depth, sum_ppfd = ppfd.add)
tube    = readRDS("matsuda_tube_experiment.rds") %>% ungroup() %>% 
  rename(mean_temp = ave.t,
         sum_ppfd = ppfd.add)

dset = bind_rows(natural = natural,
                 tube = tube,
                .id = "experiment")

pcadset =  dset %>% select(mean_temp, sum_ppfd) %>% PCA(graph = FALSE)

dset2 = bind_cols(dset, as_tibble(get_pca_ind(pcadset)$coord))

dset2 = dset2 %>% select(experiment, sum_ppfd, mean_temp, contains("Dim"), numerator, denominator)

bmodel = bf(numerator|trials(denominator) ~ Dim.1 + Dim.2 + (0 + Dim.1 + Dim.2 | experiment)) + binomial("logit")

bout = brm(bmodel, data = dset2, chains = 4, cores = 4, seed = 2021,
                  backend = "cmdstanr",
                  thread = 4,
           control = list(adapt_delta = 0.99,
                          max_treedepth = 12))

bout %>% tidy_draws() %>% pull(treedepth__) %>% range()

pp_check(bout)

dsetpredict = dset2 %>% 
  expand(mean_temp = seq(min(mean_temp),
                         max(mean_temp), 
                         length = 25),
         sum_ppfd = seq(min(sum_ppfd),
                        max(sum_ppfd),
                        length = 3),
         denominator = median(denominator)) 

dsetpredict = bind_cols(dsetpredict , as_tibble(predict.PCA(pcadset, newdata = dsetpredict)$coor))

dsetpredict = tibble(experiment = unique(dset$experiment),
       dat = list(dsetpredict, dsetpredict)) %>% 
  unnest(dat)

dsetpredict = 
  dsetpredict %>% 
  add_linpred_draws(bout, type = "response") %>% 
  mutate(rate = .value/denominator) %>% 
  group_by(mean_temp, sum_ppfd) %>% 
  mean_hdci(rate) %>% 
  mutate(sum_ppfd = factor(sum_ppfd,
                           labels = c(0, 20, 40)))


xlabel = "平均水温 (°C)"
ylabel = "ノコギリモクの放卵率 (%)"

ggplot(dset2)  +
  geom_point(aes(x = mean_temp,
                 y = numerator/denominator)) +
  geom_line(aes(x = mean_temp,
                y = rate, 
                color = sum_ppfd),
            size = 1,
            data = dsetpredict)  +
  scale_x_continuous(xlabel, limits = c(24, 31)) +
  scale_y_continuous(ylabel) +
  scale_color_manual("光量子量",
                     values = viridis::viridis(4),
                     labels = c("0", "20", "40")) +
  ggpubr::theme_pubr() +
  theme(legend.position=c(1, 1),
        legend.justification=c(1, 1),
        legend.background=element_blank(),
        legend.key.width=unit(10, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.spacing.y = unit(1, "mm"))

wh = aseries(5);wh
DPI = 300
ggsave("kyogikai/eggs.png", width = wh[2], height = wh[1], dpi = DPI, units = "mm")


ggplot(dset2 %>% select(-sum_ppfd))  +
  geom_ribbon(aes(x = mean_temp,
                ymin = .lower,
                ymax = .upper, 
                fill = sum_ppfd),
            alpha = 0.,
            data = dsetpredict)  +
  geom_line(aes(x = mean_temp,
                y = rate, 
                color = sum_ppfd),
            size = 1,
            data = dsetpredict)  +
  geom_point(aes(x = mean_temp,
                 y = numerator/denominator)) +
  scale_x_continuous(xlabel, limits = c(24, 31)) +
  scale_y_continuous(ylabel) +
  scale_color_manual("光量子量",
                     values = viridis::viridis(4)) +
  scale_fill_manual("光量子量",
                     values = viridis::viridis(4)) +
  ggpubr::theme_pubr() +
  facet_rep_grid(rows = vars(sum_ppfd)) +
  theme(legend.position=c(1, 1),
        legend.justification=c(1, 1),
        legend.background=element_blank(),
        legend.key.width=unit(10, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.spacing.y = unit(1, "mm"),
        strip.background = element_blank(),
        strip.text = element_blank())

wh = aseries(5);wh
DPI = 300
ggsave("kyogikai/eggs_temp_facet.png", width = wh[1], height = wh[1], dpi = DPI, units = "mm")










dsetpredict = dset2 %>% 
  expand(mean_temp = seq(min(mean_temp),
                         max(mean_temp), 
                         length = 3),
         sum_ppfd = seq(min(sum_ppfd),
                        max(sum_ppfd),
                        length = 25),
         denominator = median(denominator)) 

dsetpredict = bind_cols(dsetpredict , as_tibble(predict.PCA(pcadset, newdata = dsetpredict)$coor))

dsetpredict = tibble(experiment = unique(dset$experiment),
                     dat = list(dsetpredict, dsetpredict)) %>% 
  unnest(dat)

dsetpredict = 
  dsetpredict %>% 
  add_linpred_draws(bout, type = "response") %>% 
  mutate(rate = .value/denominator) %>% 
  group_by(mean_temp, sum_ppfd) %>% 
  mean_hdci(rate)


xlabel = "光量子量~(mol~m^{-2}~d^{-1})"
ylabel = "ノコギリモクの放卵率 (%)"

ggplot(dset2)  +
  geom_point(aes(x = sum_ppfd,
                 y = numerator/denominator)) +
  geom_line(aes(x = sum_ppfd,
                y = rate, 
                color = as.factor(mean_temp)),
            size = 1,
            data = dsetpredict)  +
  scale_x_continuous(parse(text = xlabel), limits = c(0, 40)) +
  scale_y_continuous(ylabel) +
  scale_color_manual("平均水温",
                     values = viridis::viridis(4),
                     labels = c("25", "28", "31")) +
  ggpubr::theme_pubr() +
  theme(legend.position=c(1, 1),
        legend.justification=c(1, 1),
        legend.background=element_blank(),
        legend.key.width=unit(10, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.spacing.y = unit(1, "mm"))

wh = aseries(5);wh
DPI = 300
ggsave("kyogikai/eggs_ppfd.png", width = wh[2], height = wh[1], dpi = DPI, units = "mm")






dsetpredict %>% ungroup() %>% 
  group_by(sum_ppfd) %>% 
  summarise(min = min(rate),
            max = max(rate))


































