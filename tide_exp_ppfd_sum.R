# tide exp
# 2021/01/16
# yuhei matsuda


#read package----

library(tidyverse)
library(lubridate)
library(gnnlab)
library(rstan)
library(rstanarm)
library(bayesplot)
library(tidybayes)

# read data----
tide_exp = read_csv(file = "~/Lab_Data/matsuday/egg_release/egg_data/release_day/tide_experiment/tide_exp2020.csv")
# set system ----
Sys.setlocale("LC_TIME","en_US.UTF-8") # ここでやると、ファイルの時刻データに問題がでる（午前・午後）


tide_exp_cont =
  tide_exp %>% 
  select(-tube_number) %>%  
  rename(datetime = date) %>% 
  mutate(datetime = ymd_hm(datetime),
         position = as.factor(position)) %>%
  mutate(egg = ifelse(is.na(egg),"Y",egg))%>%
  filter(!str_detect(egg,"D"))

tide_exp_rem =
  tide_exp %>% 
  select(-tube_number) %>%  
  rename(datetime = date) %>% 
  mutate(datetime = ymd_hm(datetime),
         position = as.factor(position)) %>%
  filter(!str_detect(egg,"D|NA")) 


## make graph for tide_exp row data----

tide_exp_cont %>% 
  filter(str_detect(datetime,"08-03|09-03|04|21|22")) %>% 
  ggplot()+
  geom_bar(aes(x = position,fill = egg),position = "dodge")

tide_exp_rem %>% 
  filter(str_detect(datetime,"08-03|09-03|04|21|22")) %>% 
  ggplot()+
  geom_bar(aes(x = position,fill = egg),position = "dodge")

## statistic analyse----

### model with logger data----

# tidedata_2020は使える

### temp data----
fnames_all = dir(path = "~/Lab_Data/kawatea/Oxygen/",pattern = "arikawagaramo_(0m|1m|surface)_200[78]",full.names = TRUE)

exp_temp =
  tibble(fnames = fnames_all) %>%
  mutate(data = map(fnames, function(x){
    read_onset(x)
  })) %>% 
  unnest(data) 

exp_temp=
  exp_temp %>% 
  separate(col = fnames,into = c("a","b","c","d","position","date"),sep = "_") %>% 
  select(position,datetime,temperature) %>% 
  mutate(position = recode(position,
                           "1m" ="middle")) %>% 
  mutate(position = factor(position,labels = c("surface","middle","0m"),levels = c("surface","middle","0m")))



exp_temp =
  exp_temp %>% 
  mutate(date = date(datetime)) %>%
  group_by(date,position) %>% 
  summarise(max.t = max(temperature),
            min.t = min(temperature),
            ave.t = mean(temperature)) %>% 
  mutate(dif.t = max.t-min.t)



###egg data----

tide_exp_rem=
  tide_exp_rem %>% 
  mutate(position = factor(position,labels = c("surface","middle","0m"),levels = c("surface","2m","bottom"))) %>% 
  mutate(egg = ifelse(egg=="N",0,1))


tide_exp_date =
  tide_exp_rem %>% 
  mutate(date = date(datetime)) %>% 
  group_by(date,position) %>% 
  summarise(rate = sum(egg)/length(egg))

###tide data ----

tidedata_2020 = read_csv("~/Lab_Data/matsuday/arikawa2020.csv")
tidedata2020_date =
  tidedata_2020 %>% 
  mutate(Date = as.Date(datetime))


tidedata2020_date =
  tidedata2020_date %>% 
  group_by(Date) %>% 
  summarise(depth_max = max(depth),
            depth_min = min(depth)) %>% 
  mutate(diff_depth = depth_max - depth_min)

tidedata2020_date=
tidedata2020_date %>% 
  rename("date" = Date)



### ppfd data----
tide_light_bottom =  dir("~/Lab_Data/kawatea/Light/",pattern = "arikawagaramo_0m_200[78]",full.names = TRUE)
tide_light_surface = dir("~/Lab_Data/kawatea/Light/",pattern = "arikawaamamo_surface_200[78]",full.names = TRUE)

bottom_light =
  tibble(tide_light_bottom) %>% 
  mutate(data = map(tide_light_bottom,read_odyssey)) %>%
  select(data) %>%
  unnest(data) %>% 
  mutate(ppfd = ppfd*0.150)

surface_light =
  tibble(tide_light_surface) %>% 
  mutate(data = map(tide_light_surface,read_odyssey)) %>% 
  select(data) %>%
  unnest(data) %>% 
  mutate(ppfd = ifelse(str_detect(datetime,"-07-"),
                       ppfd*0.1107191371679865,
                       ppfd*0.1079841270990315))



###microstation（ランバートの式からｋを計算）----

fnames_micro = dir("~/Lab_Data/kawatea/Microstation/",pattern = "arikawaamamo_surface_200[789]",full.names = TRUE)


exp_micro =
  tibble(fnames = fnames_micro) %>%
  mutate(data = map(fnames, function(x){
    read_onset(x)
  })) %>% 
  unnest(data) 

exp_micro =
exp_micro %>% 
  rename(above.ppfd = ppfd) %>% 
  select(datetime,above.ppfd) %>% 
  mutate(datetime = floor_date(datetime,unit = "minutes"))

tide_arikawa_light =
bottom_light %>% 
  rename(b.ppfd = ppfd) %>% 
  left_join(surface_light) %>% 
  rename(s.ppfd = ppfd)


tide_arikawa_light =
tide_arikawa_light %>% 
  left_join(exp_micro)


cal_k =
  tide_arikawa_light %>% 
  right_join(tidedata_2020) %>%
  filter(above.ppfd>=s.ppfd) %>% 
  pivot_longer(cols = c(b.ppfd,s.ppfd,above.ppfd),names_to = "position2",values_to = "ppfd") %>% 
  mutate(position = ifelse(position2 == "b.ppfd",depth,ifelse(position2 == "s.ppfd",0.1,0))) 

#nlsでkの推定

model = function (I0,k,z){
  log(I/I0) =  k * z + inter 
  I =  I0 * exp(k * z + inter) 
}

cal_k2 = cal_k %>%
  filter(date == as.Date("2020-08-01",tz = "JST")) %>%
  filter(position != 0) %>%
  group_by(datetime, H) %>%
  mutate(I0 = max(ppfd),
         I_I0 = ppfd/I0) %>%
  drop_na(I_I0) %>% 
  nest() %>%
  mutate(lmout = map(data, function(x){
    lm(I_I0~position, data = x)
  })) %>% 
  mutate(k = map_dbl(lmout, function(x){coef(x)[[2]]}),
         intercept = map_dbl(lmout, function(x){coef(x)[[1]]}))

cal_k2=
cal_k2 %>% 
  unnest(cols = "data") %>%
  ungroup() %>% 
  mutate(I_2m = I0*exp(k*2)) 
  # ggplot()+
  # geom_line(aes(x = datetime, y = I_2m, color = "2m"))+
  # geom_line(aes(x = datetime, y = ppfd, color = position2))


k =mean(cal_k2$k)# 減少係数

tide_arikawa_light=
tide_arikawa_light %>%
  mutate(I_2m = s.ppfd*exp(k*2)) %>% 
  # right_join(cal_k2) %>% 
  select(datetime,b.ppfd,s.ppfd,I_2m)



light_tide_date_2020 =
  tide_arikawa_light %>% 
  mutate(Date = as.Date(datetime)) %>%
  group_by(Date) %>% 
  summarise("0m" =     sum(b.ppfd),
            surface =  sum(s.ppfd),
            "middle" = sum(I_2m)) %>% 
  ungroup() %>% 
  mutate(across(c(`0m`, surface, middle), 
                function(x) {10 * 60 * x / 1000000}))  # mol/m2/day に変換

light_tide_date_2020 = light_tide_date_2020 %>% 
  pivot_longer(-Date,names_to = "position", values_to = "ppfd.add") %>% # 
  rename(date = Date)

### 統合----

tide_exp_date=
  tide_exp_date %>% 
  left_join(exp_temp,by = c("date", "position")) %>% 
  left_join(tidedata2020_date,by = "date") %>% 
  left_join(light_tide_date_2020) %>% 
  mutate(diff_depth = ifelse(str_detect(position, "mid|sur") ,0,diff_depth))

tide_exp_date=
  tide_exp_date %>% 
  mutate(frac =as.character(MASS::fractions(rate))) %>% 
  mutate(numerator = str_extract(frac, "^[0-9]+")) %>%
  mutate(denominator = str_extract(frac, "/[0-9]+")) %>%
  mutate(denominator = str_remove(denominator, "/")) %>% 
  mutate(denominator = ifelse(is.na(denominator), "1", denominator)) %>% 
  mutate(numerator = as.numeric(numerator),
         denominator = as.numeric(denominator)) %>% 
  mutate(trial = rep(1,n())) 


tide_exp_date =#一通り完成
  tide_exp_date %>% 
  mutate(d.depth = factor(ifelse(diff_depth >0,1,0),levels = c(0,1),labels = c("Static","Dynamic")))

saveRDS(tide_exp_date, "matsuda_tube_experiment.rds")

tide_exp_date_zscore =
  tide_exp_date %>% 
  ungroup() %>% 
  mutate(ave.t = scale(ave.t)[,1],
         ppfd.add = scale(ppfd.add)[,1],
         diff_depth = scale(diff_depth)[,1],
         year = year(date)) %>% 
  mutate(year = as.factor(year))

# PCA----

library(FactoMineR)
library(factoextra)

pca_result =
  tide_exp_date_zscore %>% 
  ungroup() %>% 
  select(ave.t, diff_depth, ppfd.add) %>% 
  PCA()
pca_result$var$coord

pca_result$eig




circleFun <- function(center = c(0,0), diameter = 1, npoints = 100) {
  # Calculate a unit circle for the PCA plots.
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(tibble(x = xx, y = yy))
}

pca_result %>% summary()
X = get_pca_var(pca_result)$coord %>% as_tibble(rownames = "variable")

pcascore = pca_result$eig[,2]

xlab = sprintf("PC1 (%2.1f %%)", pcascore[1])
ylab = sprintf("PC2 (%2.1f %%)", pcascore[2])

X %>% 
  mutate(variable = recode(variable,
                           ave.t = "Temperature",
                           diff_depth = "Depth",
                           ppfd.add = "PPFD")) %>% 
  ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2, color = variable),
               arrow = arrow(angle = 20, length = unit(2.5, "mm"), type = "closed"),
               show.legend = F) +
  geom_segment(aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2, color = variable)) +
  geom_path(aes(x = x, y = y), data = circleFun(diameter = 2))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_equal() +
  scale_x_continuous(xlab) +
  scale_y_continuous(ylab) +
  scale_color_manual(values = viridis::viridis(4)) +
  ggrepel::geom_label_repel(aes(x = x, y = y, label = l),
                            data = . %>%
                              mutate(l = str_extract(variable, "D|P|T"),
                                     x = 0.5*Dim.1,
                                     y = 0.5*Dim.2),
                            direction = "both",
                            nudge_x = c(-0.2, 0.2, 0.0),
                            nudge_y = c( 0.1, 0.2, 0.2),
                            seed = 2020) +
  guides(color = guide_legend(override.aes = list(size = 1),
                              label.position="left")) +
  # annotate(x = Inf, y = Inf, label = "A", vjust = 1, hjust = 1, geom = "text") +
  # annotate(geom = "rect", xmin = -1.0, xmax = 0.0, ymin = -1, ymax = -0.1, fill = "white", color = "white")+
  ggpubr::theme_pubr(FONTSIZE) +
  theme(legend.position=c(0.02, 0.02),
        legend.justification=c(0.,0.00),
        legend.title= element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.key.width=unit(10, "mm"))

ggsave("~/Egg_release/plots_final/tide_pca_ppfd_sum_1_2.png",
       width = 0.8*HEIGHT, 
       height = 0.8*HEIGHT,
       units = "mm",dpi = DPI)

# ggsave("~/Egg_release/plots_seminar_210118/tide_pca_ppfd_sum_1_2.png",width = 150,height = 150,units = "mm",dpi = 500)

xlab = sprintf("PC2 (%2.1f %%)", pcascore[2])
ylab = sprintf("PC3 (%2.1f %%)", pcascore[3])

X %>% 
  mutate(variable = recode(variable,
                           ave.t = "Temperature",
                           diff_depth = "Depth",
                           ppfd.add = "PPFD")) %>% 
  ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = Dim.2, yend = Dim.3, color = variable),
               arrow = arrow(angle = 20, length = unit(2.5, "mm"), type = "closed"),
               show.legend = F) +
  geom_segment(aes(x = 0, y = 0, xend = Dim.2, yend = Dim.3, color = variable)) +
  geom_path(aes(x = x, y = y), data = circleFun(diameter = 2))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_equal() +
  scale_x_continuous(xlab) +
  scale_y_continuous(ylab) +
  scale_color_manual(values = viridis::viridis(4)) +
  ggrepel::geom_label_repel(aes(x = x, y = y, label = l),
                            data = . %>%
                              mutate(l = str_extract(variable, "D|P|T"),
                                     x = 0.5*Dim.2,
                                     y = 0.5*Dim.3),
                            direction = "both",
                            nudge_x = c(-0.2, 0.2, 0.0),
                            nudge_y = c( 0.1, 0.2, 0.2),
                            seed = 2020) +
  guides(color = guide_legend(override.aes = list(size = 1),
                              label.position="left")) +
  # annotate(x = Inf, y = Inf, label = "A", vjust = 1, hjust = 1, geom = "text") +
  # annotate(geom = "rect", xmin = -1.0, xmax = 0.0, ymin = -1, ymax = -0.1, fill = "white", color = "white")+
  ggpubr::theme_pubr(FONTSIZE) +
  theme(legend.position=c(0.02, 0.02),
        legend.justification=c(0.,0.00),
        legend.title= element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.key.width=unit(10, "mm"))

ggsave("~/Egg_release/plots_final/tide_pca_ppfd_sum_2_3.png",
       width  = 0.8*HEIGHT,
       height = 0.8*HEIGHT, units = "mm",dpi = DPI)
# ggsave("~/Egg_release/plots_seminar_210118/tide_pca_ppfd_sum_2_3.png",width = 150,height = 150,units = "mm",dpi = 500)


xlab = sprintf("PC1 (%2.1f %%)", pcascore[1])
ylab = sprintf("PC3 (%2.1f %%)", pcascore[3])


X %>% 
  mutate(variable = recode(variable,
                           ave.t = "Temperature",
                           diff_depth = "Depth",
                           ppfd.add = "PPFD")) %>% 
  ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = Dim.1, yend = Dim.3, color = variable),
               arrow = arrow(angle = 20, length = unit(2.5, "mm"), type = "closed"),
               show.legend = F) +
  geom_segment(aes(x = 0, y = 0, xend = Dim.1, yend = Dim.3, color = variable)) +
  geom_path(aes(x = x, y = y), data = circleFun(diameter = 2))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_equal() +
  scale_x_continuous(xlab) +
  scale_y_continuous(ylab) +
  scale_color_manual(values = viridis::viridis(4)) +
  ggrepel::geom_label_repel(aes(x = x, y = y, label = l),
                            data = . %>%
                              mutate(l = str_extract(variable, "D|P|T"),
                                     x = 0.5*Dim.1,
                                     y = 0.5*Dim.3),
                            direction = "both",
                            nudge_x = c(-0.2, 0.2, 0.0),
                            nudge_y = c( 0.1, 0.2, 0.2),
                            seed = 2020) +
  guides(color = guide_legend(override.aes = list(size = 1),
                              label.position="left")) +
  # annotate(x = Inf, y = Inf, label = "A", vjust = 1, hjust = 1, geom = "text") +
  # annotate(geom = "rect", xmin = -1.0, xmax = 0.0, ymin = -1, ymax = -0.1, fill = "white", color = "white")+
  ggpubr::theme_pubr(FONTSIZE) +
  theme(legend.position=c(0.02, 0.02),
        legend.justification=c(0.,0.00),
        legend.title= element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.key.width=unit(10, "mm"))

ggsave("~/Egg_release/plots_final/tide_pca_ppfd_sum_1_3.png",
       width  = 0.8*HEIGHT,
       height = 0.8*HEIGHT,
       units = "mm",dpi = DPI)
# ggsave("~/Egg_release/plots_seminar_210118/tide_pca_ppfd_sum_2_3.png",width = 150,height = 150,units = "mm",dpi = 500)




# make models----

comps = c("PC1","PC2","PC3") %>% 
  as.tibble()

pca_result$eig %>% 
  as.tibble() %>% 
  bind_cols(comps) %>% 
  select("percentage of variance",value) %>% 
  rename(variance = "percentage of variance") %>%
  mutate(value = as.factor(value),
         variance = variance) %>% 
  ggplot()+
  geom_col(aes(x = value, 
               y = variance))+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black")
  )+
  scale_y_continuous("Proportion of variance (%)") +
  scale_x_discrete("Component")+
  ggpubr::theme_pubr(FONTSIZE+10)

ggsave("~/Egg_release/plots_final/comp_tide_ratio.png", width = WIDTH, height = HEIGHT, unit = "mm", dpi = DPI)

pca_result$ind$coord %>% 
  ggplot()+
  geom_point(aes(x = Dim.1,y = Dim.2))

pca_result$var$coord#寄与率（表の参考に）


library(rstanarm)
library(tidybayes)

pca_score = 
  pca_result$ind$coord %>% 
  as.tibble()


tide_exp_date_zscore =
  tide_exp_date_zscore %>% 
  bind_cols(pca_score)


tide_testmodel = stan_glm(cbind(numerator, denominator - numerator) ~ Dim.1 + Dim.2 + Dim.3,
                         family = binomial("logit"),data = tide_exp_date_zscore,
                         cores = 4, chains = 4, seed = 2020)


summary(tide_testmodel)
xlabel = "Model coefficients"
ylabel = "Coefficient value"
title = "Parameter's coefficient value"
cfs = brms::posterior_summary(tide_testmodel, probs = c(0.025, 0.10, 0.90, 0.975))

as_tibble(cfs, rownames = "cfs") %>% 
  filter(str_detect(cfs, "Dim")) %>% 
  mutate(cfs = recode(cfs, Dim.1 = "PCA1", 
                      Dim.2 = "PCA2",
                      Dim.3 = "PCA3")) %>% 
  ggplot() +
  geom_errorbar(aes(x = cfs,
                    ymin = Q10,
                    ymax = Q90), size = 1,
                width = 0) +
  geom_errorbar(aes(x = cfs, ymin = Q2.5, ymax = Q97.5), 
                size = 0.5,
                width = 0) +
  geom_point(aes(x = cfs, y = Estimate), size = 3) +
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  scale_x_discrete(xlabel) +
  scale_y_continuous(ylabel,
                     limits = c(-3, 1),
                     breaks = seq(-3, 1, by = 0.5))+
  ggpubr::theme_pubr(FONTSIZE)

### model vs raw----

predict_data_d1=
  tide_exp_date_zscore %>% 
  ungroup() %>% 
  tidyr::expand(ave.t = seq(min(ave.t),
                               max(ave.t), length = 11),
                diff_depth = c(min(diff_depth), max(diff_depth)),
                ppfd.add = c(min(ppfd.add), max(ppfd.add)),
                denominator = 6,
                numerator = 0)

predicted_pca_d1 =
  predict_data_d1 %>% 
  select(ave.t, diff_depth, ppfd.add) %>% 
  predict(pca_result, newdata = .)

contrib = 
  pca_result$var$contrib %>% 
  as_tibble(rownames = "variable") %>% 
  mutate(across(contains("Dim"), ~ ./100)) %>% 
  rename(d1 = Dim.1, d2 = Dim.2, d3 = Dim.3)

pout_dim1 =
  predicted_pca_d1$coord %>% as_tibble() %>% 
  bind_cols(predict_data_d1) %>% 
  add_fitted_draws(model = tide_testmodel) %>%  
  ungroup() %>% 
  group_by(diff_depth,ave.t, ppfd.add) %>% 
  # mutate(.value = exp(.value)) %>%
  mean_hdci(.value) %>% 
  rename(fit.response = .value,
         l95.response = .lower,
         u95.response = .upper)

scaling_coef = tide_exp_date %>% ungroup()%>% summarise(across(c(ave.t, diff_depth, ppfd.add), list(m = mean,s = sd)))


xlabel = "Temperature (°C)"
ylabel = "Release rate"

pout_dim1 %>% 
  ungroup() %>% 
  mutate(ave.t = (ave.t * scaling_coef$ave.t_s) + scaling_coef$ave.t_m) %>%
  mutate(diff_depth = factor(diff_depth, labels = c("Static", "Dynamic"))) %>%
  mutate(ppfd.add = factor(ppfd.add, labels = c("Minimum", "Maximum"))) %>% 
  ggplot() + 
  geom_point(aes(x = ave.t,
                 y = rate),
             size = 4,alpha = 0.5,
             data = tide_exp_date) +
  geom_line(aes(x = ave.t,
                y = fit.response,
                linetype = ppfd.add,
                color = diff_depth),
            size = 2)+
  scale_x_continuous(xlabel, limits = c(24, 31)) +
  scale_y_continuous(ylabel) +
  scale_color_manual(values = viridis::viridis(4)) +
  guides(linetype = guide_legend(title = "Daily PPFD", 
                                 override.aes = list(size = 1),
                                 label.position="left",
                                 label.hjust = 0,
                                 order = 2),
         color = guide_legend(title = "ΔDepth",
                              override.aes = list(size = 1),
                              label.position="left",
                              label.hjust = 0,
                              order = 1)) +
  ggpubr::theme_pubr(FONTSIZE+10) +
  theme(legend.position=c(1, 1),
        legend.justification=c(1, 1),
        legend.background=element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.key.width=unit(10, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.spacing.y = unit(1, "mm"))

ggsave("~/Egg_release/plots_final/tide_model_ppfd_sum_fit.png",
       width = WIDTH, height = HEIGHT, unit = "mm", dpi = DPI)

# ggsave("~/Egg_release/plots_seminar_210118/tide_model_ppfd_sum_fit.png", width = 250, height = 150, unit = "mm", dpi = 100)

### 環境毎に変える----

intercept = brms::posterior_samples(tide_testmodel) %>% as_tibble() %>% select(contains("Int"))

cfs = brms::posterior_samples(tide_testmodel) %>% as_tibble() %>% select(contains("Dim")) %>% 
  rename_with(~c("PCA1", "PCA2", "PCA3"),
              everything())

S = pca_result$var$cor %>%
  as_tibble(rownames = "variable") %>%
  rename_with(~c("PCA1", "PCA2", "PCA3"), -variable) %>%
  pivot_longer(-variable, names_to = "cfs") %>%
  mutate(signs = sign(value)) %>% select(-value)


Z = pca_result$var$contrib %>% 
  as_tibble(rownames = "variable") %>% 
  rename_with(~c("PCA1", "PCA2", "PCA3"), -variable) %>% 
  pivot_longer(-variable, names_to = "cfs") 


Y = cfs %>% gather(key = "cfs", value = "coef")

SZ = full_join(Z, S, by = c("cfs", "variable"))

YZ = full_join(Y,SZ, by = c("cfs"))

YZ = YZ %>% mutate(coef = signs * coef * value / 100 * (2 / sqrt(2))) %>% 
  select(cfs, variable, coef)

YZb = YZ %>% group_by(cfs, variable) %>% 
  summarise_at(vars(coef), list(mean = mean, 
                                Q2.5 = ~quantile(., 0.025),
                                Q10 = ~quantile(., 0.10),
                                Q90 = ~quantile(., 0.90),
                                Q97.5= ~ quantile(., 0.975)))

ggplot(YZb) +
  geom_errorbar(aes(x = variable,
                    ymin = Q10,
                    ymax = Q90), size = 1,
                width = 0) +
  geom_errorbar(aes(x = variable, ymin = Q2.5, ymax = Q97.5), 
                size = 0.5,
                width = 0) +
  geom_point(aes(x = variable, y = mean), size = 3) +
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  scale_x_discrete(xlabel) +
  scale_y_continuous(ylabel)+
  ggpubr::theme_pubr()+
  ggtitle(title) +
  facet_grid(rows = vars(cfs))


ggplot(YZb) +
  geom_errorbar(aes(x = cfs,
                    ymin = Q10,
                    ymax = Q90), size = 1,
                width = 0) +
  geom_errorbar(aes(x = cfs, ymin = Q2.5, ymax = Q97.5), 
                size = 0.5,
                width = 0) +
  geom_point(aes(x = cfs, y = mean), size = 3) +
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  scale_x_discrete(xlabel) +
  scale_y_continuous(ylabel)+
  ggpubr::theme_pubr()+
  ggtitle(title) +
  facet_grid(rows = vars(variable))


YZc = YZ %>% 
  pivot_wider(names_from = cfs,
              values_from = coef) %>% 
  unnest(everything()) %>% 
  mutate(coef = PCA1 + PCA2 + PCA3) %>% 
  group_by(variable) %>% 
  
  summarise_at(vars(coef), list(mean = mean, 
                                Q2.5 = ~quantile(., 0.025),
                                Q10 = ~quantile(., 0.10),
                                Q90 = ~quantile(., 0.90),
                                Q97.5= ~ quantile(., 0.975)))





YZc=
  YZc %>% 
  mutate(variable = factor(variable,
                           levels = c("ppfd.add", "ave.t", "diff_depth"),
                           labels = c("Daily PPFD","Mean Temperature","ΔDepth")))


xlabel = "Environmental parameter"
ylabel = "Parameter value"
title = "Parameter's coefficient value"
ggplot(YZc) +
  geom_errorbar(aes(x = variable,
                    ymin = Q10,
                    ymax = Q90), size = 1,
                width = 0) +
  geom_errorbar(aes(x = variable, ymin = Q2.5, ymax = Q97.5), 
                size = 0.5,
                width = 0) +
  geom_point(aes(x = variable, y = mean), size = 3) +
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  scale_x_discrete(xlabel) +
  scale_y_continuous(ylabel, 
                     limits = c(-3, 2),
                     breaks = seq(-3, 2, by = 1))+
  ggpubr::theme_pubr(FONTSIZE+5) 
ggsave(filename = "~/Egg_release/plots_final/tide_parameter_ppfd_sum.png",
       width = WIDTH,height = HEIGHT,dpi = DPI,units = "mm")
# ggsave(filename = "~/Egg_release/plots_seminar_210118/tide_parameter_ppfd_sum.png",width = 200,height = 150,dpi = 150,units = "mm")
YZc %>% mutate(importance = abs(mean)/sum(abs(mean)))
