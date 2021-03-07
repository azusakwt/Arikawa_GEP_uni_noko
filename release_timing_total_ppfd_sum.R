# 2021/01/13
# Yuhei Matsuda
# analys total data of release date

# read package----

library(tidyverse)
library(lubridate)
library(gnnlab)
library(lemon)
library(ggfortify)
library(rstan)
library(bayesplot)
library(tidybayes)
library(zoo)

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


# 2019 ----

## read data----

Egg_2019 = read_csv("~/Lab_Data/matsuday/egg_release/egg_data/release_day/卵放出日時表_2019.csv")
light_file_2019 = dir("~/Lab_Data/kawatea/Light/", pattern = "arikawagaramo_0m_190[78]", full.names = TRUE)
cnames = c("No", "date","time","rawvalue","calibrated_light")
water_temp_file_2019 = dir("~/Lab_Data/kawatea/Oxygen/", pattern = "arikawagaramo_0m_190[789]", full.names = TRUE)
dnames = c("No","YYYY/MM/DD","hh:mm:ss","Day","Velo[cm/s]","Dir[Deg]","Vel_EW[cm/s]","Vel_NS[cm/s]","Temp","Comp_A","Comp_B","Vel X[cm/s]","Vel Y[cm/s]","Battery[Volt]")
tide_file = dir("/home/matsuday/Egg_release/tide_data/",pattern = "2019",full.names = TRUE)
tidedata2019 = read_csv(tide_file)
weather_data_2019 = read_csv("~/Lab_Data/weather/201701_20208_Arikawa_JMA_Data.csv")

## light data edit----
light_data_2019=
  tibble(light_file_2019) %>%
  mutate(data = map(light_file_2019, read_odyssey)) %>%
  unnest(data) %>%
  filter(datetime > ymd_h("2019-08-01 0") & datetime < ymd_h("2019-09-05 0")) %>%  
  mutate(calib_value = ppfd*0.09236339709340602) %>% #光ロガー91の補正係数
  mutate(datetime = floor_date(datetime,unit = "days")) %>% 
  group_by(datetime) %>%
  summarise(daily_ppfd = sum(calib_value) * 60 * 10 / 1000000,
            ppfd.max = max(calib_value))

## water temp read and edit----
water_temp_file_2019=
  tibble(water_temp_file_2019) %>%
  mutate(data = map(water_temp_file_2019, read_onset)) %>%
  unnest(data) %>%
  mutate(datetime = floor_date(datetime,unit = "days")) %>%
  group_by(datetime) %>% 
  summarise(temp.max = max(temperature),
            temp.min = min(temperature),
            temp.ave = mean(temperature)) %>% 
  filter(datetime >= ymd_h("2019-07-31 0") & datetime < ymd_h("2019-09-05 0")) %>% 
  # mutate(temp.ave = na.approx(temp.ave))  %>% 
  mutate(across(contains("temp"), na.approx)) %>% 
  select(datetime,temp.max,temp.min,temp.ave) %>% print(n = 100)

## tide data from predict----

tidedata2019_date = 
  tidedata2019 %>%#env_sheetに合わせて編集 
  mutate(datetime = floor_date(datetime,unit = "day")) %>% 
  group_by(datetime) %>% 
  summarize(mean_depth = mean(depth),
            max_depth = max(depth),
            min_depth = min(depth))

## all environment data----
do = pi/180

all_env_data_2019=
  light_data_2019%>%
  left_join(tidedata2019_date, by = c("datetime")) %>%
  left_join(water_temp_file_2019 ,by = c("datetime")) %>%
  left_join(weather_data_2019) %>% 
  select(datetime,
         daily_ppfd,
         ppfd.max,
         mean_depth,
         max_depth,
         min_depth,
         temp.ave,
         temp.max,
         temp.min,
         temperature_air,
         wind,
         gust,
         wind_direction) %>% 
  mutate(wind_direction = factor(wind_direction,
                                 levels = c("N", "NNE","NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW", "静穏"),
                                 labels = c(seq(0, 360-(45/2), by = 45/2), NA))) %>% 
  mutate(wind_direction = as.numeric(as.character(wind_direction))) %>% 
  mutate(wind_direction = cos(do * wind_direction)) %>%
  mutate(wind_direction = round(wind_direction, digits = 3)) 

##egg data edit ----

tagname = Egg_2019 %>% slice(1)
Egg_2019 = Egg_2019 %>% slice(-1)
tagname = tagname %>% select(starts_with("X")) %>% as.matrix() %>% as.vector()
Egg_2019 = Egg_2019 %>% select(-`indivisual number`) 
colnames(Egg_2019) = c("ymd", tagname, "remarks")

Egg_2019 = Egg_2019 %>% mutate(across(-c(ymd, remarks), as.character))

release_rate_2019 = 
  Egg_2019 %>% #確認個体あたりの放出個体の割合バージョン 
  pivot_longer(-c(ymd, remarks), values_to = "release") %>% 
  select(-remarks, -name) %>% 
  rename(datetime = ymd) %>% 
  mutate(release = factor(release,
                          levels = c("Y","N"),
                          labels = c("Y","N"))) %>% 
  group_by(datetime) %>% 
  nest()%>% 
  mutate(frequency_table = map(data,function(x){
    table(x) %>% as_tibble() %>% 
      rename(release = x,frequency = n) 
  })) %>% 
  select(datetime,frequency_table) %>% 
  unnest(frequency_table) %>% 
  spread(key = release,value = frequency) %>% 
  mutate(check_indi = sum(N,Y)) %>% 
  mutate(Y_rate = (Y/check_indi)*100) %>% 
  mutate(Y_rate = ifelse(is.na(Y_rate),0,Y_rate))

release_rate_2019 =
  release_rate_2019 %>%
  ungroup() %>% 
  mutate(datetime = parse_date_time(datetime,"ymd!")) %>% 
  left_join(all_env_data_2019,by = "datetime")

release_rate_2019 =# 割合でみたほうが見た個体数が違っても比較できるので、こちらを2020とがっちゃんこ
  release_rate_2019 %>% 
  select(-c(N,Y,check_indi)) %>% 
  rename(rate = Y_rate)

data2019 = 
  release_rate_2019 %>% 
  rename(Date = datetime,
         ppfd.add = daily_ppfd) %>% 
  mutate(dif_depth = max_depth - min_depth,
         rate = rate / 100) %>% 
  mutate(Date = ymd(Date)) %>% 
  select(-mean_depth)

# 2020----
## read data----
indivi_info_2020 = read_csv(file = "~/Lab_Data/matsuday/egg_release/egg_data/individual2020/indivi_info2020.csv")
release_data_2020 = read_csv(file = "~/Lab_Data/matsuday/egg_release/egg_data/release_day/egg_release2020.csv")
autosampler_data_2020 = read_csv(file = "~/Lab_Data/matsuday/egg_release/egg_data/autosampler/autosampler2020.csv")
term_loggers_2020 = readxl::read_xlsx("~/Lab_Data/kawatea/調査期間.xlsx")
tidedata_2020 = read_csv("/home/matsuday/Egg_release/tide_data/arikawa2020.csv")

## term info ----
# ロガーと個体の観測期間＋台風により一部ロガーデータが使えない期間がある為
term_loggers_2020=
  term_loggers_2020 %>% 
  filter(location=="arikawagaramo") %>% 
  select(start_date,end_date)

term_loggers_2020 =
  term_loggers_2020 %>% 
  mutate(Date = map2(start_date,end_date,function(x,y){
    seq(x,y,by = "1 day") %>% 
      as.Date()
  })) %>% 
  unnest(Date) %>% 
  select(Date) 

## individual size info ----
#　個体情報の確認
indivi_info_2020 =
  indivi_info_2020 %>% 
  rename(tags = 個体番号,height = 高さ,width = 幅) %>% 
  summarise(ave_height = mean(height),
            sd_height = sd(height),
            ave_width = mean(width),
            sd_width = sd(width))

## release date----


release_data_2020=
  release_data_2020 %>% 
  rename_with(~c(1:12), -c("date", "sun", "checktime")) %>% 
  select(!c("5","6","10","12")) %>%
  unite(col = datetime, c(date, checktime), sep = " ") %>% 
  mutate(datetime = ymd_hms(datetime)) %>% 
  drop_na(datetime) %>% 
  pivot_longer(cols = 3:10, names_to = "ID", values_to = "value") %>% 
  mutate(value = recode(value, Y = 1, N = 0)) %>%
  group_nest(datetime) %>% 
  mutate(sum = map_dbl(data, function(x){
    x %>% summarise(S = sum(value, na.rm = T)) %>% pull(S)
  })) %>%
  mutate(number = map_dbl(data, function(x){
    x %>% drop_na() %>% nrow()
  })) %>%
  mutate(rate = sum/number)

release_data_2020 =#ロガーと合わせられるように、10分毎のデータに
  release_data_2020 %>%
  mutate(datetime = round_date(datetime,unit = "10 mins"))
## weather data edit ----
#　ここからロガーデータ、気象庁データをまとめていく
#　基本的に10分毎

weather_data_arikawa = read_csv("/home/matsuday/Egg_release/tide_data/arikawa_wheather_202007-09.csv")
weather_data_fukue   = read_csv("/home/matsuday/Egg_release/tide_data/fukue_wheather_202007-09.csv")

weather_data = weather_data_arikawa %>% 
  left_join(weather_data_fukue) %>% 
  select(!contains("gust"))#完成形

## depth data----


tidedata2020_date = tidedata_2020 %>% 
  mutate(Date = as.Date(datetime))

tidedata2020_date = tidedata2020_date %>% 
  group_by(Date) %>% 
  summarise(depth_max = max(depth),
            depth_min = min(depth)) %>% 
  mutate(diff_depth = depth_max - depth_min)


tidedata2020_10min = tidedata_2020 %>% select(datetime,depth) 

## logger temp----

fnames = dir(path = "~/Lab_Data/kawatea/Oxygen/",pattern = "arikawagaramo_0m_200[78]",full.names = TRUE)
temp_data_2020 = tibble(fnames = fnames) %>%
  mutate(data = map(fnames, function(x){
    read_onset(x)
  })) %>% 
  unnest(data) 

## light logger ----
light_data =  dir("~/Lab_Data/kawatea/Light/",pattern = "arikawagaramo_0m_200[78]",full.names = TRUE)

lightdata_10min_2020 = tibble(light_data) %>% 
  mutate(data = map(light_data,read_odyssey)) %>%#warning について先生と相談（何の警報なのか）
  select(data) %>%
  unnest(data) %>% 
  mutate(ppfd = ppfd*0.150)

lightdata_date_2020= lightdata_10min_2020 %>% 
  mutate(Date = as.Date(datetime)) %>%
  group_by(Date) %>% 
  summarise(ppfd.add = sum(ppfd)*60*10/1000000,
            ppfd.max = max(ppfd))

## all env data ----

env_data_10min_2020 = # 10分毎の環境データ
  tidedata2020_10min %>% 
  left_join(lightdata_10min_2020) %>% 
  left_join(temp_data_2020) %>% 
  left_join(weather_data)

env_data_chectime_2020 =# 卵放出データと合わせた。（10分毎）
  env_data_10min_2020 %>% 
  full_join(release_data_2020,by = "datetime")


# 一日毎の卵放出と環境データ ----

day_egg_env_2020 =env_data_chectime_2020 %>% 
  mutate(Date = as.Date(datetime)) %>% 
  group_by(Date) %>% 
  summarise(max_depth = max(depth),
            min_depth = min(depth),
            dif_depth = max_depth-min_depth,
            max_ppfd = max(ppfd),
            ppfd.add = sum(ppfd) * 60 * 10 / 1000000,
            day_eggrate = max(rate,na.rm = T),
            max_temp = max(temperature),
            min_temp = min(temperature),
            ave_temp = mean(temperature)) 

test = day_egg_env_2020 %>% drop_na()
final_data_2020 = test %>% filter(is.finite(day_eggrate))

data2020 = final_data_2020
moonage2020 = c(11.4,12.4,13.4,14.4,15.4,16.4,17.4,18.4,19.4,20.4,21.4,22.4,23.4,24.4,25.4,26.4,27.4,28.4,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)

data2020 = data2020 %>%  mutate(moonage = moonage2020, .after = Date) %>% 
  rename(rate = day_eggrate)

data2020 = data2020 %>% 
  rename(ppfd.max = max_ppfd,
         temp.max = max_temp,
         temp.min = min_temp,
         temp.ave = ave_temp)


# total data----

data_2y = bind_rows(data2019,data2020)

data_2y =  data_2y %>% 
  mutate(frac =as.character(MASS::fractions(rate))) %>% 
  mutate(numerator = str_extract(frac, "^[0-9]+")) %>%
  mutate(denominator = str_extract(frac, "/[0-9]+")) %>%
  mutate(denominator = str_remove(denominator, "/")) %>% 
  mutate(denominator = ifelse(is.na(denominator), 1, denominator)) %>% 
  mutate(numerator = as.numeric(numerator),
         denominator = as.numeric(denominator))

data_2y =
  data_2y %>% 
  select(-c(wind,wind_direction,moonage,gust,temperature_air))

data_2y_zscore =
  data_2y %>% print(n = 80) %>% 
  mutate(dif_temp = temp.max - temp.min) %>% 
  mutate(dif_temp = scale(dif_temp)[,1]) %>% 
  mutate(temp.ave = scale(temp.ave)[,1],
         temp.max = scale(temp.max)[,1],
         temp.min = scale(temp.min)[,1],
         ppfd.add = scale(ppfd.add)[,1],
         ppfd.max = scale(ppfd.max)[,1],
         max_depth = scale(max_depth)[,1],
         min_depth = scale(max_depth)[,1],
         dif_depth = scale(dif_depth)[,1],
         year = year(Date)) %>% 
  mutate(year = as.factor(year))

write_rds(data_2y, file = "matsuda_natural_experiment.rds")
################################################################################

library(FactoMineR)
library(factoextra)


pca_result =
  data_2y_zscore %>% 
  select(temp.ave, dif_depth, ppfd.add) %>% 
  PCA()


pca_result$var$coord %>% 
  as.tibble(rownames = "variable") %>%
  mutate(variable = factor(variable,levels = c("temp.ave","dif_depth","ppfd.add"))) %>% 
  pivot_longer(-variable,names_to = "Dim",values_to = "coord") %>% 
  ggplot()+
  geom_bar(aes(x = Dim,y = coord,fill =variable),stat = "identity",position = "dodge")


## stanglm -----
library(rstanarm)
library(tidybayes)

pca_score = 
  pca_result$ind$coord %>% 
  as.tibble()


data_2y_zscore =
  data_2y_zscore %>% 
  bind_cols(pca_score)

# pca_testmodel1 = glm(cbind(numerator, denominator - numerator) ~ Dim.1 + Dim.2 + Dim.3,family = quasibinomial,data = data_2y_zscore)

pca_testmodel = stan_glm(cbind(numerator, denominator - numerator) ~ Dim.1 + Dim.2 + Dim.3,
                         family = binomial("logit"),data = data_2y_zscore,
                         cores = 4, chains = 4, seed = 2020)

pca_result %>% summary()
pca_testmodel %>% summary()

pca_result$var$contrib



# model 診断
# y = pca_testmodel$data %>% pull(rate)
# yrep = posterior_predict(pca_testmodel, nsamples = 5)
# ppc_dens_overlay(y,yrep)

summary_table <- function(result){
  
  tbl <- as.data.frame(summary(result)[[12]])
  tbl <- tbl[2:nrow(tbl), ]
  colnames(tbl) <- c("Estimate", "SE", "t-value", "p-value")
  write.csv(tbl, "result.csv")
  
  return(summary(result))
  
}
#############################################################################
library(brms)



summary(pca_testmodel)
xlabel = "Model coefficients"
ylabel = "Coefficient value"
title = "Parameter's coefficient value"
cfs = posterior_summary(pca_testmodel, probs = c(0.025, 0.10, 0.90, 0.975))
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
                     limits = c(-1,1),
                     breaks = seq(-1,1, by = 0.5)) +
  ggpubr::theme_pubr()

# ggsave("~/Egg_release/plots_final/parameter.png",width = WIDTH,height = HEIGHT,units = "mm",dpi = DPI)

#環境要因毎の影響

intercept = tidy_draws(pca_testmodel) %>% select(contains("Int"))


cfs = tidy_draws(pca_testmodel) %>% 
  select(contains("Dim")) %>% 
  rename_with(~c("PCA1", "PCA2", "PCA3"),
              everything())

S = 
  pca_result$var$cor %>% 
  as_tibble(rownames = "variable") %>% 
  rename_with(~c("PCA1", "PCA2", "PCA3"), -variable) %>% 
  pivot_longer(-variable, names_to = "cfs") %>% 
  mutate(signs = sign(value)) %>% 
  rename(value2 = value)

PCA

Z = pca_result$var$contrib %>% 
  as_tibble(rownames = "variable") %>% 
  rename_with(~c("PCA1", "PCA2", "PCA3"), -variable) %>% 
  pivot_longer(-variable, names_to = "cfs") 


Y = cfs %>% gather(key = "cfs", value = "coef")

SZ = full_join(Z, S, by = c("cfs", "variable"))

YZ = full_join(Y,SZ, by = c("cfs"))


YZ = YZ %>% mutate(coef = signs * coef * value / 100) %>% select(cfs, variable, coef)

# YZ = YZ %>% mutate(coef = signs * coef * value / 100 *(2 / sqrt(2))) %>% select(cfs, variable, coef)

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

YZc =YZc %>% 
  mutate(variable = recode(variable,
                           dif_depth="ΔDepth",
                           ppfd.add="Daily PPFD",
                           temp.ave = "Mean Temperature"))


xlabel = "Environmental parameter"
ylabel = "Parameter value"
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
  scale_y_continuous(ylabel, limits = c(-1, 1),
                     breaks = seq(-1, 1, by = 0.5))+
  ggpubr::theme_pubr(FONTSIZE +5) 


# model の予測----
data_2y_zscore =
  data_2y_zscore %>%
  mutate(tide = ifelse(str_detect(Date,"2019-08-0[12]"),1,0),.after = rate) %>%
  mutate(tide = ifelse(str_detect(Date,"2019-08-1[4-7]"),1,tide)) %>%
  mutate(tide = ifelse(str_detect(Date,"2019-08-29"),1,tide)) %>%
  mutate(tide = ifelse(str_detect(Date,"2019-08-3[01]"),1,tide)) %>%
  mutate(tide = ifelse(str_detect(Date,"2020-08-0[3-6]"),1,tide)) %>%
  mutate(tide = ifelse(str_detect(Date,"2020-08-1[89]"),1,tide)) %>%
  mutate(tide = ifelse(str_detect(Date,"2020-08-20"),1,tide)) %>%
  mutate(tide = ifelse(str_detect(Date,"2020-09-0[1-4]"),1,tide)) %>%
  mutate(tide = factor(tide,labels = c("flood tide","others"),levels = c(1,0)))


# model vs raw----

predict_data_d1=
  data_2y_zscore %>% 
  ungroup() %>% 
  tidyr::expand(temp.ave = seq(min(temp.ave),
                               max(temp.ave), length = 11),
                dif_depth = c(min(dif_depth), median(dif_depth), max(dif_depth)),
                ppfd.add = c(min(ppfd.add), max(ppfd.add)),
                denominator = 6,
                numerator = 0)

predicted_pca_d1 =
  predict_data_d1 %>% 
  select(temp.ave, dif_depth, ppfd.add) %>% 
  predict(pca_result, newdata = .)

contrib = 
  pca_result$var$contrib %>% 
  as_tibble(rownames = "variable") %>% 
  mutate(across(contains("Dim"), ~ ./100)) %>% 
  rename(d1 = Dim.1, d2 = Dim.2, d3 = Dim.3)

pout_dim1 = 
  predicted_pca_d1$coord %>% as_tibble() %>% 
  bind_cols(predict_data_d1) %>% 
  add_fitted_draws(model = pca_testmodel) %>% 
  ungroup() %>% 
  group_by(dif_depth, temp.ave, ppfd.add) %>% 
  # mutate(.value = exp(.value)) %>% 
  mean_hdci(.value) %>% 
  rename(fit.response = .value,
         l95.response = .lower,
         u95.response = .upper)

scaling_coef = data_2y %>% summarise(across(c(temp.ave, dif_depth, ppfd.add), list(m = mean,s = sd)))

xlabel = "平均水温 (°C)"
ylabel = "ノコギリモクの放卵率 (%)"

pout_dim1 %>% ungroup() %>% 
  mutate(temp.ave = (temp.ave * scaling_coef$temp.ave_s) + scaling_coef$temp.ave_m) %>% 
  mutate(dif_depth = factor(dif_depth, labels = c("最小", "中央値", "最大"))) %>%
  mutate(ppfd.add = factor(ppfd.add,   labels = c("最低", 　　　　　"最高"))) %>% 
  ggplot() + 
  geom_point(aes(x = temp.ave,
                 y = 100*rate),
             size = 2, alpha = 0.5,
             data = data_2y) +
  geom_line(aes(x = temp.ave,
                y = 100*fit.response,
                color = dif_depth,
                linetype = ppfd.add),
            size = 2)+
  scale_x_continuous(xlabel, limits = c(24, 31)) +
  scale_y_continuous(ylabel) +
  scale_color_manual(values = viridis::viridis(4)) +
  guides(color = guide_legend(title = "水深差", 
                              override.aes = list(size = 1),
                               label.position="left",
                              label.hjust = 0),
         linetype = guide_legend(title = "光量子量", 
                                 override.aes = list(size = 1),
                              label.position="left",
                              label.hjust = 0)) +
  ggpubr::theme_pubr() +
  theme(legend.position=c(1, 1),
        legend.justification=c(1, 1),
        legend.background=element_blank(),
        legend.key.width=unit(10, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.spacing.y = unit(1, "mm"))




wh = aseries(5);wh
DPI = 300
ggsave("kyogikai/eggs.png", width = wh[1], height = wh[2], dpi = DPI, units = "mm")

