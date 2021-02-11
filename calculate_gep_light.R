# Calculate GEP
# Greg Nishihara
# 2021/01/10
#
# GEPの算出はここに移した

library(tidyverse)
library(lubridate)
library(gnnlab)
library(readxl)

#データの読み込み---------------------
oxygen = tibble(fnames = dir("~/Lab_Data/kawatea/Oxygen/", full = TRUE),
                type = str_extract(fnames, "xlxs|csv")) %>% 
  mutate(data = map(fnames, read_onset))

oxygen %>%tail() 

oxygen = oxygen %>% mutate(fnames = basename(fnames)) %>% 
  separate(fnames, into = c("type", "id", "location", "position", "date")) %>% 
  filter(!str_detect(position, "calibra"))

#調査期間
kikan = read_xlsx("~/Lab_Data/kawatea/調査期間.xlsx")

kikan = kikan %>% mutate(Date = map2(start_date, end_date, function(x, y){
  seq(x, y, by = "1 day") %>% as.Date
})) %>% unnest(Date)

oxygen = oxygen %>% unnest(data) %>% 
  select(location, position, datetime, oxygen = mgl, temperature) %>% 
  mutate(Date = as.Date(datetime))
oxygen %>% filter(location == "arikawaamamo") %>% arrange(Date) %>% tail()

oxy_all = left_join(kikan, oxygen, by = c("location", "Date")) %>% 
  drop_na(oxygen) %>% 
  select(location, datetime, position, oxygen,temperature, Date)
###############################
oxy = oxy_all %>% mutate(month = month(Date),
                         year = year(Date)) %>% 
  mutate(H = hour(datetime) + minute(datetime)/60)


#一次生産量の計算----------------------------
calculate_rate=function(data){
  out=mgcv::gam(oxygen~s(H, k = 20),data=data)
  y1=predict(out)
  eps=1e-6
  y2=predict(out,newdata=data_frame(H=data$H+eps))
  data %>% mutate(rate=(y2-y1)/eps)
}

#データが足りない日は削除
#１日あたり144
oxy_all = oxy_all %>%
  distinct(location, position, datetime, .keep_all = TRUE) %>% 
  group_nest(location,position, Date) %>%
  mutate(n = map_dbl(data, function(X) {
    X %>% nrow()
  })) %>%
  filter(near(n, 144)) %>% unnest(data) %>% 
  select(-n)

# 作成した関数を立てはめる-----
oxygen_rate = oxy_all %>%
  filter(str_detect(position, "0m|1m|2m|surface")) %>%
  mutate(oxygen = ifelse(oxygen <0,NA,oxygen)) %>%
  drop_na(oxygen) %>%
  mutate(H = hour(datetime)+minute(datetime)/60) %>%
  group_by(Date, location, position) %>%
  nest() %>% ungroup() %>%
  mutate(rate=map(data,calculate_rate)) %>%
  select(-data) %>%
  unnest() %>%
  select(Date, location, position, rate, datetime, oxygen,temperature)

kyphosus = 15
siganus = 19

oxygen_rate %>% 
  filter(str_detect(location, "garamo")) %>% 
  ungroup() %>% 
  group_by(Date, datetime) %>% 
  summarise(temperature = mean(temperature)) %>% 
  group_by(Date) %>% 
  summarise_at(vars(temperature),
            list(daily_sd = sd, daily_mean = mean, daily_min = min, daily_max = max, n = ~length(.)), na.rm = TRUE) %>% 
  write_csv("temperature.csv")


oxygen_rate %>% 
  filter(str_detect(location, "garamo")) %>% 
  filter(!str_detect(position, "surface")) %>% 
  ungroup() %>% 
  group_by(Date, datetime) %>% 
  summarise(temperature = mean(temperature)) %>% 
  group_by(Date) %>% 
  summarise_at(vars(temperature),
               list(daily_sd = sd, 
                    daily_mean = mean, 
                    daily_min = min, daily_max = max, n = ~length(.)), na.rm = TRUE) %>% 
  write_csv("temperature_no-surface.csv")


write_csv(oxygen_rate,"oxygen_rate.csv")

######################################################################
light = tibble(fnames = dir("~/Lab_Data/kawatea/Light/", full = TRUE),
               type = str_extract(fnames, "xlxs|csv")) %>% 
  mutate(data = map(fnames, read_odyssey))

light = light %>% mutate(fnames = basename(fnames)) %>% 
  separate(fnames, into = c("type", "id", "location", "position", "date"))

light = light %>% unnest(data) %>% 
  select(location, position, datetime, ppfd) %>% 
  mutate(Date = as.Date(datetime))

light = left_join(kikan, light, by = c("location", "Date")) %>% 
  drop_na(ppfd) %>% 
  select(location, position, Date, datetime, ppfd)

light_all = light %>% mutate(year = year(Date), month = month(Date),
                             H = hour(datetime)+ minute(datetime)/60)

write_csv(light_all, "light.csv")
