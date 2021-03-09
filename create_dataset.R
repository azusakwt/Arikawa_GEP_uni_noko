# ファイルが多いので、ここでデータをまとめる
# Greg
# 2021 Jan 14
# Calculate GEP
# Greg Nishihara
# 2021/01/14
#

library(tidyverse)
library(lubridate)
library(gnnlab)
library(readxl)
library(furrr)
plan(multisession, workers = 20) # furrr に必要

Sys.setlocale("LC_TIME", "en_US.UTF-8") # アメリカ英語に設定

KAWATE = "~/Lab_Data/kawatea/"

#データの読み込み---------------------
#調査期間
kikan = str_glue("{KAWATE}/調査期間.xlsx") %>% read_xlsx()
kikan = kikan %>% mutate(Date = map2(start_date, end_date, function(x, y){
  seq(x, y, by = "1 day") %>% as_date()
})) %>% unnest(Date) %>% select(-c(comment, remarks))

process_cem = function(x) {
  dset = read_alec(x)
  dset %>%
    mutate(datetime = floor_date(datetime, "10 minutes")) %>%
    group_by(datetime) %>%
    summarise(across(c(speed, dir, ew, ns, temperature),
                     mean))
}

cem = tibble(fnames = dir(str_glue("{KAWATE}/CEM"), full = TRUE, pattern = "arik")) %>%
  mutate(data = future_map(fnames, process_cem)) %>% unnest(data)

cem = cem %>% group_by(datetime) %>% 
  summarise(across(c(speed, dir, ew, ns, temperature),
                   mean))

cku = tibble(fnames = dir(str_glue("{KAWATE}/CKU"), full = TRUE, pattern = "arik")) %>%
  mutate(data = future_map(fnames, read_alec)) %>% unnest(data)

cku = cku %>% group_by(datetime) %>% 
  summarise(across(c(chla, turbidity, temperature),
                   mean))
# Don't use the raw data from the depth logger. Use the calculated depths
# depth = tibble(fnames = dir(str_glue("{KAWATE}/Depth"), full = TRUE)) %>% mutate(data = future_map(fnames, read_onset)) %>% unnest(data)

depth = read_csv("tidaldata.csv") %>% select(datetime, depth) # From predict_tides.R

tmp = read_csv("~/Lab_Data/weather/201701_202012_Arikawa_JMA_Data.csv")
weather = read_csv("~/Lab_Data/weather/202101_202102_Arikawa_JMA_Data.csv")
weather = bind_rows(tmp, weather) %>% distinct()

microstation = tibble(fnames = dir(str_glue("{KAWATE}/Microstation"), full = TRUE, pattern = "arik")) %>%
  mutate(data = future_map(fnames, read_onset)) %>% unnest(data) %>% 
  mutate(datetime = floor_date(datetime, "minutes")) %>% 
  select(-fnames) %>% distinct()

oxygen = tibble(fnames = dir(str_glue("{KAWATE}/Oxygen/"), full = TRUE), 
                type = str_extract(fnames, "xlxs|csv")) %>%
  mutate(data = future_map(fnames, read_onset))

oxygen = oxygen %>% mutate(fnames = basename(fnames)) %>%
  separate(fnames, into = c("type", "id", "location", "position", "date")) %>%
  filter(!str_detect(position, "calibra"))

oxygen = oxygen %>%
  unnest(data) %>%
  select(location, position, datetime, oxygen = mgl, temperature) %>%
  mutate(Date = as_date(datetime))

oxy_all = left_join(kikan, oxygen, by = c("location", "Date")) %>%
  drop_na(oxygen) %>%
  select(location, datetime, position, oxygen,　temperature, Date)

light = tibble(fnames =dir(str_glue("{KAWATE}/Light/"), full = TRUE), 
               type = str_extract(fnames, "xlxs|csv")) %>% 
  mutate(data = map(fnames, read_odyssey))

light = light %>% mutate(fnames = basename(fnames)) %>% 
  separate(fnames, into = c("type", "id", "location", "position", "date"))

light = light %>% unnest(data) %>% 
  select(location,id, position, datetime, ppfd) %>% 
  mutate(Date = as.Date(datetime))

light = left_join(kikan, light, by = c("location", "Date")) %>% 
  drop_na(ppfd) %>% 
  select(location,id, position, Date, datetime, ppfd)

light_0 = light %>% filter(position == "0m") %>% 
  select(-position)

light_1 = light_0 %>% 
  distinct(Date,location,datetime,.keep_all = TRUE) %>% 
  mutate(light_group = ifelse(ppfd>0,TRUE,FALSE)) %>% 
  select(location, datetime, ppfd.water = ppfd, light_group)
# Join all of the data

tmp = full_join(
  cem %>% select(datetime, speed, dir, ew, ns, temperature.cem = temperature),
  cku %>% select(datetime, chla, turbidity, temperature.cku = temperature),
  by = "datetime")
tmp = full_join(tmp, depth, by = "datetime")

tmp = full_join(
  tmp,
  microstation %>% select(datetime, wind, gust, ppfd.microstation = ppfd),
  by = "datetime"
)

envdata = full_join(
  tmp, 
  weather %>% select(datetime, rain, temperature.jma = temperature_air, 
                     wind.jma = wind, gust.jma = gust),
  by = "datetime"
)

oxydata = oxy_all %>% 
  filter(str_detect(location, "arik")) %>%
  filter(str_detect(position, "surface|0m|1m")) %>% 
  select(-Date) %>% 
  group_nest(location) %>% 
  mutate(data = map(data, function(X) {
    X %>% distinct() %>% 
      pivot_wider(id_cols = c(datetime),
                      names_from = position, 
                      values_from = c(oxygen, temperature)) %>% print()
  })) %>% unnest(data)

alldata = left_join(oxydata, envdata, by = "datetime")
alldata = left_join(alldata, light_1, by = c("datetime", "location"))

# alldata %>% select(datetime, matches("ppfd"))
################################################################################  
# Use PCA to impute missing values.
################################################################################

library(FactoMineR)
library(missMDA)

X = alldata %>% select(-c(location, datetime, starts_with("oxygen")))

# Xpca = PCA(X, graph = FALSE)
# Xpca %>% summary()
Xadj = imputePCA(X, ncp = 10) # takes some time to run. 最初の 10 Dim.x を使う

alldata_adj = bind_cols(
  alldata %>% select(c(location, datetime, starts_with("oxygen"))),
  Xadj$completeObs %>% as_tibble()
) %>% arrange(location, datetime)

alldata_adj %>% 
  group_nest(location,
             date = as_date(datetime)) %>% 
  mutate(hasna = map_dbl(data, function(X) {
    X %>% drop_na(oxygen_0m, oxygen_surface) %>% nrow()
  })) %>% 
  filter(near(hasna, 0)) %>% print() %>% 
  slice(1) %>% unnest()

# Keep only full data
alldata_adj = alldata_adj %>% 
  group_nest(location,
             date = as_date(datetime)) %>% 
  mutate(hasna = map_dbl(data, function(X) {
    X %>% drop_na(starts_with("oxygen")) %>% nrow()
  })) %>% 
  filter(near(hasna, 144)) %>% 
  select(location, date, data) %>% unnest(data)

###############################
alldata_adj = alldata_adj %>% mutate(month = month(datetime),
                         year = year(datetime)) %>%
  mutate(H = hour(datetime) + minute(datetime)/60)


alldata_adj %>% 
  write_csv(file = "fulldataset.csv")

