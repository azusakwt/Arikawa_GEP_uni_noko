# Arikawa Bay, Nagasaki Prefecture, Japan
# Greg Nishihara
# 2020 Jan 24
# From https://github.com/gnishihara/predict_tides

# Required packages --------------------------------------------------
library(tidyverse)
library(lubridate)
library(gnnlab)
library(oce)
library(furrr)
plan(multisession, workers = 20)

# library(googledrive)

# Load data ----------------------------------------------------------
# datalink = "https://docs.google.com/spreadsheets/d/1eefPkgzwWX4-oNsIqrQ-qzEdgnpGICYwL-F2DRxZDck/edit?usp=sharing"
# drive_download(datalink, path = "datafile.csv", type = "csv")
# arikawa = read_csv("datafile.csv")
Sys.setlocale("LC_TIME", "en_US.UTF-8") # アメリカ英語に設定

KAWATE = "~/Lab_Data/kawatea/"

kikan = str_glue("{KAWATE}/調査期間.xlsx") %>% readxl::read_xlsx()
kikan = kikan %>% mutate(Date = map2(start_date, end_date, function(x, y){
  seq(x, y, by = "1 day") %>% as_date()
})) %>% unnest(Date) %>% select(-c(comment, remarks)) %>% 
  mutate(location = str_extract(location, "amamo|garamo"))

# Depth data 
arikawa = 
  tibble(fnames = dir("~/Lab_Data/kawatea/Depth/", pattern = "amamo|garamo", full = T)) %>% 
  mutate(location = str_extract(fnames, "amamo|garamo")) %>% 
  mutate(survey = seq_along(location)) %>% 
  mutate(data = future_map(fnames, read_onset)) %>% unnest(data) %>% 
  select(-fnames) %>% 
  mutate(datetime = floor_date(datetime, "min"))

  # Weather data (atmosphereic pressure)
# 1 hpa = 0.1 kpa
tmp = read_csv("~/Lab_Data/weather/202101_202102_Fukue_JMA_Data.csv")
weather = read_csv("~/Lab_Data/weather/201701_202012_Fukue_JMA_Data.csv")
weather = bind_rows(weather, tmp) %>% distinct()

# Combine and correct depth pressure data
gravity = 9.806614
all = full_join(arikawa, weather, by = c("datetime"))

all = all %>% mutate(value = (kpa - 0.1 * hpa)/gravity) %>% mutate(Date = as_date(datetime))
all = left_join(kikan, all, by = c("location","Date")) %>% select(location, survey, Date, datetime, value) %>% drop_na(value)

latitude = 32.988152   # location of depth logger

all = all %>% group_by(location, survey) %>% 
  mutate(value = value - mean(value)) %>% 
  mutate(value = value + 5) %>% 
  group_by(datetime) %>% 
  summarise(value = mean(value))


# Estimate the tidal harmonics and calculate the tidal levels ---------
arik = as.sealevel(elevation = all$value, time = all$datetime)
arikawa_tides = tidem(arik, latitude = latitude)

# Predict tidal level ------------------------------------------------
start = ymd_h("2017-01-01 0")
end   = ymd_h("2021-03-01 0")

datetime = seq(start - hours(1), end + hours(1), by = "10 mins")
depth = predict.tidem(arikawa_tides, datetime)

tidaldata = tibble(datetime, depth) %>%
  mutate(date = as_date(datetime),
         H = hour(datetime) + minute(datetime)/60)

## Plot of the tidal levels
datebreaks = seq(start, end, by = "1 week")
tidaldata %>% filter(between(datetime, start, end - minutes(1))) %>% 
  ggplot(aes(x = datetime, y = depth)) +
  geom_line()  +
  scale_x_datetime("",
                   date_minor_breaks="days",
                   breaks=datebreaks) +
  scale_y_continuous("Depth (m)")

# Prepare the low/high tide chart ------------------------------------
tidalchart =
  tidaldata %>% mutate(dy = depth - lag(depth)) %>%
  mutate(chk = (sign(dy) == sign(lag(dy)))) %>%
  filter(!chk) %>%
  filter(between(datetime, start, end))

tidalchart %>% select(datetime, depth) %>%
  mutate(state = ifelse(depth < lag(depth), "Low", "High")) %>%
  mutate(state = coalesce(c("Low", rep(NA, n()-1)), state)) %>%
  mutate(chk = state == lead(state))


tidaldata %>% write_csv("tidaldata.csv")
