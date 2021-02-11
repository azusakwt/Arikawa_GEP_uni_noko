# Arikawa Bay, Nagasaki Prefecture, Japan
# Greg Nishihara
# 2020 Jan 24
# From https://github.com/gnishihara/predict_tides

# Required packages --------------------------------------------------
library(tidyverse)
library(lubridate)
library(oce)
library(googledrive)

# Load data ----------------------------------------------------------
datalink = "https://docs.google.com/spreadsheets/d/1eefPkgzwWX4-oNsIqrQ-qzEdgnpGICYwL-F2DRxZDck/edit?usp=sharing"
drive_download(datalink, path = "datafile.csv", type = "csv")
arikawa = read_csv("datafile.csv")
latitude = 32.988152   # location of depth logger

# Estimate the tidal harmonics and calcuate the tidal levels ---------
arik = as.sealevel(elevation = arikawa$value, time = arikawa$datetime)
arikawa_tides = tidem(arik, latitude = latitude)

# Predict tidal level ------------------------------------------------
start = ymd_h("2017-01-01 0")
end   = ymd_h("2021-01-01 0")

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
