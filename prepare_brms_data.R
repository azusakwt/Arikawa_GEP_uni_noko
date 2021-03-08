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


##################################################################

masstransfer = function(windspeed, temperature, salinity, oxygen) {
  calc_k600 = function(windspeed) {
    # Crusius and Wanninkhof 2003 L&O 48
    U10 = 1.22 * windspeed # m/s
    0.228*U10^2.2 + 0.168 # cm/h
  }
  k600 = calc_k600(windspeed)
  SCoxygen = marelac::gas_schmidt(temperature, species = "O2")
  a = ifelse(windspeed < 3, 2/3, 1/2)
  kx = k600 * (600/SCoxygen)^a # cm / h
  o2sat = marelac::gas_O2sat(salinity, temperature, method="Weiss")
  kx/100 * (o2sat - oxygen) # g / m2 / h
}

try_calc_mt = possibly(masstransfer, otherwise = NA)

# NEP の計算 ----

calculate_rate=function(data, k = 10){
  data = data %>% mutate(H = hour(datetime) + minute(datetime) / 60)
  out=mgcv::gam(oxygen~s(H, k = k, bs = "cs"), data=data)
  y1=predict(out)
  eps=1e-6
  y2=predict(out,newdata = tibble(H=data$H+eps))
  data %>% mutate(rate=(y2-y1)/eps)
}
try_calc_rate = possibly(calculate_rate, otherwise = NA)

# alldata.rds をアップデートするなら、
# alldata.rds をさきに削除してください。
 
if(!file.exists("alldata.rds")) {
  
  alldata = read_csv("fulldataset.csv")
  
  alldata = alldata %>% 
    select(-temperature.jma) %>% 
    mutate(temperature = 0.25 * (temperature_0m + 
                                  temperature_1m +
                                  temperature.cem + 
                                  temperature.cku),
           .before = oxygen_0m) %>%  
    rename(surfacetemperature = temperature_surface) %>% 
    select(-c(temperature_0m, 
              temperature_1m,
              temperature.cem, 
              temperature.cku)) %>% 
    pivot_longer(cols = (starts_with("oxygen")),
                 values_to = "oxygen")
  
  ################################################################################
  
  alldata = alldata %>% 
    group_nest(location, name, date) %>% 
    mutate(data = future_map(data, try_calc_rate, k = 10))
  
  alldata = alldata %>% unnest(data)
  
  alldata = alldata %>% 
    mutate(name = str_remove(name, "oxygen_")) %>% 
    pivot_wider(names_from = name,
                names_glue = "{.value}_{name}",
                values_from = c(oxygen, rate)) 
  
  alldata = alldata %>% 
    mutate(mt = try_calc_mt(wind, surfacetemperature, 
                            salinity = 32, oxygen_surface), 
           .before = temperature)
  
  saveRDS(alldata, file = "alldata.rds")
} else {
  alldata = readRDS("alldata.rds")
}

se = function(x, na.rm=T) {
  # Standard Error
  n = length(x) - sum(is.na(x))
  sd(x, na.rm = na.rm) / sqrt(n - 1)
}
mad = function(x, na.rm = T) {
  # Median Absolute Deviation
  m = mean(x,na.rm = na.rm)
  xbar = x[!is.na(x)] - m
  median(abs(xbar))
}


# truncate speed, chl-a, turbidity, wind, gust, and ppfd.microstation data
# since the imputePCA() might set them to negative values.
alldata = alldata %>% 
  mutate(across(c(speed, chla, turbidity, wind, gust, ppfd.microstation, surfacetemperature),
                ~ ifelse(. < 0, 0, .)))

alldata_daily = alldata %>% 
  relocate(ppfd.microstation, .before = temperature) %>% 
  mutate(depth = ifelse(str_detect(location, "garamo"),
                        depth + 3,
                        depth - 1)) %>% 
  mutate(rate = (rate_0m + rate_1m + rate_surface)/2*depth,
         ppfd.microstation = ppfd.microstation - 1.2,
         .after = mt) %>% 
  mutate(year = year(date),
         location = factor(location)) %>% ungroup() %>% 
  group_by(location, year, date) %>% 
  summarise(
    PPFD = sum(ppfd.microstation)*60*10 / 1000000,
    dDEPTH = max(depth) - min(depth),
    MT = sum(mt)/6,
    NEP = sum(rate - mt)/6,
    RP = sum((rate - mt) * (near(ppfd.microstation, 0)))/6,
    day_hours = sum((ppfd.microstation > 0))/6,
    night_hours = sum(near(ppfd.microstation, 0))/6) %>% 
  mutate(RP = RP / night_hours * 24) %>% 
  mutate(GEP = NEP - RP, 
         .after= RP) %>% 
  ungroup() %>% 
  filter(GEP > 0) %>% 
  mutate(state = ifelse(date < as_date("2019-07-01"), "Vegetated", "Desertified"))


alldata_daily_env = alldata %>% 
  mutate(year = year(date),
         location = factor(location)) %>% ungroup() %>% 
  group_by(location, year, date) %>% 
  summarise(across(c(temperature, speed, chla, turbidity, wind, rain, surfacetemperature),
                   list(mean = mean,
                        sd = sd,
                        median = median,
                        mad = mad,
                        min = min, 
                        max = max
                   ))) %>% print()


alldata_daily = full_join(alldata_daily, alldata_daily_env, by = c("location", "year", "date")) %>% 　drop_na(GEP)




clrs = viridis::cividis(3)
alldata_daily %>% select(-c(year, date, location, state)) %>% 
  select(!matches("mad|sd|min|max|hours|mean")) %>% 
  GGally::ggcorr(low = clrs[1],
                 mid = clrs[2],
                 high = clrs[3],
                 label = T)
################################################################################
library(brms)
library(tidybayes)
library(multidplyr)

# Desertified state contains no April and May data.
dset = 
  alldata_daily %>%
  select(location:PPFD,
         GEP,
         RP,
         matches("mad"),
         matches("mean")) %>% 
  mutate (RP = -RP) %>% 
  filter(GEP > 0) %>% 
  filter(RP > 0) %>% 
  mutate(state = ifelse(date < as_date("2019-07-01"), "Vegetated", "Desertified")) %>% 
  mutate(location = ifelse(str_detect(location, "ama"),
                           "Z", "S"))

dset = 
  dset %>% ungroup() %>% 
  mutate(contiguous = as.numeric(date), .after = date) %>% 
  mutate(contiguous = c(NA, diff(contiguous))) %>% 
  mutate(contiguous = ifelse(is.na(contiguous), 1, contiguous)) %>% 
  mutate(contiguous = ifelse(near(contiguous, 1), 0, 1)) %>% 
  group_by(location, state) %>% 
  mutate(contiguous = cumsum(contiguous)) %>% 
  mutate(contiguous = ifelse(str_detect(location, "Z"),
                             contiguous + 1, 
                             contiguous))

dset = dset %>% rename_with(~str_replace(., "_", "."))
dset = dset %>% rename(TEMP = temperature.mean)
dset = dset %>% mutate(month = month(date))

dset = dset %>% mutate(month.abb = factor(month, 
                                      levels = c(6:12, 1:5),
                                      labels = month.abb[c(6:12, 1:5)])) %>% 
  mutate(month = as.numeric(month.abb)) 

dset = dset %>% mutate(state = factor(state))


dset %>% 
  relocate(TEMP, state, .before = GEP) %>% 
  ungroup() %>% 
  group_by(location) %>% 
  summarise(mean_ppfd = mean(PPFD),
         sd_ppfd = sd(PPFD),
         mean_temp = mean(TEMP),
         sd_temp = sd(TEMP))

dset %>% ungroup() %>% 
  relocate(TEMP, .before = GEP) %>% 
  group_by(location) %>% 
  mutate(mean_ppfd = mean(PPFD),
         sd_ppfd = sd(PPFD),
         mean_temp = mean(TEMP),
         sd_temp = sd(TEMP), .before = PPFD) %>% ungroup() %>% 
  mutate(PPFD = (PPFD - mean_ppfd)/ sd_ppfd, .after = PPFD) %>% 
  mutate(TEMP = (TEMP - mean_temp)/ sd_temp) %>% 
  saveRDS("prepared_brms_data.rds")
