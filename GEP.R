# GEP
# Azusa Kawate
# 2021/01/09

#####################################
#package-------------
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
caluculate_rate=function(data){
  out=mgcv::gam(oxygen~s(H,bs="cs"),data=data)
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
  mutate(rate=map(data,caluculate_rate)) %>%
  select(-data) %>%
  unnest() %>%
  select(Date, location, position, rate, datetime, oxygen,temperature)

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

###############################################################
# make figure---------------
oxygen_rate = read_csv("oxygen_rate.csv")
light_all = read_csv("light.csv")
oxy = oxygen_rate %>% mutate(month = month(Date),
                             year = year(Date)) %>% 
  mutate(H = hour(datetime) + minute(datetime)/60)

oxy %>% filter(year == 2019, month == 4,
               location == "arikawagaramo") %>%
  ggplot()+
  geom_point(aes(x = H, y = oxygen, color = position), size = 0.5)+
  facet_wrap(Date ~ .)

##################################################################
light_0 = light_all %>% filter(position == "0m") %>% 
  select(-position)

light_1 = light_0 %>% 
  mutate(time= hour(datetime)) %>% 
  filter(time == 15) %>% 
  mutate(plus = ifelse(ppfd == 0,9,0)) %>% filter(plus ==9) %>% 
  mutate(datetime = datetime+hours(plus)) 

light_1 = light_0 %>% 
  drop_na(ppfd) %>% 
  distinct(Date,location,datetime,.keep_all = TRUE) %>% 
  mutate(light_group = ifelse(ppfd>0,TRUE,FALSE))

# NEP の計算 ----
product = inner_join(oxygen_rate,light_1, by = c("location", "Date", "datetime"))

product_all = 
  product %>% 
  filter(!str_detect(position, "surface")) %>% 
  group_by(Date, location, position) %>% 
  summarise(NEP = sum(rate * 1),
            RP = sum(rate * (!light_group)),
            day_hours = sum(light_group)/6,
            night_hours = sum(!light_group)/6) %>% 
  mutate(RP = RP / night_hours * 24) %>% 
  mutate(GEP = NEP - RP, 
         .after= RP) %>% 
  filter(GEP > 0) %>% 
  group_by(Date, location) %>% 
  summarise(across(c(NEP, RP, GEP), mean))

##################################################################
product_all = product_all %>% 
  mutate(year = year(Date),
         month = month(Date),
         Time = floor_date(Date, "month")) %>% 
  ungroup()

#########################################################################
ylabel = expression("GEP"~(mg~O[2]~m^{-3}~day^{-1}))
product_all %>%filter(location == "arikawagaramo") %>% 
  # mutate(RP = -RP) %>% 
  # filter(RP>0) %>%
  ggplot()+
  geom_boxplot(aes(x = month, y = -RP, group = month, fill = location),
               outlier.colour = NA)+
  geom_smooth(aes(x = month, y = -RP), 
              method = "gam", formula = y~ s(x, k = 7, bs = "cs"))+
  scale_x_continuous("", 
                     minor_breaks = 1:12,
                     breaks = c(1, 6, 12),
                     labels = month.abb[c(1, 6, 12)]) +
  scale_y_continuous(ylabel)+
  facet_wrap("year", nrow=1, switch = "x")+
  theme_classic(base_size = 16)+
  theme(axis.title.x = element_text(size = rel(0.8)),
        axis.title.y = element_text(size = rel(0.8), angle = 90),
        axis.line = element_line(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.placement = "outside") 

ggsave("figure/GEP_4year_garamo.png",
       width = 200, height = 115, units = "mm", dpi = 600)

product_all %>% 
  ggplot()+
  geom_boxplot(aes(x = month, y = GEP, group = month, fill = "location"),
               outlier.colour = NA)+
  geom_smooth(aes(x = month, y = GEP),
              method = "gam",
              formula = y~s(x,k = 10, bs = "cs"))+
  scale_x_continuous("", 
                     minor_breaks = 1:12,
                     breaks = c(1,5, 8, 12),
                     labels = month.abb[c(1, 5, 8, 12)]) +
  scale_y_continuous(ylabel)+
  theme_classic(base_size = 16)+
  theme(axis.title.x = element_text(size = rel(0.8)),
        axis.title.y = element_text(size = rel(0.8), angle = 90),
        axis.line = element_line(),
        legend.position = "none") 

WIDTH = 297/2
HEIGHT = 210/2
ggsave("figure/GEP_season_garamo.png", 
       width = WIDTH, height = HEIGHT, units = "mm", dpi = 600)
