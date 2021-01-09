# GEP
# Azusa Kawate, Greg
# 2021/01/09, 01/10

#####################################
#package-------------
library(tidyverse)
library(lubridate)
library(gnnlab)
library(readxl)


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
  geom_boxplot(aes(x = month, y = GEP, group = month, fill = location),
               outlier.colour = NA)+
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
