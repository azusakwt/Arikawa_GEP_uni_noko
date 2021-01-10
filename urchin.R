# sea urchin
# Azusa Kawate
# 2021/01/09

# package--------------
library(tidyverse)
library(lubridate)
library(gnnlab)
library(readxl)

# data----------------
urchin = read_csv("~/Lab_Data/quadrat_data/garamo_quadrats/data/urchin_counts_csv.csv") %>% 
  filter(Site == "arikawa")

urchin = urchin %>% 
  mutate(date = parse_date_time(date, "%m/%d/%y", locale = "ja_JP.utf-8"))

# figure-----------

y =data_frame(year = c(2018:2020)) %>% 
  mutate(year = as.numeric(year))
m = data_frame(month = c(1:12)) %>% 
  mutate(month = as.numeric(month))
urchin_name = data_frame(urchin = c("dsetosum","hcrassispina"))

m_y = full_join(y,m, by = character()) %>% mutate(location = "Sargassum") %>% 
  full_join(urchin_name, by = character())
m_y = m_y[-1:-6,]

urchin_day = urchin %>% 
  mutate(year = year(date), month = month(date)) %>% 
  group_by(date, urchin, year, month) %>% 
  summarise(individual = sum(frequency)) %>% 
  right_join(m_y, by = c( "year", "month", "urchin")) %>% 
  arrange(year,month)%>% print(n = Inf)

urchin_day_mean =urchin_day %>% ungroup() %>%
  group_by(year, month, urchin) %>% 
  summarise(individual = mean(individual))%>% print(n = Inf)



WIDTH = 297/2
HEIGHT = 210/2
CLRS = as.vector(palette.colors(palette = "Okabe-Ito"))
ylabel = expression("Counts"~("indv."~m^{-2}))

write_csv(urchin_day_mean, "urchins.csv")
urchin_day_mean %>% 
  mutate(urchin = recode(urchin, 
                         "hcrassispina" = "Heliocidaris crassispina",
                         "dsetosum" = "Diadema setosum")) %>% 
  group_by(year, urchin) %>% 
  summarise_at(vars(individual),
               list(sd = sd,
                    mean = mean),
               na.rm = T) %>% 
ggplot(aes(x = year, 
           y = mean, color = urchin))+
  geom_point()+
  geom_line()+
  scale_x_continuous("",
                     limits = c(2018,2021))+
  scale_y_continuous(ylabel,
                     limits = c(0, 60))+
  scale_color_manual(values = CLRS[c(2,3)], aesthetics = c("fill", "color")) +
  theme_classic(base_size = 15)+
  theme(legend.position = c(0.2,0.9),
        legend.title = element_blank())
  
  ggsave("figure/urchin_count.png",
         width = WIDTH, height = HEIGHT, units = "mm", dpi = 600)

