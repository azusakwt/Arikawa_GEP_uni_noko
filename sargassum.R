# ドムさんコドラートデータ

# package---------------
library(tidyverse)
library(lubridate)
library(readxl)

#図のサイズ
source("size.R")

# data---------------
garamo = read_csv("~/Lab_Data/quadrat_data/garamo_quadrats/data/sargassum_kelp_counts_csv.csv")

garamo = garamo %>%  mutate(date = parse_date_time(date, "%d/%m/%y", locale = "ja_JP.utf-8")) %>%
  select(date, location, life_history, quadrat, individuals)

## Join month and year ---------------
## Calculate sum, mean and sd of frequencies
garamo =
  garamo %>% mutate(life_history = factor(life_history,
                                          levels = c("juvenile", "adult_imm", "adult_mat",
                                                     "adult_post", "kelp"),
                                          labels = c("juvenile", "adult_imm", "adult_mat",
                                                     "adult_post", "kelp")))

garamo = garamo %>% mutate(location = recode(location, "Arikawa" = "Sargassum"))

y =data_frame(year = c(2018:2020)) %>% 
  mutate(year = as.numeric(year))
m = data_frame(month = c(1:12)) %>% 
  mutate(month = as.numeric(month))

m_y = full_join(y,m, by = character()) %>% mutate(location = "Sargassum")  
m_y = m_y[-1:-3,]

garamo = garamo %>%
  filter(!location == "kelp") %>% 
  group_by(location, date) %>% 
  summarise(count = sum(individuals, .keep_all = TRUE)) %>% 
  mutate(year = year(date), month = month(date)) %>% 
  right_join(m_y, by = c("year", "month", "location")) %>%
  arrange(year,month) %>% 
  print(n = Inf) 

# adult_imm：成熟個体
# adult_mat：生殖起床がある個体
# adult_post：付着器のみ
# kelp：ノコギリモク以外の海藻
xlabel = ""
garamo %>% 
  ggplot(aes(x = year+(month-1)/12, y = count))+
  geom_point()+
  geom_line()+
  scale_x_date(xlabel,
               limits = c(2018, 2021))+
  scale_y_continuous(ylabel)+
  theme_classic(base_size = 16)

# ggsave("figure/surgassum_count.png",
#        width = WIDTH, height = HEIGHT, units = "mm", dpi = 600)

garamo %>% filter(!life_history == "kelp") %>% 
  group_by(location, date, life_history) %>% 
  summarise(individuals = sum(individuals)) %>%
  ggplot(aes(x = date, y = individuals))+
  geom_point()+
  geom_line() +
  # scale_x_continuous(labels = month)+
  facet_grid(rows = "life_history", switch="y") + 
  scale_x_date(xlabel,limits = c(ymd("2018-01-01"), ymd("2021-01-01")))+
  scale_y_continuous(ylabel, limits = c(0,120),
                     breaks = c(0,30,60,90,120))+
  theme_classic()+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        legend.position = "none")

# ggsave("figure/surgassum_count_lifehistory.png",
#        width = WIDTH, height = HEIGHT, units = "mm", dpi = 600)

garamo %>% 
  group_by(date, location, month) %>% 
  summarise(count = sum(individuals)) %>% 
  ggplot()+
  geom_boxplot(aes(x = month, y = count, group = month))

garamo_mean = garamo %>% 
  group_by(location, month) %>% 
  summarise(count = mean(individuals))

garamo_mean %>% ungroup() %>% group_by(location) %>% 
  summarise(mean = mean(count))

garamo_mean %>% filter(count>5.98)

garamo = garamo %>%
  mutate(location = if_else(str_detect(life_history, "kelp"), "kelp", location)) %>% 
  group_by(location, year, month,date) %>% 
  summarise(individuals = mean(individuals))

garamo %>% print(n = Inf)

write_csv(garamo, "sargassum.csv")
# figure---------------------------
WIDTH = 297/2
HEIGHT = 210/2
CLRS = as.vector(palette.colors(palette = "Okabe-Ito"))
ylabel = expression("Counts"~("indv."~m^{-2}))

garamo  %>% 
  ggplot(aes(x =  year+(month-1)/12, y = count, color = "Sargassum macrocarpum"))+
  geom_point()+
  geom_line()+
  scale_x_continuous(xlabel,
                     limits = c(2018,2021))+
  scale_y_continuous(ylabel, limits = c(0,250))+
  scale_color_manual(values = CLRS[2], aesthetics = c("color")) +
  theme_classic(base_size = 15)+
  theme(legend.position = c(0.75,0.9),
        legend.background = element_blank(),
        legend.title = element_blank())

ggsave("figure/garamo.png", width = WIDTH, height = HEIGHT, units = "mm", dpi = 600)
##############################################################################
