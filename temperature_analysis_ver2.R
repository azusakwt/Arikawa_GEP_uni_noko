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
library(ggrepel)

###############################################################
# make figure---------------
urchins = read_csv("urchins.csv")
sargassum = read_csv("sargassum.csv")
# temperature = read_csv("temperature.csv")
temperature = read_csv("temperature_no-surface.csv")
light = read_csv("light.csv")

light = light %>% filter(str_detect(location, "garamo")) %>% 
  mutate(x = ppfd > 10) %>% 
  group_by(position, Date) %>% 
  summarise(x = sum(x)) %>% 
  mutate(x = x * 10 / 60) %>% 
  mutate(year = year(Date),
         month = month(Date)) %>% 
  group_by(position, year, month) %>% 
  summarise(x = mean(x))
  
light %>% filter(month %in% c(7,8)) %>% 
  mutate(month = factor(month,
                        levels = 1:12,
                        labels = month.abb[1:12])) %>% 
  ggplot() + 
  geom_point(aes(x = year, y = x, 
                 color = month,
                 shape = month))+
  geom_line(aes(x = year, y = x, 
               color = month,
               shape = month)) +
  scale_color_manual(values = viridis::viridis(5))
  

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

urchins = urchins %>% rename(type = urchin)
sargassum = sargassum %>% mutate(type = "Sargassum macrocarpum")

daily_t = temperature %>% 
  select(Date, daily_mean, daily_min, daily_max) %>% 
  mutate(year = year(Date),
         month = month(Date)) %>% 
  group_by(month, year) %>% 
  summarise(across(c(daily_mean,
                     daily_min,
                     daily_max),
                   mean)) 

## Bayesian Analysis ###########################################################
library(brms)
library(tidybayes)
library(bayesplot)
library(loo)

dset = expand_grid(year = 2017:2020,
                   month = 1:12) %>% full_join(daily_t) %>% 
  mutate(year = factor(year, levels = 2017:2020))

bout = brm(daily_mean ~ 
             s(month, 
               bs = "cs", by = year),
           data = dset,
           chains = 4, 
           cores = 4,
           seed = 2020,
           iter = 5000,
           control = list(adapt_delta = 0.99,
                          max_treedepth = 12),
           save_pars = save_pars(all = TRUE))

loo(bout, 
    newdata = dset %>% drop_na(), 
    cores = 4, 
    moment_match = TRUE,
    reloo = TRUE)

summary(bout)

pred = bout$data %>% 
  expand(year, 
         month = 1:12) %>% 
  add_linpred_draws(bout)

pred2 = pred %>%
  group_by(year, month) %>%
  mean_hdci(.prediction = .value) %>%
  mutate(year = factor(year))

dset %>% drop_na() %>% 
  mutate(kyphosus = daily_mean > 15,
         siganus = daily_mean > 19) %>% 
  filter(!siganus)


p1 = ggplot() +
  geom_ribbon(aes(x = month,
                  ymin = .lower,
                  ymax = .upper, fill = year),
              alpha = 0.2,
              data = pred2) +
  geom_line(aes(x = month, y = .prediction, color = year),
            data = pred2) +
  geom_point(aes(x = month, y = daily_mean, 
                 color = year,
                 shape = year),
             data = dset, 
             alpha = 0.5) +
  scale_x_continuous("Month",
                     breaks = 1:12,
                     labels = str_sub(month.abb[1:12], end = 1L)) +
  scale_y_continuous("Temperature (°C)",
                     limits = c(0, 30)) +
  scale_color_manual(values = viridis::viridis(5))+
  scale_fill_manual(values  = viridis::viridis(5)) +
  scale_shape_manual(values = c(16, 17, 15, 1, 2, 0)) +
  theme_pubr() +
  theme(legend.position = c(0.5,0),
        legend.justification = c(0.5,0),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.width = unit(20, "pt"))


pred3 = pred2 %>% 
  filter(month %in% c(1,2, 7, 8)) %>% 
  mutate(season = ifelse(!between(month, 5, 9) , "Winter","Summer"),
         month = factor(month,
                        levels = c(12, 1:11),
                        labels = month.abb[c(12, 1:11)])) %>% 
  mutate(season = factor(season,
                         levels = c("Winter", "Summer"))) %>% 
  mutate(year = as.numeric(as.character(year)))


cfs = pred3 %>% filter(str_detect(month, "Jan|Feb")) %>% 
  mutate(year = as.numeric(as.character(year)) - 2017) %>% 
  filter(year > 0) %>% 
  lm(.prediction ~ year, data = .) %>% 
  coef()

CFS = tibble(year = 2018.5,
             label = sprintf("+ %2.2f °C year\\textsuperscript{-1}", cfs[2])) %>% 
  mutate(.prediction = cfs[1] + cfs[2]*(year-2017))


p2 = 
  ggplot(pred3 %>% filter(str_detect(season, "Winter"))) +
  geom_smooth(aes(x = year, 
                  y = .prediction),
              color = "black",
              linetype = "dashed",
              method = "lm",
              formula = y~x, 
              se=F) +
  geom_errorbar(aes(x = year, 
                    y = .prediction,
                    ymin = .lower, ymax = .upper,
                    color = month),
                width = 0,
                alpha = 0.5,
                position = position_dodge(0.2)) +
  geom_point(aes(x = year, y = .prediction,
                 color = month,
                 shape = month),
             alpha = 0.5,
             position = position_dodge(0.2)) +
  geom_label_repel(aes(x = year,
                       y = .prediction,
                       label = label),
                   ylim = 5,
                   data = CFS,
                   size = 3) +
  scale_x_continuous("Year (Winter)") +
  scale_y_continuous("Temperature (°C)", 
                     limits = c(5, 20),
                     breaks = seq(5, 20, by=5)) +
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values  = viridis::viridis(4)) +
  scale_shape_manual(values = c(16, 17, 15, 1, 2, 0)) +
  theme_pubr() +
  theme(legend.position = c(0.5,1),
        legend.justification = c(0.5,1.0),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.margin = margin(t = 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.width = unit(10, "pt")) ;p2


p3 = ggplot(pred3 %>% filter(!str_detect(season, "Winter"))) +
  geom_errorbar(aes(x = year, y = .prediction,
                    ymin = .lower, ymax = .upper,
                    color = month),
                width = 0,
                alpha = 0.5,
                position = position_dodge(0.2)) +
  geom_point(aes(x = year, y = .prediction,
                 color = month,
                 shape = month),
             alpha = 0.5,
             position = position_dodge(0.2)) +
  scale_x_continuous("Year (Summer)") +
  scale_y_continuous("Temperature (°C)", 
                     limits = c(20, 30),
                     breaks = seq(20, 30, by=5)) +
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values  = viridis::viridis(4)) +
  scale_shape_manual(values = c(16, 17, 15, 1, 2, 0)) +
  theme_pubr() +
  theme(legend.position = c(0.5,1),
        legend.justification = c(0.5,1.0),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.margin = margin(t = 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.width = unit(10, "pt")) 


ggarrange(p1, ggarrange(p2,p3, ncol =2), ncol = 1)



library(tikzDevice)
options(
  tikzXelatexPackages = c(
    "\\usepackage{xeCJK}\n",
    "\\usepackage{xunicode}\n",
    "\\usepackage{amsmath,amssymb, xfrac}\n",
    "\\usepackage{unicode-math}\n",
    "\\defaultfontfeatures{Ligatures=TeX,Scale=MatchLowercase}\n",
    "\\usepackage{microtype}\n",
    "\\UseMicrotypeSet[protrusion]{basicmath}\n",
    "\\usepackage{tikz, pgf, graphics, xcolor}\n",
    "\\usepackage[active,tightpage,xetex]{preview}\n",
    "\\PreviewEnvironment{pgfpicture}\n",
    "\\setlength\\PreviewBorder{0pt}\n",
    "\\setmainfont[Path=/usr/share/fonts/Noto/]{NotoSerif-Regular.ttf}[ItalicFont = NotoSerif-Italic.ttf,
    BoldFont = NotoSerif-Bold.ttf, BoldItalicFont = Notoserif-BoldItalic.ttf]\n",
    "\\setsansfont[Path=/usr/share/fonts/Noto/]{NotoSans-Regular.ttf}[ItalicFont = NotoSans-Italic.ttf,
    BoldFont = NotoSans-Bold.ttf, BoldItalicFont = NotoSans-BoldItalic.ttf]\n",
    "\\setmonofont[Mapping=tex-ansi]{Source Code Pro}\n",
    "\\setmathfont{XITS Math}\n",
    "\\setCJKmainfont[]{Noto Serif CJK JP}\n",
    "\\setCJKsansfont[]{Noto Sans CJK JP}\n",
    "\\usepackage{textgreek}\n",
    "\\usetikzlibrary{matrix, tikzmark, fit, shapes,}\n",
    "\\usetikzlibrary{arrows, arrows.meta, shapes.geometric, positioning}\n",
    "\\usetikzlibrary{calc, intersections,through}\n",
    "\\usetikzlibrary{decorations.text, decorations.markings}\n"
  ),
  tikzDefaultEngine="xetex",
  tikzMetricPackages = c("\\usetikzpackage{calc}"))



if(!dir.exists("figure")) {
  dir.create("figure")
} else {
  # 論文に記載した場合，図の幅は　80 mm　か 160 mm.
  wh = 100
  ht = 100
  pt = 11
  DPI = 800
  ROOT = rprojroot::find_rstudio_root_file()
  
  filename = "Figure_temperature_gam.tex"
  filename = paste0(ROOT, "/figure/", filename)
  tikz(filename, 
       width = wh/25.4, 
       height = ht/25.4, 
       pointsize = pt, standAlone = TRUE)
  print(
    ggarrange(p1, ggarrange(p2,p3, ncol =2), ncol = 1)
  ) 
  dev.off() # Never forget to run dev.off() after tikz()
  tinytex::xelatex(filename) # compile tex to pdf
  img = magick::image_read(str_replace(filename, "tex", "pdf"), density = 300)
  magick::image_write(img, str_replace(filename, "tex", "png"), format = "png")
}




