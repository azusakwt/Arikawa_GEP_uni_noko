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
oxygen_rate = read_csv("oxygen_rate.csv")
light_all = read_csv("light.csv")
urchins = read_csv("urchins.csv")
sargassum = read_csv("sargassum.csv")

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


urchins = urchins %>% rename(type = urchin)
sargassum = sargassum %>% mutate(type = "Sargassum macrocarpum")

dset = full_join(urchins, sargassum)


oxy = oxygen_rate %>% mutate(month = month(Date),
                             year = year(Date)) %>% 
  mutate(H = hour(datetime) + minute(datetime)/60)

##################################################################
light_0 = light_all %>% filter(position == "0m") %>% 
  select(-position)

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
gep = product_all %>% filter(str_detect(location, "garamo")) %>% 
  mutate(location = factor(location,
                           levels = c("arikawagaramo"),
                           labels = c("Large Macroalgal Ecosystem")))

library(rstanarm)
library(tidybayes)

sout = stan_glm(GEP ~ year, data = gep, 
                family = Gamma("log"),
                cores = 4, chains = 4)

sout %>% summary()

pred = gep %>% 
  expand(year) %>% 
  add_linpred_draws(sout) %>% 
  group_by(year) %>% 
  mean_hdci(.value)


#########################################################################

xlabel = "Year"
ylabel = expression("GEP"~(mg~O[2]~m^{-3}~day^{-1}))
ylabel = "GEP (mg O\\textsubscript{2}~m\\textsuperscript{-3}~d\\textsuperscript{-1})"
clr = viridis::viridis(4)

gep_summary  = gep %>% 
  group_by(year) %>% 
  summarise_at(vars(GEP), 
               list(sd = sd, mean = mean, n = length)) %>% 
  mutate(se = sd / sqrt(n - 1)) %>% 
  mutate(type = "GEP")

pred = pred %>% 
  mutate(type = "Model expectation \\& 95\\% HDCI")
size = 0.2
p1 = ggplot()+
  geom_ribbon(aes(x = year, 
                  ymin = .lower, 
                  ymax = .upper,
                  fill = type),
              alpha = 0.2,
              data = pred) +
  geom_ribbon(aes(x = year, 
                  ymin = 0, 
                  ymax = 0,
                  fill = type),
              alpha = 0.2,
              data = gep_summary) +
  geom_line(aes(x = year, 
                y = .value,
                color = type),
            data = pred) +
  geom_errorbar(aes(x     = year, 
                    ymin  = mean - 1.96*se,
                    ymax  = mean + 1.96*se,
                    color = type),
                width = 0,
                show.legend  = F,
                data = gep_summary) +
  geom_point(aes(x     = year, 
                 y     = mean,
                 color = type),
             data = gep_summary) +
  scale_x_continuous(xlabel) +
  scale_y_continuous(ylabel, 
                     limits = c(0,20)) +
  scale_fill_manual(name = "x", 
                    values  = clr[1:2],
                    guide = guide_legend(override.aes = list(
                      shape = c(19,NA),
                      size= c(2,1),
                      fill = c(NA, clr[2]),
                      linetype = c(0,1)
                    ))) +
  scale_color_manual(name = "x", values = clr) +
  theme_pubr() +
  theme(legend.position =      c(0, 0.5),
        legend.justification = c(0, 1),
        legend.background =    element_blank(),
        legend.title =         element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) 


ylabel2 = expression("Density"~(ind.~m^{-2}))
ylabel2 = "Density~(ind.~m\\textsuperscript{-2})"

p2 = dset %>% 
  mutate(type = str_glue("\\textit{{{type}}}")) %>% 
  ggplot() +
  geom_line(aes( x = year, 
                 y = mean,
                 color = type)) +
  geom_errorbar(aes( x = year, 
                     color = type,
                     ymin = mean - 1.96 * se,
                     ymax = mean + 1.96 * se),
                width = 0) +
  geom_point(aes( x = year, 
                  y = mean,
                  color = type,
                  shape = type)) +
  scale_x_continuous(xlabel,
                     limits = c(2017, 2020)) +
  scale_y_continuous(ylabel2,
                     limits = c(0, 200)) +
  scale_color_manual(values = clr) +
  scale_fill_manual(values  = clr) +
  theme_pubr() +
  theme(legend.position =      c(0, 0.5),
        legend.justification = c(0, 0.5),
        legend.background =    element_blank(),
        legend.title =         element_blank()) 


cowplot::plot_grid(p1,p2,
                   ncol = 1,
                   align = "vh",
                   axis = "tblr",
                   labels = c("(A)", "(B)"),
                   label_x = 1,
                   hjust = 1,
                   vjust = 1)



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
  
  filename = "Figure_greg.tex"
  filename = paste0(ROOT, "/figure/", filename)
  tikz(filename, 
       width = wh/25.4, 
       height = ht/25.4, 
       pointsize = pt, standAlone = TRUE)
  print(
    cowplot::plot_grid(p1,p2,
                       ncol = 1,
                       align = "vh",
                       axis = "tblr",
                       labels = c("(A)", "(B)"),
                       label_x = 1,
                       label_size = 11,
                       label_fontface = "plain",
                       hjust = 1,
                       vjust = 1)
  ) 
  dev.off() # Never forget to run dev.off() after tikz()
  tinytex::xelatex(filename) # compile tex to pdf
  img = magick::image_read(str_replace(filename, "tex", "pdf"), density = 600)
  magick::image_write(img, str_replace(filename, "tex", "png"), format = "png")
}

