library(tidyverse);library(ggplot2)
# 1. 读取数据
load("D:/taihu2000_2019/smapc.mat.orig.Rdata")
# 2. 分类与样点顺序
font <- "Times New Roman"
taxa <- c("Bacillariophyta", "Chlorophyta", "Cyanophyta", "Chrysophyta",
          "Cryptophyta", "Euglenophyta", "Pyrrophyta", "Xanthophyta",
          "Protozoa", "Rotifera", "Cladocera", "Copepoda")
upper_sites <- c("THL00", "THL01", "THL03", "THL04")
lower_sites <- c("THL05", "THL06", "THL07", "THL08")
sites <- c(upper_sites, lower_sites)
year_breaks <- as.Date(paste0(c(2007, 2009, 2012, 2014, 2017, 2019), "-01-01"))
site_cols <- c(THL00 = "#0072B2", THL01 = "#E69F00",THL03 = "#009E73", THL04 = "#CC3377",
               THL05 = "#D55E00", THL06 = "#7733AA",THL07 = "#0099AA", THL08 = "#222222")
# 4. 提取12类浮游生物生物量
dat <- map2_dfr(smapc.mat.orig$env.int.mat, smapc.mat.orig$Site, function(d, s) {
  d %>% select(Date, all_of(taxa)) %>% mutate(Site = s, Date = as.Date(Date))
}) %>%
  pivot_longer(all_of(taxa), names_to = "Taxon", values_to = "Value")

stopifnot(setequal(unique(dat$Site), sites), !anyNA(dat$Date),
          !anyDuplicated(dat[c("Site", "Taxon", "Date")]),
          all(dat$Value[is.finite(dat$Value)] >= 0))
# 5. 各类群×样点按自身最大值缩放；缺失月份保持断线
dat <- dat %>%
  mutate(Site = factor(Site, levels = sites), Taxon = factor(Taxon, levels = taxa)) %>%
  group_by(Site, Taxon) %>%
  complete(Date = seq.Date(min(Date), max(Date), by = "month")) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(Value = if_else(is.finite(Value), Value, NA_real_),
         Maximum = if (all(is.na(Value))) NA_real_ else max(Value, na.rm = TRUE),
         Scaled = Value / if_else(Maximum > 0, Maximum, 1),
         Plot_value = if_else(Site %in% lower_sites, -Scaled, Scaled),
         Segment = cumsum(is.na(Value))) %>%
  ungroup()
# 6. 镜像时间序列：3行×4列，6个年份刻度
fig <- ggplot(dat, aes(Date, Plot_value, colour = Site, group = interaction(Site, Segment))) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = .5) +
  geom_line(linewidth = .4, alpha = 1, lineend = "round", na.rm = TRUE) +
  facet_wrap(~ Taxon, nrow = 3, ncol = 4) +
  scale_colour_manual(name = NULL, values = site_cols, breaks = sites, drop = FALSE) +
  scale_x_date(breaks = year_breaks, date_labels = "%Y",
               expand = expansion(mult = c(.01, .025))) +
  scale_y_continuous(limits = c(-1.05, 1.05), breaks = c(-1, -.5, 0, .5, 1),
                     labels = c("1.0", "0.5", "0", "0.5", "1.0"),
                     expand = expansion(mult = 0)) +
  labs(x = "Year", y = "Scaled biomass") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = .5))) +
  theme_classic(base_family = font, base_size = 12) +
  theme(
    aspect.ratio = .71,
    strip.background = element_blank(),
    strip.text = element_text(size = 12,  hjust = 0, margin = margin(b = 4)),
    axis.title = element_text(size = 12),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line = element_line(colour = "grey40", linewidth = .3),
    axis.ticks = element_line(colour = "grey40", linewidth = .3),
    axis.ticks.length = grid::unit(1, "mm"),
    panel.spacing.x = grid::unit(6, "mm"),
    panel.spacing.y = grid::unit(4, "mm"),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.key.width = grid::unit(10, "mm"),
    legend.spacing.x = grid::unit(3, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(5, 8, 4, 4))
print(fig)
ggsave(filename = "D:/相互作用/biomass.svg", plot = fig,
       width = 14, height = 9, units = "in", bg = "white")
