rm(list = ls())
library(tidyverse);library(patchwork)
setwd("D:/taihu2000_2019")
load("smap.RData")
set.seed(1234)
# ---------- 1.---------
phy_vars <- c("Chlorophyta", "Bacillariophyta", "Chrysophyta", "Cyanophyta",
              "Euglenophyta", "Cryptophyta", "Pyrrophyta", "Xanthophyta")
zoo_vars <- c("Protozoa", "Rotifera", "Cladocera", "Copepoda")
biotic_vars <- c(phy_vars, zoo_vars)
interaction_levels <- c("Competition/Inhibition","Facilitation/Promotion")
site_order <- c("THL00", "THL01", "THL03", "THL04",
                "THL05", "THL06", "THL07", "THL08")
group_of <- function(x) {
  case_when(
    x %in% phy_vars ~ "Phytoplankton",
    x %in% zoo_vars ~ "Zooplankton",
    TRUE ~ NA_character_)}
# ---------- 2.    S-map        ----------
smap_detail <- map_dfr(
  seq_len(nrow(smapc.mat.orig0)),
  \(i) {
    x <- smapc.mat.orig0$int.ccm.smap[[i]]
    if (is.null(x) || nrow(x) == 0) {
      return(tibble())}
    x %>%
      as_tibble() %>%
      mutate(Site = smapc.mat.orig0$Site[i]) %>%
      separate(
        ccmLink,into = c("from", "to"),
        sep = "->",remove = FALSE)}) %>%
  group_by(Site, from, to, ccmLink) %>%
  summarise(
    MedianIS = median(MedianIS, na.rm = TRUE),
    Perc.Pos.int = median(Perc.Pos.int, na.rm = TRUE),
    Willcox.conv = first(Willcox.conv),
    .groups = "drop")
# ---------- 3.    Best E   Library size ----------
interactions <- linkM.std0000 %>%
  filter(
    Period == "2007-2019",
    Convergence,
    Cause %in% biotic_vars,
    Effect %in% biotic_vars,
    group_of(Cause) == group_of(Effect),
    Cause != Effect,
    !is.na(E),
    !is.na(lib_size),
    !is.na(Std_L_strength),
    Std_L_strength > 0) %>%
  transmute(
    Site,
    from = Cause,
    to = Effect,
    Group = group_of(Cause),
    Best_E = E,
    Library_size = lib_size,
    Std_CCM_strength = Std_L_strength,
    ccmLink = paste(Cause, Effect, sep = "->")) %>%
  inner_join(
    smap_detail,
    by = c("Site", "from", "to", "ccmLink")) %>%
  filter(
    Willcox.conv,
    !is.na(MedianIS),
    !is.na(Perc.Pos.int)) %>%
  mutate(
    Interaction = case_when(
      MedianIS > 0 & Perc.Pos.int > 50 ~
        "Facilitation/Promotion",
      MedianIS < 0 & Perc.Pos.int < 50 ~
        "Competition/Inhibition",
      TRUE ~ NA_character_),
    Site = factor(
      Site,
      levels = rev(site_order)),
    Group = factor(
      Group,
      levels = c("Phytoplankton", "Zooplankton")),
    Interaction = factor(
      Interaction,
      levels = interaction_levels)) %>%
  filter(
    !is.na(Interaction),
    !is.na(Site))
# Frequency is calculated within each site and plankton group.
parameter_data <- interactions %>%
  select(Site, Group, Best_E, Library_size) %>%
  pivot_longer(c(Best_E, Library_size), names_to = "Parameter", values_to = "Value") %>%
  filter(is.finite(Value))
frequency_data <- parameter_data %>%
  count(Parameter, Site, Group, Value, name = "Count") %>%
  group_by(Parameter, Site, Group) %>%
  mutate(Frequency = Count / sum(Count)) %>%
  ungroup()
frequency_max <- max(frequency_data$Frequency)

plot_quantile_ribbon <- function(parameter, title, x_breaks) {
  dat <- parameter_data %>% filter(Parameter == parameter)
  stats <- dat %>% group_by(Site, Group) %>%
    summarise(Q10 = quantile(Value, .1), Q25 = quantile(Value, .25),
              Q40 = quantile(Value, .4), Median = median(Value),
              Q60 = quantile(Value, .6), Q75 = quantile(Value, .75),
              Q90 = quantile(Value, .9), .groups = "drop") %>%
    mutate(Y = as.numeric(Site),
           Dir = if_else(Group == "Phytoplankton", 1, -1))
  bands <- bind_rows(
    stats %>% transmute(Site, Group, Y, Dir, Lo = Q10, Hi = Q90, Depth = .33, Alpha = .2),
    stats %>% transmute(Site, Group, Y, Dir, Lo = Q25, Hi = Q75, Depth = .24, Alpha = .45),
    stats %>% transmute(Site, Group, Y, Dir, Lo = Q40, Hi = Q60, Depth = .14, Alpha = .9)
  )
  connectors <- stats %>% select(Site, Group, Median, Y) %>%
    pivot_wider(names_from = Group, values_from = Median)
  ggplot() +
    geom_hline(yintercept = seq_along(site_order), colour = "#E5E8EB", linewidth = .35) +
    geom_rect(data = bands,
              aes(xmin = Lo, xmax = Hi, ymin = pmin(Y, Y + Dir * Depth),
                  ymax = pmax(Y, Y + Dir * Depth), fill = Group, alpha = Alpha),
              colour = NA) +
    geom_segment(data = bands %>% filter(Lo == Hi),
                 aes(x = Lo, xend = Hi, y = Y, yend = Y + Dir * Depth,
                     colour = Group, alpha = Alpha), linewidth = 1.2, show.legend = FALSE) +
    geom_segment(data = connectors,
                 aes(x = Phytoplankton, xend = Zooplankton, y = Y + .12, yend = Y - .12),
                 colour = "#53616A", linewidth = .45, linetype = "dotted") +
    geom_segment(data = stats,
                 aes(x = Median, xend = Median, y = Y, yend = Y + Dir * .34),
                 colour = "#17252C", linewidth = .8) +
    scale_fill_manual(values = c(Phytoplankton = "#465C78", Zooplankton = "#B87763"), name = NULL) +
    scale_colour_manual(values = c(Phytoplankton = "#465C78", Zooplankton = "#B87763"), guide = "none") +
    scale_alpha_identity() +
    scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(site_order), labels = rev(site_order), expand = c(0, 0)) +
    coord_cartesian(xlim = range(dat$Value) + c(-.55, .55), ylim = c(.5, length(site_order) + .5)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_family = "serif", base_size = 11.5) +
    theme(text = element_text(face = "bold", size = 11.5),
          axis.text = element_text(face = "bold", size = 11.5, colour = "#30383C"),
          plot.title = element_text(size = 11.5, margin = margin(b = 12)),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.x = element_line(linewidth = .45, colour = "#59656B"),
          legend.text = element_text(face = "bold", size = 11.5),
          aspect.ratio = 1, plot.margin = margin(8, 12, 8, 8))
}
p1 <- plot_quantile_ribbon("Best_E", "(a)  Best embedding dimension (E)", seq(2, 12, 2))
p2 <- plot_quantile_ribbon("Library_size", "(b)  Library size (L)", seq(120, 155, 5))
final_plot <- wrap_plots(p1, p2, ncol = 2, guides = "collect") & theme(legend.position = "bottom")
final_plot
