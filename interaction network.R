rm(list = ls())
library(tidyverse);library(ggnewscale)
library(patchwork)
#
setwd("D:/taihu2000_2019")
load("smap.RData");set.seed(1234)
# ---------- 1. 数据与分组 ----------
phy_vars <- c(
  "Chlorophyta", "Bacillariophyta", "Chrysophyta", "Cyanophyta",
  "Euglenophyta", "Cryptophyta", "Pyrrophyta", "Xanthophyta")
zoo_vars <- c("Protozoa", "Rotifera", "Cladocera", "Copepoda")
biotic_vars <- c(phy_vars, zoo_vars)
interaction_levels <- c("Competition/Inhibition", "Facilitation/Promotion")
group_of <- function(x) case_when(
  x %in% phy_vars ~ "Phytoplankton",
  x %in% zoo_vars ~ "Zooplankton",
  TRUE ~ NA_character_)
# smap direction
smap_detail <- map_dfr(seq_len(nrow(smapc.mat.orig0)), ~
                         smapc.mat.orig0$int.ccm.smap[[.x]] %>%
                         as_tibble() %>%
                         mutate(Site = smapc.mat.orig0$Site[.x]) %>%
                         separate(ccmLink, c("from", "to"), sep = "->", remove = FALSE) %>%
                         group_by(Site, from, to, ccmLink) %>%
                         summarise(
                           MedianIS = median(MedianIS, na.rm = TRUE),
                           Perc.Pos.int = median(Perc.Pos.int, na.rm = TRUE),
                           Willcox.conv = first(Willcox.conv), .groups = "drop"))
# ccm rho
interactions <- linkM.std0000 %>%
  filter(Period == "2007-2019", Convergence,
         Cause %in% biotic_vars, Effect %in% biotic_vars,
         group_of(Cause) == group_of(Effect), Cause != Effect,
         !is.na(Std_L_strength), Std_L_strength > 0) %>%
  transmute(Site, from = Cause, to = Effect, Group = group_of(from),
            Std_CCM_strength = Std_L_strength,
            ccmLink = paste(from, to, sep = "->")) %>%
  inner_join(smap_detail, by = c("Site", "from", "to", "ccmLink")) %>%
  filter(!is.na(MedianIS)) %>%
  mutate(Interaction = case_when(
    Willcox.conv & MedianIS > 0 & Perc.Pos.int > 50 ~"Facilitation/Promotion",
    Willcox.conv & MedianIS < 0 & Perc.Pos.int < 50 ~"Competition/Inhibition"),
    Interaction = factor(Interaction, levels = interaction_levels)) %>%
  filter(!is.na(Interaction))
all_sites <- sort(unique(interactions$Site))
write.csv(interactions,"D:/相互作用/interaction.csv")
#
# ---------- 2. 配色、字体与样点 ----------
interaction_cols <- c("Competition/Inhibition" = "#B64342",
                      "Facilitation/Promotion" = "#0F4D92")
interaction_labels <- c("Competition/Inhibition" = "Competition",
                        "Facilitation/Promotion" = "Facilitation")
zone_cols <- c("Bay" = "#3C5488", "Estuary" = "#D98C6B","Lake centre" = "#4F8F78")
taxon_cols <- c(Chlorophyta = "#8BCF8B", Bacillariophyta = "#4DBBD5",
                Chrysophyta = "#D8B45A", Cyanophyta = "#3775BA",
                Euglenophyta = "#E6C35C", Cryptophyta = "#E9A6A1",
                Pyrrophyta = "#B64342", Xanthophyta = "#D9D2A3",
                Protozoa = "#767676", Rotifera = "#D8BCC8",
                Cladocera = "#6FAFAE", Copepoda = "#C77C52")
plot_font <- "Times New Roman";text_col <- "#272727"
# base font
theme_set(theme_gray(base_family = plot_font))
site_meta <- tibble(Site = all_sites) %>%
  mutate(Zone = case_when(
    Site %in% c("THL00", "THL01", "THL03", "THL04", "THL05", "THL06") ~ "Bay",
    Site == "THL07" ~ "Estuary",Site == "THL08" ~ "Lake centre"),
    Zone = factor(Zone, levels = names(zone_cols))) %>%
  arrange(Zone, Site) %>%
  mutate(Site_panel = factor(Site, levels = Site))
# ---------- 3. summary ----------
fmt_mean_sd <- function(x, s) {
  ifelse(is.na(x), NA_character_, sprintf("%.3f +/- %.3f", x, s))}
# rho
rho_summary <- function(data, ...) {
  data %>%group_by(...) %>%
    summarise(n = n(), mean_rho = mean(Std_CCM_strength),
              sd_rho = sd(Std_CCM_strength),
              median_rho = median(Std_CCM_strength), .groups = "drop")}
# summary interaction
stat_data <- interactions %>%
  left_join(dplyr::select(site_meta, Site, Zone), by = "Site") %>%
  mutate(
    Group = factor(Group, levels = c("Phytoplankton", "Zooplankton")),
    Interaction = factor(Interaction, levels = interaction_levels))

quality_check <- stat_data %>%
  summarise(
    rows = n(),
    unique_directed_links = n_distinct(Site, from, to),
    duplicated_directed_links = rows - unique_directed_links,
    missing_zone = sum(is.na(Zone)),
    missing_strength = sum(is.na(Std_CCM_strength)),
    min_rho = min(Std_CCM_strength, na.rm = TRUE),
    max_rho = max(Std_CCM_strength, na.rm = TRUE))
#
overall_network <- rho_summary(stat_data, Group) %>%
  mutate(mean_sd = fmt_mean_sd(mean_rho, sd_rho))
overall_interaction <- rho_summary(stat_data, Group, Interaction) %>%
  group_by(Group) %>%
  mutate(total = sum(n), pct = 100 * n / total,
         mean_sd = fmt_mean_sd(mean_rho, sd_rho)) %>%
  ungroup()
overall_wilcox <- stat_data %>%
  group_by(Group) %>%group_modify(~ {
    z <- wilcox.test(Std_CCM_strength ~ Interaction, data = .x,
                     exact = FALSE, conf.int = TRUE)
    tibble(
      n_competition = sum(.x$Interaction == interaction_levels[1]),
      n_facilitation = sum(.x$Interaction == interaction_levels[2]),
      W = unname(z$statistic), p = z$p.value,
      significance = case_when(
        z$p.value < 0.001 ~ "***", z$p.value < 0.01 ~ "**",
        z$p.value < 0.05 ~ "*", TRUE ~ "ns"))}) %>%
  ungroup()
#
zone_sites <- site_meta %>%
  distinct(Zone, Site) %>%
  count(Zone, name = "n_sites")
zone_network <- rho_summary(stat_data, Group, Zone) %>%
  left_join(zone_sites, by = "Zone") %>%
  mutate(
    n_per_site = n / n_sites,
    mean_sd = fmt_mean_sd(mean_rho, sd_rho))
zone_interaction <- rho_summary(stat_data, Group, Zone, Interaction) %>%
  complete(Group, Zone, Interaction, fill = list(n = 0L)) %>%
  left_join(zone_sites, by = "Zone") %>%
  group_by(Group, Zone) %>%
  mutate(zone_total = sum(n), pct = 100 * n / zone_total,
         n_per_site = n / n_sites, mean_sd = fmt_mean_sd(mean_rho, sd_rho)) %>%
  ungroup()
site_network <- rho_summary(stat_data, Group, Zone, Site) %>%
  mutate(mean_sd = fmt_mean_sd(mean_rho, sd_rho))
site_interaction <- tidyr::crossing(
  distinct(stat_data, Group, Zone, Site),
  Interaction = factor(interaction_levels, levels = interaction_levels)) %>%
  left_join(
    rho_summary(stat_data, Group, Zone, Site, Interaction),
    by = c("Group", "Zone", "Site", "Interaction")) %>%
  mutate(n = replace_na(n, 0L)) %>%
  group_by(Group, Zone, Site) %>%
  mutate(
    site_total = sum(n), pct = 100 * n / site_total,
    mean_sd = fmt_mean_sd(mean_rho, sd_rho)) %>%
  ungroup()
zone_result_table <- zone_interaction %>%
  mutate(Interaction = recode(
    as.character(Interaction),
    "Competition/Inhibition" = "Competition",
    "Facilitation/Promotion" = "Facilitation")) %>%
  dplyr::select(
    Group, Zone, Interaction, n, pct, n_per_site,
    mean_rho, sd_rho) %>%
  pivot_wider(
    names_from = Interaction,
    values_from = c(n, pct, n_per_site, mean_rho, sd_rho),
    names_glue = "{.value}_{Interaction}") %>%
  left_join(
    zone_network %>%
      dplyr::select( Group, Zone, total_n = n, n_sites,
              total_n_per_site = n_per_site,
              overall_mean_rho = mean_rho, overall_sd_rho = sd_rho),
    by = c("Group", "Zone"))

stat_results <- list(
  quality_check = quality_check,
  overall_network = overall_network,
  overall_interaction = overall_interaction,
  overall_wilcox = overall_wilcox,
  zone_network = zone_network,
  zone_interaction = zone_interaction,
  site_network = site_network,
  site_interaction = site_interaction,
  zone_result_table = zone_result_table)
# ---------- 4. 网络图数据 ----------
make_layout <- function(species, group, centre_x, radius) {
  tibble(name = species, Group = group, node_id = seq_along(species)) %>%
    mutate(angle = pi / 2 - 2 * pi * (node_id - 1) / length(species),
           x = centre_x + radius * cos(angle),
           y = radius * sin(angle))}

node_layout <- bind_rows(
  make_layout(phy_vars, "Phytoplankton", -0.72, 0.63),
  make_layout(zoo_vars, "Zooplankton", 0.78, 0.43))
node_xy <- dplyr::select(node_layout, name, x, y)

node_degree <- bind_rows(
  transmute(interactions, Site, name = from),
  transmute(interactions, Site, name = to)) %>%
  count(Site, name, name = "Degree")

plot_nodes <- expand_grid(Site = all_sites, name = biotic_vars) %>%
  left_join(node_layout, by = "name") %>%
  left_join(node_degree, by = c("Site", "name")) %>%
  left_join(site_meta, by = "Site") %>%
  mutate(Degree = replace_na(Degree, 0))

plot_edges <- interactions %>%
  left_join(
    rename(node_xy, from = name, x_from = x, y_from = y), by = "from") %>%
  left_join(
    rename(node_xy, to = name, x_to = x, y_to = y), by = "to") %>%
  left_join(site_meta, by = "Site") %>%
  mutate(pair_from = pmin(from, to), pair_to = pmax(from, to)) %>%
  group_by(Site, Group, pair_from, pair_to) %>%
  mutate(
    Bidirectional = n_distinct(paste(from, to, sep = "->")) > 1,
    Edge_shape = case_when(
      !Bidirectional ~ "Straight", from < to ~ "Positive",
      TRUE ~ "Negative")) %>%
  ungroup() %>%
  mutate(
    distance = sqrt((x_to - x_from)^2 + (y_to - y_from)^2),
    x_start = x_from + 0.055 * (x_to - x_from) / distance,
    y_start = y_from + 0.055 * (y_to - y_from) / distance,
    x_end = x_to - 0.080 * (x_to - x_from) / distance,
    y_end = y_to - 0.080 * (y_to - y_from) / distance,
    draw_order = match(Interaction, interaction_levels)) %>%
  arrange(draw_order, Std_CCM_strength)

edge_check <- interactions %>%
  count(Site, name = "input_n") %>%
  left_join(count(plot_edges, Site, name = "draw_n"), by = "Site") %>%
  mutate(
    draw_n = replace_na(draw_n, 0L), difference = input_n - draw_n,
    all_drawn = input_n == draw_n)

group_labels <- tidyr::crossing(
  Site = all_sites, Group = c("Phytoplankton", "Zooplankton")) %>%
  mutate(
    x = if_else(Group == "Phytoplankton", -0.72, 0.78),
    y = 0.82, label = Group) %>%
  left_join(site_meta, by = "Site")

site_stats <- interactions %>%
  group_by(Site) %>%
  summarise(
    total = n(),
    competition = 100 * mean(Interaction == interaction_levels[1]),
    facilitation = 100 * mean(Interaction == interaction_levels[2]),
    .groups = "drop") %>%
  mutate(label = sprintf(
    "n = %d  |  Comp. %.0f%%  |  Facil. %.0f%%",
    total, competition, facilitation)) %>%
  left_join(site_meta, by = "Site")

edge_mapping <- aes(
  x = x_start, y = y_start, xend = x_end, yend = y_end,
  colour = Interaction, linewidth = Std_CCM_strength,
  alpha = Interaction)
edge_arrow <- arrow(length = unit(0.75, "mm"), type = "closed")

edge_layer <- function(data, curvature = NULL) {
  args <- list(
    data = data, mapping = edge_mapping, lineend = "round",
    arrow = edge_arrow,
    show.legend = c(colour = TRUE, linewidth = TRUE, alpha = FALSE))
  if (is.null(curvature)) do.call(geom_segment, args) else
    do.call(geom_curve, c(args, list(curvature = curvature, angle = 90)))}
# ---------- 5. 网络图 ----------
p_network <- ggplot() +
  geom_segment(data = site_meta,
               aes(x = -1.43, xend = 1.43, y = 1.03, yend = 1.03, colour = Zone),
               linewidth = 2.2, lineend = "butt") +
  scale_colour_manual(
    values = zone_cols, name = "Zone",
    guide = guide_legend(order = 1, ncol = 1, byrow = TRUE,
                         override.aes = list(linewidth = 2.4, alpha = 1))) +
  ggnewscale::new_scale_colour() +
  edge_layer(filter(plot_edges, Edge_shape == "Straight")) +
  edge_layer(filter(plot_edges, Edge_shape == "Positive"), 0.26) +
  edge_layer(filter(plot_edges, Edge_shape == "Negative"), -0.26) +
  scale_colour_manual(values = interaction_cols, breaks = interaction_levels,
                      labels = interaction_labels, name = "Interaction", drop = FALSE) +
  scale_alpha_manual(
    values = set_names(rep(0.86, 2), interaction_levels), guide = "none") +
  scale_linewidth_continuous(
    range = c(0.20, 0.95), breaks = c(0.2, 0.4, 0.6, 0.8),
    limits = c(0, max(plot_edges$Std_CCM_strength, na.rm = TRUE)),
    name = "Std.CCM rho") +
  geom_point(data = plot_nodes, aes(x, y, fill = name), shape = 21, size = 2.1,
             colour = "#D8D8D8", stroke = 0.25, alpha = 0.18,show.legend = FALSE) +
  geom_point(data = filter(plot_nodes, Degree > 0),
             aes(x, y, size = Degree, fill = name), shape = 21,
             colour = "#4D4D4D", stroke = 0.32, key_glyph = "point",
             show.legend = c(size = FALSE, fill = TRUE)) +
  scale_fill_manual(values = taxon_cols, breaks = biotic_vars, name = "Taxon") +
  scale_size_continuous(range = c(2.7, 5.2), guide = "none") +
  geom_text(data = group_labels, aes(x, y, label = label),
            family = plot_font, fontface = "bold", size = 4, colour = text_col) +
  geom_text(data = site_stats, aes(x = 0, y = -0.98, label = label),
            family = plot_font, fontface = "bold", size = 4, colour = text_col) +
  facet_wrap(~Site_panel, ncol = 4, drop = FALSE) +
  guides(colour = guide_legend(order = 2, ncol = 1, byrow = TRUE,
                               override.aes = list(linewidth = 0.9, alpha = 1)),
         linewidth = guide_legend(order = 3, ncol = 1, byrow = TRUE,
                                  override.aes = list(colour = "#606060", alpha = 1) ),
         fill = guide_legend(order = 4, ncol = 1, byrow = TRUE,
                             override.aes = list(shape = 21, size = 3, stroke = 0.32,
                                                 colour = "#4D4D4D", alpha = 1))) +
  coord_fixed(xlim = c(-1.48, 1.48), ylim = c(-1.08, 1.08),
              clip = "off", expand = FALSE) +
  theme_void(base_size = 12, base_family = plot_font) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        strip.background = element_rect(),
        strip.text = element_text(face = "bold", size = 12, colour = text_col,margin = margin(3.5, 2, 3.5, 2)),
        panel.spacing = unit(0.24, "cm"),
        legend.key = element_rect(fill = "white", colour = NA),
        legend.background = element_blank(),
        legend.position = "right", legend.direction = "vertical",
        legend.box = "vertical", legend.box.just = "left",
        legend.title = element_text(face = "bold", size = 12, colour = text_col),
        legend.text = element_text(size = 12, colour = text_col),
        legend.key.height = unit(3.2, "mm"),
        legend.key.width = unit(5.5, "mm"),
        legend.spacing.y = unit(0.6, "mm"),
        legend.box.spacing = unit(1.5, "mm"),
        legend.margin = margin(0, 0, 0, 3),
        plot.margin = margin(5))
p_network
# ---------- 6. 统计图 ----------
theme_stats <- theme_classic(base_size = 8, base_family = plot_font) +
  theme(text = element_text(family = plot_font),
    axis.line = element_line(linewidth = 0.35, colour = text_col),
    axis.ticks = element_line(linewidth = 0.35, colour = text_col),
    axis.title = element_text(family = plot_font, size = 8, colour = text_col),
    axis.text = element_text(family = plot_font, size = 8, colour = text_col),
    strip.background = element_blank(),
    strip.text = element_text(family = plot_font, size = 8, face = "bold", colour = text_col),
    legend.key = element_blank(),
    legend.title = element_text(),
    legend.text = element_text(),
    legend.key.height = unit(3.5, "mm"),
    panel.grid = element_blank(),
    plot.tag = element_text(family = plot_font, size = 8, face = "bold"),
    plot.margin = margin(4, 5, 4, 4))

p_composition <- ggplot(
  zone_interaction, aes(Zone, pct, fill = Interaction)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
  facet_wrap(~Group, nrow = 1) +
  scale_fill_manual(
    values = interaction_cols, breaks = interaction_levels,
    labels = interaction_labels, name = "Interaction", drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 25),
    labels = scales::label_percent(scale = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion of directed links") +
  theme_stats +theme(axis.text.x = element_text(angle = 25, hjust = 1))

p_site_links <- ggplot(site_network, aes(Zone, n, colour = Zone)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.48,
               linewidth = 0.45, show.legend = FALSE) +
  geom_point(aes(fill = Zone), shape = 21, size = 2.4, stroke = 0.35,
             position = position_jitter(width = 0.075, height = 0),
             show.legend = FALSE) +
  facet_wrap(~Group, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = zone_cols) +
  scale_fill_manual(values = zone_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.14))) +
  labs(x = NULL, y = "Directed links per site") +
  theme_stats +theme(axis.text.x = element_text(angle = 25, hjust = 1))

wilcox_labels <- overall_wilcox %>%
  left_join(stat_data %>%group_by(Group) %>%
              summarise(y = max(Std_CCM_strength) + 0.065, .groups = "drop"),
            by = "Group") %>%
  mutate(x = 1.5,label = if_else(
    p < 0.001, "Wilcoxon p < 0.001", sprintf("Wilcoxon p = %.3f", p)))

p_strength <- ggplot(
  stat_data, aes(Interaction, Std_CCM_strength, fill = Interaction)) +
  geom_violin(width = 0.78, scale = "width", trim = TRUE,
              alpha = 0.22, colour = NA, show.legend = FALSE) +
  geom_boxplot(
    width = 0.22, linewidth = 0.38, outlier.shape = NA,
    colour = text_col, alpha = 0.92, show.legend = FALSE) +
  geom_jitter(
    width = 0.10, size = 1, alpha = 0.5,
    colour = text_col, show.legend = FALSE) +
  geom_text(
    data = wilcox_labels, aes(x, y, label = label),
    inherit.aes = FALSE, size = 2.55, family = plot_font) +
  facet_wrap(~Group, nrow = 1) +
  scale_fill_manual(values = interaction_cols) +
  scale_x_discrete(labels = interaction_labels) +
  scale_y_continuous(limits = c(0, max(wilcox_labels$y) + 0.035),
                     expand = expansion(mult = c(0.01, 0))) +
  labs(x = NULL, y = "Standardized CCM rho") +
  theme_stats + theme(axis.text.x = element_text(angle = 18, hjust = 1))

p_zone_strength <- zone_interaction %>%
  filter(n > 0) %>%
  ggplot(aes(Zone, mean_rho, colour = Interaction)) +
  geom_errorbar(
    aes(ymin = pmax(0, mean_rho - sd_rho),
        ymax = pmin(1, mean_rho + sd_rho)),
    width = 0.12, linewidth = 0.45,
    position = position_dodge(width = 0.34)) +
  geom_point(
    aes(fill = Interaction), shape = 21, size = 2.35, stroke = 0.4,
    position = position_dodge(width = 0.34)) +
  facet_wrap(~Group, nrow = 1) +
  scale_colour_manual(values = interaction_cols, guide = "none") +
  scale_fill_manual(values = interaction_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 0.85), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Mean standardized CCM rho (mean +/- SD)") +
  theme_stats +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p_statistics <- (p_composition | p_site_links) /
  (p_strength | p_zone_strength) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(text = element_text(family = plot_font),legend.position = "right")
# ---------- 7. 输出 ----------
output_dir <- "D:/相互作用/CCM_network_outputs"
dir.create(output_dir, showWarnings = FALSE)
iwalk(stat_results, ~
        readr::write_csv(.x, file.path(output_dir, paste0(.y, ".csv"))))
write_csv(stat_data, file.path(output_dir, "figure_source_data.csv"))
write_csv(edge_check, file.path(output_dir, "edge_check.csv"))
walk(list(quality_check, overall_network, overall_interaction, overall_wilcox,
          zone_network, zone_interaction, edge_check),print)

p_network
p_statistics
write.csv(interactions,"D:/相互作用/interaction.csv")






# 1. 计算每个样点内浮游植物/浮游动物的列宽。
mosaic_group <- stat_data %>%
  count(Site, Group, name = "group_n") %>%
  group_by(Site) %>%
  arrange(Group, .by_group = TRUE) %>%
  mutate(
    site_n = sum(group_n),
    group_share = group_n / site_n,
    xmin = lag(cumsum(group_share), default = 0),
    xmax = cumsum(group_share),
    xmid = (xmin + xmax) / 2,
    Group_label = as.character(Group)) %>%
  ungroup() %>%
  left_join(dplyr::select(site_meta, Site, Site_panel), by = "Site")

# 2. 在每个生物组内部按竞争/协同比例纵向切分。
mosaic_data <- stat_data %>%
  group_by(Site, Group, Interaction) %>%
  summarise(
    N = n(),
    MedianIS_mean = mean(MedianIS, na.rm = TRUE),
    .groups = "drop") %>%
  left_join(
    dplyr::select(
      mosaic_group, Site, Group, Site_panel, site_n,
      group_n, group_share, xmin, xmax, xmid),
    by = c("Site", "Group")) %>%
  group_by(Site, Group) %>%
  arrange(Interaction, .by_group = TRUE) %>%
  mutate(
    within_group_prop = N / group_n,
    ymin = lag(cumsum(within_group_prop), default = 0),
    ymax = cumsum(within_group_prop),
    ymid = (ymin + ymax) / 2,
    tile_area = (xmax - xmin) * (ymax - ymin),
    expected_area = N / site_n) %>%
  ungroup()

# 面积核验：tile_area 应等于 n / site_n（仅允许浮点误差）。
mosaic_area_check <- mosaic_data %>%
  summarise(max_area_error = max(abs(tile_area - expected_area), na.rm = TRUE))

smap_fill_limit <- max(abs(mosaic_data$MedianIS_mean), na.rm = TRUE)
if (!is.finite(smap_fill_limit) || smap_fill_limit == 0) smap_fill_limit <- 1

# 用绝对系数的第 75 百分位作为色阶转折点，避免少数极端值使主体颜色发白。
smap_abs_values <- abs(mosaic_data$MedianIS_mean)
smap_colour_knot <- unname(quantile(
  smap_abs_values[is.finite(smap_abs_values) & smap_abs_values > 0],
  probs = 0.75, na.rm = TRUE))
if (!is.finite(smap_colour_knot) || smap_colour_knot <= 0) {
  smap_colour_knot <- 0.35 * smap_fill_limit
}
smap_colour_knot <- min(
  max(smap_colour_knot, 0.08 * smap_fill_limit),
  0.65 * smap_fill_limit)

p_smap_mosaic <- ggplot() +
  geom_rect(
    data = mosaic_data,
    aes(
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      fill = MedianIS_mean),
    colour = "#FCFBF8", linewidth = 0.70) +
  geom_text(
    data = mosaic_group,
    aes(x = xmid, y = 0, label = Group_label),
    inherit.aes = FALSE,
    family = plot_font, fontface = "bold", size = 2.85,
    colour = "#424242", vjust = 1.65) +
  scale_fill_gradientn(
    colours = c(
      "#963F46", "#D7A0A1", "#F1EEE8", "#A8C6C3", "#2E6F73"),
    values = scales::rescale(c(
      -smap_fill_limit, -smap_colour_knot, 0,
      smap_colour_knot, smap_fill_limit)),
    limits = c(-smap_fill_limit, smap_fill_limit),
    oob = scales::squish,
    name = "Mean S-map coefficient") +
  facet_wrap(~Site_panel, nrow = 2, ncol = 4, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0)) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = 0)) +
  coord_cartesian(
    xlim = c(0, 1), ylim = c(0, 1),
    expand = FALSE, clip = "off") +
  labs(
    title = "S-map competition–facilitation mosaic",
    subtitle = paste0(
      "Width = plankton-group frequency; height = interaction share; ",
      "fill = mean MedianIS"),
    x = NULL,
    y = "Proportion within plankton group") +
  guides(fill = guide_colourbar(
    title.position = "top", title.hjust = 0.5,
    barwidth = unit(48, "mm"), barheight = unit(3.2, "mm"),
    frame.colour = "#B8B3AA", ticks.colour = "#77736C")) +
  theme_minimal(base_size = 10, base_family = plot_font) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5,
                              colour = "#292929", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 9.2, hjust = 0.5, colour = "#68645E",
                                 margin = margin(b = 6)),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "#FCFBF8", colour = NA),
    panel.border = element_rect(colour = "#C8C3B8", fill = NA, linewidth = 0.45),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10.8, colour = "#343434",
                              margin = margin(2, 2, 2, 2)),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.y = element_text(face = "bold", size = 9.8, colour = "#424242"),
    axis.title.x = element_blank(),
    panel.spacing.x = unit(0.35, "cm"),
    panel.spacing.y = unit(1.9, "lines"),
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.title = element_text(face = "bold", size = 9.8, colour = "#343434"),
    legend.text = element_text(size = 8.8, colour = "#56524D"),
    plot.margin = margin(5, 20, 18, 15))

p_smap_mosaic

# 如需导出：
# ggsave("Smap_interaction_mosaic.pdf", p_smap_mosaic,
#        width = 12, height = 7.2, units = "in", device = cairo_pdf)
# ggsave("Smap_interaction_mosaic.png", p_smap_mosaic,
#        width = 12, height = 7.2, units = "in", dpi = 600, bg = "white")
rm(list = ls())
library(tidyverse);library(sf);library(cowplot)
#
setwd("D:/taihu2000_2019")
load("smap.RData")

# ---------- 1. 数据与分组 ----------
phy_vars <- c(
  "Chlorophyta", "Bacillariophyta", "Chrysophyta", "Cyanophyta",
  "Euglenophyta", "Cryptophyta", "Pyrrophyta", "Xanthophyta")
zoo_vars <- c("Protozoa", "Rotifera", "Cladocera", "Copepoda")
biotic_vars <- c(phy_vars, zoo_vars)
all_sites <- c("THL00", "THL01", "THL03", "THL04",
               "THL05", "THL06", "THL07", "THL08")
interaction_levels <- c("Competition/Inhibition", "Facilitation/Promotion")

group_of <- function(x) case_when(
  x %in% phy_vars ~ "Phytoplankton",
  x %in% zoo_vars ~ "Zooplankton",
  TRUE ~ NA_character_)

smap_detail <- map_dfr(seq_len(nrow(smapc.mat.orig0)), ~
                         smapc.mat.orig0$int.ccm.smap[[.x]] %>%
                         as_tibble() %>%
                         mutate(Site = smapc.mat.orig0$Site[.x]) %>%
                         separate(ccmLink, c("from", "to"), sep = "->", remove = FALSE) %>%
                         group_by(Site, from, to, ccmLink) %>%
                         summarise(
                           MedianIS = median(MedianIS, na.rm = TRUE),
                           Perc.Pos.int = median(Perc.Pos.int, na.rm = TRUE),
                           Willcox.conv = first(Willcox.conv), .groups = "drop"))

interactions <- linkM.std0000 %>%
  filter(Period == "2007-2019", Convergence,
         Cause %in% biotic_vars, Effect %in% biotic_vars,
         group_of(Cause) == group_of(Effect), Cause != Effect,
         !is.na(Std_L_strength), Std_L_strength > 0) %>%
  transmute(Site, from = Cause, to = Effect, Group = group_of(from),
            Std_CCM_strength = Std_L_strength,
            ccmLink = paste(from, to, sep = "->")) %>%
  inner_join(smap_detail, by = c("Site", "from", "to", "ccmLink")) %>%
  filter(!is.na(MedianIS)) %>%
  mutate(
    Interaction = case_when(
      Willcox.conv & MedianIS > 0 & Perc.Pos.int > 50 ~ "Facilitation/Promotion",
      Willcox.conv & MedianIS < 0 & Perc.Pos.int < 50 ~ "Competition/Inhibition"),
    Interaction = factor(Interaction, levels = interaction_levels)) %>%
  filter(!is.na(Interaction))

# ---------- 2. 字体、配色与尺度 ----------
font <- "Times New Roman";font_size <- 14
text_size <- font_size/ggplot2::.pt
ink <- "#252525"
width_mm <- 240;height_mm <- 150
map.dir <- "C:/Users/A/Taihu_sampling_map"

cols <- c(
  "#DEB99B", "#5ECC6D", "#5DAFD9", "#7ED1E4",
  "#EA9527", "#F16E1D", "#6E4821", "#A4B423",
  "#C094DF", "#DC95D8", "#326530", "#50C0C9")

paper.cols <- c(
  red = "#DD6B6B", blue = "#527EA6", purple = "#7E74A8")

node.fill <- setNames(cols, biotic_vars)
interaction_cols <- setNames(
  unname(paper.cols[c("red", "blue")]), interaction_levels)

edge_range <- c(0.08, 0.85)
node_range <- c(1.4, 5.0)
rho.max <- max(interactions$Std_CCM_strength)

theme_set(theme_void(base_family = font, base_size = font_size))
cowplot::set_null_device("agg")
small.theme <- theme(
  text = element_text(family = font, size = font_size, colour = ink),
  legend.position = "none", plot.margin = margin(0),
  plot.background = element_rect(fill = "white", colour = NA))

# ---------- 3. 网络布局 ----------
make_layout <- function(species, group, centre_x, radius){
  tibble(name = species, Group = group, node_id = seq_along(species)) %>%
    mutate(angle = pi/2 - 2*pi*(node_id-1)/length(species),
           x = centre_x + radius*cos(angle), y = radius*sin(angle))}

node_layout <- bind_rows(
  make_layout(phy_vars, "Phytoplankton", -0.70, 0.60),
  make_layout(zoo_vars, "Zooplankton", 0.88, 0.44))
node_xy <- dplyr::select(node_layout, name, x, y)

node_degree <- bind_rows(
  transmute(interactions, Site, name = from),
  transmute(interactions, Site, name = to)) %>%
  count(Site, name, name = "Degree")

plot_nodes <- crossing(Site = all_sites, name = biotic_vars) %>%
  left_join(node_layout, by = "name") %>%
  left_join(node_degree, by = c("Site", "name")) %>%
  mutate(Degree = replace_na(Degree, 0))

degree.max <- max(plot_nodes$Degree)
degree_size <- function(x){
  sqrt(scales::rescale(x, to = node_range^2, from = c(0, degree.max)))}
plot_nodes <- plot_nodes %>% mutate(NodeSize = degree_size(Degree))

plot_edges <- interactions %>%
  left_join(rename(node_xy, from = name, x_from = x, y_from = y), by = "from") %>%
  left_join(rename(node_xy, to = name, x_to = x, y_to = y), by = "to") %>%
  mutate(pair_from = pmin(from, to), pair_to = pmax(from, to)) %>%
  group_by(Site, Group, pair_from, pair_to) %>%
  mutate(
    Bidirectional = n() > 1,
    Edge_shape = case_when(
      !Bidirectional ~ "Straight", from < to ~ "Positive", TRUE ~ "Negative")) %>%
  ungroup() %>%
  mutate(
    distance = sqrt((x_to-x_from)^2 + (y_to-y_from)^2),
    x_start = x_from + 0.15*(x_to-x_from)/distance,
    y_start = y_from + 0.15*(y_to-y_from)/distance,
    x_end = x_to - 0.18*(x_to-x_from)/distance,
    y_end = y_to - 0.18*(y_to-y_from)/distance) %>%
  arrange(Std_CCM_strength)

edge.mapping <- aes(
  x = x_start, y = y_start, xend = x_end, yend = y_end,
  colour = Interaction, linewidth = Std_CCM_strength)

edge_layer <- function(data, curvature = NULL){
  args <- list(
    data = data, mapping = edge.mapping, lineend = "round", alpha = 1,
    arrow = grid::arrow(length = grid::unit(0.55, "mm"), type = "closed"),
    show.legend = FALSE)
  if(is.null(curvature)) do.call(geom_segment, args) else
    do.call(geom_curve, c(args, list(curvature = curvature)))}

make_site_network <- function(site){
  edges <- filter(plot_edges, Site == site)
  nodes <- filter(plot_nodes, Site == site)
  ggplot() +
    edge_layer(filter(edges, Edge_shape == "Straight")) +
    edge_layer(filter(edges, Edge_shape == "Positive"), 0.23) +
    edge_layer(filter(edges, Edge_shape == "Negative"), -0.23) +
    geom_point(data = nodes, aes(x, y, size = NodeSize, fill = name),
               shape = 21, colour = "white", stroke = 0.25, alpha = 1) +
    geom_point(data = filter(nodes, Degree == 0), aes(x, y, size = NodeSize),
               shape = 21, fill = "white", colour = "#C6CBCE", stroke = 0.23) +
    scale_colour_manual(values = interaction_cols) +
    scale_linewidth_continuous(range = edge_range, limits = c(0, rho.max)) +
    scale_fill_manual(values = node.fill) + scale_size_identity() +
    coord_fixed(xlim = c(-1.57, 1.60), ylim = c(-0.82, 0.84),
                expand = FALSE, clip = "off") +
    labs(title = site) + small.theme +
    theme(plot.title = element_text(
      family = font, size = font_size, face = "bold",
      hjust = 0.5, margin = margin(b = 2)))}

network.plots <- setNames(lapply(all_sites, make_site_network), all_sites)

drawn.n <- sapply(network.plots, function(p){
  sum(vapply(p$layers, function(z){
    if(is.data.frame(z$data) && "ccmLink" %in% names(z$data)) nrow(z$data) else 0L
  }, integer(1)))})
stopifnot(all(drawn.n == sapply(all_sites, function(z){
  sum(interactions$Site == z)})))

# ---------- 4. 中央采样点图 ----------
sites <- read_csv(file.path(map.dir, "Taihu_sites.csv"), show_col_types = FALSE)
lake <- st_read(file.path(map.dir, "data/taihu.gmt"), quiet = TRUE) %>%
  filter(id == "relation/1126533") %>% st_make_valid() %>% st_transform(4326)
stopifnot(nrow(lake) == 1)

# 在投影坐标中构建 10 km 比例尺
bar.xy <- st_sfc(st_point(c(119.93, 30.91)), crs = 4326) %>%
  st_transform(32651) %>% st_coordinates()
bar <- tibble(x = bar.xy[1,1] + c(0, 10000), y = bar.xy[1,2]) %>%
  st_as_sf(coords = c("x", "y"), crs = 32651) %>% st_transform(4326) %>%
  st_coordinates() %>% as_tibble()

p.map <- ggplot() +
  geom_sf(data = lake, fill = "grey", colour = "black",
          linewidth = 0.5, alpha = 0.5) +
  geom_point(data = sites, aes(Longitude, Latitude), shape = 21,
             size = 3, fill = unname(paper.cols["red"]),
             colour = "white", stroke = 0.20) +
  annotate("text", x = 120.095, y = 31.455, label = "Meiliang\nBay",
           family = font, size = text_size,
           colour = unname(paper.cols["red"]), lineheight = 0.95) +
  annotate("segment", x = bar$X[1], xend = bar$X[2],
           y = bar$Y[1], yend = bar$Y[2], linewidth = 0.35) +
  geom_segment(data = bar, aes(x = X, xend = X, y = Y, yend = Y + 0.009),
               linewidth = 0.35) +
  annotate("text", x = mean(bar$X), y = bar$Y[1] - 0.033, label = "10 km",
           family = font, size = text_size) +
  coord_sf(crs = 32651, default_crs = st_crs(4326),
           xlim = c(119.85, 120.69), ylim = c(30.84, 31.64), expand = FALSE) +
  small.theme

# ---------- 5. 图块位置与连接线 ----------
boxes <- tribble(
  ~Site, ~x, ~y, ~w, ~h, ~anchor.x, ~anchor.y,
  "THL00", 0.050, 0.735, 0.232, 0.245, 0.175, 0.735,
  "THL01", 0.260, 0.735, 0.232, 0.245, 0.420, 0.735,
  "THL03", 0.505, 0.735, 0.232, 0.245, 0.615, 0.735,
  "THL04", 0.718, 0.735, 0.232, 0.245, 0.810, 0.735,
  "THL05", 0.050, 0.490, 0.232, 0.235, 0.282, 0.600,
  "THL06", 0.050, 0.255, 0.232, 0.235, 0.282, 0.365,
  "THL07", 0.718, 0.490, 0.232, 0.235, 0.718, 0.600,
  "THL08", 0.718, 0.255, 0.232, 0.235, 0.718, 0.365)

map.ranges <- ggplot_build(p.map)$layout$panel_params[[1]]
map.x <- 0.355;map.y <- 0.230;map.w <- 0.290
map.h <- map.w*width_mm/height_mm*
  diff(map.ranges$y_range)/diff(map.ranges$x_range)

station.xy <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(32651) %>% st_coordinates()
leaders <- sites %>%
  mutate(
    x = map.x + map.w*(station.xy[,1]-map.ranges$x_range[1])/diff(map.ranges$x_range),
    y = map.y + map.h*(station.xy[,2]-map.ranges$y_range[1])/diff(map.ranges$y_range)) %>%
  left_join(dplyr::select(boxes, Site, anchor.x, anchor.y), by = "Site")

leaders.direct <- filter(leaders, Site != "THL06")
leader.06 <- filter(leaders, Site == "THL06")
leader.06.path <- tibble(
  x = c(leader.06$x, 0.374, leader.06$anchor.x),
  y = c(leader.06$y, leader.06$y, leader.06$anchor.y))

# ---------- 6. 类群图例 ----------
taxa.key <- tibble(name = biotic_vars, TaxonID = seq_along(biotic_vars)) %>%
  mutate(
    x = if_else(TaxonID <= 8, 0.02 + ((TaxonID-1) %% 4)*0.17,
                0.71 + ((TaxonID-9) %% 2)*0.16),
    y = if_else(TaxonID %in% c(1:4, 9:10), 0.45, 0.12))

p.taxa <- ggplot(taxa.key) +
  geom_point(aes(x, y, fill = name),
             shape = 21, size = 3, colour = "white",
             stroke = 0.25, alpha = 1) +
  geom_text(aes(x = x+0.021, y = y, label = name), hjust = 0,
            family = font, size = text_size, colour = ink) +
  annotate("text", x = 0.02, y = 0.86, label = "Phytoplankton", hjust = 0,
           family = font, fontface = "bold", size = text_size,
           colour = unname(paper.cols["blue"])) +
  annotate("text", x = 0.71, y = 0.86, label = "Zooplankton", hjust = 0,
           family = font, fontface = "bold", size = text_size,
           colour = unname(paper.cols["purple"])) +
  scale_fill_manual(values = node.fill) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1),
                  expand = FALSE, clip = "off") +
  small.theme

# ---------- 7. 公共图例 ----------
rho.key <- tibble(rho = c(0.2, 0.4, 0.6, 0.8),
                  x = c(0.345, 0.430, 0.515, 0.600)) %>%
  mutate(width = scales::rescale(rho, to = edge_range, from = c(0, rho.max)))
degree.key <- tibble(degree = c(0, 5, 10, degree.max),
                     x = c(0.755, 0.818, 0.881, 0.944)) %>%
  mutate(size = degree_size(degree))

p.key <- ggplot() +
  annotate("segment", x = 0.012, xend = 0.042, y = 0.72, yend = 0.72,
           colour = unname(interaction_cols[1]), linewidth = 0.36,
           arrow = grid::arrow(length = grid::unit(0.75, "mm"), type = "closed")) +
  annotate("text", x = 0.055, y = 0.72, label = "Negative net effect", hjust = 0,
           family = font, size = text_size) +
  annotate("segment", x = 0.012, xend = 0.042, y = 0.25, yend = 0.25,
           colour = unname(interaction_cols[2]), linewidth = 0.36,
           arrow = grid::arrow(length = grid::unit(0.75, "mm"), type = "closed")) +
  annotate("text", x = 0.055, y = 0.25, label = "Positive net effect", hjust = 0,
           family = font, size = text_size) +
  annotate("text", x = 0.345, y = 0.88, label = "Standardised CCM \u03c1", hjust = 0,
           family = font, size = text_size) +
  geom_segment(data = rho.key,
               aes(x, xend = x+0.027, y = 0.36, yend = 0.36, linewidth = width),
               colour = "#707E85") +
  geom_text(data = rho.key, aes(x = x+0.035, y = 0.36, label = rho),
            family = font, size = text_size, hjust = 0) +
  annotate("text", x = 0.735, y = 0.88, label = "Total degree", hjust = 0,
           family = font, size = text_size) +
  geom_point(data = degree.key, aes(x, y = 0.36, size = size),
             shape = 21, fill = "#7D8B92", colour = "white", stroke = 0.25) +
  geom_point(data = filter(degree.key, degree == 0), aes(x, y = 0.36, size = size),
             shape = 21, fill = "white", colour = "#C6CBCE", stroke = 0.23) +
  geom_text(data = degree.key, aes(x = x+0.014, y = 0.36, label = degree),
            family = font, size = text_size, hjust = 0) +
  scale_size_identity() + scale_linewidth_identity() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
  small.theme

# ---------- 8. 指北针 ----------
compass.points <- map_dfr(0:3, function(k){
  angle <- k*pi/2
  tibble(x = c(0, -0.17, 0, 0, 0.17, 0),
         y = c(0, 0.25, 1, 0, 0.25, 1),
         part = rep(c("dark", "light"), each = 3)) %>%
    mutate(x.rot = x*cos(angle)-y*sin(angle),
           y.rot = x*sin(angle)+y*cos(angle), id = paste(k, part))})
compass.ring <- tibble(angle = seq(0, 2*pi, length.out = 121)) %>%
  mutate(x = 0.62*cos(angle), y = 0.62*sin(angle))
compass.labels <- tribble(
  ~x, ~y, ~label,
  0, 1.34, "N", 1.34, 0, "E",
  0, -1.34, "S", -1.34, 0, "W")

p.compass <- ggplot() +
  geom_path(data = compass.ring, aes(x, y), linewidth = 0.5, colour = ink) +
  geom_polygon(data = compass.points,
               aes(x.rot, y.rot, group = id, fill = part),
               colour = ink, linewidth = 0.5) +
  geom_text(data = compass.labels, aes(x, y, label = label),
            family = font, size = text_size, colour = ink) +
  scale_fill_manual(values = c(dark = ink, light = "white")) +
  coord_fixed(xlim = c(-1.70, 1.70), ylim = c(-1.70, 1.70),
              expand = FALSE, clip = "off") +
  small.theme +
  theme(plot.background = element_rect(fill = NA, colour = NA))

# ---------- 9. 拼图 ----------
final.plot <- ggdraw() +
  draw_plot(p.map, x = map.x, y = map.y, width = map.w, height = map.h) +
  geom_segment(data = leaders.direct,
               aes(x, y, xend = anchor.x, yend = anchor.y),
               inherit.aes = FALSE, linewidth = 1, colour = "#C0C8CD") +
  geom_path(data = leader.06.path, aes(x, y),
            inherit.aes = FALSE, linewidth = 1, colour = "#C0C8CD") +
  draw_plot(p.compass, x = 0.666, y = 0.607,
            width = 0.070, height = 0.070*width_mm/height_mm) +
  draw_plot(p.key, x = 0.025, y = 0.140, width = 0.950, height = 0.100) +
  draw_plot(p.taxa, x = 0.025, y = 0.007, width = 0.950, height = 0.125)

for(i in seq_len(nrow(boxes))){
  box <- boxes[i, ]
  final.plot <- final.plot +
    draw_plot(network.plots[[box$Site]],
              x = box$x, y = box$y, width = box$w, height = box$h)}

final.plot <- final.plot + theme(aspect.ratio = height_mm/width_mm)
print(final.plot)
