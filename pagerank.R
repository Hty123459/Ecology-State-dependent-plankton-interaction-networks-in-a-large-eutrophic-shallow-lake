rm(list = ls())

library(tidyverse)
library(igraph)
library(patchwork)

filter <- dplyr::filter
summarise <- dplyr::summarise
crossing <- tidyr::crossing

setwd("D:/taihu2000_2019")
load("smap.RData")
set.seed(1234)

# ---------- 1. Variables ----------
phy_vars <- c(
  "Chlorophyta", "Bacillariophyta", "Chrysophyta", "Cyanophyta",
  "Euglenophyta", "Cryptophyta", "Pyrrophyta", "Xanthophyta"
)

zoo_vars <- c("Protozoa", "Rotifera", "Cladocera", "Copepoda")
biotic_vars <- c(phy_vars, zoo_vars)

interaction_levels <- c(
  "Competition/Inhibition",
  "Facilitation/Promotion"
)

group_of <- function(x) {
  case_when(
    x %in% phy_vars ~ "Phytoplankton",
    x %in% zoo_vars ~ "Zooplankton",
    TRUE ~ NA_character_
  )
}

# ---------- 2. Extract intragroup biotic interactions ----------
smap_detail <- map_dfr(seq_len(nrow(smapc.mat.orig0)), ~
                         smapc.mat.orig0$int.ccm.smap[[.x]] %>%
                         as_tibble() %>%
                         mutate(Site = smapc.mat.orig0$Site[.x]) %>%
                         separate(ccmLink, c("from", "to"), sep = "->", remove = FALSE) %>%
                         group_by(Site, from, to, ccmLink) %>%
                         summarise(
                           MedianIS = median(MedianIS, na.rm = TRUE),
                           Perc.Pos.int = median(Perc.Pos.int, na.rm = TRUE),
                           Willcox.conv = first(Willcox.conv),
                           .groups = "drop"
                         )
)

interactions <- linkM.std0000 %>%
  filter(
    Period == "2007-2019",
    Convergence,
    Cause %in% biotic_vars,
    Effect %in% biotic_vars,
    group_of(Cause) == group_of(Effect),
    Cause != Effect,
    !is.na(Std_L_strength),
    Std_L_strength > 0
  ) %>%
  transmute(
    Site,
    from = Cause,
    to = Effect,
    Group = group_of(from),
    Std_CCM_strength = Std_L_strength,
    ccmLink = paste(from, to, sep = "->")
  ) %>%
  inner_join(
    smap_detail,
    by = c("Site", "from", "to", "ccmLink")
  ) %>%
  filter(!is.na(MedianIS)) %>%
  mutate(
    Interaction = case_when(
      Willcox.conv & MedianIS > 0 & Perc.Pos.int > 50 ~
        "Facilitation/Promotion",
      Willcox.conv & MedianIS < 0 & Perc.Pos.int < 50 ~
        "Competition/Inhibition"
    ),
    Interaction = factor(Interaction, levels = interaction_levels)
  ) %>%
  filter(!is.na(Interaction))

site_order <- c("THL00", "THL01", "THL03", "THL04", "THL05", "THL06", "THL07", "THL08")
all_sites <- site_order[site_order %in% unique(interactions$Site)]

site_meta <- tibble(Site = all_sites) %>%
  mutate(
    Zone = case_when(
      Site %in% c("THL00", "THL01", "THL03", "THL04", "THL05", "THL06") ~ "Bay",
      Site == "THL07" ~ "Estuary",
      Site == "THL08" ~ "Lake centre",
      TRUE ~ NA_character_
    ),
    Zone = factor(Zone, levels = c("Bay", "Estuary", "Lake centre")),
    Site_panel = factor(Site, levels = site_order)
  )

# ---------- 3. PageRank calculation ----------
make_group_graph <- function(d, node_set) {
  adj <- matrix(
    0,
    nrow = length(node_set),
    ncol = length(node_set),
    dimnames = list(node_set, node_set)
  )
  
  if (nrow(d) > 0) {
    for (i in seq_len(nrow(d))) {
      adj[d$from[i], d$to[i]] <- d$Std_CCM_strength[i]
    }
  }
  
  list(
    graph = graph_from_adjacency_matrix(
      adj,
      mode = "directed",
      weighted = TRUE,
      diag = FALSE
    ),
    adj = adj
  )
}

pagerank_data <- crossing(
  Site = all_sites,
  Group = c("Phytoplankton", "Zooplankton")
) %>%
  group_by(Site, Group) %>%
  group_modify(~ {
    node_set <- if (.y$Group == "Phytoplankton") phy_vars else zoo_vars
    
    d <- interactions %>%
      filter(Site == .y$Site, Group == .y$Group)
    
    graph_obj <- make_group_graph(d, node_set)
    g <- graph_obj$graph
    adj <- graph_obj$adj
    
    pr <- page_rank(g, weights = E(g)$weight)$vector
    
    g_rev <- graph_from_adjacency_matrix(
      t(adj),
      mode = "directed",
      weighted = TRUE,
      diag = FALSE
    )
    
    cr <- page_rank(g_rev, weights = E(g_rev)$weight)$vector
    
    tibble(
      name = names(pr),
      PageRank = as.numeric(pr),
      CheiRank = as.numeric(cr),
      Integrated_Centrality = sqrt(PageRank * CheiRank),
      PageRank_rank = rank(-PageRank, ties.method = "min"),
      CheiRank_rank = rank(-CheiRank, ties.method = "min"),
      Integrated_rank = rank(-Integrated_Centrality, ties.method = "min")
    )
  }) %>%
  ungroup()

pagerank_summary <- pagerank_data %>%
  group_by(Site, Group) %>%
  summarise(
    mean_PageRank = mean(PageRank, na.rm = TRUE),
    mean_CheiRank = mean(CheiRank, na.rm = TRUE),
    mean_Integrated_Centrality = mean(Integrated_Centrality, na.rm = TRUE),
    top_PageRank_nodes = paste(
      name[order(PageRank, decreasing = TRUE)][seq_len(min(3, n()))],
      collapse = ", "
    ),
    top_CheiRank_nodes = paste(
      name[order(CheiRank, decreasing = TRUE)][seq_len(min(3, n()))],
      collapse = ", "
    ),
    top_Integrated_nodes = paste(
      name[order(Integrated_Centrality, decreasing = TRUE)][seq_len(min(3, n()))],
      collapse = ", "
    ),
    .groups = "drop"
  )

from_centrality <- pagerank_data %>%
  dplyr::select(
    Site, Group, from = name,
    from_PageRank = PageRank,
    from_CheiRank = CheiRank,
    from_Integrated_Centrality = Integrated_Centrality,
    from_PageRank_rank = PageRank_rank,
    from_CheiRank_rank = CheiRank_rank,
    from_Integrated_rank = Integrated_rank
  )

to_centrality <- pagerank_data %>%
  dplyr::select(
    Site, Group, to = name,
    to_PageRank = PageRank,
    to_CheiRank = CheiRank,
    to_Integrated_Centrality = Integrated_Centrality,
    to_PageRank_rank = PageRank_rank,
    to_CheiRank_rank = CheiRank_rank,
    to_Integrated_rank = Integrated_rank
  )

interactions_with_pagerank <- interactions %>%
  left_join(from_centrality, by = c("Site", "Group", "from")) %>%
  left_join(to_centrality, by = c("Site", "Group", "to")) %>%
  mutate(
    Edge_Integrated_Centrality =
      sqrt(from_Integrated_Centrality * to_Integrated_Centrality)
  )

# ---------- 4. Combined PageRank-CheiRank scatter plot ----------
taxon_cols <- c(
  Chlorophyta = "#267449",
  Bacillariophyta = "#2C8390",
  Chrysophyta = "#D1A34C",
  Cyanophyta = "#285B88",
  Euglenophyta = "#D9915B",
  Cryptophyta = "#896D9D",
  Pyrrophyta = "#B75B58",
  Xanthophyta = "#8A9655",
  Protozoa = "#456F79",
  Rotifera = "#D9915B",
  Cladocera = "#267449",
  Copepoda = "#285B88"
)

plot_font <- "serif"
text_col <- "#111111"

pagerank_plot_data <- pagerank_data %>%
  left_join(
    dplyr::select(site_meta, Site, Zone, Site_panel),
    by = "Site"
  ) %>%
  mutate(
    Group = factor(Group, levels = c("Phytoplankton", "Zooplankton")),
    name = factor(name, levels = biotic_vars),
    Site_panel = factor(Site_panel, levels = site_order)
  ) %>%
  arrange(desc(Integrated_Centrality))

make_pagerank_facet <- function(d) {
  
  x_lim <- range(d$PageRank, na.rm = TRUE)
  y_lim <- range(d$CheiRank, na.rm = TRUE)
  x_pad <- max(diff(x_lim) * 0.12, 0.012)
  y_pad <- max(diff(y_lim) * 0.12, 0.012)
  size_limits <- range(d$Integrated_Centrality, na.rm = TRUE)
  
  ggplot(d, aes(PageRank, CheiRank)) +
    geom_point(
      aes(size = Integrated_Centrality, fill = name, shape = Group),
      colour = "white",
      stroke = 0.25,
      alpha = 0.94
    ) +
    facet_wrap(~Site_panel, ncol = 4) +
    scale_shape_manual(
      values = c(Phytoplankton = 21, Zooplankton = 23),
      guide = "none"
    ) +
    scale_fill_manual(
      values = taxon_cols,
      breaks = biotic_vars,
      limits = biotic_vars,
      drop = FALSE,
      name = "Taxon"
    ) +
    scale_size_continuous(
      range = c(2.7, 6.0),
      limits = size_limits,
      guide = "none"
    ) +
    coord_cartesian(
      xlim = c(x_lim[1] - x_pad, x_lim[2] + x_pad),
      ylim = c(y_lim[1] - y_pad, y_lim[2] + y_pad),
      clip = "off"
    ) +
    labs(
      x = "PageRank",
      y = "CheiRank",
      title = "",
      subtitle = ""
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        ncol = 6,
        byrow = TRUE,
        override.aes = list(
          shape = c(rep(21, length(phy_vars)), rep(23, length(zoo_vars))),
          size = 4.5,
          colour = "white",
          alpha = 1
        )
      )
    ) +
    theme_classic(base_size = 13, base_family = plot_font) +
    theme(
      plot.title = element_text(
        face = "plain",
        size = 16,
        hjust = 0,
        colour = text_col,
        margin = margin(b = 12)
      ),
      plot.subtitle = element_text(
        size = 12,
        colour = text_col,
        margin = margin(b = 9)
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "plain",
        size = 12.5,
        colour = text_col,
        margin = margin(5, 2, 5, 2)
      ),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 0.5, colour = text_col),
      axis.ticks = element_line(linewidth = 0.50, colour = text_col),
      axis.ticks.length = unit(1.7, "mm"),
      axis.title = element_text(size = 13, colour = text_col),
      axis.text = element_text(size = 11, colour = text_col),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "left",
      legend.direction = "horizontal",
      legend.title = element_text(face = "plain", size = 11.5, colour = text_col),
      legend.text = element_text(size = 10, colour = text_col),
      legend.key = element_blank(),
      legend.key.height = unit(5.0, "mm"),
      legend.spacing.x = unit(2.0, "mm"),
      legend.spacing.y = unit(1.3, "mm"),
      legend.box.spacing = unit(3.0, "mm"),
      panel.spacing = unit(0.55, "cm"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(9, 12, 8, 9)
    )
}

p_pagerank_combined <- make_pagerank_facet(pagerank_plot_data)

write.csv(pagerank_plot_data,"D:/相互作用/pagerank.csv")
# ---------- 5. Output objects ----------
pagerank_data
pagerank_summary
interactions_with_pagerank

p_pagerank_combined


save(linkM.std0000,smapc.mat.orig0,file = "D:/相互作用/test.Rdata")
