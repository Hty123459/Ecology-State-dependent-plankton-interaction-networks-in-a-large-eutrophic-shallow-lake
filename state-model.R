library(tidyverse)

# 1. 字体、符号与分位线样式
font <- "Times New Roman"
fs <- 12
ts <- fs / ggplot2::.pt
type_order <- c(1, 3, 5, 7, 2, 4, 6, 8)
time_shapes <- c("t-1" = 16, "t" = 17, "t+1" = 15)
q_alpha <- c("0.05" = .25, "0.50" = 0.65, "0.95" = 1)
q_width <- c("0.05" = .45, "0.50" = 0.65, "0.95" = 1)
q_colour <- "#C40000"

# 2. 概念示意参数：三条分位线均为直线
pars <- tribble(
  ~type, ~tau, ~a, ~b,
  1, "0.05", .22, .55,    1, "0.50", .42, .70,    1, "0.95", .62, .85,
  2, "0.05", 1.25, -.85,  2, "0.50", 1.38, -.65,  2, "0.95", 1.48, -.45,
  3, "0.05", -1.40, .45,  3, "0.50", -1.20, .60,  3, "0.95", -1.00, .75,
  4, "0.05", -.68, -.55,  4, "0.50", -.48, -.70,  4, "0.95", -.28, -.85,
  5, "0.05", -1.35, 2.20, 5, "0.50", -1.10, 2.20, 5, "0.95", -.85, 2.20,
  6, "0.05", .85, -2.20,  6, "0.50", 1.10, -2.20, 6, "0.95", 1.35, -2.20,
  7, "0.05", .30, .65,    7, "0.50", .72, .35,    7, "0.95", 1.38, -.25,
  8, "0.05", -1.38, .25,  8, "0.50", -.72, -.35,  8, "0.95", -.30, -.65
) %>%
  mutate(x = .10, xend = .84, y = a + b * x, yend = a + b * xend,
         tau = factor(tau, levels = c("0.05", "0.95", "0.50"))) %>%
  arrange(type, tau)

# 3. 示意散点，不代表实际拟合结果
set.seed(20260818)
dots <- map_dfr(1:8, function(i) {
  x <- sort(runif(42, .10, .84))
  m <- filter(pars, type == i, tau == "0.50")
  tibble(type = i, x = x, y = m$a + m$b * x + rnorm(42, 0, if_else(i %in% 7:8, .18, .15)),
         time = factor(sample(rep(c("t", "t-1", "t+1"), each = 14)), levels = names(time_shapes)))
})

# 4. 压缩无数据区域，顶部预留 Type 标签空间
bounds <- map_dfr(1:8, function(i) {
  ids <- if (i %% 2 == 1) c(i, i + 1) else c(i - 1, i)
  p <- filter(pars, type %in% ids)
  extent <- max(abs(c(filter(dots, type %in% ids)$y, p$y, p$yend)))
  yl <- if (i %in% c(1, 2, 7)) c(-.025, 1.06) * extent else
    if (i %in% c(3, 4, 8)) c(-1.06, .025) * extent else c(-1.06, 1.06) * extent
  yl[2] <- yl[2] + .18 * diff(yl)
  tibble(type = i, ymin = yl[1], ymax = yl[2], label_y = yl[2] - .035 * diff(yl))
})

# 5. 固定面板顺序，State 仅标注于第一行
pars <- pars %>% mutate(type = factor(type, levels = type_order))
dots <- dots %>% mutate(type = factor(type, levels = type_order))
bounds <- bounds %>% mutate(label = paste("Type", type), type = factor(type, levels = type_order))
state_labels <- c(`1` = "State 1", `3` = "State 2", `5` = "State 3", `7` = "State 4",
                  `2` = "", `4` = "", `6` = "", `8` = "")

# 6. 黑色散点、红色实线；用透明度和线宽突出中位线
fig <- ggplot() +
  geom_blank(data = bounds, aes(x = .055, y = ymin)) +
  geom_blank(data = bounds, aes(x = .885, y = ymax)) +
  geom_hline(yintercept = 0, colour = "#BFC3C6", linewidth = .28, linetype = "22") +
  geom_point(data = dots, aes(x, y, shape = time), colour = "black", size = 1.65) +
  geom_segment(data = pars, aes(x = x, xend = xend, y = y, yend = yend, alpha = tau, linewidth = tau),
               colour = q_colour, linetype = "solid", lineend = "round") +
  geom_text(data = bounds, aes(x = .08, y = label_y, label = label),
            hjust = 0, vjust = 1, family = font, size = ts, colour = "black") +
  facet_wrap(~ type, ncol = 4, scales = "free_y", labeller = labeller(type = state_labels)) +
  scale_shape_manual(name = "Time", values = time_shapes, breaks = names(time_shapes)) +
  scale_alpha_manual(name = "Quantile", values = q_alpha, breaks = names(q_alpha)) +
  scale_linewidth_manual(values = q_width, guide = "none") +
  scale_x_continuous(breaks = c(.10, .84), labels = c("Low", "High"),
                     limits = c(.055, .885), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = 0, expand = expansion(mult = 0)) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 2, colour = "black")),
    alpha = guide_legend(order = 2, override.aes = list(
      colour = q_colour, linewidth = unname(q_width), linetype = "solid"))
  ) +
  labs(x = "Environmental gradient", y = "S-map coefficient") +
  theme_classic(base_family = font, base_size = fs) +
  theme(
    text = element_text(family = font, size = fs, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = fs, face = "bold", margin = margin(t = 0, b = 5)),
    axis.title = element_text(size = fs),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text = element_text(size = fs, colour = "black"),
    axis.line = element_line(colour = "#92999E", linewidth = .25),
    axis.ticks = element_line(colour = "#92999E", linewidth = .25),
    axis.ticks.length = grid::unit(.8, "mm"),
    panel.spacing.x = grid::unit(5, "mm"),
    panel.spacing.y = grid::unit(2, "mm"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = fs),
    legend.text = element_text(size = fs),
    legend.key.width = grid::unit(10, "mm"),
    legend.spacing.x = grid::unit(3, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(t = 3),
    plot.margin = margin(6, 8, 4, 4)
  )

# 7. 仅显示，不保存
print(fig)