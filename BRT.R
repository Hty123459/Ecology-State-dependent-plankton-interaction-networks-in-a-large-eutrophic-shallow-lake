library(dismo); library(gbm); library(tidyverse)
library(quantreg); library(openxlsx); library(zoo)

set.seed(1234)
setwd("D:/taihu2000_2019")
#source("brt_function.R")
select <- dplyr::select; filter <- dplyr::filter

years <- c(2007, 2019 + 11 / 12)
link_file <- "D:/相互作用/CCM_network_outputs/figure_source_data.csv"
output_dir <- "D:/相互作用/BRT_QR_grouped_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load("brt.orig.Rdata")
smap <- smapc.mat.orig

# 21个候选变量在每个BRT中共同参与竞争
drivers <- c(
  "Chlorophyta", "Bacillariophyta", "Chrysophyta", "Cyanophyta",
  "Euglenophyta", "Cryptophyta", "Pyrrophyta", "Xanthophyta",
  "Protozoa", "Rotifera", "Cladocera", "Copepoda",
  "ZP", "WT", "Wind_speed", "Precip", "CN", "CP", "TP", "TN", "NP")

links <- read_csv(link_file, show_col_types = FALSE) %>%
  transmute(
    Site, ccmLink, Cause = from, Effect = to,
    From_type = Group, To_type = Group, Interaction, MedianIS,
    IS.direction05 = if_else(MedianIS > 0, "Positive", "Negative"),
    Network = Group) %>%
  filter(Network %in% c("Phytoplankton", "Zooplankton")) %>%
  distinct()

stopifnot(
  nrow(links) == 377,
  all(links$Cause != links$Effect),
  !anyDuplicated(links[c("Site", "ccmLink")]))

analyse_network <- function(network) {
  link_data <- links %>% filter(Network == network) %>% select(-Network)
  brt <- run_brt(link_data, smap, drivers, years)
  brt$models <- NULL
  tag <- function(x) bind_cols(Network = rep(network, nrow(x)), x)
  brt[c("results", "important", "skipped")] <-
    map(brt[c("results", "important", "skipped")], tag)
  qr <- brt$important %>%
    fit_qr(smap, years) %>%
    classify_qr(smap, years, dominance = .95) %>%
    clean_qr(keep = 1:9) %>%
    tag()
  list(brt = brt, qr = qr, expected = nrow(link_data))}

write_network <- function(x, network) {
  brt <- x$brt; qr <- x$qr
  qr_tail <- qr %>% filter(both_tail_sig)
  qr_three <- qr %>% filter(all_three_sig)
  run_check <- tibble(
    Network = network,
    expected_links = x$expected,
    successful_brt = nrow(brt$results),
    skipped_brt = nrow(brt$skipped),
    accounted_links = successful_brt + skipped_brt,
    complete = accounted_links == expected_links,
    candidate_predictors = length(drivers),
    importance_threshold_pct = round(100 / length(drivers), 4),
    important_predictor_rows = nrow(brt$important),
    qr_rows = nrow(qr),
    qr_tail_significant_rows = nrow(qr_tail),
    qr_all_three_significant_rows = nrow(qr_three))
  tabs <- list(
    README = tibble(
      Item = c("Response network", "Model scope", "Candidate predictors",
               "BRT importance rule", "Time window", "Minimum observations", "Random seed"),
      Value = c(network,
        "All 21 predictors entered jointly; no predictor-group stratification",
        paste(drivers, collapse = ", "),
        sprintf("Relative influence > 100/21 = %.4f%%", 100 / length(drivers)),
        "2007-01 to 2019-12", "20", "1234")),
    BRT_Results = brt$results,
    BRT_Important = brt$important,
    QR_All = qr,
    QR_TailSig = qr_tail,
    QR_AllThreeSig = qr_three,
    Skipped = brt$skipped,
    Run_Check = run_check,
    Summary = qr_summary(qr) %>% mutate(Network = network, .before = 1),
    Variables = tibble(Predictor = drivers))
  write.xlsx(tabs,
    file.path(output_dir, paste0("BRT_QR_", network, ".xlsx")),
    overwrite = TRUE)}

c("Phytoplankton", "Zooplankton") %>%
  set_names() %>%
  map(analyse_network) %>%
  iwalk(write_network)

