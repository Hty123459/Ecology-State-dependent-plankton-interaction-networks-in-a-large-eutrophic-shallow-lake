rm(list = ls())

library(tidyverse)
library(gbm)
library(quantreg)
library(openxlsx)
library(zoo)
library(patchwork)
library(showtext)

set.seed(42)
setwd("D:/taihu2000_2019/")
source("brt_function.R")
showtext_auto(FALSE)

# Settings

years <- c(2007,2019 + 11/12)
input_dir <- "C:/Users/A/brt_qr_grouped_results"
files <- c(
  Phytoplankton = file.path(input_dir,"BRT_QR_Phytoplankton.xlsx"),
  Zooplankton = file.path(input_dir,"BRT_QR_Zooplankton.xlsx"))

data_env <- new.env()
load("brt.orig.Rdata",envir = data_env)
smap <- data_env$smapc.mat.orig
rm(data_env)

exclude <- c("Euglenophyta","Pyrrophyta","Xanthophyta")
fill_cols <- c(Upper = "#e5b5bc",Lower = "#b0c4de")
line_cols <- c(Upper = "#ce1126",Lower = "#145b9d")
time_breaks <- c(2007,2011,2015,2019)

fixed_y_label <- function(x)
  formatC(x,format = "fg",digits = 3,width = 6)

base_theme <- theme_classic(
  base_size = 10,base_family = "Times New Roman")+
  theme(
    aspect.ratio = 0.618,
    axis.text = element_text(size = 9,color = "black"),
    axis.title = element_text(size = 10,color = "black"),
    axis.line = element_line(linewidth = 0.4),
    axis.line.x = element_blank(),
    axis.ticks = element_line(linewidth = 0.35),
    axis.ticks.length = grid::unit(2,"pt"),
    plot.title = element_text(size = 10,face = "bold",hjust = 0),
    legend.title = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8),
    legend.key.width = grid::unit(12,"pt"),
    legend.key.height = grid::unit(8,"pt"),
    legend.margin = margin(0,0,0,0),
    legend.box.spacing = grid::unit(2,"pt"),
    plot.margin = margin(6,8,6,6))

compact_theme <- theme(
  plot.margin = margin(1,1,1,1),
  plot.title = element_text(
    size = 9,face = "bold",hjust = 0,margin = margin(b = 2)),
  plot.title.position = "plot",
  axis.title = element_text(size = 9),
  axis.text = element_text(size = 8,color = "black"),
  legend.text = element_text(size = 7.5),
  legend.title = element_blank(),
  legend.key.width = grid::unit(9,"pt"),
  legend.key.height = grid::unit(7,"pt"),
  legend.spacing.x = grid::unit(1,"pt"),
  legend.margin = margin(0,0,0,0),
  legend.box.margin = margin(0,0,0,0),
  legend.box.spacing = grid::unit(0,"pt"))

# Quantile regression and BRT

make_qr_panel <- function(row,env){
  driver <- row$env_var
  link_used <- find_link(
    row$ccmLink_used,row$Cause,row$Effect,env)
  
  stopifnot(
    !is.na(link_used),
    all(c(driver,link_used)%in%names(env)))
  
  d <- tibble(
    x = env[[driver]],
    y = env[[link_used]])%>%
    filter(is.finite(x),is.finite(y))
  
  stopifnot(nrow(d) > 3,n_distinct(d$x) > 1)
  
  x_grid <- tibble(
    x = seq(min(d$x),max(d$x),length.out = 200))
  
  qr_lines <- map_dfr(c(0.05,0.50,0.95),function(tau){
    fit <- quantreg::rq(y~x,data = d,tau = tau)
    tibble(
      x = x_grid$x,
      y = as.numeric(predict(fit,newdata = x_grid)),
      Quantile = sprintf("%.2f",tau))})
  
  set.seed(42)
  pd <- fit_curve(d)
  yr <- range(d$y)
  pr <- range(pd$y)
  
  if(diff(yr) == 0) yr <- yr + c(-1,1)
  if(diff(pr) == 0) pr <- pr + c(-1,1)
  
  pd$y_plot <- (pd$y-pr[1])/diff(pr)*diff(yr)+yr[1]
  
  brt_axis <- local({
    y0 <- yr[1]
    p0 <- pr[1]
    dy <- diff(yr)
    dp <- diff(pr)
    function(z) (z-y0)/dy*dp+p0})
  
  link_label <- paste(row$Cause,"\u2192",row$Effect)
  
  ggplot()+
    geom_hline(
      yintercept = 0,linetype = "dashed",
      color = "grey55",linewidth = 0.3)+
    geom_point(
      data = d,aes(x,y),color = "black",
      shape = 8,size = 0.8,alpha = 0.4)+
    geom_line(
      data = qr_lines,
      aes(x,y,group = Quantile,alpha = Quantile),
      color = line_cols[["Upper"]],linewidth = 0.6)+
    geom_line(
      data = pd,aes(x,y_plot),
      color = line_cols[["Lower"]],linewidth = 0.7)+
    geom_rug(
      data = d,aes(x = x),sides = "b",
      color = "black",alpha = 0.12,
      length = grid::unit(0.02,"npc"))+
    scale_alpha_manual(
      name = "Quantile",
      values = c("0.05" = 0.30,"0.50" = 0.65,"0.95" = 1),
      breaks = c("0.05","0.50","0.95"))+
    scale_x_continuous(
      name = chem_label(driver),labels = sci_label)+
    scale_y_continuous(
      name = link_label,labels = fixed_y_label,
      sec.axis = sec_axis(
        brt_axis,name = "Predicted response",
        labels = fixed_y_label))+
    guides(alpha = guide_legend(
      override.aes = list(color = line_cols[["Upper"]])))+
    base_theme+
    theme(
      aspect.ratio = 1,
      plot.title = element_text(size = 8,face = "bold"),
      axis.text = element_text(size = 7,color = "black"),
      axis.text.y.left = element_text(
        family = "Times New Roman",margin = margin(r = 2)),
      axis.text.y.right = element_text(
        family = "Times New Roman",margin = margin(l = 2)),
      axis.title = element_text(size = 8),
      axis.title.y.left = element_text(
        size = 8,angle = 90,margin = margin(r = 4)),
      axis.title.y.right = element_text(
        size = 8,angle = 270,margin = margin(l = 4)),
      axis.line.x = element_line(linewidth = 0.35),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      plot.margin = margin(2,1,2,1))}

# Time-series panels

make_triplet <- function(row,env){
  driver <- row$env_var
  link_label <- paste(row$Cause,"\u2192",row$Effect)
  link_used <- find_link(
    row$ccmLink_used,row$Cause,row$Effect,env)
  
  stopifnot(
    !is.na(link_used),
    all(c(row$Cause,row$Effect,driver,link_used)%in%names(env)))
  
  ts <- tibble(
    Time = env$Time,
    Cause = env[[row$Cause]],
    Effect = env[[row$Effect]],
    Driver = env[[driver]],
    IS = env[[link_used]])%>%
    arrange(Time)
  
  is_max <- max(abs(ts$IS),na.rm = TRUE)
  stopifnot(is.finite(is_max),is_max > 0)
  
  ts <- ts%>%
    mutate(
      across(c(Cause,Effect,Driver),~fill_na(scale01(.x))),
      IS = fill_na(IS/is_max),
      Effect_plot = -Effect)
  
  component_labels <- c(
    Upper = chem_label(row$Cause),
    Lower = chem_label(row$Effect))
  
  driver_labels <- c(
    Upper = chem_label(driver),
    Lower = link_label)
  
  p1 <- ggplot(ts,aes(Time))+
    geom_ribbon(
      aes(ymin = 0,ymax = Cause,fill = "Upper"),
      alpha = 0.9)+
    geom_ribbon(
      aes(ymin = Effect_plot,ymax = 0,fill = "Lower"),
      alpha = 0.9)+
    geom_line(
      aes(y = Cause,color = "Upper"),linewidth = 0.55)+
    geom_line(
      aes(y = Effect_plot,color = "Lower"),linewidth = 0.55)+
    geom_hline(
      yintercept = 0,color = "#8B3038",linewidth = 0.45)+
    scale_fill_manual(
      values = fill_cols,breaks = c("Upper","Lower"),
      labels = component_labels,guide = "none")+
    scale_color_manual(
      values = line_cols,breaks = c("Upper","Lower"),
      labels = component_labels)+
    scale_y_continuous(
      name = "Scaled value",
      breaks = c(-1,-0.5,0,0.5,1),
      labels = c("1.0","0.5","0","0.5","1.0"),
      expand = expansion(mult = 0))+
    scale_x_continuous(
      name = "Year",breaks = time_breaks,
      expand = expansion(mult = 0))+
    coord_cartesian(
      xlim = years,ylim = c(-1.05,1.05))+
    labs(title = paste0(
      row$Clean_Type," | ",row$Site," | n = ",row$n_obs))+
    base_theme
  
  p2 <- ggplot(ts,aes(Time))+
    geom_ribbon(
      aes(ymin = 0,ymax = Driver,fill = "Upper"),
      alpha = 0.65)+
    geom_ribbon(
      aes(
        ymin = pmin(IS,0),
        ymax = pmax(IS,0),
        fill = "Lower"),
      alpha = 0.7)+
    geom_line(
      aes(y = Driver,color = "Upper"),linewidth = 0.55)+
    geom_line(
      aes(y = IS,color = "Lower"),linewidth = 0.65)+
    geom_hline(
      yintercept = 0,color = "#8B3038",linewidth = 0.45)+
    scale_fill_manual(
      values = fill_cols,breaks = c("Upper","Lower"),
      labels = driver_labels,guide = "none")+
    scale_color_manual(
      values = line_cols,breaks = c("Upper","Lower"),
      labels = driver_labels)+
    scale_y_continuous(
      name = "Scaled value",
      breaks = c(-1,-0.5,0,0.5,1),
      expand = expansion(mult = 0))+
    scale_x_continuous(
      name = "Year",breaks = time_breaks,
      expand = expansion(mult = 0))+
    coord_cartesian(
      xlim = years,ylim = c(-1.05,1.05))+
    labs(title = "")+
    base_theme
  
  p3 <- make_qr_panel(row,env)+
    labs(title = paste0(
      "State ",row$Regime," | ",row$Clean_Type))+
    guides(alpha = "none")
  
  p1|p2|p3}

# Select examples and display

process_file <- function(network,path){
  message("\nProcessing ",network)
  
  qr_plot <- openxlsx::read.xlsx(
    path,sheet = "QR_AllThreeSig")%>%
    mutate(
      across(where(is.character),trimws),
      Clean_Type = factor(
        Clean_Type,levels = paste("Type",1:8)))%>%
    filter(all_three_sig,!is.na(Clean_Type))
  
  examples_preferred <- qr_plot%>%
    filter(
      !Cause%in%exclude,
      !Effect%in%exclude,
      !is.na(n_obs),
      !is.na(env_rel_inf))%>%
    arrange(
      Clean_Type,desc(n_obs),desc(env_rel_inf))%>%
    group_by(Clean_Type)%>%
    slice_head(n = 1)%>%
    ungroup()
  
  missing_types <- setdiff(
    levels(droplevels(qr_plot$Clean_Type)),
    examples_preferred$Clean_Type)
  
  examples <- bind_rows(
    examples_preferred,
    qr_plot%>%
      filter(
        Clean_Type%in%missing_types,
        !is.na(n_obs),
        !is.na(env_rel_inf))%>%
      arrange(
        Clean_Type,desc(n_obs),desc(env_rel_inf))%>%
      group_by(Clean_Type)%>%
      slice_head(n = 1)%>%
      ungroup())%>%
    arrange(Clean_Type)
  
  examples%>%
    select(
      Clean_Type,Site,ccmLink,
      env_var,n_obs,env_rel_inf)%>%
    print(n = Inf)
  
  sites <- unique(qr_plot$Site)
  env_list <- setNames(
    map(sites,~prep_env(
      smap$env.int.mat[[match(.x,smap$Site)]],years)),
    sites)
  
  triplets <- map(
    seq_len(nrow(examples)),
    ~make_triplet(
      examples[.x,],
      env_list[[examples$Site[.x]]]))
  
  pages <- split(
    seq_along(triplets),
    ceiling(seq_along(triplets)/4))%>%
    map(function(idx){
      panels <- map(idx,function(i){
        row <- examples[i,]
        label <- paste0(
          row$Site," | State ",row$Regime,
          " | ",row$Clean_Type," | n = ",row$n_obs)
        
        list(
          triplets[[i]][[1]]+
            labs(title = label)+
            compact_theme+
            theme(aspect.ratio = 0.618),
          
          triplets[[i]][[2]]+
            labs(title = NULL,y = NULL)+
            compact_theme+
            theme(aspect.ratio = 0.618),
          
          triplets[[i]][[3]]+
            labs(title = NULL)+
            compact_theme+
            theme(
              aspect.ratio = 1,
              legend.title = element_text(size = 7.5),
              axis.title.y = element_text(size = 8),
              axis.title.y.left = element_text(
                size = 8,angle = 90,margin = margin(r = 4)),
              axis.title.y.right = element_text(
                size = 8,angle = 270,margin = margin(l = 4))))})%>%
        flatten()
      
      panels <- c(
        panels,
        rep(list(plot_spacer()),12-length(panels)))
      
      wrap_plots(
        panels,ncol = 3,nrow = 4,
        heights = rep(1,4))})
  
  walk(pages,print)
  invisible(NULL)}

walk2(names(files),files,process_file)