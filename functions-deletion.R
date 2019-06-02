library(dplyr)
library(ordinal)
library(ggplot2)
library(grid)

cnvPlot <- function(Genetics_ranges)
{
  Genetics_ranges %>%
    mutate(Patient.ID = Patient.ID %>% as.character) ->
    Genetics_ranges

  Genetics_ranges %>%
    filter(Result.type == "coordinates") %>%
    ggplot() +
    aes(x = Patient.ID) +
    scale_x_discrete(labels = NULL, limits = unique(Genetics_ranges$Patient.ID[order(desc(Genetics_ranges$min))])) +
    geom_linerange(size = 1, aes(color = Gain_Loss, ymin = Start/1e6, ymax = End/1e6)) +
    scale_color_manual(values = c(Gain = "blue", Loss = "red")) +
    geom_rect(ymin = 50674641/1e6,
              ymax = 50733211/1e6,
              xmin = 1,
              xmax = Genetics_ranges %>% distinct(Patient.ID) %>% nrow,
              alpha = .01,
              fill = "grey") +
    geom_point(data = Genetics_ranges %>% filter(Result.type == "mutation"), aes(y = Start/1e6, shape = Result.type), alpha = .4) +
    scale_shape_manual(values = c(mutation = 8)) +
    ylab("Chromosomic coordinates (Mb)") +
    xlab(NULL) +
    theme_minimal() +
    theme(axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = c(0.1, 0.1),
          legend.title = element_blank()) +
    coord_flip() +
    annotate("segment", x = nrow(Genetics_ranges %>% distinct(Patient.ID)), xend = 1, y = 52, yend = 52, arrow = arrow(type = "closed", angle = 20)) +
    annotate("text", x = nrow(Genetics_ranges %>% distinct(Patient.ID))/2, y = 51.95, label = "Patients", angle = 90, hjust = .5, vjust = 0) +
    annotate("text", x = nrow(Genetics_ranges %>% distinct(Patient.ID))/2, y = 52.05, label = "Decreasing deletion size", angle = 90, hjust = .5, vjust = 1)
}

delAnalysis <- function(var, data)
{
  print(var)
  x <- data[[var]]

  data %>%
    mutate(size = 50818468 - min) -> data

  if (is.numeric(x))
  {
    # Run full model only if multiple genders
    if (data$Gender[!is.na(x)] %>% factor %>% nlevels > 1)
      model <- lm(x ~ size + Gender, data = data)
    else
      model <- lm(x ~ size, data = data)
    
    p <- summary(model)$coefficients["size","Pr(>|t|)"]
    icl <- ((confint(model) %>% exp) ** 1e6)["size", 1]
    or <- (coef(model)[["size"]] %>% exp) ** 1e6
    icu <- ((confint(model) %>% exp) ** 1e6)["size", 2]
  }
  else if (is.ordered(x))
  {
    model <- clm(x ~ scale(size) + Gender + scale(Age), data = data)
    p <- summary(model)$coefficients["scale(size)","Pr(>|z|)"]
    icl <- confint(model)["scale(size)", 1] %>% exp
    or <- coef(model)[["scale(size)"]] %>% exp
    icu <- confint(model)["scale(size)", 2] %>% exp
  }
  else
  {
    model <- glm(x ~ size + Gender + Age, family = binomial, data = data)
    p <- summary(model)$coefficients["size","Pr(>|z|)"]
    icl <- ((confint(model) %>% exp) ** 1e6)["size", 1]
    or <- (coef(model)[["size"]] %>% exp) ** 1e6
    icu <- ((confint(model) %>% exp) ** 1e6)["size", 2]
  }

  data.frame(Variable = var,
             p = p,
             icl = icl,
             or = or,
             icu = icu,
             stringsAsFactors = F,
             check.names = F)
}

delPlotGroup <- function(plotdata, delplot)
{
  plotdata %>%
    arrange(p.adj) %>%
    mutate(text = text %>% fct_inorder) %>%
    group_split(text) %>%
    map(delPlotOne) %>%
    wrap_plots(nrow = 1) -> groupplot

  plotdata$text %>%
    str_length %>%
    max -> max_text_length

  delplot +
    ggtitle(unique(plotdata$Group)) +
    theme(plot.title = element_text(size = 36, margin = margin(b = max_text_length * 7, t = 10, l = 10))) -> delplot

  delplot +
    groupplot +
    plot_layout(widths = c(1, 5))
}

delPlotOne <- function(plotdata)
{
  plotdata$p.adj %>%
    unique %>%
    prettyNum(digits = 2) %>%
    str_c("p = ", .) -> pvalue

  plotdata %>%
    distinct(text, .keep_all = T) %>%
    mutate_at(vars(or, icl, icu), ~prettyNum(., digits = 2)) %>%
    transmute(OR = str_c("OR = ", or, " [", icl, "-", icu, "]")) %>%
    pull -> OR

  str_c(pvalue, OR, sep = "\n") -> xtitle

  plotdata$text %>%
    unique -> titre

  plotdata %>%
    mutate(Patient.ID = Patient.ID %>% as.character) %>%

    ggplot() +
    aes(x = Patient.ID, ymin = 0, ymax = 1, color = value) +
    labs(title = titre,
         x = NULL,
         y = xtitle) +
    geom_linerange(size = 1) +
    scale_x_discrete(labels = NULL, limits = unique(plotdata$Patient.ID[order(desc(plotdata$min))])) +
    scale_y_continuous(labels = NULL) +
    scale_color_grey(start = .9, end = .2, na.value = "white") +
    theme(plot.title = element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.title.x = element_text(size = 9),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0, 0, 1, 0), "lines")) +
    coord_flip()
}
