library(dplyr)
library(ordinal)

cnvPlot <- function(Genetics_ranges)
{
  Genetics_ranges %>%
    mutate(End = ifelse(Gain_Loss == "Mutation", End + 5e4, End)) %>%
    ggplot(aes(x = Patient.ID, ymin = Start/1e6, ymax = End/1e6)) +
    geom_linerange(size = 1, aes(color = Gain_Loss)) +
    scale_color_manual(values = c(Gain = "blue", Loss = "red", Mutation = "darkred")) +
    scale_x_discrete(labels = NULL, limits = unique(Genetics_ranges$Patient.ID[order(desc(Genetics_ranges$min))])) +
    theme(axis.ticks.y = element_blank(), axis.line.y = element_blank(), legend.position = c(0.1,0.1), legend.title = element_blank()) +
    ylab("Chromosomic coordinates") +
    xlab("Patients") +
    geom_rect(ymin = 50674641/1e6, ymax = 50733212/1e6, xmin = 1, xmax = nrow(Genetics_ranges), alpha = .01, fill = "grey") +
    coord_flip()
}

delAnalysis <- function(x, size)
{
  if (is.numeric(x))
  {
    model <- lm(x ~ size + data$Gender)
    p <- summary(model)$coefficients["size","Pr(>|t|)"]
    icl <- ((confint(model) %>% exp) ** 1e6)["size", 1]
    or <- (coef(model)[["size"]] %>% exp) ** 1e6
    icu <- ((confint(model) %>% exp) ** 1e6)["size", 2]
  }
  else if (is.ordered(x))
  {
    model <- clm(x ~ scale(size) + data$Gender + scale(data$Age))
    p <- summary(model)$coefficients["scale(size)","Pr(>|z|)"]
    icl <- confint(model)["scale(size)", 1] %>% exp
    or <- coef(model)[["scale(size)"]] %>% exp
    icu <- confint(model)["scale(size)", 2] %>% exp
  }
  else
  {
    model <- glm(x ~ size + data$Gender + data$Age, family = binomial)
    p <- summary(model)$coefficients["size","Pr(>|z|)"]
    icl <- ((confint(model) %>% exp) ** 1e6)["size", 1]
    or <- (coef(model)[["size"]] %>% exp) ** 1e6
    icu <- ((confint(model) %>% exp) ** 1e6)["size", 2]
  }

  c(p = p,icl = icl,or = or,icu = icu)
}

delPlot <- function(var, data, results_ranges)
{
  data %>%
    select(Patient.ID, var = one_of(var)) %>%
    mutate(var = var %>% ordered) %>%
    ggplot(aes(x = Patient.ID, ymin = 0, ymax = 1, color = var)) +
    geom_linerange(size = 1) +
    scale_x_discrete(labels = NULL, limits = data$Patient.ID[order(desc(data$min))]) +
    scale_y_continuous(labels = NULL) +
    scale_color_grey(start = .9, end = .2, na.value = "white") +
    theme(axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "none") +
    xlab(NULL) +
    ylab(paste0("p = ", results_ranges$p.adj[results_ranges$Variable == var] %>% format(digits = 4))) +
    ggtitle(results_ranges$text[results_ranges$Variable == var]) +
    theme(plot.title = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 9)) +
    theme(axis.title.x = element_text(angle = -45, size = 9, vjust = 0)) +
    theme(plot.margin = grid::unit(c(0,0,1,0), "cm")) +
    coord_flip()
}
