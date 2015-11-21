library(dplyr)
library(ordinal)

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
    icl <- confint(model)["scale(size)", 1]
    or <- coef(model)[["scale(size)"]]
    icu <- confint(model)["scale(size)", 2]
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