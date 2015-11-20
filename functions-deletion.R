library(epicalc)
library(dplyr)
library(ggplot2)

delAnalysis <- function(genetics_ranges, data, depvar)
{
  genetics_ranges <- filter(genetics_ranges, Chr_Gene == "22", Gain_Loss == "Loss" | Result.type == "mutation")
  genetics_ranges <- mutate(genetics_ranges, size = End - Start + 1)

  patients <- group_by(genetics_ranges, Patient.ID) %>% summarize(size = sum(size)) %>% arrange(desc(size))
  patients <- add_rownames(patients, var = "y")

  df <- merge(data, patients)

  if (is.numeric(df[[depvar]]))
  {
    model <- lm(df[[depvar]] ~ df$size + df$Gender)
    table <- regress.display(model)$table[-1,]
  } else if (is.ordered(df[[depvar]]))
  {
    df2 <<- df
    depvar <<- depvar
    model <- polr(df2[[depvar]] ~ scale(df2$size) + df2$Gender + scale(df2$Age))
    table <- ordinal.or.display(model)
  } else
  {
    if (depvar == "Is.the.patient's.menstrual.cycle.regular?_ currently")
      model <- glm(df[[depvar]] ~ df$size + df$Age, family = binomial)
    else
      model <- glm(df[[depvar]] ~ df$size + df$Gender + df$Age, family = binomial)
    table <- logistic.display(model)$table
  }

  or  <- as.numeric(table[1, 1])
  icl <- as.numeric(table[1, 2])
  icu <- as.numeric(table[1, 3])
  p   <- as.numeric(table[1, 4])

  sd <- sd(patients$size)

  if (is.ordered(df[[depvar]]))
  {
    or <- (or ** (1 / sd)) ** 1e6
    icl <- (icl ** (1 / sd)) ** 1e6
    icu <- (icu ** (1 / sd)) ** 1e6
  }
  else
  {
    or <- or ** 1e6
    icl <- icl ** 1e6
    icu <- icu ** 1e6
  }

  c(p,icl,or,icu)
}