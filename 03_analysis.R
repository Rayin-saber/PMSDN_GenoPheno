library(tidyverse)
library(desctable)
library(patchwork)

source("functions-deletion.R")

readRDS("data.rds") %>%
  list2env(envir = globalenv())

article <- list()
article$nb_pat$demo <- nrow(Demographics)
article$nb_pat$clin <- nrow(Clinical)
article$nb_pat$dev <- nrow(Developmental)

# Create main data frame -------------------------------------------------------
Clinical %>%
  full_join(Developmental) %>%
  inner_join(Demographics) -> data

rbind(Clin_vars, Dev_vars) -> vars
article$PRO$sel <- nrow(vars)
article$nb_pat$pheno <- nrow(data)

rm(Clinical, Developmental)
rm(Clin_vars, Dev_vars)

article$nb_pat$all_gen <- Genetics_ranges %>% distinct(Patient.ID) %>% nrow
article$genetics$all <- nrow(Genetics_ranges)

# Manage Genetics_ranges -------------------------------------------------------
Genetics_ranges %>%
  mutate(Gain_Loss = ifelse(Result.type == "mutation", "Mutation", Gain_Loss)) -> Genetics_ranges

article$genetics$X22$del <- Genetics_ranges %>% filter(Gain_Loss == "Loss") %>% nrow
article$genetics$X22$dup <- Genetics_ranges %>% filter(Gain_Loss == "Gain") %>% nrow
article$genetics$X22$mut <- Genetics_ranges %>% filter(Gain_Loss == "Mutation") %>% nrow
article$nb_pat$mut <- Genetics_ranges %>% filter(Gain_Loss == "Mutation") %>% distinct(Patient.ID) %>% nrow

# Create the del extent var ----------------------------------------------------
create_del_extent_var <- function(ranges, first = F)
{
  if (!first)
    ranges %>% select(-min) -> ranges
  
  ranges %>%
    full_join(ranges %>%
              filter(Gain_Loss != "Gain") %>%
              group_by(Patient.ID) %>%
              summarise(min = min(Start)))
}

Genetics_ranges %>%
  create_del_extent_var(T) %>%
  arrange(desc(min)) -> Genetics_ranges

# First plot with all patients and every genetic alteration --------------------
Genetics_ranges %>% cnvPlot -> article$CNV_plot
# save_plot(plot = CNV_plot, filename = "CNV_plot.svg", device = svglite, base_height = 12, base_width =  30)

# Fix patients with small gaps (size < 5% total deletion length) between two large deletions ---------------------
Genetics_ranges %>%
  filter(Gain_Loss == "Loss") %>%
  mutate(length = End - Start) %>%
  full_join(group_by(., Patient.ID) %>%
            summarise(total = sum(length))) %>%
  group_by(Patient.ID) %>%
  filter(n() > 1) %>%
  arrange(Patient.ID, Start) %>%
  mutate(gap = lead(Start) - End,
         ratio = gap / total) %>%
  ungroup %>%
  mutate(End = case_when(ratio < .05 ~ lead(End),
                         T ~ End)) %>%
  filter(ratio < .05) %>%
  select(Patient.ID, Gain_Loss, Start, End) %>%
  mutate(Result.type = "coordinates") -> merged

Genetics_ranges %>%
  anti_join(merged) %>%
  bind_rows(merged) -> Genetics_ranges

# Keep only terminal deletions/mutations ---------------------------------------
Genetics_ranges %>%
  filter(Gain_Loss != "Gain",
         End > 50674641) -> Genetics_ranges # shank3 : 50674642-50733210

article$nb_pat$term_del <- Genetics_ranges %>% distinct(Patient.ID) %>% nrow

# Keep only patients who have phenotypic data ----------------------------------
Genetics_ranges %>%
  semi_join(data) -> Genetics_ranges

# Re-compute del extent var ----------------------------------------------------
Genetics_ranges %>%
  create_del_extent_var %>%
  arrange(desc(min)) -> Genetics_ranges

# Keep only patients with geno & pheno data and build dataframe ----------------
data %>%
  inner_join(Genetics_ranges %>%
             distinct(Patient.ID, .keep_all = T) %>%
             select(Patient.ID, min)) %>%
  arrange(desc(min)) -> data

article$nb_pat$mut_sel <- Genetics_ranges %>% filter(Gain_Loss == "Mutation") %>% distinct(Patient.ID) %>% nrow
article$genetics$range <- Genetics_ranges %>% filter(Gain_Loss == "Loss") %>% mutate(size = End - min) %>% .$size %>% summary

# Final plot of the deletions for the included patients ------------------------
Genetics_ranges %>% cnvPlot -> article$DEL_plot

rm(Genetics_ranges)

article$nb_pat$included <- nrow(data)

# Table 1 : demographics and comparison
Demographics %>%
  mutate(included = Patient.ID %in% data$Patient.ID,
         Race = ifelse(Race == "", NA, Race) %>% factor,
         included = included %>% factor) -> Demographics

Demographics %>%
  ungroup() %>%
  select(Age, Gender, Race, included) %>%
  group_by(included) %>%
  desctable -> article$table1

# PRO selection
for (var in grep("(communication|motor)\\.milestones", names(data), value = T))
{
  levels <- c("0 - 3 months", "4 - 7 months",	"8 - 11 months",	"12 - 18 months",	"19 - 24 months",	"25 - 30 months",	"31 - 36 months", "4 - 5 years", "6 - 7 years",	"8 - 9 years", "10 years +")
  unclassed <- unclass(data[[var]])
  lvl_low <- 1
  while (sum(summary(data[[var]])[1:lvl_low]) < 6 & lvl_low < 11)
    lvl_low <- lvl_low + 1

  lvl_high <- 11
  while (sum(summary(data[[var]])[lvl_high:11]) < 6 & lvl_high > 1)
    lvl_high <- lvl_high - 1

  if (lvl_high - lvl_low < 2)
    next

  unclassed[unclassed < lvl_low] <- lvl_low
  unclassed[unclassed > lvl_high] <- lvl_high

  if (lvl_low > 1)
    levels[lvl_low] <- gsub("\\d+ - (\\d+ \\w+)", "\\1 -", levels[lvl_low], perl = T)
  if (lvl_high < 11)
    levels[lvl_high] <- gsub("(\\d+) - \\d+ (\\w+)", "\\1 \\2 +", levels[lvl_high], perl = T)
  data[[var]] <- ordered(unclassed, levels = lvl_low:lvl_high, labels = levels[lvl_low:lvl_high])
}
rm(unclassed, lvl_low, lvl_high, levels)

for (var in vars$Variable)
  if (!is.numeric(data[[var]]))
    if (any(summary(data[[var]])[1:nlevels(data[[var]])] < 6) | grepl("_Other", var))
      data[[var]] <- NULL

vars %>%
  filter(Variable %in% names(data)) -> vars

# for (var in vars$Variable)
#   for (lvl in levels(data[[var]]))
#     if (!is.na(lvl))
#       vars[vars$Variable == var, lvl] <- summary(factor(data[[var]]))[lvl]

article$PRO$incl <- vars %>% nrow
article$completion <- sapply(data[vars$Variable], function (x) {is.na(x)[is.na(x) == F] %>% length/nrow(data)}) %>% summary

# Association analysis ---------------------------------------------------------
data %>%
  select(-Patient.ID, -Birthdate, -Gender, -Age, -Age_months, -Race, -Country, -min) %>%
  names %>%
  map(delAnalysis, data) %>%
  bind_rows %>%
  mutate(p.adj = p.adjust(p, "fdr")) %>%
  arrange(p.adj) %>%
  left_join(vars) -> results_ranges

results_ranges %>%
  mutate_at(vars(-p, -p.adj, -or), ~ prettyNum(., digits = 3)) %>%
  mutate(IC_OR = str_c("[", icl, "-", icu, "]", sep = " ")) %>%
  select(Variable, text,  p, p.adj, or, IC_OR, Group) -> article$results_ranges

data %>%
  arrange(desc(min)) %>%
  mutate_if(is.factor, . %>% ordered %>% as.numeric) %>%
  gather(var, value, -Patient.ID, -Birthdate, -Gender, -Race, -Country, -Age, -Age_months, -min) %>%
  mutate(value = value %>% ordered) %>%
  left_join(results_ranges, by = c("var" = "Variable")) -> dataplot

# Plots-------------------------------------------------------------------------
dataplot %>%
  filter(Group == "Renal conditions") %>%
  delPlotGroup(article$DEL_plot) -> kidney

dataplot %>%
  group_split(Group) %>%
  map(delPlotGroup, article$DEL_plot) %>%
  setNames(dataplot$Group %>% unique %>% sort) -> article$plots

article$plots %>%
  map2(str_c("resultats/", names(.), ".png"),
       ~ ggsave(plot = .x,
                filename = .y,
                width = 15,
                height = 10))

article$results_ranges %>%
  select(Group, Var = text, p, p.adj, or, IC_OR, fullname = Variable) %>%
  arrange(p.adj) %>%
  write_excel_csv("resultats/resultats.csv")

saveRDS(article, "article.rds")

# ROC curves--------------------------------------------------------------------
# Keep only the maximum extent of deletion for each patient
Genetics_ranges %>%
  group_by(Patient.ID) %>%
  summarize(Start = min(Start)) -> Start
rocvars <- filter(results_ranges, Range_p_corrected < .05) %>% select(Group, Variable)
data$Patient.ID <- as.numeric(data$Patient.ID)
datarocs <- merge(Start,data[c("Patient.ID",rocvars$Variable)])
for (var in rocvars$Variable)
  if (nlevels(datarocs[[var]]) > 2 | nlevels(datarocs[[var]]) == 0)
    datarocs[var] <- NULL

rocvars <- rocvars[rocvars$Variable %in% names(datarocs)[-(1:2)],]

for (var in rocvars$Variable)
{
  print(var)
  dataroc <- select(datarocs, one_of("Start", var))
  thres_bt <- unlist(mclapply(1:10000, mc.cores = 8, function(x)
  {
    dataroc_sample <- sample_n(dataroc, nrow(dataroc), replace = T)
    roc <- tryCatch({roc(dataroc_sample[[var]], dataroc_sample$Start, percent = T)}, error = function(e){roc(dataroc[[var]], dataroc$Start, percent = T)})

    youden <- roc$sensitivities + roc$specificities - 100
    roc$thresholds[which(youden == max(youden))]
  }))

  roc <- roc(dataroc[[var]], dataroc$Start, percent = T)
  youden <- roc$sensitivities + roc$specificities - 100
  thres <- roc$thresholds[which(youden == max(youden))]

  hist(thres_bt, main = var, freq = F, xlim = c(40e6, 51e6))
  lines(density(thres_bt, adjust = 2))
  abline(v = mean(thres_bt), col = "orange")
  abline(v = thres, col = "darkorange", lty = 2)
  abline(v = quantile(thres_bt, probs = .025), col = "green")
  abline(v = quantile(thres_bt, probs = .975), col = "red")


  rocvars$min[rocvars$Variable == var] <- max(dataroc$Start[dataroc$Start < thres])
  rocvars$max[rocvars$Variable == var] <- min(dataroc$Start[dataroc$Start > thres])
  rocvars$auc[rocvars$Variable == var] <- as.numeric(auc(roc))
  rocvars$thres[rocvars$Variable == var] <- thres
  rocvars$icl[rocvars$Variable == var] <- quantile(thres_bt, probs = .025)
  rocvars$icu[rocvars$Variable == var] <- quantile(thres_bt, probs = .975)
}

# rocvars <- rocvars[rocvars$auc > 70,]
write.csv(rocvars, "rocvars.csv")
rocvars <- read.csv("rocvars.csv")
rocvars <- merge(rocvars, group_by(rocvars, Group) %>% summarise(gmin = min(min)))
rocvars <- merge(rocvars, group_by(rocvars, Group) %>% summarise(gmax = max(max)))
rocvars <- arrange(rocvars, gmin, Group,min)
rocvars$Variable <- factor(rocvars$Variable, levels = rocvars$Variable)
rocvars$Group <- factor(rocvars$Group, levels = rev(unique(rocvars$Group)))

ggplot(rocvars) +
  geom_linerange(aes(x = Variable, ymin = min / 1000000, ymax = max / 1000000, colour = Group), size = 8) +
  geom_linerange(aes(x = Variable, ymin = icl / 1e6, ymax = icu / 1e6, colour = Group), size = 8, alpha = .5) +
  #geom_linerange(aes(x=Variable, ymin = gmin/1000000, ymax = gmax/1000000, colour = Group), size = 8, alpha = .5) +
  geom_text(aes(x = Variable, y = gmax / 1000000, label = srt_c(as.integer(auc), "%")), hjust = 0) +
  scale_y_continuous(name = "Position") +
  coord_flip()
