library(stringr)
library(EMCluster)
library(dplyr)
library(ggplot2)
library(magrittr)
library(cowplot)
library(pROC)
library(parallel)
library(svglite)
library(DT)
library(grid)


source("functions-deletion.R")

load("data.Rda")

# Create main data frame -------------------------------------------------------
Demographics %>%
  full_join(Clinical) %>%
  full_join(Developmental) %>%
  ungroup -> data

rbind(Clin_vars, Dev_vars) -> vars

rm(Demographics, Clinical, Developmental)
rm(Clin_vars, Dev_vars)

# Clean variable names ---------------------------------------------------------
vars %>%
  filter(Group == "Self-Help.Skills") %>%
  mutate(text = Variable %>% gsub("_.*", "", ., perl = T)) %>%
  bind_rows(vars %>%
            filter(Group != "Self-Help.Skills") %>%
            mutate(text = Variable %>% gsub(".*?_(.*?)_.*", "\\1", ., perl = T)) %>%
            mutate(text = text %>% gsub("_ currently", "", ., perl = T)) %>%
            mutate(text = text %>% sub(".*?_", "", ., perl = T)) %>%
            mutate(Group = Group %>% factor)) %>%
  mutate(text = text %>% sub("\\(.*?\\)$", "", ., perl = T)) %>%
  mutate(text = text %>% gsub("\\.", " ", ., perl = T))

# Manage Genetics_ranges -------------------------------------------------------
Genetics_ranges %<>%
  filter(Chr_Gene == "22") %>%
  mutate(Patient.ID = as.character(Patient.ID)) %>%
  mutate(Gain_Loss = ifelse(Result.type == "mutation", "Mutation", Gain_Loss)) %>%
  select(-Result.type, -Chr_Gene)

# Create the del extend var ----------------------------------------------------
Genetics_ranges %<>%
  full_join(Genetics_ranges %>%
            filter(Gain_Loss != "Gain") %>%
            group_by(Patient.ID) %>%
            summarize(min = min(Start))) %>%
  arrange(desc(min))

# First plot with all patients and every genetic alteration --------------------
Genetics_ranges %>% cnvPlot -> CNV_plot

# Keep only terminal deletions/mutations ---------------------------------------
Genetics_ranges %<>% filter(Gain_Loss != "Gain", End > 50500000) # shank3 = 50674641

# And only patients who have phenotypic data ----------------------------------
Genetics_ranges %<>% semi_join(data)

# Re-compute del extent var ----------------------------------------------------
Genetics_ranges %<>% select(-min)
Genetics_ranges %<>%
  full_join(Genetics_ranges %>%
            group_by(Patient.ID) %>%
            summarize(min = min(Start))) %>%
  arrange(desc(min))

# Final plot of the deletions for the included patients ------------------------
Genetics_ranges %>% cnvPlot -> DEL_plot

# Keep only patients with geno & pheno data and build dataframe ----------------
data %<>%
  inner_join(Genetics_ranges %>%
             distinct(Patient.ID) %>%
             select(Patient.ID, min)) %>%
  arrange(desc(min))

# Association analysis ---------------------------------------------------------
data %>%
  select(-Patient.ID, -Birthdate, -Gender, -Ancestral.Background, -Country, -min, -`Is.the.patient's.menstrual.cycle.regular?_ currently`) %>%
  sapply(delAnalysis, 50818468 - data$min) %>% # end of chr22
  t %>%
    data.frame %>%
    add_rownames("Variable") %>%
    mutate(p.adj = p.adjust(p, "fdr")) %>%
    arrange(p.adj) %>%
    left_join(vars) -> results_ranges

write.csv2(results_ranges, "results_ranges.csv")

# Plots-------------------------------------------------------------------------
for (group in unique(results_ranges$Group))
{
  print(group)
  results_ranges %>%
    filter(Group == group) %>%
    arrange(p.adj) %>%
    .[["Variable"]] %>%
    lapply(delPlot, data, results_ranges) %>%
    plot_grid(DEL_plot, plotlist = ., align = "h", nrow = 1, rel_widths = c(2, rep(1, length(.)))) %>%
    save_plot(filename = str_c("plots_test/",group, ".svg"), device = svglite, base_height = 12, base_width =  4 + 2 * length(.))
}

# ROC curves--------------------------------------------------------------------
# Keep only the maximum extent of deletion for each patient
Start <- Genetics_ranges %>% group_by(Patient.ID) %>% summarize(Start = min(Start))
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
