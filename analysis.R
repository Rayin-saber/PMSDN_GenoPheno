library(dplyr)
library(ggplot2)
library(magrittr)
library(cowplot)
library(pROC)
library(parallel)
library(RSvgDevice)
library(DT)

source("functions-deletion.R")

load("data.Rda")

# Create main data frame -------------------------------------------------------
full_join(Demographics, Clinical) %>%
  full_join(Developmental) %>%
  ungroup -> data

vars <- rbind(Clin_vars, Dev_vars)

rm(Demographics, Clinical, Developmental)
rm(Clin_vars, Dev_vars)

# Manage Genetics_ranges -------------------------------------------------------
Genetics_ranges %<>%
  filter(Chr_Gene == "22") %>%
  mutate(Patient.ID = as.character(Patient.ID)) %>%
  mutate(Gain_Loss = ifelse(Result.type == "mutation", "Mutation", Gain_Loss)) %>%
  select(-Result.type, -Chr_Gene)

# Create the del extend var ----------------------------------------------------
Genetics_ranges %<>%
  full_join(
    Genetics_ranges %>%
      filter(Gain_Loss != "Gain") %>%
      group_by(Patient.ID) %>%
      summarize(min = min(Start))
  ) %>%
  arrange(desc(min))

# First plot with all patients and every genetic alteration --------------------
Genetics_ranges %>%
  mutate(End = ifelse(Gain_Loss == "Mutation", End + 5e4, End)) %>%
  ggplot(aes(x = Patient.ID, ymin = Start/1e6, ymax = End/1e6)) +
  geom_linerange(size = 1, aes(color = Gain_Loss)) +
  scale_color_manual(values = c(Gain = "blue", Loss = "red", mutation = "darkred")) +
  scale_x_discrete(labels = NULL, limits = unique(Genetics_ranges$Patient.ID[order(desc(Genetics_ranges$min))])) +
  theme(axis.ticks.y = element_blank(), legend.position = "left") +
  ylab("Chromosomic position") +
  xlab("Patients") +
  geom_rect(ymin = 50674641/1e6, ymax = 50733212/1e6, xmin = 1, xmax = nrow(Genetics_ranges), alpha = .01, fill = "grey") +
  coord_flip() -> CNV_plot

# Keep only terminal deletions/mutations ---------------------------------------
Genetics_ranges %<>% filter(Gain_Loss != "Gain", End > 50500000) # shank3 = 50674641

# And only patients who have phenotypic data ----------------------------------
Genetics_ranges %<>% semi_join(data)

# Re-compute del extent var ----------------------------------------------------
Genetics_ranges %<>% select(-min)
Genetics_ranges %<>%
  full_join(
    Genetics_ranges %>%
      group_by(Patient.ID) %>%
      summarize(min = min(Start))
  ) %>%
  arrange(desc(min))

# Final plot of the deletions for the included patients ------------------------
Genetics_ranges %>%
  mutate(End = ifelse(Gain_Loss == "Mutation", End + 5e4, End)) %>%
  ggplot(aes(x = Patient.ID, ymin = Start/1e6, ymax = End/1e6)) +
  geom_linerange(size = 1, aes(color = Gain_Loss)) +
  scale_color_manual(values = c(Loss = "red", Mutation = "darkred")) +
  scale_x_discrete(labels = NULL, limits = unique(Genetics_ranges$Patient.ID[order(desc(Genetics_ranges$min))])) +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank(), legend.position = "left", legend.title = element_blank()) +
  ylab("Chromosomic position") +
  xlab("Patients") +
  geom_rect(ymin = 50674641/1e6, ymax = 50733212/1e6, xmin = 0.5, xmax = nrow(Genetics_ranges), alpha = .01, fill = "grey") +
  coord_flip() -> DEL_plot

# Keep only patients with geno & pheno data and build dataframe ----------------
data %<>%
  inner_join(
    Genetics_ranges %>%
      distinct(Patient.ID) %>%
      select(Patient.ID, min)
  ) %>%
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

# Plot for one var (test)-------------------------------------------------------
data %>%
  ggplot(aes(x = Patient.ID, ymin = 0, ymax = 1, color = data[["Has.the.patient.been.diagnosed.with.any.of.the.following.kidney.conditions?_Vesicoureteral Reflux (Renal reflux)_Conditions"]])) +#`Has.the.patient.been.diagnosed.with.any.of.the.following.kidney.conditions?_Vesicoureteral Reflux (Renal reflux)_Conditions`)) +
  geom_linerange(size = 1) +
  scale_x_discrete(labels = NULL, limits = data$Patient.ID[order(desc(data$min))]) +
  scale_y_continuous(labels = NULL) +
  scale_color_grey(start = .9, end = .2, na.value = "white") +
  theme(axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "none") +
  xlab("") +
  ggtitle("Vesico-ureteral Reflux") +
  theme(plot.title = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 10)) +
  coord_flip() -> g

plot_grid(DEL_plot, g, g, g, g, align = "h", nrow = 1, rel_widths = c(4, 1, 1, 1, 1)) %>%
  save_plot(filename = "test.svg", device = devSVG, nrow = 1, base_height = 7, base_width = 12)

# Plots-------------------------------------------------------------------------
dir.create("delplotsrangesvg")
for (var in vars$Variable[1:100])
  delPlotRange(Genetics_ranges, data, results_ranges, var, noOutput = T, bnw = T)

for (var in vars$Variable[vars$Group == "Renal.-.Kidney"])
  delPlotRange(Genetics_ranges, data, results_ranges, var, noOutput = T, bnw = T)

for (var in vars$Variable[vars$Group == "Gross.Motor.Development"])
  delPlotRange(Genetics_ranges, data, results_ranges, var, noOutput = T, bnw = T)

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
  geom_linerange(aes(x = Variable, ymin = min/1000000, ymax = max/1000000, colour = Group), size = 8) +
  geom_linerange(aes(x = Variable, ymin = icl/1e6, ymax = icu/1e6, colour = Group), size = 8, alpha = .5) +
  #geom_linerange(aes(x=Variable, ymin = gmin/1000000, ymax = gmax/1000000, colour = Group), size = 8, alpha = .5) +
  geom_text(aes(x = Variable, y = gmax/1000000, label = paste0(as.integer(auc), "%")), hjust=0) +
  scale_y_continuous(name = "Position") +
  coord_flip()