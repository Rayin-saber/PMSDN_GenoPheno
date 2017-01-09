library(cosmosR)
library(readr)
library(rvest)
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

article <- list()
article$nb_pat$demo <- nrow(Demographics)
article$nb_pat$clin <- nrow(Clinical)
article$nb_pat$dev <- nrow(Developmental)

# Create main data frame -------------------------------------------------------
Demographics %>%
  inner_join(Clinical) %>%
  full_join(Developmental) %>%
  ungroup -> data

rbind(Clin_vars, Dev_vars) -> vars
article$PRO$sel <- nrow(vars)
article$nb_pat$pheno <- nrow(data)

rm(Clinical, Developmental)
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
  mutate(text = text %>% gsub("\\.", " ", ., perl = T)) %>%
  mutate(text = text %>% sub("With PMS, some children develop understandable verbal speech while others do not  If the patient is over age 2, please choose the response that most closely matches your child's abilities today", "Verbal speech", ., fixed = T)) %>%
  mutate(text = text %>% sub("Has the patient ever exhibited symptoms of or been diagnosed with ADD/ADHD or any other significant attention or hyperactivity issues warranting intervention?", "Symptoms/diagnostic of ADD/ADHD", ., fixed = T)) %>%
  mutate(text = text %>% sub("Has less post-operative pain than would be expected", "Less post-op pain than expected", ., fixed = T)) %>%
  mutate(text = text %>% sub("Needs a great deal of post-operative pain medication", "Great need of post-op pain med", ., fixed = T)) %>%
  mutate(text = text %>% sub("Needs less post-operative pain medication than would be expected", "Less post-op pain med than expected", ., fixed = T)) %>%
  mutate(text = text %>% sub("What is the highest number of ear infections the patient has had in a given year?", "Highest number of ear infections in a year", ., fixed = T)) %>%
  mutate(text = text %>% sub("Has the patient had ear infections after age three, which required recurrent courses of antibiotics?", "Ear infections age>3, requiring antibiotics", ., fixed = T)) %>%
  mutate(text = text %>% sub("If the patient ever had a metabolic panel test, was it abnormal?", "Abnormal metabolic panel test", ., fixed = T)) %>%
  mutate(text = text %>% sub("Excessive irritability especially following mealtimes", "Excessive irritability following mealtimes", ., fixed = T)) %>%
  mutate(text = text %>% sub("Two or more serious sinus infections within one year", "2+ serious sinus infections within a year", ., fixed = T)) %>%
  mutate(text = text %>% sub("Two or more months on antibiotics with little effect", "2+ months on antibiotics with little effect", ., fixed = T)) %>%
  mutate(text = text %>% sub("Persistent thrush in mouth or fungal infection on skin", "Persistent thrush in mouth/fungal infection on skin", ., fixed = T)) %>%
  mutate(text = text %>% sub("Need for intravenous antibiotics to clear infections", "Need for IV antibiotics to clear infections", ., fixed = T)) %>%
  mutate(text = text %>% sub("If not, has the delay resulted from immune system concerns?", "Delay resulted from immune system concerns", ., fixed = T)) %>%
  mutate(text = text %>% sub("Has the patient ever had a severe adverse reaction to any vaccines?", "Severe adverse reaction to any vaccines", ., fixed = T)) %>%
  mutate(text = text %>% sub("Bloating is apparent at all times but varies in severity", "Bloating apparent at all times but varies in severity", ., fixed = T)) %>%
  mutate(text = text %>% sub("Concerns about the mother's ability to deliver vaginally", "Concerns about the mother's ability to deliver vaginally", ., fixed = T)) %>%
  mutate(text = text %>% sub("Difficulty going back to sleep after nighttime awakening", "Difficulty going back to sleep after nighttime awakening", ., fixed = T)) %>%
  mutate(text = text %>% sub("Have difficulty sustaining attention in tasks or play activities", "Difficulty sustaining attention in tasks", ., fixed = T)) %>%
  mutate(text = text %>% sub("Leave seat in classroom or in other situations in which remaining seated is expected", "Leave seat when remaining seated is expected", ., fixed = T)) %>%
  mutate(text = text %>% sub("Run around or climb excessively in situations in which it is inappropriate", "Run around or climb excessively when inappropriate", ., fixed = T)) %>%
  mutate(text = text %>% sub("Have difficulty playing or engaging in leisure activities quietly", "Difficulty playing quietly", ., fixed = T)) %>%
  mutate(text = text %>% sub("Has the patient ever exhibited excessive mouthing/chewing not associated with teething?", "Excessive chewing not associated with teething", ., fixed = T)) %>%
  mutate(text = text %>% sub("Demonstrate understanding of the use of familiar objects", "Understanding of the use of familiar objects", ., fixed = T)) %>%
  mutate(text = text %>% sub("Act out recent experiences using gestures, body languages, words", "Act out recent experiences using gestures, etc.", ., fixed = T)) %>%
  mutate(text = text %>% sub("Observe actions of others and imitate them at a later time", "Observe actions and imitate them later", ., fixed = T)) %>%
  mutate(text = text %>% sub("Has the patient ever experienced regression in cognitive abilities?", "Regression in cognitive abilities", ., fixed = T)) %>%
  mutate(text = text %>% sub("What is the patient's ability to follow directions?", "Ability to follow directions", ., fixed = T)) %>%
  mutate(text = text %>% sub("Look at, reach and grasp objects placed at a distance", "Look at, reach and grasp distant objects", ., fixed = T)) %>%
  mutate(text = text %>% sub("Does the patient have difficulties with motor planning or motor learning?", "Difficulties with motor planning or learning", ., fixed = T)) %>%
  mutate(text = text %>% sub("Does the patient have difficulties maintaining balance?", "Difficulties maintaining balance", ., fixed = T)) %>%
  mutate(text = text %>% sub("If the patient has difficulties maintaining balance, does the inability to balance depend on sight?", "Inability to balance dependant on sight", ., fixed = T)) %>%
  mutate(text = text %>% sub("Mixed hypo- and hyper-sensitivity to touch throughout the body", "Mixed hypo- and hyper-sensitivity to touch", ., fixed = T)) %>%
  mutate(text = text %>% sub("What is the level of patient's response to movement?", "Response to movement", ., fixed = T)) %>%
  mutate(text = text %>% sub("Initiate social interaction by smiling, moving arms and/or vocalizing", "Initiate social interaction", ., fixed = T)) %>%
  mutate(text = text %>% sub("Has the patient ever had poor eye contact?", "Poor eye contact", ., fixed = T)) %>%
  mutate(text = text %>% sub("Cling to caregivers/familiar adult in presence of a stranger", "Cling to familiar adult in presence of a stranger", ., fixed = T)) -> vars

article$nb_pat$all_gen <- Genetics_ranges %>% distinct(Patient.ID) %>% nrow
article$genetics$all <- nrow(Genetics_ranges)
article$genetics$Xmes <- Genetics_ranges$Chr_Gene %>% factor %>% levels

# Manage Genetics_ranges -------------------------------------------------------
Genetics_ranges %<>%
  filter(Chr_Gene == "22") %>%
  mutate(Patient.ID = as.character(Patient.ID)) %>%
  mutate(Gain_Loss = ifelse(Result.type == "mutation", "Mutation", Gain_Loss)) %>%
  select(-Result.type, -Chr_Gene)

article$genetics$X22$del <- Genetics_ranges %>% filter(Gain_Loss == "Loss") %>% nrow
article$genetics$X22$dup <- Genetics_ranges %>% filter(Gain_Loss == "Gain") %>% nrow
article$genetics$X22$mut <- Genetics_ranges %>% filter(Gain_Loss == "Mutation") %>% nrow
article$nb_pat$mut <- Genetics_ranges %>% filter(Gain_Loss == "Mutation") %>% distinct(Patient.ID) %>% nrow

# Create the del extend var ----------------------------------------------------
Genetics_ranges %<>%
  full_join(Genetics_ranges %>%
            filter(Gain_Loss != "Gain") %>%
            group_by(Patient.ID) %>%
            summarize(min = min(Start))) %>%
  arrange(desc(min))

# First plot with all patients and every genetic alteration --------------------
Genetics_ranges %>% cnvPlot -> article$CNV_plot
# save_plot(plot = CNV_plot, filename = "CNV_plot.svg", device = svglite, base_height = 12, base_width =  30)

# Keep only terminal deletions/mutations ---------------------------------------
Genetics_ranges %<>% filter(Gain_Loss != "Gain", End > 50674641) # shank3 : 50674641

article$nb_pat$term_del <- Genetics_ranges %>% distinct(Patient.ID) %>% nrow

# Keep only patients who have phenotypic data ----------------------------------
Genetics_ranges %<>% semi_join(data)

# Re-compute del extent var ----------------------------------------------------
Genetics_ranges %<>% select(-min)
Genetics_ranges %<>%
  full_join(Genetics_ranges %>%
            group_by(Patient.ID) %>%
            summarize(min = min(Start))) %>%
  arrange(desc(min))

# Keep only patients with geno & pheno data and build dataframe ----------------
data %<>%
  inner_join(Genetics_ranges %>%
             distinct(Patient.ID, .keep_all = T) %>%
             select(Patient.ID, min)) %>%
  arrange(desc(min))


article$nb_pat$mut_sel <- Genetics_ranges %>% filter(Gain_Loss == "Mutation") %>% distinct(Patient.ID) %>% nrow
article$genetics$range <- Genetics_ranges %>% filter(Gain_Loss == "Loss") %>% mutate(size = End - min) %>% .$size %>% summary

# Final plot of the deletions for the included patients ------------------------
Genetics_ranges %>% cnvPlot -> article$DEL_plot

rm(Genetics_ranges)

article$nb_pat$included <- nrow(data)

# Table 1 : demographics and comparison
Demographics %<>%
  mutate(included = Patient.ID %in% data$Patient.ID)
Demographics$Race[Demographics$Race == ""] <- NA
Demographics$included %<>% factor
Demographics$Race %<>% factor

Demographics %>% select(Age, Gender, Race, included) %>% desc_groupe("included")
read_file("desc_groupe_included.html") %>% read_html %>% html_table %>% .[[1]] -> article$table1

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

vars %<>% filter(Variable %in% names(data))

# for (var in vars$Variable)
#   for (lvl in levels(data[[var]]))
#     if (!is.na(lvl))
#       vars[vars$Variable == var, lvl] <- summary(factor(data[[var]]))[lvl]

article$PRO$incl <- vars %>% nrow
article$completion <- sapply(data[vars$Variable], function (x) {is.na(x)[is.na(x) == F] %>% length/nrow(data)}) %>% summary

# Association analysis ---------------------------------------------------------
data %>%
  select(-Patient.ID, -Birthdate, -Gender, -Age, -Age_months, -Race, -Country, -min, -`Is.the.patient's.menstrual.cycle.regular?_ currently`) %>%
  sapply(delAnalysis, 50818468 - data$min) %>% # end of chr22
  t %>%
    data.frame %>%
    add_rownames("Variable") %>%
    mutate(p.adj = p.adjust(p, "fdr")) %>%
    arrange(p.adj) %>%
    left_join(vars) -> results_ranges

results_ranges %>%
  mutate_each(funs(prettyNum(., digits = 3)), -p, -p.adj, -or) %>%
  mutate(IC_OR = str_c("[", icl, "-", icu, "]", sep = " ")) %>%
  select(Variable, text,  p, p.adj, or, IC_OR, Group) -> article$results_ranges
# write.csv2("results_ranges.csv", row.names = F)

save(data,article,results_ranges,vars, file = "results.Rdata")
load("results.Rdata")

library(tidyr)
library(purrr)
library(magrittr)
library(gtable)

data %<>%
  mutate(Patient.ID = Patient.ID %>% as.numeric,
         limit = order(min)) %>%
  dmap_if(is.factor, . %>% ordered %>% as.numeric) %>%
  gather(var, value, -Patient.ID, -Birthdate, -Gender, -Race, -Country, -Age, -Age_months, -min, -limit) %>%
  mutate(value = value %>% ordered) %>%
  left_join(results_ranges, by = c("var" = "Variable"))

# Plots-------------------------------------------------------------------------
data %>%
  filter(Group == "Renal.-.Kidney") %>%
  arrange(p.adj) %>%
  mutate(text = text %>% factor(levels = text %>% unique),
         p.adj = p.adj %>% prettyNum(digits = 3)) %>%
  ggplot(aes(x = limit, ymin = 0, ymax = 1, color = value)) +
  facet_grid(. ~ text) +
  geom_linerange(size = 1) +
  geom_text(color = "black", angle = -45, aes(x = -12, y = 0.5, label = paste0("p = ", p.adj))) +
  scale_y_continuous(labels = NULL) +
  scale_color_grey(start = .9, end = .2, na.value = "white") +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        strip.text = element_text(vjust = 0, hjust = 0, angle = 75),
        strip.background = element_blank())+
coord_flip()

unique(results_ranges$Group) %>%
  sort %>%
  sapply(simplify = F, function(group)
{
  print(group)
  results_ranges %>%
    filter(Group == group) %>%
    arrange(p.adj) %>%
    .[["Variable"]] %>%
    lapply(delPlot, data, results_ranges) %>%
    plot_grid(article$DEL_plot, plotlist = ., align = "hv", nrow = 1, rel_widths = c(6, rep(1, length(.)))) #%>%
    # save_plot(filename = str_c("plots/",group, ".png"), base_height = 12) #, base_width =  6 + .75 * length(.$layers), limitsize = F)
    # save_plot(filename = str_c("plots/",group, ".svg"), device = svglite, base_height = 12) #, base_width = 6 + .75 * length(.$layers))
}) -> article$plots

save(article, file = "article.Rdata")

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
