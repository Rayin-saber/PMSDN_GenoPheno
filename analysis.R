source("process.R")
source("functions-deletion.R")

export_date = as.Date("2015-03-20")

Demographics      <- processDemographics(         noOutput = T)
Demographics$Gender               <- factor(Demographics$Gender)
Demographics$Ancestral.Background <- factor(Demographics$Ancestral.Background)
Demographics$Country              <- factor(Demographics$Country)

Clinical          <- processFile("Clinical",      noOutput = T)
Developmental     <- processFile("Developmental", noOutput = T)

Genetics <- read.csv("Genetics.csv", stringsAsFactors = F)
Genetics <- filter(Genetics,Genetic.Status == "Results Verified")

patients <- intersect(unique(Genetics$Patient.ID),Clinical$Patient.ID)
Clinical <- Clinical[Clinical$Patient.ID %in% patients,]
Genetics <- Genetics[Genetics$Patient.ID %in% patients,]
Developmental <- Developmental[Developmental$Patient.ID %in% patients,]

# ==== Phenotypic variable selection ====
# Clinical
# Preselection
Clin_vars <- read.csv2("clin_vars.csv",      stringsAsFactors = F)
Clinical <- Clinical[c("Patient.ID",Clin_vars$Variable)]

# Group creation
Clin_vars_groups <- read.csv2("clin_vars_groups.csv", stringsAsFactors = F)
Clinical_groups <- select(Clinical, Patient.ID)
for (group in Clin_vars_groups$Variable)
{
  Clinical_groups[group] <- apply(select(Clinical, starts_with(group)), 1, function(x)
    {
      if (any(x == "Yes"))
        "Yes"
      else
        x[1]
    })
}

for (var in Clin_vars_groups$Variable)
{
  Clinical_groups[[var]][Clinical_groups[[var]] == "Not applicable"] <- NA
  for (lvl in levels(factor(Clinical_groups[[var]], exclude = "")))
    if (!is.na(lvl))
      Clin_vars_groups[Clin_vars_groups$Variable == var, lvl] <- summary(factor(Clinical_groups[[var]]))[lvl]
}
Clin_vars_groups <- filter(Clin_vars_groups, Yes > 5, No > 5)
Clinical_groups <- Clinical_groups[c("Patient.ID", Clin_vars_groups$Variable)]

# Variable and group selection
for (var in Clin_vars$Variable)
{
  Clinical[[var]][Clinical[[var]] == "Not applicable"] <- NA
  for (lvl in levels(factor(Clinical[[var]], exclude = "")))
    if (!is.na(lvl))
      Clin_vars[Clin_vars$Variable == var, lvl] <- summary(factor(Clinical[[var]]))[lvl]
}
write.csv2(Clin_vars, "clin_vars_summary.csv")
# Manual step of selecting variables
Clin_vars <- read.csv2("clin_vars_select.csv", stringsAsFactors = F)
Clinical <- Clinical[c("Patient.ID",Clin_vars$Variable)]
for (var in Clin_vars$Variable)
{
  Clinical[[var]][Clinical[[var]] == "Not applicable"] <- NA
  for (lvl in levels(factor(Clinical[[var]], exclude = "")))
    if (!is.na(lvl))
      Clin_vars[Clin_vars$Variable == var, lvl] <- summary(factor(Clinical[[var]]))[lvl]
}

# Developmental
# Preselection
Dev_vars  <- read.csv2("dev_vars.csv",       stringsAsFactors = F)
Developmental <- Developmental[c("Patient.ID", Dev_vars$Variable)]

# Group creation
Dev_vars_groups <- read.csv2("dev_vars_groups.csv", stringsAsFactors = F)
Developmental_groups <- select(Developmental, Patient.ID)
for (group in Dev_vars_groups$Variable)
{
  Developmental_groups[group] <- apply(select(Developmental, starts_with(group)), 1, function(x)
  {
    if (any(x == "Yes"))
      "Yes"
    else
      x[1]
  })
}

for (var in Dev_vars_groups$Variable)
{
  Developmental_groups[[var]][Developmental_groups[[var]] == "Not applicable" | Developmental_groups[[var]] == "N/A" | Developmental_groups[[var]] == "Not applicable (has severe motor impairment)"] <- NA
  for (lvl in levels(factor(Developmental_groups[[var]], exclude = "")))
    if (!is.na(lvl))
      Dev_vars_groups[Dev_vars_groups$Variable == var, lvl] <- summary(factor(Developmental_groups[[var]]))[lvl]
}
Dev_vars_groups <- filter(Dev_vars_groups, Yes > 5, No > 5)
Developmental_groups <- Developmental_groups[c("Patient.ID", Dev_vars_groups$Variable)]


# Variable and group selection
for (var in Dev_vars$Variable)
{
  Developmental[[var]][Developmental[[var]] == "Not applicable" | Developmental[[var]] == "N/A" | Developmental[[var]] == "Not applicable (has severe motor impairment)"] <- NA
  for (lvl in levels(factor(Developmental[[var]], exclude = "")))
    if (!is.na(lvl))
      Dev_vars[Dev_vars$Variable == var, lvl] <- summary(factor(Developmental[[var]]))[lvl]
}
write.csv2(Dev_vars, "dev_vars_summary.csv")
# Manual selection
Dev_vars <- read.csv2("dev_vars_select.csv", stringsAsFactors = F)
Developmental <- Developmental[c("Patient.ID", Dev_vars$Variable)]
for (var in Dev_vars$Variable)
{
  Developmental[[var]][Developmental[[var]] == "Not applicable" | Developmental[[var]] == "N/A" | Developmental[[var]] == "Not applicable (has severe motor impairment)"] <- NA
  for (lvl in levels(factor(Developmental[[var]], exclude = "")))
    if (!is.na(lvl))
      Dev_vars[Dev_vars$Variable == var, lvl] <- summary(factor(Developmental[[var]]))[lvl]
}

# ==== Genetics ====
Genetics_ranges    <- processRanges(Genetics)
Genetics_genes     <- processGenes(Genetics)
Genetics_genes_bin <- processGenes(Genetics, bin = T)
Genetics_pathways  <- processPathways(Genetics, Genetics_genes)

# Genes modified in at least 5 patients
Genes <- data.frame(t(apply(Genetics_genes_bin[-(1:2)], 2, function(x){summary(factor(x))[1:2]})))
Genes <- add_rownames(Genes, var = "Gene")
Genes <- filter(Genes, TRUE. > 5, FALSE. > 5)
Genetics_genes <- select(Genetics_genes, Patient.ID, one_of(Genes$Gene))
Genetics_genes_bin <- select(Genetics_genes_bin, Patient.ID, one_of(Genes$Gene))

# Pathways modified in at least 5 patients
Genetics_pathways_bin <- data.frame(apply(Genetics_pathways, 2, function(x){x > 0}))
names(Genetics_pathways_bin) <- names(Genetics_pathways)
Genetics_pathways_bin$Patient.ID <- Genetics_pathways$Patient.ID
Pathways <- data.frame(t(apply(Genetics_pathways_bin[-1], 2, function(x){summary(factor(x))})))
Pathways <- add_rownames(Pathways, var = "Pathway")
Pathways <- filter(Pathways, TRUE. > 5, FALSE. > 5)
Genetics_pathways <- select(Genetics_pathways, Patient.ID, one_of(Pathways$Pathway))
Genetics_pathways_bin <- select(Genetics_pathways_bin, Patient.ID, one_of(Pathways$Pathway))

# ==== Analyses ====
# Deletions
CNVplot(Genetics_ranges)
dir.create("delplots")
for (var in Clin_vars$Variable)
  Clin_vars$Range_p[Clin_vars$Variable == var] <- delPlot(Genetics_ranges, Clinical[c("Patient.ID", var)], noOutput = F)

for (var in Clin_vars_groups$Variable)
  Clin_vars_groups$Range_p[Clin_vars_groups$Variable == var] <- delPlot(Genetics_ranges, Clinical_groups[c("Patient.ID", var)], noOutput = F)

for (var in Dev_vars$Variable)
  Dev_vars$Range_p[Dev_vars$Variable == var] <- delPlot(Genetics_ranges, Developmental[c("Patient.ID", var)], noOutput = F)

for (var in Dev_vars_groups$Variable)
  Dev_vars_groups$Range_p[Dev_vars_groups$Variable == var] <- delPlot(Genetics_ranges, Developmental_groups[c("Patient.ID", var)], noOutput = F)

results_ranges <- rbind(Clin_vars[c("Group","Variable","Range_p")], Clin_vars_groups[c("Group","Variable","Range_p")], Dev_vars[c("Group","Variable","Range_p")], Dev_vars_groups[c("Group","Variable","Range_p")])
results_ranges$Range_p_corrected <- p.adjust(results$Range_p, "fdr")
results_ranges <- arrange(results, Range_p_corrected)
write.csv2(results_ranges, "results_ranges.csv")

# Genes
results_genes <- rbind(Clin_vars[1:2], Clin_vars_groups[1:2], Dev_vars[1:2], Dev_vars_groups[1:2])
for (var in Clin_vars$Variable)
{
  print(var)
  for (gene in Genes$Gene)
  {
    results_genes[results_genes$Variable == var, gene] <- fisher.test(factor(Clinical[[var]], exclude = c("", NA)), factor(Genetics_genes_bin[[gene]]), workspace = 1e8)$p.value
  }
}

for (var in Clin_vars_groups$Variable)
{
  print(var)
  for (gene in Genes$Gene)
  {
    results_genes[results_genes$Variable == var, gene] <- fisher.test(factor(Clinical_groups[[var]], exclude = c("", NA)), factor(Genetics_genes_bin[[gene]]), workspace = 1e8)$p.value
  }
}

Genetics_genes_bin <- filter(Genetics_genes_bin, Patient.ID %in% Developmental$Patient.ID)
for (var in Dev_vars$Variable)
{
  print(var)
  for (gene in Genes$Gene)
  {
    results_genes[results_genes$Variable == var, gene] <- fisher.test(factor(Developmental[[var]], exclude = c("", NA)), factor(Genetics_genes_bin[[gene]]), workspace = 1e8)$p.value
  }
}

for (var in Dev_vars_groups$Variable)
{
  print(var)
  for (gene in Genes$Gene)
  {
    results_genes[results_genes$Variable == var, gene] <- fisher.test(factor(Developmental_groups[[var]], exclude = c("", NA)), factor(Genetics_genes_bin[[gene]]), workspace = 1e8)$p.value
  }
}
results_genes_mat <- apply(results_genes[-(1:2)], 2, p.adjust, method = "bonferroni")
#results_genes_mat <- t(apply(results_genes_mat, 1, p.adjust, method = "fdr"))
rownames(results_genes_mat) <- results_genes$Variable
results_genes_mat <- results_genes_mat[which(apply(results_genes_mat, 1, function(x){min(x) < .05}) == TRUE), which(apply(results_genes_mat, 2, function(x){min(x) < .05}) == TRUE)]
results_genes_mat <- -log(results_genes_mat)
heatmap(results_genes_mat, col = colorRampPalette(c("white","blue"))(100))

# Pathways
results_pathways <- rbind(Clin_vars[1:2], Clin_vars_groups[1:2], Dev_vars[1:2], Dev_vars_groups[1:2])
for (var in Clin_vars$Variable)
{
  print(var)
  for (pathway in Pathways$Pathway)
  {
    results_pathways[results_pathways$Variable == var, pathway] <- fisher.test(factor(Clinical[[var]], exclude = c("", NA)), factor(Genetics_pathways_bin[[pathway]]), workspace = 1e8)$p.value
  }
}

for (var in Clin_vars_groups$Variable)
{
  print(var)
  for (pathway in Pathways$Pathway)
  {
    results_pathways[results_pathways$Variable == var, pathway] <- fisher.test(factor(Clinical_groups[[var]], exclude = c("", NA)), factor(Genetics_pathways_bin[[pathway]]), workspace = 1e8)$p.value
  }
}

Genetics_pathways_bin <- filter(Genetics_pathways_bin, Patient.ID %in% Developmental$Patient.ID)
for (var in Dev_vars$Variable)
{
  print(var)
  for (pathway in Pathways$Pathway)
  {
    results_pathways[results_pathways$Variable == var, pathway] <- fisher.test(factor(Developmental[[var]], exclude = c("", NA)), factor(Genetics_pathways_bin[[pathway]]), workspace = 1e8)$p.value
  }
}

for (var in Dev_vars_groups$Variable)
{
  print(var)
  for (pathway in Pathways$Pathway)
  {
    results_pathways[results_pathways$Variable == var, pathway] <- fisher.test(factor(Developmental_groups[[var]], exclude = c("", NA)), factor(Genetics_pathways_bin[[pathway]]), workspace = 1e8)$p.value
  }
}

results_pathways_mat <- apply(results_pathways[-(1:2)], 2, p.adjust, method = "bonferroni")
#results_pathways_mat <- t(apply(results_pathways_mat, 1, p.adjust, method = "fdr"))
rownames(results_pathways_mat) <- results_pathways$Variable
results_pathways_mat <- results_pathways_mat[which(apply(results_pathways_mat, 1, function(x){min(x) < .05}) == TRUE), which(apply(results_pathways_mat, 2, function(x){min(x) < .05}) == TRUE)]
results_pathways_mat <- -log(results_pathways_mat)
heatmap(results_pathways_mat, col = colorRampPalette(c("white","blue"))(100))
