source("functions-deletion.R")

load("Phenotypes.Rda")
load("Genetics.Rda")

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

Genetics_genes_bin_ <- Genetics_genes_bin
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
results_genes_mat <- apply(results_genes[-(1:2)], 2, p.adjust, method = "fdr")
rownames(results_genes_mat) <- results_genes$Variable
results_genes_mat <- results_genes_mat[which(apply(results_genes_mat, 1, function(x){min(x) < .05}) == TRUE), which(apply(results_genes_mat, 2, function(x){min(x) < .05}) == TRUE)]
results_genes_mat <- -log(results_genes_mat)
heatmap(results_genes_mat, Colv = NA, col = colorRampPalette(c("white","blue"))(100))

# ==== Pathways ====
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
