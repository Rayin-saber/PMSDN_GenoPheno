source("functions-deletion.R")
library(glmnet)
library(doMC)

registerDoMC(cores = 8)

load("Phenotypes.Rda")
load("Genetics.Rda")

# Deletions
CNVplot(Genetics_ranges)

data <- merge(Demographics, Clinical,             all = T)
data <- merge(data,         Clinical,             all = T)
data <- merge(data,         Clinical_groups,      all = T)
data <- merge(data,         Developmental,        all = T)
data <- merge(data,         Developmental_groups, all = T)

rm(Demographics, Clinical, Clinical_groups, Developmental, Developmental_groups)

vars <- rbind(Clin_vars, Clin_vars_groups, Dev_vars, Dev_vars_groups)

rm(Clin_vars, Clin_vars_groups, Dev_vars, Dev_vars_groups)

results_ranges <- vars
results_ranges$Range_p <- NA
results_ranges$ICl <- NA
results_ranges$OR <- NA
results_ranges$ICu <- NA

for (var in vars$Variable)
  results_ranges[results_ranges$Variable == var, c("Range_p", "ICl", "OR", "ICu")] <- delAnalysis(Genetics_ranges, data, var)

results_ranges$Range_p_corrected <- p.adjust(results_ranges$Range_p, "fdr")
results_ranges <- arrange(results_ranges, Range_p_corrected)
write.csv2(results_ranges, "results_ranges.csv")

# Genes
results_lasso <- vars
results_lasso[Genes$Gene] <- 0

for (var in vars$Variable)
{
  print(var)
  df <- merge(data[c("Patient.ID", var)], Genetics_genes_bin)
  df <- df[complete.cases(df),]
  x <- as.matrix(df[,-(1:2)])
  multi <- F
  if (is.numeric(df[[var]]))
  {
    cv.fit <- cv.glmnet(x, df[[var]], alpha = .9, parallel = T)
  } else if (is.ordered(df[[var]]))
  {
    cv.fit <- cv.glmnet(x, df[[var]], alpha = .9, parallel = T, family = "multinomial", type.multinomial = "grouped")
    multi <- T
  } else
  {
    cv.fit <- cv.glmnet(x, df[[var]], alpha = .9, parallel = T, family = "binomial")
  }
  plot(cv.fit)
  if (multi)
  {
    lasso <- as.matrix(coef(cv.fit, s = "lambda.min")[[nlevels(df[[var]])]])
  } else
  {
    lasso <- as.matrix(coef(cv.fit, s = "lambda.min"))
  }
  results_lasso[results_lasso$Variable == var, -(1:2)] <- lasso[-1]
}


results_genes_p   <- vars
results_genes_ICl <- vars
results_genes_ICu <- vars
results_genes_OR  <- vars

for (var in vars$Variable)
{
  print(var)
  df <- merge(data, Genetics_genes_bin)
  for (gene in Genes$Gene)
  {
    if (is.numeric(df[[var]]))
    {
      model <- lm(df[[var]] ~ df[[gene]] + df$Gender)
      table <- regress.display(model)$table
    } else if (is.ordered(df[[var]]))
    {
      tryCatch({
        model <- polr(df[[var]] ~ df[[gene]] + df$Gender + scale(df$Age))
        table <- ordinal.or.display(model)
      }, error = function(e) {
        table <<- NA
      })
    } else
    {
      if (names(df[var]) == "Is.the.patient's.menstrual.cycle.regular?_ currently")
        model <- glm(df[[var]] ~ df[[gene]] + df$Age, family = binomial)
      else
        model <- glm(df[[var]] ~ df[[gene]] + df$Gender + df$Age, family = binomial)
      table <- logistic.display(model)$table
    }
    tryCatch({
      results_genes_OR[results_genes_OR$Variable == var, gene]   <- table["df[[gene]]TRUE", 1]
      results_genes_ICl[results_genes_ICl$Variable == var, gene] <- table["df[[gene]]TRUE", 2]
      results_genes_ICu[results_genes_ICu$Variable == var, gene] <- table["df[[gene]]TRUE", 3]
      results_genes_p[results_genes_p$Variable == var, gene]     <- table["df[[gene]]TRUE", 4]
    }, error = function(e) {
      results_genes_OR[results_genes_OR$Variable == var, gene]   <<- NA
      results_genes_ICl[results_genes_ICl$Variable == var, gene] <<- NA
      results_genes_ICu[results_genes_ICu$Variable == var, gene] <<- NA
      results_genes_p[results_genes_p$Variable == var, gene]     <<- NA
    })
  }
}

dir.create("delplots")
for (var in vars$Variable)
  delPlotGenes(Genetics_ranges, data, results_ranges, results_genes_p, results_lasso, Genes, var, noOutput = F)

dir.create("delplotsrangesvg")
for (var in vars$Variable)
  delPlotRange(Genetics_ranges, data, results_ranges, var, noOutput = T)

results_genes_mat <- t(apply(results_genes_p[-(1:2)], 1, p.adjust, method = "fdr"))
results_genes_mat_or <- apply(results_genes_OR[-(1:2)], 2, as.numeric)
rownames(results_genes_mat) <- results_genes_p$Variable
rownames(results_genes_mat_or) <- results_genes_OR$Variable
results_genes_mat_or <- results_genes_mat_or[which(apply(results_genes_mat, 1, function(x){min(x) < .05}) == TRUE), which(apply(results_genes_mat, 2, function(x){min(x) < .05}) == TRUE)]
results_genes_mat <- results_genes_mat[which(apply(results_genes_mat, 1, function(x){min(x) < .05}) == TRUE), which(apply(results_genes_mat, 2, function(x){min(x) < .05}) == TRUE)]
results_genes_mat <- -log(results_genes_mat, base = 10)
results_genes_mat_or <- abs(log(results_genes_mat_or, base = 10))
rowInd <- heatmap(results_genes_mat, margins = c(7,20), Colv = NA, col = colorRampPalette(c("white","blue"))(100))$rowInd
heatmap(results_genes_mat_or[rowInd,], margins = c(7,20), Colv = NA, Rowv = NA, col = colorRampPalette(c("white","blue"))(100))

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
