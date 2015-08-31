source("functions-deletion.R")
library(glmnet)
library(doMC)

registerDoMC(cores = 8)

load("Phenotypes.Rda")
load("Genetics.Rda")
Genetics_ranges <- Genetics
rm(Genetics)

# Deletions----------------------------------------------------------------------------------
data <- merge(Demographics, Clinical,             all = T)
data <- merge(data,         Clinical,             all = T)
#data <- merge(data,         Clinical_groups,      all = T)
data <- merge(data,         Developmental,        all = T)
#data <- merge(data,         Developmental_groups, all = T)

rm(Demographics, Clinical, Developmental)

vars <- rbind(Clin_vars, Dev_vars)

rm(Clin_vars, Dev_vars)

# First plot of the deletions
CNVplot(Genetics_ranges)

# Keep only terminal deletions/mutations on chr22
Genetics_ranges <- filter(Genetics_ranges, Chr_Gene == "22", Gain_Loss != "Gain", End > 50000000)
data <- filter(data, Patient.ID %in% Genetics_ranges$Patient.ID)

# Second plot of the deletions
CNVplot(Genetics_ranges)

results_ranges <- vars
results_ranges$Range_p <- NA
results_ranges$ICl <- NA
results_ranges$OR <- NA
results_ranges$ICu <- NA

for (var in vars$Variable)
  results_ranges[results_ranges$Variable == var, c("Range_p", "ICl", "OR", "ICu")] <- delAnalysis(Genetics_ranges, data, var)

rm(df2, depvar)

results_ranges$Range_p_corrected <- p.adjust(results_ranges$Range_p, "fdr")
results_ranges <- arrange(results_ranges, Range_p_corrected)
write.csv2(results_ranges, "results_ranges.csv")

# Plots--------------------------------------------------------------------------
dir.create("delplotsrangesvg")
for (var in vars$Variable[1:100])
  delPlotRange(Genetics_ranges, data, results_ranges, var, noOutput = T, bnw = T)

for (var in vars$Variable[vars$Group == "Renal.-.Kidney"])
  delPlotRange(Genetics_ranges, data, results_ranges, var, noOutput = T, bnw = T)

for (var in vars$Variable[vars$Group == "Gross.Motor.Development"])
  delPlotRange(Genetics_ranges, data, results_ranges, var, noOutput = T, bnw = T)

# ROC curves---------------------------------------------------------------------
# Keep only the maximum extent of deletion for each patient
library(pROC)
Start <- Genetics_ranges %>% group_by(Patient.ID) %>% summarize(Start = min(Start))
vars2 <- results_ranges$Variable[results_ranges$Range_p_corrected < .05]
data$Patient.ID <- as.numeric(data$Patient.ID)
dataroc <- merge(Start,data[c("Patient.ID",vars2)])
for (var in vars2)
  if (nlevels(dataroc[[var]]) > 2)
    dataroc[var] <- NULL

for (var in names(dataroc[-(1:2)]))
  plot(roc(dataroc[[var]], dataroc$Start, percent = T), print.thres = "best", print.auc = T, main = var)

# Lasso---------------------------------------------------------------------------------------
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

rm(lasso, x, cv.fit, multi, df)

# GWAS----------------------------------------------------------------------------------
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

rm(table, model)

dir.create("delplots")
for (var in vars$Variable)
  delPlotGenes(Genetics_ranges, data, results_ranges, results_genes_p, results_lasso, Genes, var, noOutput = F)


results_genes_mat <- t(apply(results_genes_p[-(1:2)], 1, p.adjust, method = "bonferroni"))
#results_genes_mat_or <- apply(results_genes_OR[-(1:2)], 2, as.numeric)
rownames(results_genes_mat) <- results_genes_p$Variable
#rownames(results_genes_mat_or) <- results_genes_OR$Variable
#results_genes_mat_or <- results_genes_mat_or[which(apply(results_genes_mat, 1, function(x){min(x) < .05}) == TRUE), which(apply(results_genes_mat, 2, function(x){min(x) < .05}) == TRUE)]
results_genes_mat <- results_genes_mat[which(apply(results_genes_mat, 1, function(x){min(x) < .05}) == TRUE), which(apply(results_genes_mat, 2, function(x){min(x) < .05}) == TRUE)]
results_genes_mat <- -log(results_genes_mat, base = 10)
#results_genes_mat_or <- abs(log(results_genes_mat_or, base = 10))
rowInd <- heatmap(results_genes_mat, margins = c(7,20), Colv = NA, col = colorRampPalette(c("white","blue"))(100))$rowInd
#heatmap(results_genes_mat_or[rowInd,], margins = c(7,20), Colv = NA, Rowv = NA, col = colorRampPalette(c("white","blue"))(100))

# Misc-----------------------------------------------------------------------
lasso <- as.matrix(results_lasso[-(1:2)])
rownames(lasso) <- results_lasso$Variable
lasso[lasso != 0] <- 1
lasso <- lasso * results_genes_mat
lasso <- 10 ** (-lasso)

lasso <- lasso[apply(lasso, 1, function(x) {any(x < .05, na.rm = T)}), apply(lasso, 2, function(x) {any(x < .05, na.rm = T)})]

library(RCurl)

pubmed <- Genes["Gene"]
for (gene in Genes$Gene)
{
  esearch <- as.character(getURLContent(paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmax=300&term=",gene)))
  esearch <- gsub("\\n", "", esearch, perl = T)
  count <- as.numeric(sub("^.*<eSearchResult><Count>(\\d+)<\\/Count>.*$","\\1", esearch, perl = T))
  if (count == 0)
    next
  esearch <- gsub("^.*<IdList><Id>", "", esearch, perl = T)
  esearch <- gsub("<\\/Id><\\/IdList>.*$", "", esearch, perl = T)
  esearch <- gsub("<\\/Id><Id>", ",", esearch, perl = T)
  pubmed$pmid[pubmed$Gene == gene] <- esearch

  efetch <- as.character(getURLContent(paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id=",esearch)))
}

# ==== Pathways ====
# results_pathways <- rbind(Clin_vars[1:2], Clin_vars_groups[1:2], Dev_vars[1:2], Dev_vars_groups[1:2])
# for (var in Clin_vars$Variable)
# {
#   print(var)
#   for (pathway in Pathways$Pathway)
#   {
#     results_pathways[results_pathways$Variable == var, pathway] <- fisher.test(factor(Clinical[[var]], exclude = c("", NA)), factor(Genetics_pathways_bin[[pathway]]), workspace = 1e8)$p.value
#   }
# }
#
# for (var in Clin_vars_groups$Variable)
# {
#   print(var)
#   for (pathway in Pathways$Pathway)
#   {
#     results_pathways[results_pathways$Variable == var, pathway] <- fisher.test(factor(Clinical_groups[[var]], exclude = c("", NA)), factor(Genetics_pathways_bin[[pathway]]), workspace = 1e8)$p.value
#   }
# }
#
# Genetics_pathways_bin <- filter(Genetics_pathways_bin, Patient.ID %in% Developmental$Patient.ID)
# for (var in Dev_vars$Variable)
# {
#   print(var)
#   for (pathway in Pathways$Pathway)
#   {
#     results_pathways[results_pathways$Variable == var, pathway] <- fisher.test(factor(Developmental[[var]], exclude = c("", NA)), factor(Genetics_pathways_bin[[pathway]]), workspace = 1e8)$p.value
#   }
# }
#
# for (var in Dev_vars_groups$Variable)
# {
#   print(var)
#   for (pathway in Pathways$Pathway)
#   {
#     results_pathways[results_pathways$Variable == var, pathway] <- fisher.test(factor(Developmental_groups[[var]], exclude = c("", NA)), factor(Genetics_pathways_bin[[pathway]]), workspace = 1e8)$p.value
#   }
# }
#
# results_pathways_mat <- apply(results_pathways[-(1:2)], 2, p.adjust, method = "bonferroni")
# #results_pathways_mat <- t(apply(results_pathways_mat, 1, p.adjust, method = "fdr"))
# rownames(results_pathways_mat) <- results_pathways$Variable
# results_pathways_mat <- results_pathways_mat[which(apply(results_pathways_mat, 1, function(x){min(x) < .05}) == TRUE), which(apply(results_pathways_mat, 2, function(x){min(x) < .05}) == TRUE)]
# results_pathways_mat <- -log(results_pathways_mat)
# heatmap(results_pathways_mat, col = colorRampPalette(c("white","blue"))(100))
