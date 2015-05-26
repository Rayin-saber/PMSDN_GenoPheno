source("process.R")
source("functions-deletion.R")

export_date = as.Date("2015-03-20")

Demographics      <- processDemographics(         noOutput = T)
Demographics$Gender               <- factor(Demographics$Gender)
Demographics$Ancestral.Background <- factor(Demographics$Ancestral.Background)
Demographics$Country              <- factor(Demographics$Country)

Clinical          <- processFile("Clinical",      noOutput = T)
Developmental     <- processFile("Developmental", noOutput = T)

# ==== Phenotypic variable selection ====
# Clinical
# Preselection
Clin_vars <- read.csv2("clin_vars.csv",      stringsAsFactors = F)
Clinical <- Clinical[c("Patient.ID",Clin_vars$Variable)]

# Group creation
Clin_vars_groups <- read.csv2("clin_vars_groups.csv", stringsAsFactors = F, header = F)
Clinical_groups <- select(Clinical, Patient.ID)
for (group in Clin_vars_groups$V2)
{
  Clinical_groups[group] <- apply(select(Clinical, starts_with(group)), 1, function(x)
    {
      if (any(x == "Yes"))
        "Yes"
      else
        x[1]
    })
}

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

for (var in Clin_vars_groups$V2)
{
  Clinical_groups[[var]][Clinical_groups[[var]] == "Not applicable"] <- NA
  for (lvl in levels(factor(Clinical_groups[[var]], exclude = "")))
    if (!is.na(lvl))
      Clin_vars_groups[Clin_vars_groups$V2 == var, lvl] <- summary(factor(Clinical_groups[[var]]))[lvl]
}
Clin_vars_groups <- filter(Clin_vars_groups, Yes > 5, No > 5)
Clinical_groups <- Clinical_groups[c("Patient.ID", Clin_vars_groups$V2)]

# Developmental
# Preselection
Dev_vars  <- read.csv2("dev_vars.csv",       stringsAsFactors = F)
Developmental <- Developmental[c("Patient.ID", Dev_vars$Variable)]

# Group creation
Dev_vars_groups <- read.csv2("dev_vars_groups.csv", stringsAsFactors = F, header = F)
Developmental_groups <- select(Developmental, Patient.ID)
for (group in Dev_vars_groups$V2)
{
  Developmental_groups[group] <- apply(select(Developmental, starts_with(group)), 1, function(x)
  {
    if (any(x == "Yes"))
      "Yes"
    else
      x[1]
  })
}

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

for (var in Dev_vars_groups$V2)
{
  Developmental_groups[[var]][Developmental_groups[[var]] == "Not applicable" | Developmental_groups[[var]] == "N/A" | Developmental_groups[[var]] == "Not applicable (has severe motor impairment)"] <- NA
  for (lvl in levels(factor(Developmental_groups[[var]], exclude = "")))
    if (!is.na(lvl))
      Dev_vars_groups[Dev_vars_groups$V2 == var, lvl] <- summary(factor(Developmental_groups[[var]]))[lvl]
}
Dev_vars_groups <- filter(Dev_vars_groups, Yes > 5, No > 5)
Developmental_groups <- Developmental_groups[c("Patient.ID", Dev_vars_groups$V2)]

# ==== Genetics ====
Genetics <- read.csv("Genetics.csv", stringsAsFactors = F)
Genetics <- filter(Genetics,Genetic.Status == "Results Verified")
Genetics_ranges   <- processRanges(Genetics)
Genetics_genes    <- processGenes(Genetics)
Genetics_pathways <- processPathways(Genetics, Genetics_genes)

for (var in Clin_vars$Variable)
  Clin_vars$Range_p[Clin_vars$Variable == var] <- delPlot(Genetics_ranges, Clinical[c("Patient.ID",var)], noOutput = F)


# Analyses
Clinical <- Clinical[, grep("_Other",  names(Clinical), invert = T)]
Clinical <- Clinical[, grep("Unsure",  names(Clinical), invert = T)]
Clinical <- Clinical[, grep("_at age", names(Clinical), invert = T)]
Clinical <- as.data.frame(apply(Clinical, 2, factor, exclude = ""))

Clin_yes <- apply(Clinical, 2, function(x){summary(factor(x))["Yes"] < 5})
names(Clin_yes[Clin_yes == T])[!is.na(names(Clin_yes[Clin_yes == T]))]
summary(Clinical)

Adult <- Adult[, grep("at\\.what\\.age\\?", names(Adult), invert = T)]
Adult <- Adult[, grep("_at age",            names(Adult), invert = T)]
Adult <- as.data.frame(apply(Adult, 2, factor, exclude = ""))

Adul_yes <- apply(Adult, 2, function(x){summary(factor(x))["Yes"] < 5})
names(Adul_yes[Adul_yes == T])[!is.na(names(Adul_yes[Adul_yes == T]))]
summary(Adult)


Developmental <- Developmental[, grep("_at age", names(Developmental), invert = T)]
Developmental <- Developmental[, grep("_Other",  names(Developmental), invert = T)]
Developmental <- as.data.frame(apply(Developmental, 2, factor, exclude = ""))

Dev_yes <- apply(Developmental, 2, function(x){summary(factor(x))["Yes"] < 5})
names(Dev_yes[Dev_yes == T])[!is.na(names(Dev_yes[Dev_yes == T]))]
summary(Developmental[1:100])

# Genes modified in at least 5 patients
length(which(apply(Genetics_genes, 2, function(x){summary(factor(x))["2"] <= nrow(Genetics_genes) - 5})))

# Pathways modified in at least 5 patients
length(which(apply(Genetics_pathways, 2, function(x){summary(factor(x))["0"] <= nrow(Genetics_pathways) - 5})))

