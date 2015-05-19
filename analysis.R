source("process.R")
source("functions-deletion.R")

export_date = as.Date("2015-03-20")

Demographics      <- processDemographics(         noOutput = T)
Demographics$Gender               <- factor(Demographics$Gender)
Demographics$Ancestral.Background <- factor(Demographics$Ancestral.Background)
Demographics$Country              <- factor(Demographics$Country)

Clinical          <- processFile("Clinical",      noOutput = T)
Developmental     <- processFile("Developmental", noOutput = T)

# Variable selection
#Clinical <- Clinical[, grep("_Other",  names(Clinical), invert = T)]
#Clinical <- Clinical[, grep("Unsure",  names(Clinical), invert = T)]
#Clinical <- Clinical[, grep("_at age", names(Clinical), invert = T)]
Clin_vars <- read.csv2("clin_vars.csv",      stringsAsFactors = F)
#Clinical <- select(Clinical, Patient.ID, one_of(Clin_vars$Variable))

#Developmental <- Developmental[, gerp("_Other", names(Developmental), invert = T)]
Dev_vars  <- read.csv2("dev_vars.csv",       stringsAsFactors = F)

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

