source("process.R")

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
Demographics <- Demographics[Demographics$Patient.ID %in% patients,]

# ==== Phenotypic variable selection ====
# ==== Clinical ====
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

# Data types
Clinical$`What.was.the.patient's.gestational.age?` <- ifelse(Clinical$`What.was.the.patient's.gestational.age?` == "Over 40 weeks", 41, as.numeric(sub(" weeks", "", Clinical$`What.was.the.patient's.gestational.age?`)))
Clinical$`What.was.the.birth.weight.of.the.patient?.` <- ifelse(Clinical$`What.was.the.birth.weight.of.the.patient?.` == "More than 5.44 kg", 5.44, as.numeric(sub(" kg", "", Clinical$`What.was.the.birth.weight.of.the.patient?.`)))
Clinical$What.was.the.birth.length <- ifelse(Clinical$What.was.the.birth.length == "Greater than 65.41 centimeters", 65.41, as.numeric(Clinical$What.was.the.birth.length))
Clinical$`What.was.the.head.circumference?` <- ifelse(Clinical$`What.was.the.head.circumference?` == "Greater than 65.41 centimeters", 65.41, as.numeric(Clinical$`What.was.the.head.circumference?`))
Clinical$`What.were.the.baby's.APGAR.scores.at.1.minute.and.5.minutes?_1.minute` <- as.numeric(Clinical$`What.were.the.baby's.APGAR.scores.at.1.minute.and.5.minutes?_1.minute`)
Clinical$`What.were.the.baby's.APGAR.scores.at.1.minute.and.5.minutes?_5.minutes` <- as.numeric(Clinical$`What.were.the.baby's.APGAR.scores.at.1.minute.and.5.minutes?_5.minute`)

Clinical$`What.is.the.highest.number.of.ear.infections.the.patient.has.had.in.a.given.year?`[Clinical$`What.is.the.highest.number.of.ear.infections.the.patient.has.had.in.a.given.year?` == "7 - 9" | Clinical$`What.is.the.highest.number.of.ear.infections.the.patient.has.had.in.a.given.year?` == "10 +"] <- "7 +"
Clinical$`What.is.the.highest.number.of.ear.infections.the.patient.has.had.in.a.given.year?` <- factor(Clinical$`What.is.the.highest.number.of.ear.infections.the.patient.has.had.in.a.given.year?`, levels = c("0", "1 - 4", "4 - 6", "7 +"), ordered = T)

for (var in names(Clinical[-1]))
  if (is.character(Clinical[[var]]))
    Clinical[var] <- factor(Clinical[[var]], levels = c("No", "Yes"), exclude = c("", NA))

for (var in names(Clinical_groups[-1]))
  Clinical_groups[var] <- factor(Clinical_groups[[var]], levels = c("No", "Yes"), exclude = c("", NA))

# Summary and selection for groups
for (var in Clin_vars_groups$Variable)
  for (lvl in levels(Clinical_groups[[var]]))
    if (!is.na(lvl))
      Clin_vars_groups[Clin_vars_groups$Variable == var, lvl] <- summary(factor(Clinical_groups[[var]]))[lvl]

Clin_vars_groups <- filter(Clin_vars_groups, Yes > 5, No > 5)
Clinical_groups <- Clinical_groups[c("Patient.ID", Clin_vars_groups$Variable)]

# Summary and selection for individual vars
for (var in names(Clinical[-1]))
  if (!is.numeric(Clinical[[var]]))
    if (any(summary(Clinical[[var]])[1:nlevels(Clinical[[var]])] < 6) | grepl("_Other", var))
      Clinical[[var]] <- NULL

Clin_vars <- Clin_vars[Clin_vars$Variable %in% names(Clinical[-1]),]

for (var in Clin_vars$Variable)
  for (lvl in levels(Clinical[[var]]))
    if (!is.na(lvl))
      Clin_vars[Clin_vars$Variable == var, lvl] <- summary(factor(Clinical[[var]]))[lvl]
write.csv2(Clin_vars, "clin_vars_summary.csv")

# ==== Developmental ====
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

# Data types
for (var in names(Developmental))
  if (any(Developmental[[var]] == "Always", na.rm = T))
    Developmental[var] <- factor(Developmental[[var]], levels = c("No", "Rarely", "Sometimes", "Always"), ordered = T)


for (var in grep("(communication|motor)\\.milestones", names(Developmental), value = T))
  Developmental[var] <- factor(Developmental[[var]], levels = c("0 - 3 months", "4 - 7 months",	"8 - 11 months",	"12 - 18 months",	"19 - 24 months",	"25 - 30 months",	"31 - 36 months", "4 - 5 years", "6 - 7 years",	"8 - 9 years", "10 years +"), ordered = T)

for (var in grep("(communication|motor)\\.milestones", names(Developmental), value = T))
{
  levels <- c("0 - 3 months", "4 - 7 months",	"8 - 11 months",	"12 - 18 months",	"19 - 24 months",	"25 - 30 months",	"31 - 36 months", "4 - 5 years", "6 - 7 years",	"8 - 9 years", "10 years +")
  unclassed <- unclass(Developmental[[var]])
  lvl_low <- 1
  while (sum(summary(Developmental[[var]])[1:lvl_low]) < 6 & lvl_low < 11)
    lvl_low <- lvl_low + 1

  lvl_high <- 11
  while (sum(summary(Developmental[[var]])[lvl_high:11]) < 6 & lvl_high > 1)
    lvl_high <- lvl_high - 1

  if (lvl_high - lvl_low < 2)
    next

  unclassed[unclassed < lvl_low] <- lvl_low
  unclassed[unclassed > lvl_high] <- lvl_high

  if (lvl_low > 1)
    levels[lvl_low] <- gsub("\\d+ - (\\d+ \\w+)", "\\1 -", levels[lvl_low], perl = T)
  if (lvl_high < 11)
    levels[lvl_high] <- gsub("(\\d+) - \\d+ (\\w+)", "\\1 \\2 +", levels[lvl_high], perl = T)
  Developmental[[var]] <- ordered(unclassed, labels = levels[lvl_low:lvl_high])
}

Developmental$`How.many.words.are.used.in.a.typical.sentence?_ currently`[Developmental$`How.many.words.are.used.in.a.typical.sentence?_ currently` == "6 - 7 words" | Developmental$`How.many.words.are.used.in.a.typical.sentence?_ currently` == "8 - 10 words"] <- "6 + words"
Developmental$`How.many.words.are.used.in.a.typical.sentence?_ currently` <- factor(Developmental$`How.many.words.are.used.in.a.typical.sentence?_ currently`, levels = c("2 - 3 words", "4 - 5 words", "6 + words"), ordered = T)

Developmental$`What.is.the.patient's.ability.to.follow.directions?_ currently` <- factor(Developmental$`What.is.the.patient's.ability.to.follow.directions?_ currently`, levels = c("Does not understand", "Understands and follows one-step command (e.g., Go to chair)", "Understands and follows two-step command (e.g., Get your coat and go to door)", "Understands and follows complex multi-step commands"), ordered = T)
Developmental$`What.is.the.level.of.patient's.response.to.movement?_ currently` <- factor(Developmental$`What.is.the.level.of.patient's.response.to.movement?_ currently`, levels = c("Hypo-sensitive (e.g. seeks exaggerated, fast movement)", "Normal", "Hyper-sensitive (e.g. avoids playground equipment, prefers sedentary tasks)"), labels = c("Hypo-sensitive", "Normal", "Hyper-sensitive"), ordered = T)
Developmental$`What.is.the.level.of.patient's.response.to.sound?_ currently` <- factor(Developmental$`What.is.the.level.of.patient's.response.to.sound?_ currently`, levels = c("Hyper-sensitive (e.g. easily startles at loud noises)", "Normal", "Hypo-sensitive (e.g. does not appear to react to excessively loud noises)"), labels = c("Hypo-sensitive", "Normal", "Hyper-sensitive"), ordered = T)

Developmental$`What.is.the.patient's.pain.tolerance.level?_ currently`[Developmental$`What.is.the.level.of.patient's.response.to.sound?_ currently` == "Lower than usual" | Developmental$`What.is.the.patient's.pain.tolerance.level?_ currently` == "Normal"] <- "Normal or lower than usual"
Developmental$`What.is.the.patient's.pain.tolerance.level?_ currently` <- factor(Developmental$`What.is.the.patient's.pain.tolerance.level?_ currently`, levels = c("Normal or lower than usual", "Higher than usual"), labels = c("Normal or hypo-sensitive", "Hyper-sensitive"))

Developmental$`With.PMS,.some.children.develop.understandable.verbal.speech.while.others.do.not..If.the.patient.is.over.age.2,.please.choose.the.response.that.most.closely.matches.your.child's.abilities.today._ currently` <- factor(Developmental$`With.PMS,.some.children.develop.understandable.verbal.speech.while.others.do.not..If.the.patient.is.over.age.2,.please.choose.the.response.that.most.closely.matches.your.child's.abilities.today._ currently`, levels = c("Non-verbal – does not use sounds, uses only gestures and pointing", "Non-verbal – uses sounds to communicate meaning, such as grunts", "Non-verbal – uses word-like sounds to communicate meaning but does not state words clearly", "Verbal – uses meaningful single words or phrases that are understandable to others", "Verbal – uses meaningful words and sentences that are understandable to others"), ordered = T)

Developmental$`Is.the.patient.able.to.feed.himself/herself.independently?_1_Answer_ currently` <- factor(Developmental$`Is.the.patient.able.to.feed.himself/herself.independently?_1_Answer_ currently`, levels = c("No", "Without utensils (hand/finger feeding", "With utensils", "Both (with utensils and hand/finger feeding"), labels = c("No", "Without utensils", "With utensils", "Both"), ordered = T)

for (var in names(Developmental[-1]))
  if (is.character(Developmental[[var]]))
    Developmental[var] <- factor(Developmental[[var]], levels = c("No", "Yes"), exclude = c("", NA))

for (var in names(Developmental_groups[-1]))
  if (is.character(Developmental_groups[[var]]))
    Developmental_groups[var] <- factor(Developmental_groups[[var]], levels = c("No", "Yes"), exclude = c("", NA))

for (var in Dev_vars_groups$Variable)
  for (lvl in levels(Developmental_groups[[var]]))
    if (!is.na(lvl))
      Dev_vars_groups[Dev_vars_groups$Variable == var, lvl] <- summary(factor(Developmental_groups[[var]]))[lvl]
Dev_vars_groups <- filter(Dev_vars_groups, Yes > 5, No > 5)
Developmental_groups <- Developmental_groups[c("Patient.ID", Dev_vars_groups$Variable)]

# Variable and group selection

# Summary and selection for individual vars
for (var in names(Developmental[-1]))
  if (!is.numeric(Developmental[[var]]))
    if (any(summary(Developmental[[var]])[1:nlevels(Developmental[[var]])] < 6) | grepl("_Other", var))
      Developmental[[var]] <- NULL

    Dev_vars <- Dev_vars[Dev_vars$Variable %in% names(Developmental[-1]),]

for (var in Dev_vars$Variable)
  for (lvl in levels(Developmental[[var]]))
    if (!is.na(lvl))
      Dev_vars[Dev_vars$Variable == var, lvl] <- summary(factor(Developmental[[var]]))[lvl]
write.csv2(Dev_vars, "dev_vars_summary.csv")


Clin_vars <- select(Clin_vars, Group, Variable)
Clin_vars_groups <- select(Clin_vars_groups, Group, Variable)
Dev_vars <- select(Dev_vars, Group, Variable)
Dev_vars_groups <- select(Dev_vars_groups, Group, Variable)
save(Clin_vars, Clin_vars_groups, Dev_vars, Dev_vars_groups, Clinical, Clinical_groups, Developmental, Developmental_groups, Demographics, file = "Phenotypes.Rda")

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

# Reorder genes based on position on chr22
hg38 <- read.delim("refGene.txt.hg38", stringsAsFactors = F)
hg38 <- filter(hg38, name2 %in% Genes$Gene, chrom == "chr22") %>% select(chrom, txStart, txEnd, name2) %>% arrange(name2) %>% distinct %>% group_by(name2) %>% summarise(txStart = min(txStart), txEnd = max(txEnd))
Genes <- merge(Genes, hg38, by.x = "Gene", by.y = "name2") %>% arrange(txStart)
Genetics_genes_bin <- select(Genetics_genes_bin, Patient.ID, one_of(Genes$Gene))

Genetics_genes_bin$Patient.ID    <- as.character(Genetics_genes_bin$Patient.ID)
Genetics_pathways$Patient.ID     <- as.character(Genetics_pathways$Patient.ID)
Genetics_pathways_bin$Patient.ID <- as.character(Genetics_pathways_bin$Patient.ID)
Genetics_ranges$Patient.ID       <- as.character(Genetics_ranges$Patient.ID)
Genetics_ranges$Start            <- as.numeric(Genetics_ranges$Start)
Genetics_ranges$End              <- as.numeric(Genetics_ranges$End)

save(Genes, Genetics_genes_bin, Genetics_pathways, Genetics_pathways_bin, Genetics_ranges, file = "Genetics.Rda")
