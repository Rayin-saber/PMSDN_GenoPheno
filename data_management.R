source("functions-process.R")

export_date <- as.Date("2015-12-04")

# ==== Demographics ====
Demographics      <- processDemographics(         noOutput = T)

Demographics$Gender               <- factor(Demographics$Gender)
Demographics$Ancestral.Background <- factor(Demographics$Ancestral.Background)
Demographics$Country              <- factor(Demographics$Country)

# ==== Clinical ====
Clinical          <- processFile("Clinical",      noOutput = T)

# Preselection
Clin_vars <- read.csv2("clin_vars.csv",      stringsAsFactors = F)
Clinical <- Clinical[c("Patient.ID",Clin_vars$Variable)]

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
  {
    Clinical[var] <- gsub(" (by imputation)", "", Clinical[[var]], fixed = T)
    Clinical[var] <- factor(Clinical[[var]], levels = c("No", "Yes"), exclude = c("", NA))
  }

# ==== Developmental ====
Developmental     <- processFile("Developmental", noOutput = T)

# Preselection
Dev_vars  <- read.csv2("dev_vars.csv",       stringsAsFactors = F)
Developmental <- Developmental[c("Patient.ID", Dev_vars$Variable)]

# Data types
for (var in names(Developmental))
  if (any(Developmental[[var]] == "Always", na.rm = T))
    Developmental[var] <- factor(Developmental[[var]], levels = c("No", "Rarely", "Sometimes", "Always"), ordered = T)


for (var in grep("(communication|motor)\\.milestones", names(Developmental), value = T))
  Developmental[var] <- factor(Developmental[[var]], levels = c("0 - 3 months", "4 - 7 months",	"8 - 11 months",	"12 - 18 months",	"19 - 24 months",	"25 - 30 months",	"31 - 36 months", "4 - 5 years", "6 - 7 years",	"8 - 9 years", "10 years +"), ordered = T)


# Developmental$`How.many.words.are.used.in.a.typical.sentence?_ currently`[Developmental$`How.many.words.are.used.in.a.typical.sentence?_ currently` == "6 - 7 words" | Developmental$`How.many.words.are.used.in.a.typical.sentence?_ currently` == "8 - 10 words"] <- "6 + words"
Developmental$`How.many.words.are.used.in.a.typical.sentence?_ currently` <- factor(Developmental$`How.many.words.are.used.in.a.typical.sentence?_ currently`, levels = c("2 - 3 words", "4 - 5 words", "6 - 7 words", "8 - 10 words"), ordered = T)

Developmental$`What.is.the.patient's.ability.to.follow.directions?_ currently` <- factor(Developmental$`What.is.the.patient's.ability.to.follow.directions?_ currently`, levels = c("Does not understand", "Understands and follows one-step command (e.g., Go to chair)", "Understands and follows two-step command (e.g., Get your coat and go to door)", "Understands and follows complex multi-step commands"), ordered = T)
Developmental$`What.is.the.level.of.patient's.response.to.movement?_ currently` <- factor(Developmental$`What.is.the.level.of.patient's.response.to.movement?_ currently`, levels = c("Hypo-sensitive (e.g. seeks exaggerated, fast movement)", "Normal", "Hyper-sensitive (e.g. avoids playground equipment, prefers sedentary tasks)"), labels = c("Hypo-sensitive", "Normal", "Hyper-sensitive"), ordered = T)
Developmental$`What.is.the.level.of.patient's.response.to.sound?_ currently` <- factor(Developmental$`What.is.the.level.of.patient's.response.to.sound?_ currently`, levels = c("Hyper-sensitive (e.g. easily startles at loud noises)", "Normal", "Hypo-sensitive (e.g. does not appear to react to excessively loud noises)"), labels = c("Hypo-sensitive", "Normal", "Hyper-sensitive"), ordered = T)

# Developmental$`What.is.the.patient's.pain.tolerance.level?_ currently`[Developmental$`What.is.the.level.of.patient's.response.to.sound?_ currently` == "Lower than usual" | Developmental$`What.is.the.patient's.pain.tolerance.level?_ currently` == "Normal"] <- "Normal or lower than usual"
Developmental$`What.is.the.patient's.pain.tolerance.level?_ currently` <- factor(Developmental$`What.is.the.patient's.pain.tolerance.level?_ currently`, levels = c("Lower than usual", "Normal", "Higher than usual"), labels = c("Hypo-sensitive", "Normal", "Hyper-sensitive"), ordered = T)

Developmental$`With.PMS,.some.children.develop.understandable.verbal.speech.while.others.do.not..If.the.patient.is.over.age.2,.please.choose.the.response.that.most.closely.matches.your.child's.abilities.today._ currently` <- factor(Developmental$`With.PMS,.some.children.develop.understandable.verbal.speech.while.others.do.not..If.the.patient.is.over.age.2,.please.choose.the.response.that.most.closely.matches.your.child's.abilities.today._ currently`, levels = c("Non-verbal – does not use sounds, uses only gestures and pointing", "Non-verbal – uses sounds to communicate meaning, such as grunts", "Non-verbal – uses word-like sounds to communicate meaning but does not state words clearly", "Verbal – uses meaningful single words or phrases that are understandable to others", "Verbal – uses meaningful words and sentences that are understandable to others"), ordered = T)

Developmental$`Is.the.patient.able.to.feed.himself/herself.independently?_1_Answer_ currently` <- factor(Developmental$`Is.the.patient.able.to.feed.himself/herself.independently?_1_Answer_ currently`, levels = c("No", "Without utensils (hand/finger feeding", "With utensils", "Both (with utensils and hand/finger feeding"), labels = c("No", "Without utensils", "With utensils", "Both"), ordered = T)

for (var in names(Developmental[-1]))
  if (is.character(Developmental[[var]]))
  {
    Developmental[var] <- gsub(" (by imputation)", "", Developmental[[var]], fixed = T)
    Developmental[var] <- factor(Developmental[[var]], levels = c("No", "Yes"), exclude = c("", NA))
  }

# ==== Genetics ====
Genetics_ranges <- processGenetics(noOutput = T)

Genetics_ranges$Start <- as.numeric(Genetics_ranges$Start)
Genetics_ranges$End   <- as.numeric(Genetics_ranges$End)

# ==== Export ====
Clin_vars <- select(Clin_vars, Group, Variable)
Dev_vars <- select(Dev_vars, Group, Variable)

save(Clin_vars, Dev_vars, Clinical, Developmental, Demographics, Genetics_ranges, file = "data.Rda")
