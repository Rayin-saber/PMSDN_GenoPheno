library(tidyverse)

# ==== Demographics ====
read_csv("data/demographics.csv") %>%
  mutate_at(vars(Gender, Race, Country), factor) -> Demographics

# ==== Clinical ====
# Variables
read_csv2("clin_vars.csv") -> Clin_vars

# Data
read_csv("data/clinical.csv") %>%
# Variable selection
  select(Patient.ID, one_of(Clin_vars$Variable)) %>%
# Data types
  mutate(`What.was.the.patient's.gestational.age?`                                            = ifelse(`What.was.the.patient's.gestational.age?` == "Over 40 weeks", 41, as.numeric(sub(" weeks", "", `What.was.the.patient's.gestational.age?`))),
         `What.was.the.birth.weight.of.the.patient?.`                                         = ifelse(`What.was.the.birth.weight.of.the.patient?.` == "More than 5.44 kg", 5.44, as.numeric(sub(" kg", "", `What.was.the.birth.weight.of.the.patient?.`))),
         What.was.the.birth.length                                                            = ifelse(What.was.the.birth.length == "Greater than 65.41 centimeters", 65.41, as.numeric(What.was.the.birth.length)),
         `What.was.the.head.circumference?`                                                   = ifelse(`What.was.the.head.circumference?` == "Greater than 65.41 centimeters", 65.41, as.numeric(`What.was.the.head.circumference?`)),
         `What.is.the.highest.number.of.ear.infections.the.patient.has.had.in.a.given.year?`  = ifelse(`What.is.the.highest.number.of.ear.infections.the.patient.has.had.in.a.given.year?` == "7 - 9" | `What.is.the.highest.number.of.ear.infections.the.patient.has.had.in.a.given.year?` == "10 +", "7 +", `What.is.the.highest.number.of.ear.infections.the.patient.has.had.in.a.given.year?`)) %>%
  mutate_at(vars(`What.were.the.baby's.APGAR.scores.at.1.minute.and.5.minutes?_1.minute.APGAR.score`,
                 `What.were.the.baby's.APGAR.scores.at.1.minute.and.5.minutes?_5.minutes.APGAR.score`),
            as.numeric) %>%
  mutate_at(vars(`What.is.the.highest.number.of.ear.infections.the.patient.has.had.in.a.given.year?`),
            ~ ordered(. , levels = c("0", "1 - 4", "4 - 6", "7 +"))) %>%
  mutate_at(vars(starts_with("Which.of.the.following.medications")),
            ~ ordered(., levels = c("No reduction in seizures", "Some reduction (50%)", "Significant reduction (&gt;50%)", "Complete seizure control"))) %>%
  mutate_at(vars(contains("ketogenic")),
            ~ ordered(., levels = c("No effect", "Possibly effective", "Very effective"))) %>%
  # mutate_at(vars(contains("diagnosed.with.diabetes")),
            # ~ factor(., levels = c("Type I", "Type II", "Not applicable"))) %>%
  # Discard "by imputation", and coerce all remaining variables to Yes/No
  mutate_if(is.character,
            ~ str_replace_all(., fixed(" (by imputation)"), "") %>%
              factor(levels = c("No", "Yes"))) -> Clinical

# ==== Developmental ====
# Variables
read_csv2("dev_vars.csv") -> Dev_vars

read_csv("data/developmental.csv") %>%
  select(Patient.ID, one_of(Dev_vars$Variable)) %>%
# Data types
  mutate_if(~ str_detect(., "Always") %>% any(na.rm = T),
            ~ ordered(., levels = c("No", "Rarely", "Sometimes", "Always"))) %>%
  mutate_at(vars(starts_with("At.what.age.did")),
            ~ str_replace(., " years old", "") %>%
              as.numeric) %>%
  mutate_at(vars(matches("(communication|motor)\\.milestones")),
            ~ ordered(., levels = c("0 - 3 months", "4 - 7 months",	"8 - 11 months",	"12 - 18 months",	"19 - 24 months",	"25 - 30 months",	"31 - 36 months", "4 - 5 years", "6 - 7 years",	"8 - 9 years", "10 years +"))) %>%
  mutate_at(vars(`How.many.words.are.used.in.a.typical.sentence?_ currently`),
            ~ ordered(., levels = c("2 - 3 words", "4 - 5 words", "6 - 7 words", "8 - 10 words"))) %>%
  mutate_at(vars(`What.is.the.patient's.ability.to.follow.directions?_ currently`),
            ~ str_replace(., " \\(.*?\\)", "") %>%
              ordered(levels = c("Does not understand", "Understands and follows one-step command", "Understands and follows two-step command", "Understands and follows complex multi-step commands"))) %>%
  mutate_at(vars(`What.is.the.level.of.patient's.response.to.movement?_ currently`),
            ~ str_replace(., " \\(.*?\\)", "") %>%
              ordered(levels = c("Hypo-sensitive", "Normal", "Hyper-sensitive"))) %>%
  mutate_at(vars(`What.is.the.level.of.patient's.response.to.sound?_ currently`),
            ~ str_replace(., " \\(.*?\\)", "") %>%
              ordered(levels = c("Hyper-sensitive", "Normal", "Hypo-sensitive"))) %>%
  mutate_at(vars(`What.is.the.patient's.pain.tolerance.level?_ currently`),
            ~ ordered(., levels = c("Lower than usual", "Normal", "Higher than usual"))) %>%
  mutate_at(vars(`Is.the.patient.able.to.feed.himself/herself.independently?_1_Answer_ currently`),
            ~ ordered(., levels = c("No", "Without utensils (hand/finger feeding)", "With utensils", "Both (with utensils and hand/finder feeding)"), labels = c("No", "Without utensils", "With utensils", "Both"))) %>%
  mutate_at(vars(`How.long.were.special.feeds.required?`),
            ~ ordered(., levels = c("Up to 1 month", "1 - 3 months", "4 - 6 months", "7 - 11 months", "1 - 3 years", "4 + years"))) %>%
  mutate_at(vars(`With.PMS,.some.children.develop.understandable.verbal.speech.while.others.do.not..If.the.patient.is.over.age.2,.please.choose.the.response.that.most.closely.matches.your.child's.abilities.today._ currently`),
            ~ ordered(., levels = c("Non-verbal &ndash; does not use sounds, uses only gestures and pointing",
                                    "Non-verbal &ndash; uses sounds to communicate meaning, such as grunts",
                                    "Non-verbal &ndash; uses word-like sounds to communicate meaning but does not state words clearly",
                                    "Verbal &ndash; uses meaningful single words or phrases that are understandable to others",
                                    "Verbal &ndash; uses meaningful words and sentences that are understandable to others"))) %>%
  # Discard "by imputation", and coerce all remaining variables to Yes/No
  mutate_if(is.character,
            ~ str_replace_all(., fixed(" (by imputation)"), "") %>%
              factor(levels = c("No", "Yes"))) -> Developmental

# ==== Genetics ====
read_csv("data/genetic_hg38.csv") %>%
  select(Patient.ID,
         Result.type,
         Gain_Loss = Array.Type,
         Start = Array.Start.latest,
         End = Array.End.latest,
         seq_pos = SequencingGenomicChange) %>%
  mutate(seq_pos = seq_pos %>%
                     str_replace_all("[ ,]", "") %>%
                     str_replace("^chr22\\(GRC[Hh]37\\):", "") %>%
                     str_replace("^g\\.", "") %>%
                     str_replace("^chr22:", "") %>%
                     str_extract("^\\d+") %>%
                     as.numeric) %>%
  mutate(seq_pos = seq_pos - 51113070 + 50674642) %>% # convert hg19 -> hg38 for mutation positions
  mutate(Start = ifelse(Result.type == "Sequencing", seq_pos, Start),
         End = ifelse(Result.type == "Sequencing", 50733210, End)) %>%
  filter(! Start %>% is.na) %>%
  select(-seq_pos) %>%
  mutate(Result.type = Result.type %>% fct_recode("coordinates" = "Array",
                                                  "mutation" = "Sequencing")) -> Genetics_ranges

# ==== Export ====
saveRDS(list(Clin_vars = Clin_vars,
             Dev_vars = Dev_vars,
             Clinical = Clinical,
             Developmental = Developmental,
             Demographics = Demographics,
             Genetics_ranges = Genetics_ranges),
        file = "data.rds")
