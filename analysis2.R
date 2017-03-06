library(tidyverse)
library(stringr)
library(magrittr)
library(cowplot)
library(svglite)
library(grid)
library(ordinal)

load("data.Rda")

# Create main data frame -------------------------------------------------------
Demographics %>%
  inner_join(Clinical) %>%
  full_join(Developmental) %>%
  ungroup -> data

rbind(Clin_vars, Dev_vars) -> vars

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

# Manage Genetics_ranges -------------------------------------------------------
Genetics_ranges %<>%
  filter(Chr_Gene == "22",
         Gain_Loss == "Loss",
         End > 50674641) %>% # shank3 : 50674641
  mutate(Patient.ID = as.character(Patient.ID)) %>%
  mutate(Gain_Loss = ifelse(Result.type == "mutation", "Mutation", Gain_Loss)) %>%
  select(-Result.type, -Chr_Gene)

# Keep only patients who have phenotypic data ----------------------------------
Genetics_ranges %<>% semi_join(data)

# Re-compute del extent var ----------------------------------------------------
Genetics_ranges %<>%
  full_join(Genetics_ranges %>%
            group_by(Patient.ID) %>%
            summarize(min = min(Start), End = max(End))) %>%
  arrange(desc(min))

# Keep only patients with geno & pheno data and build dataframe ----------------
data %<>%
  inner_join(Genetics_ranges %>%
             distinct(Patient.ID, .keep_all = T) %>%
             select(Patient.ID, min, End)) %>%
  arrange(desc(min))


rm(Genetics_ranges)

# Analyses --------------------------------

data %>%
  select(Patient.ID, 
         verbal = `With.PMS,.some.children.develop.understandable.verbal.speech.while.others.do.not..If.the.patient.is.over.age.2,.please.choose.the.response.that.most.closely.matches.your.child's.abilities.today._ currently`, 
         walk = `How.old.was.the.patient.when.the.following.gross.motor.milestones.were.met?_Walk.unassisted`,
         seizure = `Does.the.patient.have.seizures?_ currently`,
         regression = `Has.the.patient.ever.experienced.regression.in.cognitive.abilities?`,
         min, 
         End,
         Gender, 
         Age) %>%
  mutate(size = End - min) -> data2

labels <- c(verbal = "Ability to communicate verbally",
            walk = "Age of acquisition of unassisted walk",
            seizure = "Seizures",
            regression = "Regression in cognitive abilities")

c("verbal", "walk", "seizure", "regression") %>%
  map(function(x)
      {
        data2 %>%
          select(Patient.ID, min, End, var = one_of(x)) %>%
          ggplot() +
          aes(x = Patient.ID, ymin = min/1e6, ymax = End/1e6) +
          geom_linerange(size = 1, aes(color = var)) +
          # geom_rect(ymin = 50674641/1e6, ymax = 50733212/1e6, xmin = 1, xmax = nrow(data2 %>% distinct(Patient.ID)), alpha = .01, fill = "grey") +
          scale_x_discrete(labels = NULL, limits = unique(data2$Patient.ID[order(data2$size)])) +
          scale_color_grey(start = .9, end = .2, na.value = "white", name = labels[x]) +
          ylab("Chromosomic coordinates") +
          xlab(NULL) +
          theme(axis.ticks.y = element_blank(),
                axis.line.y = element_blank(),
                legend.position = c(0, 0),
                legend.justification = c(0, 0)) +
      coord_flip() +
      annotate("segment", x = nrow(data2 %>% distinct(Patient.ID)), xend = 1, y = 52, yend = 52, arrow = arrow(type = "closed", angle = 20)) +
      annotate("text", x = nrow(data2 %>% distinct(Patient.ID))/2, y = 51.95, label = "Patients", angle = 90, hjust = .5, vjust = 0) +
      annotate("text", x = nrow(data2 %>% distinct(Patient.ID))/2, y = 52.05, label = "Decreasing deletion size", angle = 90, hjust = .5, vjust = 1) +
      annotate("rect", ymin = 50674641/1e6, ymax = 50733212/1e6, xmin = 1, xmax = nrow(data2 %>% distinct(Patient.ID)), fill = rgb(0,0,0,0), color = "black")
      }) %>%
`names<-`(c("verbal", "walk", "seizure", "regression")) -> plots

c("verbal", "walk", "seizure", "regression") %>%
  map(~ggsave(device = svglite, plot = plots[[.x]], filename = str_c(.x, ".svg")))

clm(data = data2, formula = verbal ~ scale(size) + Gender + Age) -> model_verbal
clm(data = data2, formula = walk ~ scale(size) + Gender + Age) -> model_walk
glm(data = data2, formula = seizure ~ size + Gender + Age, family = "binomial") -> model_seizure
glm(data = data2, formula = regression ~ size + Gender + Age, family = "binomial") -> model_regression
  

summary(model_verbal)$coefficients["scale(size)","Pr(>|z|)"]
confint(model_verbal)["scale(size)",] %>% exp
coef(model_verbal)[["scale(size)"]] %>% exp

summary(model_walk)$coefficients["scale(size)","Pr(>|z|)"]
confint(model_walk)["scale(size)",] %>% exp
coef(model_walk)[["scale(size)"]] %>% exp

summary(model_seizure)$coefficients["size","Pr(>|z|)"]
((confint(model_seizure) %>% exp) ** 1e6)["size",]
(coef(model_seizure)[["size"]] %>% exp) ** 1e6

summary(model_regression)$coefficients["size","Pr(>|z|)"]
((confint(model_regression) %>% exp) ** 1e6)["size",]
(coef(model_regression)[["size"]] %>% exp) ** 1e6
