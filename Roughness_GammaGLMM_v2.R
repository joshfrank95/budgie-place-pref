

# --- 0. Load libraries, set paths ####
# 1. Base / stats packages that conflict with dplyr (load FIRST)
library(MASS) # Prevents MASS::select from overwriting dplyr::select

# 2. Core Tidyverse (load AFTER base/modeling packages)
library(tidyverse) # Loads dplyr, ggplot2, tidyr, lubridate, etc.

# 3. Auxiliary plotting & modeling packages
library(car)
library(fitdistrplus)
library(glmmTMB)
library(DHARMa)
library(parameters)
library(performance)
library(insight)
library(bayesplot)
library(posterior)
library(cowplot)

# 4. Explicitly resolve key function conflicts
select <- dplyr::select
summarise <- dplyr::summarise
summarize <- dplyr::summarize
filter <- dplyr::filter


#clean space
#ls()
#rm(list=ls())

#directory = ''
#Humans or Birds
Exp_I_raw = read.csv('.../Roughness_detailed_bird.csv', header = T, sep = ',')
Exp_II = read.csv('.../Roughness_detailed_human.csv', header = T, sep = ',')
Exp_III = read.csv('.../ARI Roughness/data/Exp_III/Exp_III_Birds.csv', header = T, sep = ',')

# --- 1. Data Preparation ####
Exp_I_raw <- Exp_I_raw %>% filter(time > 0.3)
Exp_I_raw <- Exp_I_raw %>% filter(condition != "Free")
#Exp_I_raw <- Exp_I_raw %>% filter(session > 1)

#Exp_I_raw <- Exp_I_raw %>% filter(subject != "Famine") # Famine: Clearly not understood apparatus. 
# Famine session 1: only unpleasant. session 2: Only silent & unpleasant. Session 3: Only unpleasant

Exp_II <- Exp_II %>% filter(time > 0.3)
Exp_II <- Exp_II %>% filter(condition != "FREE_STATE")

Exp_III <- Exp_III %>% filter(time > 0.3)
Exp_III <- Exp_III %>% filter(condition != "Free")

Exp_I_raw$condition <- fct_relevel(Exp_I_raw$condition,"Pleasant","Unpleasant","Silence")
Exp_II$condition <- fct_relevel(Exp_II$condition,"Pleasant","Unpleasant","Silence")
Exp_III$condition <- fct_relevel(Exp_III$condition,"Pleasant","Unpleasant","Silence")


#colony_R_subjects <- c("Banach", "Mercury", "Bender", "Moonbeam", "Lemon", 
#                       "Hades", "Hobo", "Houdini", "Summer", "Famine", 
#                       "Pebbles", "Storm", "Afina", "Hilbert", "Fourier", "Pistachio")
#colony_L_subjects <- c("Kiwi", "Peanut", "Sir", "Bilbo", "Dino", "Smoothie", "Gretchen", 
#                       "Zeilinger", "Ice", "Bandit", "Dynamo", "Popeye", "Euler")
#Exp_I_raw <- Exp_I_raw %>% mutate(colony = case_when(
#  subject %in% colony_R_subjects ~ "R",
#  subject %in% colony_L_subjects ~ "L",
#  TRUE ~ NA_character_))

female_subjects <- c("Dino", "Kiwi", "Gretchen", "Lemon", "Summer", "Famine", 
                     "Pebbles", "Storm","Smoothie","Houdini","Pistachio",
                     "Dove", "Celeste", "Pineapple", "Calico", "Micro",
                     "Elisabeth"
                     )
male_subjects <- c("Peanut", "Sir", "Bilbo", "Zeilinger", "Ice", "Bandit", 
                   "Dynamo", "Popeye", "Euler", "Banach", "Mercury", "Bender", 
                   "Moonbeam", "Hades", "Hobo", "Afina", "Hilbert", "Fourier",
                   "Winter", "Pippin", "Pikachu", "Pickachu", "Wolferl")


Exp_I <- Exp_I_raw %>% mutate(sex = case_when(
  subject %in% female_subjects ~ "female",
  subject %in% male_subjects ~ "male",
  TRUE ~ NA_character_))

unique(Exp_I$subject)
unique(Exp_I$sex)
Exp_I %>%
  distinct(subject, sex) %>%  # Keep only one row per person
  count(sex)

# Quick check whether all heard everything => Yes.
# Gretchen didn't hear Silence but did hear all conditions
# that were tested here
subjects_to_keep <- Exp_I %>%
  group_by(subject) %>%
  summarise(n_cond = n_distinct(condition)) %>%
  #filter(n_cond >= 3) %>% 
  filter(n_cond < 3) %>% 
  pull(subject)

# More cleanup Exp. III
Exp_III <- Exp_III %>% mutate(sex = case_when(
  subject %in% female_subjects ~ "female",
  subject %in% male_subjects ~ "male",
  TRUE ~ NA_character_))
subjects_to_keep_III <- Exp_III %>% # Micro didn't hear everything. She only sat un Unpleasant the whole 3 sessions
  group_by(subject) %>%
  summarise(n_cond = n_distinct(condition)) %>%
  #filter(n_cond >= 3) %>% 
  filter(n_cond < 3) %>% 
  pull(subject)

Exp_III <- Exp_III %>%
  filter(subject != "Micro")

#Aggregating Data.
ag1 <-with(Exp_I, aggregate(cbind(time)~subject+session+condition+sex, FUN="sum"))
ag2 <-with(Exp_II, aggregate(cbind(time)~subject+condition+sex+session, FUN="sum"))
ag3 <-with(Exp_III, aggregate(cbind(time)~subject+condition+sex+session, FUN="sum"))


# --- 2. Descriptive Analyses ####

summary(Exp_I)
summary(Exp_II)
summary(Exp_III)

# Data distribution: We suspected a lognormal Distribution
performance::check_distribution(Exp_I$time)
performance::check_distribution(Exp_II$time)
performance::check_distribution(Exp_III$time)

# Miscellaneous

# Count number of switches between conditions => get mean. (-1 accounts for the first row, i.e., first landing)
row_counts_I <- Exp_I %>% group_by(subject, session) %>% summarise(number_of_rows_I = dplyr::n()-1)
meanSwitch_I = mean(row_counts_I$number_of_rows_I)
row_counts_II <- Exp_II %>% group_by(subject, session) %>% summarise(number_of_rows_II = dplyr::n()-1)
meanSwitch_II = mean(row_counts_II$number_of_rows_II)
row_counts_III <- Exp_III %>% group_by(subject, session) %>% summarise(number_of_rows_III = dplyr::n()-1)
meanSwitch_III = mean(row_counts_III$number_of_rows_III)
print(meanSwitch_I)
print(meanSwitch_II)
print(meanSwitch_III)

# Mean number of switches per session
# print(meanSwitch_I)
# 84.25862
# print(meanSwitch_II)
# 27.9322
# print(meanSwitch_III)
# 46.36364

# ============ =

# Time not spent on any perch

Exp_Ifree = read.csv('.../Roughness_detailed_bird.csv', header = T, sep = ',')
Exp_Ifree <- Exp_Ifree %>% aggregate(cbind(time)~subject, FUN="sum")
Exp_Ifree$time_subtracted <- (21600 - Exp_Ifree$time) # The `time_subtracted` column now holds the new values.
Exp_Ifree_mean = mean(Exp_Ifree$time_subtracted, na.rm = TRUE)
print(Exp_Ifree_mean) # Exp. I: "Free" was unfortunately not explicitly logged, thus this approach

# Experiment I (Birds)
# 4090.583 seconds

Exp_IIfree = read.csv('.../Roughness_detailed_human.csv', header = T, sep = ',')
Exp_IIfree <- Exp_IIfree %>% filter(condition == "FREE_STATE") %>% aggregate(cbind(time)~subject, FUN="sum") 
Exp_IIfree_mean = mean(Exp_IIfree$time, na.rm = TRUE)
print(Exp_IIfree_mean) # Exp II: "Free" was logged

# Experiment II (Humans)
# 289.2848 seconds

Exp_IIIfree = read.csv('.../data/Exp_III/Exp_III_Birds.csv', header = T, sep = ',')
Exp_IIIfree <- Exp_IIIfree %>% filter(condition == "Free") %>% aggregate(cbind(time)~subject, FUN="sum") 
Exp_IIIfree_mean = mean(Exp_IIIfree$time, na.rm = TRUE)
print(Exp_IIIfree_mean) # Exp II: "Free" was logged

# Experiment III (Birds)
# 2853.189 seconds

# ============ =

# Counts per condition
Trial_counts_I <- Exp_I %>%
  group_by(condition) %>%
  summarize(count = n())
print(Trial_counts_I)

# Pleasant     1586
# Unpleasant   1652
# Silence      1707

Trial_counts_II <- Exp_II %>%
  group_by(condition) %>%
  summarize(count = n())
print(Trial_counts_II)

# Pleasant     891
# Unpleasant   586
# Silence      230

Trial_counts_III <- Exp_III %>%
  group_by(condition) %>%
  summarize(count = n())
print(Trial_counts_III)

# Pleasant     587
# Unpleasant   455
# Silence      528

# --- --- 2.1 Total time per condition ####
avg_total_per_session_condition_I <- Exp_I %>%
  group_by(subject, session, condition) %>%
  summarise(subject_total_time = sum(time, na.rm = TRUE), .groups = "drop") %>%
  group_by(condition) %>%
  summarise(avg_total_time = mean(subject_total_time, na.rm = TRUE), .groups = "drop")
print(avg_total_per_session_condition_I)

# condition  avg_total_time
# Pleasant            2205.
# Unpleasant          1900.
# Silence             2391.

avg_total_per_session_condition_II <- Exp_II %>%
  group_by(subject, session, condition) %>%
  summarise(subject_total_time = sum(time, na.rm = TRUE), .groups = "drop") %>%
  group_by(condition) %>%
  summarise(avg_total_time = mean(subject_total_time, na.rm = TRUE), .groups = "drop")
print(avg_total_per_session_condition_II)

# condition  avg_total_time
# Pleasant            115. 
# Unpleasant           70.8
# Silence              22.6

avg_total_per_session_condition_III <- Exp_III %>%
  group_by(subject, session, condition) %>%
  summarise(subject_total_time = sum(time, na.rm = TRUE), .groups = "drop") %>%
  group_by(condition) %>%
  summarise(avg_total_time = mean(subject_total_time, na.rm = TRUE), .groups = "drop")
print(avg_total_per_session_condition_III)

# condition  avg_total_time
# Pleasant            2384.
# Unpleasant          2415.
# Silence             3182.

# --- --- 2.2 Mean & Median times per trial per condition #####

Trial_times_I <- Exp_I %>%
  group_by(condition) %>%
  aggregate(cbind(time)~condition, FUN="mean")
print(Trial_times_I)

# mean
#> print(Trial_times_I)
#condition     time
#1   Pleasant 75.07886
#2 Unpleasant 60.96297
#3    Silence 70.04121

# median
# condition   time
# Pleasant    7
# Unpleasant  4.725157
# Silence     8.638542

Trial_times_II <- Exp_II %>%
  group_by(condition) %>%
  aggregate(cbind(time)~condition, FUN="median")
print(Trial_times_II)

# mean
# condition     time
# Pleasant      7.582334
# Unpleasant    7.006174
# Silence       5.687775

# median
# condition     time
# Pleasant      1.403743
# Unpleasant    1.401852
# Silence       3.009669

Trial_times_III <- Exp_III %>%
  group_by(condition) %>%
  aggregate(cbind(time)~condition, FUN="median")
print(Trial_times_III)

# mean
# condition     time
# Pleasant      109.6467
# Unpleasant    123.9665
# Silence       156.7072

# median
# condition      time
# Pleasant       7.843076
# Unpleasant     10.550050
# Silence        12.063900

Exp_I$time <- as.numeric(Exp_I$time)
Exp_I$session <- factor(Exp_I$session)
Exp_II$time <- as.numeric(Exp_II$time)
Exp_II$session <- factor(Exp_II$session)
Exp_III$time <- as.numeric(Exp_III$time)
Exp_III$session <- factor(Exp_III$session)

Total_times_III <- Exp_III %>%
  group_by(condition) %>%
  aggregate(cbind(time)~condition, FUN="sum")
print(Total_times_III)

# --- --- 2.3 Plot raw data ####

library(tidyverse)

birds_raw_plot <- Exp_I %>%
  group_by(subject, condition) %>%
  summarise(
    mean_val = mean(time, na.rm = TRUE),
    sum_val  = sum(time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(species = "Bird")

humans_raw_plot <- Exp_II %>%
  group_by(subject, condition) %>%
  summarise(
    mean_val = mean(time, na.rm = TRUE),
    sum_val  = sum(time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(species = "Human")

III_raw_plot <- Exp_III %>%
  group_by(subject, condition, session) %>%
  summarise(
    mean_val = mean(time, na.rm = TRUE),
    sum_val  = sum(time, na.rm = TRUE),
    .groups = "drop"
  )

# Merge everything
#plot_raw_both <- bind_rows(birds_raw_plot, humans_raw_plot) %>%
plot_raw_III <- III_raw_plot %>%
  #mutate(sex = ifelse(sex %in% c("f", "female"), "Female", "Male")) %>%
  pivot_longer(cols = c(mean_val, sum_val), 
               names_to = "metric", 
               values_to = "value") %>%
  mutate(metric = case_when(
    metric == "mean_val" ~ "Mean Time per Trial",
    metric == "sum_val"  ~ "Total Time per Session"
  ))

# Normalised plot
#plot_norm <- plot_raw_both %>%
plot_normIII <- plot_raw_III %>%
  group_by(subject, metric, session) %>%
  mutate(
    centered_score = value - mean(value, na.rm = TRUE)
  ) %>%
  ungroup()

ggplot(plot_normIII, aes(x = condition, y = centered_score, fill = condition)) +
  geom_line(aes(group = subject), color = "gray", alpha = 0.4) +
  geom_jitter(width = 0.05, alpha = 0.4, size = 1.2) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) + 
  facet_wrap(metric ~ session, scales = "free") +
  labs(y = "Relative Time (Mean-Centered within Subject)",
       x = "Condition") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(strip.placement = "outside",
        strip.text.y = element_text(angle = 270),
        legend.position = "none")

# --- 3. Main Analyses ####
# --- --- 3.1 Main Analysis 1: Total Time (Both Species) ####

sum_bird <- Exp_I %>%
  group_by(subject, session, sex, condition) %>%
  dplyr::summarise(
    sum_time = sum(time, na.rm = TRUE),
    n_bouts = n(), # Good to keep track of how many bouts made up that mean
    .groups = "drop"
  )

sum_human <- Exp_II %>%
  group_by(subject, session, sex, condition) %>%
  dplyr::summarise(
    sum_time = sum(time, na.rm = TRUE),
    n_bouts = n(), .groups = "drop"
  )

# Add "species" columns
sum_bird$species <- "bird"
sum_human$species <- "human"

# Combine the two data frames
sum_both <- rbind(sum_bird, sum_human)

sum_both$sex[sum_both$sex %in% c("f")] <- "female"
sum_both$sex[sum_both$sex %in% c("m")] <- "male"

sum_both$species <- factor(sum_both$species)
sum_both$sex <- factor(sum_both$sex)
sum_both$session <- factor(sum_both$session)

performance::check_distribution(sum_both$sum_time)

# CHECK NEED FOR DISP
m_no_disp_totals = glmmTMB(log(sum_time) ~ condition * (session + sex + species) + (1|subject), 
                    data = sum_both, family = gaussian())
# Simulate residuals for your model WITHOUT dispersion modeling
res = simulateResiduals(m_no_disp_totals)

# Plot residuals against your main predictors (Species, Sex, etc.)
testCategorical(res, catPred = sum_both$species)
testCategorical(res, catPred = sum_both$sex)
testCategorical(res, catPred = sum_both$condition)
testCategorical(res, catPred = sum_both$session)

# Species could not be taken into account in dispformulas, because it 
# led to failure in KS test and Quantile Deviations Res vs Predicted

# MODELS

m_both_total_glm_max = glmmTMB(log(sum_time) ~ condition * (session + sex + species) + (1|subject), 
                            dispformula = ~condition+sex+session,
                            data = sum_both, 
                            family = gaussian())

m_both_total_glm_mid = glmmTMB(log(sum_time) ~ condition * (session + sex) + species + (1|subject), 
                           dispformula = ~condition+sex+session,
                           data = sum_both, 
                           family = gaussian())

m_both_total_glm_min1 = glmmTMB(log(sum_time) ~ condition * session + sex + species + (1|subject), 
                               dispformula = ~condition+sex+session,
                               data = sum_both, 
                               family = gaussian())

m_both_total_glm_min2 = glmmTMB(log(sum_time) ~ condition * sex + session + species + (1|subject), 
                                dispformula = ~condition+sex+session,
                                data = sum_both, 
                                family = gaussian())

m_both_total_glm_add = glmmTMB(log(sum_time) ~ condition + sex + session + species + (1|subject), 
                                dispformula = ~condition+sex+session,
                                data = sum_both, 
                                family = gaussian())

m_both_total_glm_cov = glmmTMB(log(sum_time) ~ sex + session + species + (1|subject), 
                               dispformula = ~condition+sex+session,
                               data = sum_both, 
                               family = gaussian())

m_both_total_glm_null = glmmTMB(log(sum_time) ~ (1|subject), 
                               dispformula = ~condition+sex+session,
                               data = sum_both, 
                               family = gaussian())

m_both_total_glm_ultra = glmmTMB(log(sum_time) ~ (condition * (sex + session))* species + (1|subject), 
                                dispformula = ~condition+sex+session,
                                data = sum_both, 
                                family = gaussian())

anova(m_both_total_glm_null, m_both_total_glm_cov, m_both_total_glm_add, 
      m_both_total_glm_min1, m_both_total_glm_min2, m_both_total_glm_max,
      m_both_total_glm_mid, m_both_total_glm_ultra)

# Double-check winning model
anova(m_both_total_glm_max, m_both_total_glm_cov, test = 'Chisq') # Winner vs Null
drop1(m_both_total_glm_max, test = 'Chisq') # Test interaction terms significant
# !!! Here, the interaction of condition:sex is n.s., BUT removal fails KS test. => KEPT.

# Check diagnostics
simulation_output <- simulateResiduals(m_both_total_glm_max); plot(simulation_output)
testResiduals(simulation_output)

summary(m_both_total_glm_max)
Anova(m_both_total_glm_max, type = 3); #summary(m1)

emmeans(m_both_total_glm_max, pairwise ~ condition|session, type = "response", adjust = "Bonferroni")
emmeans(m_both_total_glm_max, pairwise ~ condition|species, type = "response", adjust = "Bonferroni")

# --- --- 3.2 Main Analysis 2: Trials (Both Species) ####


mean_bird <- Exp_I %>%
  group_by(subject, session, sex, condition) %>%
  dplyr::summarise(
    mean_time = mean(time, na.rm = TRUE),
    n_bouts = n(), # Good to keep track of how many bouts made up that mean
    .groups = "drop"
  )

mean_human <- Exp_II %>%
  group_by(subject, session, sex, condition) %>%
  dplyr::summarise(
    mean_time = mean(time, na.rm = TRUE),
    n_bouts = n(), # Good to keep track of how many bouts made up that mean
    .groups = "drop"
  )

# Add "species" columns
mean_bird$species <- "bird"
mean_human$species <- "human"

# Combine the two data frames
mean_both <- rbind(mean_bird, mean_human)

mean_both$sex[mean_both$sex %in% c("f")] <- "female"
mean_both$sex[mean_both$sex %in% c("m")] <- "male"

mean_both$species <- factor(mean_both$species)
mean_both$sex <- factor(mean_both$sex)
mean_both$session <- factor(mean_both$session)
summary(mean_both)

performance::check_distribution(mean_both$mean_time)

# 1. Fitting a version of your best model WITHOUT the dispformula
m_no_disp = glmmTMB(log(mean_time) ~ condition * (session + sex) + species + (1|subject), 
                    data = mean_both, family = gaussian())
# Simulate residuals for your model WITHOUT dispersion modeling
res = simulateResiduals(m_no_disp)

# Plot residuals against your main predictors (Species, Sex, etc.)
testCategorical(res, catPred = mean_both$species)
testCategorical(res, catPred = mean_both$sex)
testCategorical(res, catPred = mean_both$condition)
testCategorical(res, catPred = mean_both$session)

# Model selection
m_both_mean_null = glmmTMB(log(mean_time) ~ session+sex+species+
                            (1|subject), 
                          dispformula = ~condition+species+sex,
                          data = mean_both, family = gaussian())
m_both_mean_min = glmmTMB(log(mean_time) ~ condition+session+sex+species+
                        (1|subject), 
                      dispformula = ~condition+species+sex, 
                      data = mean_both, family = gaussian())
m_both_mean = glmmTMB(log(mean_time) ~ condition*(session+sex)+species+
                         (1|subject), 
                       dispformula = ~condition+species+sex,
                       data = mean_both, family = gaussian())
m_both_mean_mid = glmmTMB(log(mean_time) ~ condition*(session+sex+species)+
                            (1|subject), 
                          dispformula = ~condition+species+sex,
                          data = mean_both, family = gaussian())
m_both_mean_max = glmmTMB(log(mean_time) ~ condition*session*sex*species+
                            (1|subject), 
                          dispformula = ~condition+species+sex,
                          data = mean_both, family = gaussian())
m_both_mean_ultra = glmmTMB(log(mean_time) ~ (condition*(session+sex))*species+
                            (1|subject), 
                          dispformula = ~condition+species+sex,
                          data = mean_both, family = gaussian())
anova(m_both_mean_null, m_both_mean_min, m_both_mean, 
      m_both_mean_mid, m_both_mean_max,
      m_both_mean_ultra)

# Double-check winning model
anova(m_both_mean_mid, m_both_mean_null, test = 'Chisq') # Winner vs Null
drop1(m_both_mean_mid, test = 'Chisq') # Test interaction terms significant

simulation_output <- simulateResiduals(m_both_mean_mid); plot(simulation_output)


summary(m_both_mean_mid)
Anova(m_both_mean_mid, type = 3); #summary(m1)


#emmeans(m_both_mean, pairwise ~ condition, type = "response")
#emmeans(m_both_mean, pairwise ~ sex, type = "response")
emmeans(m_both_mean_mid, pairwise ~ condition|session, type = "response", adjust = "Bonferroni")
emmeans(m_both_mean_mid, pairwise ~ condition|sex, type = "response", adjust = "Bonferroni")
emmeans(m_both_mean_mid, pairwise ~ condition|species, type = "response", adjust = "Bonferroni")

# --- --- 3.3 Main Analysis 3: Trial counts ####

library(emmeans)
library(car)
library(glmmTMB) #(Brooks et al., 2017; McGillycuddy et al., 2025)

Trial_counts_Ib <- Exp_I %>%
  group_by(condition, subject, session, sex) %>%
  summarize(number = n())

Trial_counts_IIb <- Exp_II %>%
  group_by(condition, subject, session, sex) %>%
  summarize(number = n())

# Add "species" columns
Trial_counts_Ib$species <- "bird"
Trial_counts_IIb$species <- "human"

# Combine the two data frames
trial_numbers <- rbind(Trial_counts_Ib, Trial_counts_IIb)

trial_numbers$sex[trial_numbers$sex %in% c("f")] <- "female"
trial_numbers$sex[trial_numbers$sex %in% c("m")] <- "male"

trial_numbers$species <- factor(trial_numbers$species)
trial_numbers$sex <- factor(trial_numbers$sex)
trial_numbers$session <- factor(trial_numbers$session)

performance::check_distribution(trial_numbers$number)

# 1. Fitting a version of your best model WITHOUT the dispformula
m_no_disp_numbers = glmmTMB(number ~ condition * (session + sex + species) + (1|subject), 
                    data = trial_numbers, family = Gamma(link='log'))
# Simulate residuals for your model WITHOUT dispersion modeling
res = simulateResiduals(m_no_disp_numbers)

# Plot residuals against your main predictors (Species, Sex, etc.)
testCategorical(res, catPred = trial_numbers$species) # SIG
testCategorical(res, catPred = trial_numbers$sex) # SIG
testCategorical(res, catPred = trial_numbers$condition) # SIG
testCategorical(res, catPred = trial_numbers$session) # NS

m1Counts = glmmTMB(number ~ condition*(sex+session+species)+(1|subject),
                 data = trial_numbers, family = Gamma(link='log'),
                 dispformula = ~species+sex+condition)
m1Counts_max = glmmTMB(number ~ (condition*(sex+session))*species+(1|subject),
                   data = trial_numbers, family = Gamma(link='log'),
                   dispformula = ~species+sex+condition)
m1Counts_min = glmmTMB(number ~ condition*species+(1|subject),
                       data = trial_numbers, family = Gamma(link='log'),
                       dispformula = ~species+sex+condition)
m1Counts_null = glmmTMB(number ~ (sex+session)*species+(1|subject),
                       data = trial_numbers, family = Gamma(link='log'),
                       dispformula = ~species+sex+condition)
m1Counts_add = glmmTMB(number ~ condition+sex+session+species+(1|subject),
                        data = trial_numbers, family = Gamma(link='log'),
                        dispformula = ~species+sex+condition)
anova(m1Counts, m1Counts_null, m1Counts_min, m1Counts_max, m1Counts_add)
Anova(m1Counts, type = 3)

emt = emmeans(m1Counts, pairwise ~ condition|species, type = "response", adjust = "Bonferroni"); summary(emt)
emt = emmeans(m1Counts, pairwise ~ condition|sex, type = "response", adjust = "Bonferroni"); summary(emt)

# --- --- Bonus: Correlations between count and time ####
# Correlations between trial count per subject and mean trial time
###
# Experiment 1
###
Corr_I_counts <- Exp_I %>%
  group_by(subject,condition) %>%
  summarize(count = n())
Corr_I_counts <- Corr_I_counts %>% 
  arrange(subject, condition)

Corr_I_times <-with(Exp_I, aggregate(cbind(time)~subject+condition, FUN="mean")) # median: bad residuals
Corr_I_times <- Corr_I_times %>% 
  arrange(subject, condition)

Corr_I_combined <- Corr_I_times %>%
  left_join(Corr_I_counts, by = c("subject", "condition"))

Corr_I_combined$subject <- factor(Corr_I_combined$subject)

#install.packages("rmcorr")
library(rmcorr)
rmc_result <- rmcorr(participant = subject, 
                     measure1 = time, 
                     measure2 = count, 
                     dataset = Corr_I_combined)
plot(rmc_result)
rmc_result

###
# Experiment 2
###

# Moderate negative correlation between count and mean time (rho = -0.6237294, p < 0.0001)

Corr_II_counts <- Exp_II %>%
  group_by(subject,condition) %>%
  summarize(count = n())
Corr_II_counts <- Corr_II_counts %>% 
  arrange(subject, condition)

Corr_II_times <-with(Exp_II, aggregate(cbind(time)~subject+condition, FUN="mean")) # median: bad residuals
Corr_II_times <- Corr_II_times %>% 
  arrange(subject, condition)

performance::check_distribution(Corr_II_counts$count)
performance::check_distribution(Corr_II_times$time)

Corr_II_combined <- Corr_II_times %>%
  left_join(Corr_II_counts, by = c("subject", "condition"))

Corr_II_combined$subject <- factor(Corr_II_combined$subject)

#install.packages("rmcorr")
library(rmcorr)
rmc_result2 <- rmcorr(participant = subject, 
                     measure1 = time, 
                     measure2 = count, 
                     dataset = Corr_II_combined)
plot(rmc_result2)
rmc_result2

# --- 4. Listener Rating Survey ####
# --- --- 4.1 Data Loading and Preparation ####

rating_data = read.csv('.../roughness_ratings_raw.csv', header = T, sep = ',')
names(data)

# 1. Clean up those long names first
# Focus on columns 29 to 64 (the ratings)
clean_data <- rating_data %>%
  rename_with(~ paste0("S", str_extract(., "\\d+"), "_Pleasantness"), contains("pleasant")) %>%
  rename_with(~ paste0("S", str_extract(., "\\d+"), "_Liking"), contains("much.do.you.like")) %>%
  rename_with(~ paste0("S", str_extract(., "\\d+"), "_Interest"), contains("interesting"))

long_data <- clean_data %>%
  pivot_longer(
    cols = starts_with("S"), 
    names_to = "temp_id", 
    values_to = "rating"
  ) %>%
  # Separate Sound ID and Rating Type
  separate(temp_id, into = c("sound", "rating_type"), sep = "_")

# Ensure rating is an ordered factor for CLMM
long_data <- long_data %>%
  mutate(rating = factor(rating, ordered = TRUE))

long_data$condition <- "Unpleasant"  # Set a default
long_data$condition[long_data$sound %in% c("S1", "S2", "S3", "S4", "S5", "S6")] <- "Pleasant"

long_data <- long_data %>% rename(subject = What.is.your.anonymized.abbreviation...e.g...f12.or.m8.)
long_data <- long_data %>%   rename(age = What.is.your.age.)
long_data <- long_data %>%   rename(musicianship = Which.of.the.following.descriptions.fit.you.as.a.musician.)
long_data <- long_data %>%   rename(sex = What.is.your.sex.)

long_data <- long_data %>%
  pivot_wider(
    names_from = rating_type, 
    values_from = rating
  )
names(long_data)
long_data = subset(long_data, select = c(age, sex, subject, 
                                         musicianship,
                                         condition, sound,
                                         Pleasantness, Liking, Interest))

long_data <- long_data %>%
  mutate(musicianship = str_remove(musicianship, " \\(.*"))

# --- --- --- Factorize

long_data$condition <- as.factor(long_data$condition)
long_data$sex <- factor(long_data$sex)
long_data$subject <- factor(long_data$subject)
long_data$sound <- factor(long_data$sound)
long_data$condition <- fct_relevel(long_data$condition, "Unpleasant","Pleasant")
long_data= droplevels(long_data[long_data$condition != '',])

library(ordinal)

long_data$Pleasantness <- factor(long_data$Pleasantness, ordered = TRUE)
long_data$Liking <- factor(long_data$Liking, ordered = TRUE)
long_data$Interest <- factor(long_data$Interest, ordered = TRUE)
# --- --- 4.2 Main Ratings Analysis ####
# --- --- -- 4.2.1 Pleasantness ####
m0 = clmm(Pleasantness ~ condition + 
            (1 + condition | subject) + 
            (1| sound),
          data = long_data)
m0s = clmm(Pleasantness ~ condition*sex + 
             (1 + condition | subject) + 
             (1| sound),
           data = long_data)
m0sAdd = clmm(Pleasantness ~ condition+sex + 
                (1 + condition | subject) + 
                (1 | sound),
              data = long_data)
m0_null = clmm(Pleasantness ~ sex + 
                 (1 + condition | subject) + 
                 (1 | sound),
               data = long_data)
anova(m0, m0s, m0sAdd, m0_null) #==> Sex (interact) did NOT improve model fit

# ==> m0sAdd wins.

subj_ranef <- ranef(m0sAdd)$subject
# Q-Q plot of subjects
qqnorm(subj_ranef$`(Intercept)`, main = "Q-Q Plot: Subject Random Effects")
qqline(subj_ranef$`(Intercept)`, col = "red")
# ==> Q-Q plot is neat. COuld be better, but does not look catastrophic.

print(car::Anova(m0sAdd, type = 3))
summary(m0sAdd)
emmeans(m0sAdd, ~ condition, adjust = 'Bonferroni', type = 'response')
print(pairs(posthoc))

# --- --- -- 4.2.2 Liking ####
m1 = clmm(Liking ~ condition + (1 + condition | subject) + (1 | sound),
          data = long_data)
m1s <- clmm(Liking ~ condition * sex + (1 + condition | subject) + (1 | sound), 
            data = long_data)
m1sAdd <- clmm(Liking ~ condition + sex + (1 + condition | subject) + (1 | sound), 
               data = long_data)
m1s_0 <- clmm(Liking ~ sex + (1 + condition | subject) + (1 | sound), 
              data = long_data)
anova(m1, m1s, m1sAdd, m1s_0)

# Note: Including Sex as a random slope in subjects led to Hessian not pos. definit.

# ==> m1sAdd wins

subj_ranef <- ranef(m1sAdd)$subject
# Q-Q plot of subjects
qqnorm(subj_ranef$`(Intercept)`, main = "Q-Q Plot: Subject Random Effects")
qqline(subj_ranef$`(Intercept)`, col = "red")
# ==> Q-Q plot looks really good

print(car::Anova(m1sAdd, type = 3))
# addiditive structure above
summary(m1s)
posthoc <- emmeans(m1sAdd, ~ condition, adjust = 'Bonferroni', type = "response")
print(pairs(posthoc))


# --- --- -- 4.2.3 Interest ####
m2 = clmm(Interest ~ condition +
            (1 + condition | subject) + 
            (1 | sound),
          data = long_data)
m2s <- clmm(Interest ~ condition * sex + 
              (1 + condition | subject) + (1 | sound), 
            data = long_data)
m2sAdd <- clmm(Interest ~ condition + sex + 
                 (1 + condition | subject) + (1 | sound), 
               data = long_data)
m2s_0 <- clmm(Interest ~ sex + 
                (1 + condition | subject) + (1 | sound), 
              data = long_data)
anova(m2, m2s, m2sAdd, m2s_0) # More complex one with sex didn't improve it

# ==> m2sAdd wins

subj_ranef <- ranef(m2sAdd)$subject
# Q-Q plot of subjects
qqnorm(subj_ranef$`(Intercept)`, main = "Q-Q Plot: Subject Random Effects")
qqline(subj_ranef$`(Intercept)`, col = "red")
# ==> Q-Q plot looks ok

print(car::Anova(m2sAdd, type = 3))
summary(m2sAdd)
posthoc <- emmeans(m2sAdd, ~ condition, adjust = 'Bonferroni', type = "response")
print(pairs(posthoc))

# --- 5. Plotting the raw (i.e., modeled) data ####

# Totals
ag1a <-with(Exp_I, aggregate(cbind(time)~subject+condition+sex, FUN="sum"))
Exp_1ag <- ag1a
Exp_2ag <- ag2
Exp_1ag <- Exp_1ag %>%
  dplyr::select(time, condition, subject, sex)
Exp_2ag <- Exp_2ag %>%
  dplyr::select(time, condition, subject, sex)
Exp_1ag$species <- "Budgerigars"
Exp_2ag$species <- "Humans"
Both_Total <- rbind(Exp_1ag, Exp_2ag)

Both_Total %>%
  group_by(species, condition) %>%
  summarise(
    mean_time = mean(time, na.rm = TRUE),
    n = n()  # This shows how many trials/observations are in each group
  )

# Trials
Exp_1 <- Exp_I
Exp_2 <- Exp_II
Exp_1 <- Exp_1 %>%
  dplyr::select(time, condition, subject)
Exp_2 <- Exp_2 %>%
  dplyr::select(time, condition, subject)
Exp_1$species <- "Budgerigars"
Exp_2$species <- "Humans"
Both_Trial <- rbind(Exp_1, Exp_2)

####
####

# First plot for total time
anno_df1 <- data.frame(
  species = "Humans", 
  start = c("Pleasant", "Pleasant", "Unpleasant"),
  end = c("Unpleasant", "Silence", "Silence"),
  y = c(130, 145, 110),
  label = c("***", "***", "***")
)
pd <- position_dodge(0.2)
plot_total <- ggplot(Both_Total, aes(condition, time, group=condition, color=condition, fill=condition))+
  theme_bw()+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  stat_summary(fun.y=mean, geom="point", size = 5, alpha=0.8)+
  theme_classic()+
  facet_wrap( ~ species, scales = "free_y")+
  scale_color_manual(values = c("Pleasant" = "#1E2F84", "Unpleasant" = "#E8B4BD", "Silence" = "#DEDAD9"))+
  scale_fill_manual(values = c("Pleasant" = "#1E2F84", "Unpleasant" = "#E8B4BD", "Silence" = "#DEDAD9"))+
  stat_summary(fun.data = 'mean_se', geom = "errorbar", width=0, alpha=0.6)+
  xlab("Condition") + ylab("Mean Total Time (s)")+
  geom_signif(
    data = anno_df1,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    manual = TRUE,
    inherit.aes = FALSE,# Crucial: stops it from looking for 'condition' in anno_df
    tip_length = 0.01
  )+
  theme(legend.position = 'none',
        text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, color = "#327264"),
        axis.title.y = element_text(size = 16, color = "#327264"),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

# Second plot for trial time
anno_df2 <- data.frame(
  species = "Humans", 
  start = c("Pleasant", "Unpleasant"),
  end = c("Silence", "Silence"),
  y = c(8.5, 8.0),
  label = c("***", "**")
)
plot_trial <- ggplot(Both_Trial, aes(condition, time, group=condition, color=condition, fill=condition))+
  theme_bw()+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  stat_summary(fun.y=mean, geom="point", size = 5, alpha=0.8)+
  theme_classic()+
  facet_wrap( ~ species, scales = "free_y")+
  scale_color_manual(values = c("Pleasant" = "#1E2F84", "Unpleasant" = "#E8B4BD", "Silence" = "#DEDAD9"))+
  scale_fill_manual(values = c("Pleasant" = "#1E2F84", "Unpleasant" = "#E8B4BD", "Silence" = "#DEDAD9"))+
  stat_summary(fun.data = 'mean_se', geom = "errorbar", width=0, alpha=0.6)+
  xlab("Condition") + ylab("Mean Time per Trial (s)")+
  geom_signif(
    data = anno_df2,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    manual = TRUE,
    inherit.aes = FALSE,# Crucial: stops it from looking for 'condition' in anno_df
    tip_length = 0.0005
  )+
  theme(legend.position = 'none',
        text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, color = "#327264"),
        axis.title.y = element_text(size = 16, color = "#327264"),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

# Combine Plots
library(cowplot)
combined_plot <- cowplot::plot_grid(plot_total, plot_trial,
                           ncol = 1,
                           labels = c("A", "B"),
                           rel_heights = c(1, 1))
combined_plot
ggsave(paste0(directory, 'RanRegCombined_Results_cowplot.pdf'),
       plot = combined_plot,
       width = 8, height = 8)

# --- Bonus: Plotting Counts vs Trials ####

library(tidyverse)

data_TimeAndCount <- tribble(
  ~species, ~condition, ~count, ~median_time,
  "Birds",  "Pleasant", 848,    7,
  "Birds",  "Unpleasant", 970,    3,
  "Birds",  "Silent",  901,    12,
  "Humans", "Pleasant", 891,    1.403743,
  "Humans", "Unpleasant", 586,    1.401852,
  "Humans", "Silent",  230,    3.009669
)

# For the order Pleasant, Unpleasant, Silent
data_TimeAndCount <- data_TimeAndCount %>%
  mutate(condition = factor(
    condition,
    levels = c("Pleasant", "Unpleasant", "Silent")
  ))

data_TimeAndCount <- data_TimeAndCount %>%
  mutate(species = factor(
    species,
    levels = c("Humans", "Birds")
  ))

# Create the scatter plot
TimeAndCount_plot <- ggplot(data_TimeAndCount,
                        aes(y = count,
                            x = median_time,
                            shape = species,
                            color = condition)) +
  # Add the points
  geom_point(size = 5) +
  # Manually set colors for each condition

  scale_color_manual(values = c(
    "Pleasant" = "#1E2F84",
    "Unpleasant" = "#E8B4BD",
    "Silent" = "#DEDAD9"
  )) +
  # Manually set shapes for each species
  scale_shape_manual(values = c(
    "Birds" = 8,  # 17 is a filled triangle
    "Humans" = 15 # 19 is a filled circle
  )) +
  # Add labels and a title
  labs(
    y = "Total Trial Count",
    x = "Median Time per Trial (s)",
    color = "Condition",
    shape = "Species"
  ) +
  theme_minimal() +
  guides(
    shape = guide_legend(order = 1), # Species legend first
    color = guide_legend(order = 2)  # Condition legend second
  ) +
  theme(
    legend.position = "right", # Position the legend at the bottom
    axis.title.y = element_text(color = "#327264"),
    axis.title.x = element_text(color = "#327264")
  )

# Print the plot
print(TimeAndCount_plot)

# --- 6. Experiment III ####

# --- --- 6.1 Total Times ####
options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly"))


sum_III <- Exp_III %>%
  group_by(subject, session, sex, condition) %>%
  dplyr::summarise(
    sum_time = sum(time, na.rm = TRUE),
    n_bouts = n(), # Good to keep track of how many bouts made up that mean
    .groups = "drop"
  )

sum_III$sex <- factor(sum_III$sex)
sum_III$session <- factor(sum_III$session)

performance::check_distribution(sum_III$sum_time) # Lognormal

# CHECK NEED FOR DISP
m_no_disp_totals_III = glmmTMB(log(sum_time) ~ condition * (session + sex) + (1|subject), 
                           data = sum_III, family = gaussian())
# Simulate residuals for your model WITHOUT dispersion modeling
res_III = simulateResiduals(m_no_disp_totals_III)

# Plot residuals against your main predictors (Species, Sex, etc.)
testCategorical(res_III, catPred = sum_III$sex) # good
testCategorical(res_III, catPred = sum_III$condition) #good
testCategorical(res_III, catPred = sum_III$session) # good

# No need for dispformula

# MODELS

m_III_total_glm_max = glmmTMB(log(sum_time) ~ condition * (session + sex) + (1|subject), 
                               data = sum_III, 
                               family = gaussian())

m_III_total_glm_mid = glmmTMB(log(sum_time) ~ condition * (session + sex) + (1|subject), 
                               data = sum_III, 
                               family = gaussian())

m_III_total_glm_min1 = glmmTMB(log(sum_time) ~ condition * session + sex + (1|subject), 
                                data = sum_III, 
                                family = gaussian())

m_III_total_glm_min2 = glmmTMB(log(sum_time) ~ condition * sex + session + (1|subject), 
                                data = sum_III, 
                                family = gaussian())

m_III_total_glm_add = glmmTMB(log(sum_time) ~ condition + sex + session + (1|subject), 
                               data = sum_III, 
                               family = gaussian())

m_III_total_glm_null = glmmTMB(log(sum_time) ~ (1|subject), 
                                data = sum_III, 
                                family = gaussian())

m_III_total_glm_cond = glmmTMB(log(sum_time) ~ condition + (1|subject), 
                               data = sum_III, 
                               family = gaussian())

anova(m_III_total_glm_null, m_III_total_glm_add, 
      m_III_total_glm_min1, m_III_total_glm_min2, m_III_total_glm_max,
      m_III_total_glm_mid, m_III_total_glm_cond)
# Well... The null model wins. So let's check the next simpler model (with 2nd lowest AIC/BIC)

# Check diagnostics
simulation_output <- simulateResiduals(m_III_total_glm_cond); plot(simulation_output)
testResiduals(simulation_output)

summary(m_III_total_glm_cond)
Anova(m_III_total_glm_cond, type = 3); #summary(m1)


# Model: sum_time ~ condition + (1|subject)
#Response: sum_time
#                 Chisq Df Pr(>Chisq)    
# (Intercept) 892.6503  1     <2e-16 ***
# condition     3.0987  2     0.2124  

# --- --- 6.2 Trial Times ####

mean_III <- Exp_III %>%
  group_by(subject, session, sex, condition) %>%
  dplyr::summarise(
    mean_time = mean(time, na.rm = TRUE),
    n_bouts = n(), # Good to keep track of how many bouts made up that mean
    .groups = "drop"
  )

# Combine the two data frames
mean_both <- rbind(mean_III)

mean_III$sex[mean_III$sex %in% c("f")] <- "female"
mean_III$sex[mean_III$sex %in% c("m")] <- "male"

mean_III$sex <- factor(mean_III$sex)
mean_III$session <- factor(mean_III$session)
summary(mean_III)

performance::check_distribution(mean_III$mean_time)

# 1. Fitting a version of your best model WITHOUT the dispformula
mIII_no_disp = glmmTMB(log(mean_time) ~ condition * (session + sex) + (1|subject), 
                    data = mean_III, family = gaussian())
# Simulate residuals for your model WITHOUT dispersion modeling
res = simulateResiduals(mIII_no_disp)

# Plot residuals against your main predictors (Species, Sex, etc.)
testCategorical(res, catPred = mean_III$sex) # Significant
testCategorical(res, catPred = mean_III$condition) # NS
testCategorical(res, catPred = mean_III$session) # NS

# Model selection
m_III_mean_null = glmmTMB(log(mean_time) ~ session+sex+
                             (1|subject), 
                           data = mean_III, family = gaussian())
m_III_mean_min = glmmTMB(log(mean_time) ~ condition+session+sex+
                            (1|subject), 
                          data = mean_III, family = gaussian())
m_III_mean_mid = glmmTMB(log(mean_time) ~ condition*(session+sex)+
                            (1|subject), 
                          data = mean_III, family = gaussian())
m_III_mean_max = glmmTMB(log(mean_time) ~ condition*session*sex+
                            (1|subject), 
                          data = mean_III, family = gaussian())
m_III_mean_sess = glmmTMB(log(mean_time) ~ condition*session+sex+
                           (1|subject), 
                         data = mean_III, family = gaussian())
m_III_mean_cond = glmmTMB(log(mean_time) ~ condition +
                            (1|subject),
                          data = mean_III, family = gaussian())
unique(mean_III$condition)
anova(m_III_mean_null, m_III_mean_min, m_III_mean_cond,
      m_III_mean_mid, m_III_mean_max, m_III_mean_sess)
# Ok, I mean, max wins, but has bad fit (residuals vs predicted).
# Also, there's just ridiculous significant terms: A three-way interaction of
# Sex, session and condition. With that few birds, that's too spurious.
# Let's then just keep the one with condition.

simulation_output <- simulateResiduals(m_III_mean_cond); plot(simulation_output)
# The model with only cond is healthy 

summary(m_III_mean_cond)
Anova(m_III_mean_cond, type = 3); #summary(m1)

# --- --- 6.3 Trial Count ####

Trial_count_III <- Exp_III %>%
  group_by(condition, subject, session, sex) %>%
  summarize(number = n())

Trial_count_III$sex[Trial_count_III$sex %in% c("f")] <- "female"
Trial_count_III$sex[Trial_count_III$sex %in% c("m")] <- "male"

Trial_count_III$sex <- factor(Trial_count_III$sex)
Trial_count_III$session <- factor(Trial_count_III$session)

performance::check_distribution(Trial_count_III$number)

# 1. Fitting a version of your best model WITHOUT the dispformula
m_no_disp_numbersIII = glmmTMB(number ~ condition * (session + sex) + (1|subject), 
                            data = Trial_count_III, family = Gamma(link='log'))
# Simulate residuals for your model WITHOUT dispersion modeling
res = simulateResiduals(m_no_disp_numbersIII)

# Plot residuals against your main predictors (Species, Sex, etc.)
testCategorical(res, catPred = Trial_count_III$sex) # SIG
testCategorical(res, catPred = Trial_count_III$condition) # NS
testCategorical(res, catPred = Trial_count_III$session) # NS

mIIICounts_max = glmmTMB(number ~ condition*(sex+session)+(1|subject),
                       data = Trial_count_III, family = Gamma(link='log'),
                       dispformula = ~sex)
mIIICounts_min = glmmTMB(number ~ condition+(1|subject),
                       data = Trial_count_III, family = Gamma(link='log'),
                       dispformula = ~sex)
mIIICounts_null = glmmTMB(number ~ (sex+session)+(1|subject),
                        data = Trial_count_III, family = Gamma(link='log'),
                        dispformula = ~sex)
mIIICounts_add = glmmTMB(number ~ condition+sex+session+(1|subject),
                       data = Trial_count_III, family = Gamma(link='log'),
                       dispformula = ~sex)
anova(mIIICounts_null, mIIICounts_min, mIIICounts_max, mIIICounts_add)
Anova(mIIICounts_min, type = 3)

# --- Miscellaneous ####

#   
#         _...._
#        /      \
#       /  X __ X .
#       (    \/  )
#      )          (
#     (  -  -  -  )
#    (  /          )
#    (/            )
#     [          ]
#-------/l\    /l\-----------.
#--------.------.------------'
#       \ ___/
#        \ /
#
#



# --- 6. Experiment III <-  Follow Up on actual Roughness with birds ####
# Testing roughness as known by humans on the birds because the original study
# we replicated/translated actually was dubious regarding whether roughness 
# was actually manipulated in a controlled manner.

# First, check out the data
# Check missing counts specifically for model variables
colSums(is.na(Exp_III[, c("time", "condition", "sex", "subject")]))
head(Exp_III)
hist(Exp_III$time, breaks = 50)

library(tidyverse)

# Aggregate to mean time per subject & condition, then plot
Exp_III %>%
  group_by(subject, condition) %>%
  summarise(mean_trial_time = mean(time, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = condition, y = mean_trial_time, fill = condition)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  labs(
    title = "Mean Time per Trial across Conditions",
    x = "Condition",
    y = "Mean Trial Duration (seconds)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Aggregate to total time per subject & condition
Exp_III %>%
  group_by(subject, condition) %>%
  summarise(total_time = sum(time, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = condition, y = total_time, fill = condition)) +
  stat_summary(fun = "mean", geom = "col", alpha = 0.6, width = 0.5) +
  geom_point(position = position_jitter(width = 0.1), size = 2.5) +
  labs(
    title = "Total Cumulative Time across Conditions",
    x = "Condition",
    y = "Total Duration (seconds)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# --- --- Mean Time per trial ####

Exp_III$sex[Exp_III$sex %in% c("f")] <- "female"
Exp_III$sex[Exp_III$sex %in% c("m")] <- "male"
Exp_III$sex <- factor(Exp_III$sex)
Exp_III$session <- factor(Exp_III$session)

performance::check_distribution(Exp_III$time) # 97% lognormal

# CHECK NEED FOR DISP
m_Exp_III_trial_nodisp = glmmTMB(time ~ condition * sex + (1|subject), 
                           data = Exp_III, family = gaussian())
# Simulate residuals for your model WITHOUT dispersion modeling
res3 = simulateResiduals(m_Exp_III_trial_nodisp)

# Plot residuals against your main predictors (Species, Sex, etc.)
testCategorical(res3, catPred = Exp_III$sex) # Levene's test sig, Need disp
testCategorical(res3, catPred = Exp_III$condition) # Levene's test NON sig, no disp necessary
# But interaction could be tested.
head(Exp_III,10)
m_Exp_III_trial = glmmTMB(time ~ condition * sex + (1+condition||subject), 
                          dispformula = ~sex,
                                 data = Exp_III, family = gaussian())
sim_ExpIII_trial <- simulateResiduals(m_Exp_III_trial); plot(sim_ExpIII_trial)

Anova(m_Exp_III_trial, type = "III")
emmeans(m_Exp_III_trial, pairwise ~ condition)

# --- --- Total time ####

performance::check_distribution(ag3$time) # 50% lognormal; 16% Gamma


m_Exp_III_total = glmmTMB(time ~ condition + sex + (1+condition||subject), 
                          data = ag3, family = Gamma(link = 'log'))
sim_ExpIII_total <- simulateResiduals(m_Exp_III_total); plot(sim_ExpIII_total)

Anova(m_Exp_III_total, type = "II")
emmeans(m_Exp_III_total, pairwise ~ condition)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ####
# End ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ####
# Leftovers / old ####
# --- --- 3.1 - BIRDS, total time ####

library(DHARMa)
library(glmmTMB)
sum_bird <- Exp_I %>%
  group_by(subject, session, sex, condition) %>%
  dplyr::summarise(
    sum_time = sum(time, na.rm = TRUE),
    n_bouts = n(), # Good to keep track of how many bouts made up that mean
    .groups = "drop"
  )

b_total_gl_1 = glmmTMB(sum_time ~ condition+sex+session+
                         (1|subject),
                       #weights = n_bouts,
                       dispformula = ~ condition,
                       data = sum_bird, family = Gamma(link = 'log'))
simulation_output <- simulateResiduals(b_total_gl_1); plot(simulation_output)

b_total_gl_2 = glmmTMB(sum_time ~ session*condition+sex+
                         (1|subject),
                       #weights = n_bouts,
                       dispformula = ~ condition,
                       data = sum_bird, family = Gamma(link = 'log'))
simulation_output <- simulateResiduals(b_total_gl_2); plot(simulation_output)

b_total_gl_3 = glmmTMB(sum_time ~ (session+sex)*condition+
                         (1|subject), # had to remove (condition|subject) due to convergence issues
                       #weights = n_bouts,
                       dispformula = ~ condition,
                       data = sum_bird, family = Gamma(link = 'log'))
simulation_output <- simulateResiduals(b_total_gl_3); plot(simulation_output)

anova(b_total_gl_1, b_total_gl_2, b_total_gl_3)
summary(b_total_gl_2)
Anova(b_total_gl_2, type = 3); #summary(m1)

emt_b_tot_sess = emmeans(b_total_gl_2, pairwise ~ condition|session, type = "response"); summary(emt_b_tot_sess)

plot_emt_b_tot_sess <- as.data.frame(emt_b_tot_sess$emmeans)
pd <- position_dodge(-0.2)
ggplot(plot_emt_b_tot_sess, aes(x = session, y = response, color = condition, group = condition)) +
  geom_line(size = 1, position = pd) +
  geom_point(size = 4, position = pd) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1, alpha = 0.3, position = pd) +
  labs(x = "Session Number", y = "Predicted Total Time (seconds)", color = "Condition") +
  theme_minimal() + scale_color_brewer(palette = "Set1")


# --- --- 3.2- HUMANS, total time ####

sum_human <- Exp_II %>%
  group_by(subject, session, sex, condition) %>%
  dplyr::summarise(
    sum_time = sum(time, na.rm = TRUE),
    n_bouts = n(), .groups = "drop"
  )
performance::check_distribution(sum_human$sum_time)

h_total_gl = glmmTMB(sum_time ~ condition*(sex+session)+
                       (1|subject), 
                     dispformula = ~condition,
                     #weights = n_bouts, #On kills disptest and outlier test
                     data = sum_human, family = gaussian())
simulation_output <- simulateResiduals(h_total_gl); plot(simulation_output)
summary(h_total_gl)
Anova(h_total_gl, type = 3); #summary(m1)

emt = emmeans(h_total_gl, pairwise ~ condition, type = "response"); summary(emt)
emt_h_tot_sess = emmeans(h_total_gl, pairwise ~ condition|session, type = "response"); summary(emt_h_tot_sess)
emt_h_tot_sex = emmeans(h_total_gl, pairwise ~ condition|sex, type = "response"); summary(emt_h_tot_sex)

# Plot estimated total time by session
plot_emt_h_tot_sess <- as.data.frame(emt_h_tot_sess$emmeans)
ggplot(plot_emt_h_tot_sess, aes(x = session, y = response, color = condition, group = condition)) +
  geom_line(size = 1, position = pd) +
  geom_point(size = 4, position = pd) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1, alpha = 0.3, position = pd) +
  labs(x = "Session Number", y = "Predicted Total Time (seconds)", color = "Condition") +
  theme_minimal() + scale_color_brewer(palette = "Set3")

# Plot estimated total time by sex
plot_emt_h_tot_sex <- as.data.frame(emt_h_tot_sex$emmeans)
ggplot(plot_emt_h_tot_sex, aes(x = sex, y = response, fill = condition)) +
  geom_col(position = position_dodge(0.9)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                position = position_dodge(0.9), width = 0.2) +
  labs(x = "Sex", y = "Predicted Total Time (seconds)", fill = "Condition") +
  theme_minimal() + scale_fill_brewer(palette = "Set3")

# --- --- 3.3 BOTH, total time - NORMALISED (worse) ####

# Add "species" columns
sum_bird$species <- "bird"
sum_human$species <- "human"

# Combine the two data frames
sum_both <- rbind(sum_bird, sum_human)

sum_both_norm <- sum_both %>%
  group_by(subject, condition) %>%
  mutate(
    centered_time = sum_time / sum(sum_time, na.rm = TRUE)
  ) %>%
  ungroup()
sum_both_norm$sex[sum_both_norm$sex %in% c("f")] <- "female"
sum_both_norm$sex[sum_both_norm$sex %in% c("m")] <- "male"

sum_both_norm$species <- factor(sum_both_norm$species)
sum_both_norm$sex <- factor(sum_both_norm$sex)
sum_both_norm$session <- factor(sum_both_norm$session)

performance::check_distribution(sum_both_norm$centered_time)

# Standard transformation to pull 0/1 slightly inward
n <- nrow(sum_both_norm)
sum_both_norm$centered_time_scaled <- (sum_both_norm$centered_time * (n - 1) + 0.5) / n

m_both_total_beta = glmmTMB(centered_time_scaled ~ condition * (session + sex) + species + (1+condition|subject), 
                            dispformula = ~condition + sex + species,
                            data = sum_both_norm, 
                            family = beta_family(link = "logit"))

# The _mid model failed outlier and dispersion test
m_both_total_beta_mid = glmmTMB(centered_time_scaled ~ condition * (session + sex + species) + (1+condition|subject), 
                                dispformula = ~condition + sex + species,
                                data = sum_both_norm, 
                                family = beta_family(link = "logit"))
# Check diagnostics
simulation_output <- simulateResiduals(m_both_total_beta); plot(simulation_output)
testResiduals(simulation_output)

summary(m_both_total_beta)
Anova(m_both_total_beta, type = 3); #summary(m1)

emmeans(m_both_total_beta, pairwise ~ condition|session, type = "response")


# --- --- 4.1 - BIRDS, per trial ####

mean_bird <- Exp_I %>%
  group_by(subject, session, sex, condition) %>%
  dplyr::summarise(
    mean_time = mean(time, na.rm = TRUE),
    n_bouts = n(), # Good to keep track of how many bouts made up that mean
    .groups = "drop"
  )


performance::check_distribution(Exp_I$time)

b_trial_gl = glmmTMB(mean_time ~ condition*session + sex + 
                       (1+condition|subject), 
                     dispformula = ~ condition,
                     weights = n_bouts,
                     data = mean_bird, family = Gamma(link='log'))
b_trial_gl1 = glmmTMB(mean_time ~ condition+session + sex + 
                        (1+condition|subject), 
                      dispformula = ~ condition,
                      weights = n_bouts,
                      data = mean_bird, family = Gamma(link='log'))
b_trial_gl2 = glmmTMB(mean_time ~ session + sex + 
                        (1+condition|subject), 
                      dispformula = ~ condition,
                      weights = n_bouts,
                      data = mean_bird, family = Gamma(link='log'))
b_trial_gl3 = glmmTMB(mean_time ~ condition*(session + sex) + 
                        (1+condition|subject), 
                      dispformula = ~ condition,
                      weights = n_bouts,
                      data = mean_bird, family = Gamma(link='log'))

simulation_output <- simulateResiduals(b_trial_gl3); plot(simulation_output)

anova(b_trial_gl, b_trial_gl1, b_trial_gl2, b_trial_gl3)
# We evaluated the model structure using Likelihood Ratio Tests (LRT) and AIC. 
# The inclusion of an interaction between condition and session significantly 
# improved model fit compared to an additive model (chi^2(4) = 351.81, 
# p < 0.0001, Delta AIC = 344).
# The interaction between sex and condition was non-significant and was thus
# not included in the final model.

Anova(b_trial_gl, type = 3); #summary(m1)

#emt = emmeans(b_trial_gl, pairwise ~ session, type = "response"); summary(emt)
emt_b_trial_sess = emmeans(b_trial_gl, pairwise ~ condition|session, type = "response"); summary(emt_b_trial_sess)

plot_emt_b_trial_sess <- as.data.frame(emt_b_trial_sess$emmeans)
pd <- position_dodge(-0.2)
ggplot(plot_emt_b_trial_sess, aes(x = session, y = response, color = condition, group = condition)) +
  geom_line(size = 1, position = pd) +
  geom_point(size = 4, position = pd) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1, alpha = 0.3, position = pd) +
  labs(x = "Session Number", y = "Predicted Mean Trial Time (seconds)", color = "Condition") +
  theme_minimal() + scale_color_brewer(palette = "Set1")

# --- --- 4.2 - HUMANS, per trial ####

mean_human <- Exp_II %>%
  group_by(subject, session, sex, condition) %>%
  dplyr::summarise(
    mean_time = mean(time, na.rm = TRUE),
    n_bouts = n(), # Good to keep track of how many bouts made up that mean
    .groups = "drop"
  )

#performance::check_distribution(Exp_II$time)

#Weighting Formula in Humans led to significant deviations in
# residuals vs predicted (DHARMa)
# Dispformula removal made residuals even better.

h_trial_gl = glmmTMB(mean_time ~ condition*session + sex + 
                       (1+condition|subject), 
                     dispformula = ~ condition,
                     #weights = n_bouts,
                     data = mean_human, family = Gamma(link='log'))
h_trial_gl1 = glmmTMB(mean_time ~ condition+session + sex + 
                        (1+condition|subject), 
                      dispformula = ~ condition,
                      #weights = n_bouts,
                      data = mean_human, family = Gamma(link='log'))
h_trial_gl2 = glmmTMB(mean_time ~ session + sex + 
                        (1+condition|subject), 
                      dispformula = ~ condition,
                      #weights = n_bouts,
                      data = mean_human, family = Gamma(link='log'))
h_trial_gl3 = glmmTMB(mean_time ~ session + sex*condition + 
                        (1+condition|subject), 
                      #dispformula = ~ condition,
                      #weights = n_bouts,
                      data = mean_human, family = Gamma(link='log'))
h_trial_gl4 = glmmTMB(mean_time ~ (session + sex)*condition + 
                        (1+condition|subject), 
                      #dispformula = ~ condition,
                      #weights = n_bouts,
                      data = mean_human, family = Gamma(link='log'))

simulation_output <- simulateResiduals(h_trial_gl4); plot(simulation_output)

#3 and 4 failed to converge
anova(h_trial_gl, h_trial_gl1, h_trial_gl2)

# We evaluated the model structure using Likelihood Ratio Tests (LRT) and AIC. 

# The interaction between sex and condition was non-significant and was thus
# not included in the final model.

Anova(h_trial_gl, type = 3); #summary(m1)

emt = emmeans(h_trial_gl, pairwise ~ condition, type = "response"); summary(emt)
emt = emmeans(h_trial_gl, pairwise ~ condition|session, type = "response"); summary(emt)

plot_data <- as.data.frame(emt$emmeans)
pd <- position_dodge(-0.2)
ggplot(plot_data, aes(x = session, y = response, color = condition, group = condition)) +
  geom_line(size = 1, position=pd) +
  geom_point(size = 4, position=pd) +
  geom_errorbar(aes(ymin = response-SE, ymax = response+SE), width = 0.1, alpha = 0.3, position = pd) +
  labs(x = "Session Number", y = "Predicted Mean Bout Time (seconds)", color = "Condition") +
  theme_minimal() + scale_color_brewer(palette = "Set1")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ####
#
