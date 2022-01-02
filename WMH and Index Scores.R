setwd("~/Dropbox (Partners HealthCare)/Dickerson lab/Projects/Active_In Progress/Flaherty White Matter Abnormalities AD/Flaherty WMH Longitudinal Cognition (Tau Study Data)")
library(tidyverse)

#import cohort, then pull years of education from RedCap registry
cohort <- read_csv('Master_WMHLongCog_2020-06-03.csv')
redcap_registry <- read_csv('RedCap Reports/DickersonMasterEnrol-AllRegistryWNACCIDs_DATA_2020-08-05_2357.csv')
redcap_registry <- redcap_registry %>% select(-birth_date)
cohort <- left_join(cohort, redcap_registry)

cohort <- cohort %>% select(subject_id, id_nacc, diagnostic_cohort, sex, years_education, scan_date)
cohort$sex <- factor(cohort$sex, levels = c('M','F'), labels = c(0,1))
cohort <- mutate(cohort, sex = case_when(sex == 0 ~ as.numeric(0),
                                         sex == 1 ~ as.numeric(1)))#for z score calculation



#baseline white matter T1 hypointensity volumes and global thickness for each subject id
aseg <- read_tsv('WMH_table.txt')
aseg_wm <- aseg %>% select(`Measure:volume`, `WM-hypointensities`, CerebralWhiteMatterVol) %>%
  mutate(corrected_wmh = (`WM-hypointensities`/CerebralWhiteMatterVol)*1000) %>%
  rename(subject_id = `Measure:volume`)
lh_aparc <- read_tsv('aparc_lh_stats.txt')
lh_aparc_total <- lh_aparc %>% select(lh.aparc.thickness, lh_MeanThickness_thickness) %>%
  rename(subject_id = lh.aparc.thickness)
rh_aparc <- read_tsv('aparc_rh_stats.txt')
rh_aparc_total <- rh_aparc %>% select(rh.aparc.thickness, rh_MeanThickness_thickness) %>%
  rename(subject_id = rh.aparc.thickness)
aparc_aseg_totals <- full_join(aseg_wm, lh_aparc_total)
aparc_aseg_totals <- full_join(aparc_aseg_totals, rh_aparc_total)
aparc_aseg_totals <- aparc_aseg_totals %>% mutate(global_MeanThickness = rowMeans(.[5:6]))

#baseline moca index scores for each subject id
moca <- read_csv('RedCap Reports/NACCDownloads-MoCAScores_DATA_2020-08-05_2332.csv')
moca <- moca %>% select(-redcap_repeat_instrument, -redcap_repeat_instance, -sex, -educ) %>% 
  rename(id_nacc = ptid) %>%
  mutate(mocatots = case_when(mocatots != 88 ~ mocatots)) %>% #only keep total scores not equal to NA
  filter(visityr >= 2017) #no visits earlier than start of Tau Study are included

moca <- left_join(cohort, moca) #pull moca scores for the cohort

moca <- moca %>% filter(!duplicated(subject_id)) %>% #select first visit (exlude subject ids duplicated earlier in the list)
  mutate(age = visityr-birthyr)

moca_full <- moca %>% filter_all(all_vars(!grepl(98,.))) %>% #exclude subjects with verbal refusals
  mutate_at(vars(mocatots : mocaorct), ~case_when(. = . %in% c(95,96,97) ~ 0, #change 'missing' scores to 0
                                                  . = . != 88 ~ .)) %>% #only keep values not equal to 88
  mutate(mocatots = mocatrai+ mocacube+ mocacloc+ mocaclon+ mocacloh+ mocanami+ mocadigi + mocalett + mocaser7 + 
           mocarepe + mocaflue+ mocaabst+ mocarecn+ mocaordt+ mocaormo+ mocaoryr+ mocaordy+ mocaorpl+ mocaorct)  %>%
  select(subject_id, redcap_event_name, mocatots)

moca <- full_join(moca, moca_full, by = c('subject_id', 'redcap_event_name')) %>%
  mutate(mocatots = coalesce(mocatots.x,mocatots.y)) %>%
  mutate(mocatots = case_when(mocatots != 88 ~ mocatots)) %>% #only keep total scores not equal to 88
  select(-mocatots.x,-mocatots.y) %>%
  group_by(id_nacc) %>% arrange(subject_id, redcap_event_name) %>%
  mutate(moca_zscore = (mocatots-(26.187 + (0.351*sex) + (-0.077*age) + (0.332*years_education)))/2.485558132)


#calculate MoCA index scores for the cohort
moca <- moca %>% ungroup() %>% #verbal refusal and not administered turned into NA, cognitive and physical problems turned into 0
  mutate_at(vars(mocarecc,mocarecr),
            ~case_when(. = . %in% c(88, 96, 95) ~ 0,
                       is.na(.) ~ 0,
                       TRUE ~ as.numeric(.))) %>% #not administered category cue and recog = 0
  mutate_at(vars(mocatrai : mocaorct),
            ~case_when(. = . %in% c(96, 95) ~ 0,
                       . != 98 ~ as.numeric(.))) 

moca <- moca %>%  #generate raw subscores
  mutate(
    moca_memory_index = (mocarecn*3) + (mocarecc*2) + mocarecr,
    moca_executive_index = mocatrai + mocacloc + `mocaclon` + `mocacloh` + `mocadigi` + `mocalett` + `mocaser7` + `mocaflue` +`mocaabst`,
    moca_attention_index = `mocaregi` + `mocadigi` + `mocalett` + `mocaser7` + `mocarepe`,
    moca_language_index = `mocanami` + `mocarepe` + `mocaflue`,
    moca_visuospatial_index = mocacube + mocacloc + `mocaclon` +`mocacloh` + `mocanami`,
    moca_orientation_index = `mocaordt` + `mocaormo` + `mocaoryr` + `mocaordy` + `mocaorpl` +mocaorct
  )

#join wmh and moca scores
wmh_and_moca <- full_join(aparc_aseg_totals, moca) %>% 
  select(subject_id:sex, birthyr, age, everything()) #rearrange columns
#manually check that scan_date is within 3 months of NACC/MoCA date (visitmo, visitday, visityr)

#####Linear Regression Analyses

#Notes from Baseline MRI Analysis (without MoCA scores)
#age and WMH correlate in this cohort (p = 0.001). 
#Adding global thickness as a covariate improves the model. 
#Correlation between WMH and Global Thickness has trendline significance (p = 0.09). 
#There is no correlation between age and global thickness.
#adding diagnosis improves the model, but we have small Ns.
#no statictically significant difference in age between dx (excluding one AD case)
#no statistically significant difference in WMH between dx (excluding one AD case)
#there is a statistically significant different in global thickness between dx (p = 0.03) (excluding one AD case)

#first without covariates
wmh_dx <- wmh_and_moca %>% filter(diagnostic_cohort != 'AD') #only look at PPA-L and PCA

wmh_age <- lm(formula = corrected_wmh ~ age, data = wmh_and_moca)
summary(wmh_age)

wmh_memory <- lm(formula = corrected_wmh ~ moca_memory_index, data = wmh_and_moca)
summary(wmh_memory) #not significant

wmh_exec <- lm(formula = corrected_wmh ~ moca_executive_index, data = wmh_and_moca)
summary(wmh_exec)

wmh_attn <- lm(formula = corrected_wmh ~ moca_attention_index, data = wmh_and_moca)
summary(wmh_attn) #trendline p = 0.08

wmh_lang <- lm(formula = corrected_wmh ~ moca_language_index, data = wmh_and_moca)
summary(wmh_lang)

wmh_vis <- lm(formula = corrected_wmh ~ moca_visuospatial_index, data = wmh_and_moca)
summary(wmh_vis)

wmh_ori <- lm(formula = corrected_wmh ~ moca_orientation_index, data = wmh_and_moca)
summary(wmh_ori)


#Does adding cohort (+) improve the model? 
#if yes, check for interactions (*)
#if interactions don't improve the model, check statistical significance of the "adding cohort" model
#if interactions improve the model, check statistical significance of interactions model.

wmh_dx <- wmh_and_moca %>% filter(diagnostic_cohort != 'AD') #only look at PPA-L and PCA

wmh_memory_cohort <- lm(formula = corrected_wmh ~ moca_memory_index + diagnostic_cohort, data = wmh_dx)
anova(wmh_memory, wmh_memory_cohort)

wmh_exec_cohort <- lm(formula = corrected_wmh ~ moca_executive_index + diagnostic_cohort, data = wmh_dx)
anova(wmh_exec, wmh_exec_cohort)

wmh_attn_cohort <- lm(formula = corrected_wmh ~ moca_attention_index + diagnostic_cohort, data = wmh_dx)
anova(wmh_attn, wmh_attn_cohort)

wmh_lang_cohort <- lm(formula = corrected_wmh ~ moca_language_index + diagnostic_cohort, data = wmh_dx)
anova(wmh_lang, wmh_lang_cohort)

wmh_vis_cohort <- lm(formula = corrected_wmh ~ moca_visuospatial_index + diagnostic_cohort, data = wmh_dx)
anova(wmh_vis, wmh_vis_cohort)

wmh_ori_cohort <- lm(formula = corrected_wmh ~ moca_orientation_index + diagnostic_cohort, data = wmh_dx)
anova(wmh_ori, wmh_ori_cohort) #trendline p=0.07
wmh_ori_cohort_interaction <- lm(formula = corrected_wmh ~ moca_orientation_index * diagnostic_cohort, data = wmh_dx)
anova(wmh_ori_cohort, wmh_ori_cohort_interaction)
summary(wmh_ori_cohort) 


#does adding age improve the model?

wmh_memory_age <- lm(formula = corrected_wmh ~ moca_memory_index + age, data = wmh_and_moca)
anova(wmh_memory, wmh_memory_age) #p = 0.002
wmh_memory_age_interaction <- lm(formula = corrected_wmh ~moca_memory_index * age, data = wmh_and_moca)
anova(wmh_memory_age, wmh_memory_age_interaction)
summary(wmh_memory_age) #wmh and memory not significant, wmh and age significant (p = 0.0019)
anova(wmh_age, wmh_memory_age)

wmh_exec_age <- lm(formula = corrected_wmh ~ moca_executive_index + age, data = wmh_and_moca)
anova(wmh_exec, wmh_exec_age)
wmh_exec_age_interaction <- lm(formula = corrected_wmh ~ moca_executive_index * age, data = wmh_and_moca)
anova(wmh_exec_age, wmh_exec_age_interaction)
summary(wmh_exec_age) #wmh and memory not significant, wmh and age significant (p = 0.007)


wmh_attn_age <- lm(formula = corrected_wmh ~ moca_attention_index + age, data = wmh_and_moca)
anova(wmh_attn, wmh_attn_age)
wmh_attn_age_interaction <- lm(formula = corrected_wmh ~ moca_attention_index * age, data = wmh_and_moca)
anova(wmh_attn_age, wmh_attn_age_interaction)
wmh_attn_age_thick_interaction <- lm(formula = corrected_wmh ~ moca_attention_index * age * global_MeanThickness, data = wmh_and_moca)
anova(wmh_attn_age_interaction, wmh_attn_age_thick_interaction)
summary(wmh_attn_age_interaction) 
#there is a main effect of attention on wmh, with an interaction effect between attention and age

library('sjPlot')
plot_model(wmh_attn_age_interaction, type = 'int')

#to plot this, bin age into 1 SD below the mean or less, around the mean, and 1 SD above the mean or more

age_sd <- c(-Inf,mean(wmh_and_moca$age),Inf)
age_sd <- round(age_sd) #rounded because base ages are calculated from years only
wmh_and_moca$age_bins <- cut(wmh_and_moca$age,
                             breaks = age_sd,
                             labels = c("ages 52-68","ages 68-83"))
wmh_and_moca <- arrange(wmh_and_moca, age_bins)

library('viridis')
ggplot(wmh_and_moca, aes(moca_attention_index, corrected_wmh, color = age)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y ~ x) +
  theme_bw() +
  scale_color_viridis()


wmh_lang_age <- lm(formula = corrected_wmh ~ moca_language_index + age, data = wmh_and_moca)
anova(wmh_lang, wmh_lang_age)
anova(wmh_age, wmh_lang_age)
wmh_lang_age_interaction <- lm(formula = corrected_wmh ~ moca_language_index * age, data = wmh_and_moca)
anova(wmh_lang_age, wmh_lang_age_interaction)
summary(wmh_lang_age) #wmh and age correlate, but not wmh and language

wmh_vis_age <- lm(formula = corrected_wmh ~ moca_visuospatial_index + age, data = wmh_and_moca)
anova(wmh_vis, wmh_vis_age)
anova(wmh_age, wmh_vis_age)
wmh_vis_age_interaction <- lm(formula = corrected_wmh ~ moca_visuospatial_index * age, data = wmh_and_moca)
anova(wmh_vis_age, wmh_vis_age_interaction)
summary(wmh_vis_age)  #wmh and age correlate, but not wmh and visuospatial

wmh_ori_age <- lm(formula = corrected_wmh ~ moca_orientation_index +age, data = wmh_and_moca)
anova(wmh_ori, wmh_ori_age)
anova(wmh_age, wmh_ori_age)
wmh_ori_age_interaction <- lm(formula = corrected_wmh ~ moca_orientation_index *age, data = wmh_and_moca)
anova(wmh_ori_age, wmh_ori_age_interaction)
summary(wmh_ori_age) #wmh and age correlate, but not wmh and visuospatial

wmh_tots <- lm(formula = corrected_wmh ~ mocatots, data = wmh_and_moca)
wmh_tots_age <- lm(formula = corrected_wmh ~ mocatots + age, data = wmh_and_moca)
anova(wmh_tots, wmh_tots_age)
anova(wmh_age, wmh_tots_age)
wmh_tots_age_interaction <- lm(formula = corrected_wmh ~ mocatots * age, data = wmh_and_moca)
anova(wmh_tots_age, wmh_tots_age_interaction)
summary(wmh_tots_age)

wmh_z <- lm(formula = corrected_wmh ~ moca_zscore, data = wmh_and_moca)
summary(wmh_z)
