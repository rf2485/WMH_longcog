setwd("~/Dropbox (Partners HealthCare)/Dickerson lab/Projects/Active_In Progress/Flaherty White Matter Abnormalities AD/Flaherty WMH Longitudinal Cognition (Tau Study Data)")
library(tidyverse)

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

#average lh and rh thickness per Jess's email from 7/20/20
aparc_aseg_totals <- aparc_aseg_totals %>% mutate(global_MeanThickness = rowMeans(.[5:6]))

#look at differences by dx and demographics
redcap_registry <- read_csv('RedCap Reports/DickersonMasterEnrol-AllRegistryWNACCIDs_DATA_2020-08-05_2357.csv')
redcap_registry <- redcap_registry %>% select(-birth_date)
cohort <- read_csv('Master_WMHLongCog_2020-06-03.csv')
cohort <- left_join(cohort, redcap_registry)
cohort <- cohort %>% select(subject_id, id_nacc, diagnostic_cohort, sex, years_education)
age <- moca %>% ungroup() %>% filter(redcap_event_name == 't1_arm_1') %>% select(subject_id, age)
cohort <- inner_join(cohort, age)

#merge with moca and demographic data
aparc_aseg_totals <- inner_join(cohort, aparc_aseg_totals)

#separate by cohort for statistical analyses
pca <- filter(aparc_aseg_totals, diagnostic_cohort == 'PCA')
ppa_l <- filter(aparc_aseg_totals, diagnostic_cohort == 'PPA_L')

thickness_age <- lm(formula = global_MeanThickness ~ age, data = aparc_aseg_totals)
thickness_wmh <- lm(formula = global_MeanThickness ~ corrected_wmh, data = aparc_aseg_totals)

wmh <- lm(formula = corrected_wmh ~ age, data = aparc_aseg_totals)
wmh_cohort <- lm(formula = corrected_wmh ~ age + diagnostic_cohort, data = aparc_aseg_totals)
anova(wmh, wmh_cohort) #wmh_cohort better

dx_data <- aparc_aseg_totals %>% filter(diagnostic_cohort != 'AD')
dx_wmh <- lm(formula = corrected_wmh ~ diagnostic_cohort, data = dx_data)
dx_thickness <- lm(formula = global_MeanThickness ~ diagnostic_cohort, data = dx_data)
dx_age <- lm(formula = age ~ diagnostic_cohort, data = dx_data)

wmh_cohort_interaction <- lm(formula = corrected_wmh ~ age * diagnostic_cohort, data = aparc_aseg_totals)
anova(wmh_cohort, wmh_cohort_interaction) #wmh_cohort better
summary(wmh_cohort)

wmh_thickness <- lm(formula = corrected_wmh ~ age + diagnostic_cohort + global_MeanThickness, data = aparc_aseg_totals)
anova(wmh_cohort, wmh_thickness) #wmh_thickness better
summary(wmh_thickness)

wmh_thickness_interaction <- lm(lm(formula = corrected_wmh ~ age * diagnostic_cohort + global_MeanThickness, data = aparc_aseg_totals))
anova(wmh_thickness, wmh_thickness_interaction) #wmh_thickness better

ggplot(aparc_aseg_totals, aes(age, aparc_aseg_totals$corrected_wmh, color = diagnostic_cohort, shape = diagnostic_cohort)) +
  geom_point() +
  geom_smooth(method="glm", se=FALSE) +
  theme_bw() +
  #scale_color_manual(values = c('EOAD' = 'purple', 'Control' = 'darkgreen')) +
  #theme(legend.position = c(0.2,0.5)) +
  xlab('Age (years)') + 
  ylab(expression(atop('WM Hypointensity Volume', paste('log(mm'^3, ')', sep = '')))) +
  #stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
  theme(legend.text = element_text(size = 26), axis.title = element_text(size = 26), title = element_text(size = 26), legend.title = element_text(size = 26, face = 'bold'))

ggplot(aparc_aseg_totals, aes(global_MeanThickness, aparc_aseg_totals$corrected_wmh, color = diagnostic_cohort, shape = diagnostic_cohort)) +
  geom_point() +
  geom_smooth(method="glm", se=FALSE) +
  theme_bw() +
  #scale_color_manual(values = c('EOAD' = 'purple', 'Control' = 'darkgreen')) +
  #theme(legend.position = c(0.2,0.5)) +
  xlab('Age (years)') + 
  ylab(expression(atop('Global Mean Thickness', paste('mm'^3, ')', sep = '')))) +
  #stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
  theme(legend.text = element_text(size = 26), axis.title = element_text(size = 26), title = element_text(size = 26), legend.title = element_text(size = 26, face = 'bold'))


