setwd("~/Dropbox (Partners HealthCare)/Dickerson lab/Projects/Active_In Progress/Flaherty White Matter Abnormalities AD/Flaherty WMH Longitudinal Cognition (Tau Study Data)")
library(tidyverse)
library(gridExtra)
library(ggpubr)

redcap_registry <- read_csv('RedCap Reports/DickersonMasterEnrol-AllRegistryWNACCIDs_DATA_2020-08-05_2357.csv')
redcap_registry <- redcap_registry %>% select(-birth_date)
cohort <- read_csv('Master_WMHLongCog_2020-06-03.csv')
cohort <- left_join(cohort, redcap_registry)
cohort <- cohort %>% select(subject_id, id_nacc, diagnostic_cohort, sex, years_education)
cohort$sex <- factor(cohort$sex, levels = c('M','F'), labels = c(0,1))
cohort <- mutate(cohort, sex = case_when(sex == 0 ~ as.numeric(0),
                                   sex == 1 ~ as.numeric(1)))#for z score calculation


moca <- read_csv('RedCap Reports/NACCDownloads-MoCAScores_DATA_2020-08-05_2332.csv')
moca <- moca %>% select(-redcap_repeat_instrument, -redcap_repeat_instance, -sex, -educ) %>% rename(id_nacc = ptid) %>%
  inner_join(cohort,.) %>% 
  mutate(mocatots = as.numeric(mocatots)) %>% 
  mutate_at(vars(mocatots:mocaorct), ~case_when(. = . != 88 ~ .)) %>%
  filter(visityr >= 2017)

####Global Scores

moca_full <- moca %>% filter_all(all_vars(!grepl(98,.))) %>% #exclude subjects with verbal refusals
  mutate_at(vars(mocatots : mocaorct), ~case_when(. = . %in% c(95,96,97) ~ 0, #change 'missing' scores to 0
                                                  . = . != 88 ~ .)) %>% #only keep values not equal to 88
  mutate(mocatots = mocatrai+ mocacube+ mocacloc+ mocaclon+ mocacloh+ mocanami+ mocadigi + mocalett + mocaser7 + 
           mocarepe + mocaflue+ mocaabst+ mocarecn+ mocaordt+ mocaormo+ mocaoryr+ mocaordy+ mocaorpl+ mocaorct)  %>%
  select(subject_id, redcap_event_name, mocatots)

moca <- full_join(moca, moca_full, by = c('subject_id', 'redcap_event_name')) %>%
  mutate(mocatots = coalesce(mocatots.x,mocatots.y),
         age = visityr-birthyr) %>%
  mutate(mocatots = case_when(mocatots != 88 ~ mocatots)) %>% #only keep total scores not equal to 88
  select(-mocatots.x,-mocatots.y) %>%
  group_by(id_nacc) %>% arrange(subject_id, redcap_event_name) %>%
  mutate(total_change_from_baseline = mocatots - mocatots[1],
         moca_zscore = (mocatots-(26.187 + (0.351*sex) + (-0.077*age) + (0.332*years_education)))/2.485558132)

#moca_t1 <- filter(moca, redcap_event_name == 't1_arm_1')

moca_t1 <- moca %>% filter(redcap_event_name == 't1_arm_1') %>% mutate(total_change_from_baseline = na_if(total_change_from_baseline,0))
moca_long <- moca %>% filter(redcap_event_name != 't1_arm_1')
moca <- rbind(moca_t1, moca_long) %>% group_by(id_nacc) %>% arrange(subject_id, redcap_event_name)
  #z score formula taken from NACC Z Score Calculator
 # mutate(moca_zscore = (mocatots-(26.187 + (0.351*sex) + (-0.077*age) + (0.332*years_education)))/2.485558132)

#####Subscores

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

moca <- moca %>% #change from baseline
  group_by(id_nacc) %>% arrange(subject_id, redcap_event_name) %>%
  mutate(
    moca_memory_change = moca_memory_index - moca_memory_index[1],
    moca_executive_change = moca_executive_index - moca_executive_index[1],
    moca_attention_change = moca_attention_index - moca_attention_index[1],
    moca_language_change = moca_language_index - moca_language_index[1],
    moca_visuospatial_change = moca_visuospatial_index - moca_visuospatial_index[1],
    moca_orientation_change = moca_orientation_index - moca_orientation_index[1] 
    )

moca_t1 <- moca %>% filter(redcap_event_name == 't1_arm_1') %>% #all Change variables turned into NA at t1_arm_1
  mutate_at(vars(ends_with("Change")), ~na_if(.,0))
moca_long <- moca %>% filter(redcap_event_name != 't1_arm_1')
moca_t2 <- moca %>% filter(redcap_event_name == 't2_arm_1')
moca <- rbind(moca_t1, moca_long) %>% group_by(id_nacc) %>% arrange(subject_id, redcap_event_name) %>% select(subject_id:visityr, birthyr, age, everything())


####Summary Statistics
means_total <- colMeans(moca[sapply(moca, is.numeric)]) %>% bind_rows()
sd_total <- sapply(moca, sd) %>% bind_rows()

timepoint_means_nas <- moca %>% group_by(redcap_event_name) %>% 
  select(redcap_event_name,mocatots, moca_zscore, moca_executive_index, moca_attention_index) %>% 
  na.omit() %>% summarise_all(list(~mean(.), ~sd(.)))
timepoint_means_long <- moca %>% group_by(redcap_event_name) %>% select(total_change_from_baseline, moca_executive_change, moca_attention_change)%>%
  na.omit() %>% summarise_all(list(~mean(.), ~sd(.)))
timepoint_means <- moca %>% group_by(redcap_event_name) %>% select(-subject_id:-mocahear) %>% 
  select(-mocatots,-total_change_from_baseline, -moca_zscore, -moca_executive_index, -moca_attention_index, -moca_executive_change, -moca_attention_change) %>%
  summarise_all(list(~mean(.), ~sd(.)))
timepoint_means <- full_join(timepoint_means, timepoint_means_nas) %>% full_join(., timepoint_means_long)

t1_dx_means_nas <- moca_t1 %>% group_by(diagnostic_cohort) %>% 
  select(mocatots, moca_zscore, moca_executive_index, moca_attention_index) %>% 
  na.omit() %>% summarise_all(list(~mean(.), ~sd(.)))
t1_dx_means <- moca_t1 %>% group_by(diagnostic_cohort) %>% select(-subject_id:-mocahear) %>% 
  select(-mocatots, -moca_zscore, -moca_executive_index, -moca_attention_index, -total_change_from_baseline, -moca_executive_change, -moca_attention_change, -moca_memory_change, 
         -moca_language_change, -moca_visuospatial_change, -moca_orientation_change) %>%
  summarise_all(list(~mean(.), ~sd(.)))
t1_dx_means <- full_join(t1_dx_means, t1_dx_means_nas)

timepoint_dx_means_nas <- moca %>% group_by(redcap_event_name, diagnostic_cohort) %>% 
  select(redcap_event_name,mocatots, moca_zscore, moca_executive_index, moca_attention_index) %>% 
  na.omit() %>% summarise_all(list(~mean(.), ~sd(.)))
timepoint_dx_means_long <- moca %>% group_by(redcap_event_name, diagnostic_cohort) %>% select(total_change_from_baseline, moca_executive_change, moca_attention_change) %>%
  na.omit() %>% summarise_all(list(~mean(.), ~sd(.)))
timepoint_dx_means <- moca %>% group_by(redcap_event_name, diagnostic_cohort) %>% select(-subject_id:-mocahear) %>% 
  select(-mocatots,-total_change_from_baseline, -moca_zscore, -moca_executive_index, -moca_attention_index, -moca_executive_change, -moca_attention_change) %>%
  summarise_all(list(~mean(.), ~sd(.)))
timepoint_dx_means <- full_join(timepoint_dx_means, timepoint_dx_means_nas) %>% full_join(., timepoint_dx_means_long)

####Plots

#total scores at t1 (histogram)
ggplot(moca_t1, aes(x = mocatots)) +
  geom_histogram(alpha = 0.8, fill = 'blue') +
  geom_vline(xintercept = mean(na.omit(moca$mocatots)), size = 1) +
  geom_vline(xintercept = mean(na.omit(moca$mocatots)) - (sd(na.omit(moca$mocatots))), linetype = 'dashed') +
  geom_vline(xintercept = mean(na.omit(moca$mocatots)) + (sd(na.omit(moca$mocatots))), linetype = 'dashed') +
  theme_bw()

#total z scores at t1 (histogram)
ggplot(moca_t1, aes(x = moca_zscore)) +
  geom_histogram(alpha = 0.8, fill = 'blue') +
  geom_vline(xintercept = mean(na.omit(moca$moca_zscore)), size = 1) +
  geom_vline(xintercept = mean(na.omit(moca$moca_zscore)) - (sd(na.omit(moca$moca_zscore))), linetype = 'dashed') +
  geom_vline(xintercept = mean(na.omit(moca$moca_zscore)) + (sd(na.omit(moca$moca_zscore))), linetype = 'dashed') +
  theme_bw()

#change scores at t2 (histogram)
ggplot(moca_t2, aes(x = total_change_from_baseline)) +
  geom_histogram(alpha = 0.8) +
  geom_vline(xintercept = mean(na.omit(moca$total_change_from_baseline)), size = 1) +
  geom_vline(xintercept = mean(na.omit(moca$total_change_from_baseline)) - (sd(na.omit(moca$total_change_from_baseline))), linetype = 'dashed') +
  geom_vline(xintercept = mean(na.omit(moca$total_change_from_baseline)) + (sd(na.omit(moca$total_change_from_baseline))), linetype = 'dashed') +
  theme_bw()

#t1 index scores (histogram)
memory_index <- ggplot(moca_t1, aes(x = moca_memory_index)) +
  geom_histogram(alpha = 0.8) +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_memory_index)), size = 1) +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_memory_index)) - (sd(na.omit(moca_t1$moca_memory_index))), linetype = 'dashed') +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_memory_index)) + (sd(na.omit(moca_t1$moca_memory_index))), linetype = 'dashed') +
  theme_bw() + theme(legend.position = 'none')

executive_index <- ggplot(moca_t1, aes(x = moca_executive_index)) +
  geom_histogram(alpha = 0.8) +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_executive_index)), size = 1) +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_executive_index)) - (sd(na.omit(moca_t1$moca_executive_index))), linetype = 'dashed') +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_executive_index)) + (sd(na.omit(moca_t1$moca_executive_index))), linetype = 'dashed') +
  theme_bw() + theme(legend.position = 'none')

attention_index <- ggplot(moca_t1, aes(x = moca_attention_index)) +
  geom_histogram(alpha = 0.8) +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_attention_index)), size = 1) +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_attention_index)) - (sd(na.omit(moca_t1$moca_attention_index))), linetype = 'dashed') +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_attention_index)) + (sd(na.omit(moca_t1$moca_attention_index))), linetype = 'dashed') +
  theme_bw()+ theme(legend.position = 'none')

language_index <- ggplot(moca_t1, aes(x = moca_language_index)) +
  geom_histogram(alpha = 0.8) +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_language_index)), size = 1) +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_language_index)) - (sd(na.omit(moca_t1$moca_language_index))), linetype = 'dashed') +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_language_index)) + (sd(na.omit(moca_t1$moca_language_index))), linetype = 'dashed') +
  theme_bw()+ theme(legend.position = 'none')

visuospatial_index <- ggplot(moca_t1, aes(x = moca_visuospatial_index)) +
  geom_histogram(alpha = 0.8) +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_visuospatial_index)), size = 1) +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_visuospatial_index)) - (sd(na.omit(moca_t1$moca_visuospatial_index))), linetype = 'dashed') +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_visuospatial_index)) + (sd(na.omit(moca_t1$moca_visuospatial_index))), linetype = 'dashed') +
  theme_bw()+ theme(legend.position = 'none')

orientation_index <- ggplot(moca_t1, aes(x = moca_orientation_index)) +
  geom_histogram(alpha = 0.8) +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_orientation_index)), size = 1) +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_orientation_index)) - (sd(na.omit(moca_t1$moca_orientation_index))), linetype = 'dashed') +
  geom_vline(xintercept = mean(na.omit(moca_t1$moca_orientation_index)) + (sd(na.omit(moca_t1$moca_orientation_index))), linetype = 'dashed') +
  theme_bw()+ theme(legend.position = 'none')

grid.arrange(memory_index,executive_index,attention_index,language_index,visuospatial_index,orientation_index, nrow = 3)

g <- arrangeGrob(memory_index,executive_index,attention_index,language_index,visuospatial_index,orientation_index, nrow = 3)
ggsave("histogram_t1_moca_index_scores.pdf", g)

#total scores by timepoint (bar chart)
ggplot(timepoint_means, aes(x=redcap_event_name, y=mocatots_mean, fill = redcap_event_name)) + 
  geom_point(data = moca, aes(x=redcap_event_name, y = mocatots)) +
  geom_col(alpha=0.5) +
  geom_errorbar(data = timepoint_means, aes(ymin = mocatots_mean-mocatots_sd, ymax = mocatots_mean+mocatots_sd)) +
  stat_compare_means(comparisons = list(c('t1_arm_1','t2_arm_1'), c('t1_arm_1','t3_arm_1')), label.y = c(35,40),
                    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
  theme_bw()

#total scores by cohort (bar chart)
ggplot(t1_dx_means, aes(x=diagnostic_cohort, y=mocatots_mean, fill = diagnostic_cohort)) + 
  geom_point(data = moca, aes(x=diagnostic_cohort, y = mocatots)) +
  geom_col(alpha=0.5) +
  geom_errorbar(data = t1_dx_means, aes(ymin = mocatots_mean-mocatots_sd, ymax = mocatots_mean+mocatots_sd)) +
  geom_signif(comparisons = list(c('PCA','PPA_L')), map_signif_level = TRUE) +
  theme_bw()

#total scores by cohort and timepoint (bar chart)
ggplot(timepoint_dx_means, aes(x=diagnostic_cohort, y=mocatots_mean, fill = redcap_event_name)) + 
  geom_col(alpha=0.5, position = 'dodge', width = 0.5) +
  geom_errorbar(aes(ymin = mocatots_mean-mocatots_sd, ymax = mocatots_mean+mocatots_sd), position = 'dodge', width = 0.5) +
  #geom_signif(comparisons = list(c('PCA','PPA_L')), map_signif_level = TRUE) +
  theme_bw()

#total z scores by timepoint (bar chart)
ggplot(timepoint_means, aes(x=redcap_event_name, y=moca_zscore_mean, fill = redcap_event_name)) + 
  geom_point(data = moca, aes(x=redcap_event_name, y = moca_zscore)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(data = timepoint_means, aes(ymin = moca_zscore_mean-moca_zscore_sd, ymax = moca_zscore_mean+moca_zscore_sd)) +
  stat_compare_means(comparisons = list(c('t1_arm_1','t2_arm_1'), c('t1_arm_1','t3_arm_1')), label.y = c(0.75, 1.5),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
  theme_bw()

#total z scores by cohort (bar chart)
ggplot(t1_dx_means, aes(x=diagnostic_cohort, y=moca_zscore_mean, fill = diagnostic_cohort)) + 
  geom_point(data = moca, aes(x=diagnostic_cohort, y = moca_zscore)) +
  geom_col(alpha=0.5) +
  geom_errorbar(data = t1_dx_means, aes(ymin = moca_zscore_mean-moca_zscore_sd, ymax = moca_zscore_mean+moca_zscore_sd)) +
  geom_signif(comparisons = list(c('PCA','PPA_L')), map_signif_level = TRUE) +
  theme_bw()

#total z scores by cohort and timepoint (bar chart)
ggplot(timepoint_dx_means, aes(x=diagnostic_cohort, y=moca_zscore_mean, fill = redcap_event_name)) + 
  geom_col(alpha=0.5, position = 'dodge', width = 0.5) +
  geom_errorbar(aes(ymin = moca_zscore_mean-moca_zscore_sd, ymax = moca_zscore_mean+moca_zscore_sd), position = 'dodge', width = 0.5) +
  #geom_signif(comparisons = list(c('PCA','PPA_L')), map_signif_level = TRUE) +
  theme_bw()

#change scores by cohort (bar chart)
ggplot(timepoint_dx_means, aes(x=diagnostic_cohort, y=total_change_from_baseline_mean, fill = redcap_event_name)) + 
  geom_col(alpha=0.5, position = 'dodge', width = 0.5) +
  geom_errorbar(aes(ymin = total_change_from_baseline_mean-total_change_from_baseline_sd, ymax = total_change_from_baseline_mean+total_change_from_baseline_sd), position = 'dodge', width = 0.5) +
  #geom_signif(comparisons = list(c('PCA','PPA_L')), map_signif_level = TRUE) +
  theme_bw()

#index scores by timepoint (boxplot)
memory_index <- ggplot(timepoint_means, aes(x=redcap_event_name, y=moca_memory_index_mean, fill = redcap_event_name)) + 
  geom_point(data = moca, aes(x=redcap_event_name, y = moca_memory_index)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(data = timepoint_means, aes(ymin = moca_memory_index_mean-moca_memory_index_sd, ymax = moca_memory_index_mean+moca_memory_index_sd)) +
  geom_signif(comparisons = list(c('t1_arm_1','t2_arm_1')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

executive_index <- ggplot(timepoint_means, aes(x=redcap_event_name, y=moca_executive_index_mean, fill = redcap_event_name)) + 
  geom_point(data = moca, aes(x=redcap_event_name, y = moca_executive_index)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(data = timepoint_means, aes(ymin = moca_executive_index_mean-moca_executive_index_sd, ymax = moca_executive_index_mean+moca_executive_index_sd)) +
  geom_signif(comparisons = list(c('t1_arm_1','t2_arm_1')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

attention_index <- ggplot(timepoint_means, aes(x=redcap_event_name, y=moca_attention_index_mean, fill = redcap_event_name)) + 
  geom_point(data = moca, aes(x=redcap_event_name, y = moca_attention_index)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(data = timepoint_means, aes(ymin = moca_attention_index_mean-moca_attention_index_sd, ymax = moca_attention_index_mean+moca_attention_index_sd)) +
  geom_signif(comparisons = list(c('t1_arm_1','t2_arm_1')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

language_index <- ggplot(timepoint_means, aes(x=redcap_event_name, y=moca_language_index_mean, fill = redcap_event_name)) + 
  geom_point(data = moca, aes(x=redcap_event_name, y = moca_language_index)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(data = timepoint_means, aes(ymin = moca_language_index_mean-moca_language_index_sd, ymax = moca_language_index_mean+moca_language_index_sd)) +
  geom_signif(comparisons = list(c('t1_arm_1','t2_arm_1')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

visuospatial_index <- ggplot(timepoint_means, aes(x=redcap_event_name, y=moca_visuospatial_index_mean, fill = redcap_event_name)) + 
  geom_point(data = moca, aes(x=redcap_event_name, y = moca_visuospatial_index)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(data = timepoint_means, aes(ymin = moca_visuospatial_index_mean-moca_visuospatial_index_sd, ymax = moca_visuospatial_index_mean+moca_visuospatial_index_sd)) +
  geom_signif(comparisons = list(c('t1_arm_1','t2_arm_1')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

orientation_index <- ggplot(timepoint_means, aes(x=redcap_event_name, y=moca_orientation_index_mean, fill = redcap_event_name)) + 
  geom_point(data = moca, aes(x=redcap_event_name, y = moca_orientation_index)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(data = timepoint_means, aes(ymin = moca_orientation_index_mean-moca_orientation_index_sd, ymax = moca_orientation_index_mean+moca_orientation_index_sd)) +
  geom_signif(comparisons = list(c('t1_arm_1','t2_arm_1')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

grid.arrange(memory_index,executive_index,attention_index,language_index,visuospatial_index,orientation_index, nrow = 3)


#index scores by cohort (boxplot)
memory_index <- ggplot(t1_dx_means, aes(x=diagnostic_cohort, y=moca_memory_index_mean, fill = diagnostic_cohort)) + 
  geom_point(data = moca_t1, aes(x=diagnostic_cohort, y = moca_memory_index)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(data = t1_dx_means, aes(ymin = moca_memory_index_mean-moca_memory_index_sd, ymax = moca_memory_index_mean+moca_memory_index_sd)) +
  geom_signif(comparisons = list(c('PCA', 'PPA_L')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

executive_index <- ggplot(t1_dx_means, aes(x=diagnostic_cohort, y=moca_executive_index_mean, fill = diagnostic_cohort)) + 
  geom_point(data = moca_t1, aes(x=diagnostic_cohort, y = moca_executive_index)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(data = t1_dx_means, aes(ymin = moca_executive_index_mean-moca_executive_index_sd, ymax = moca_executive_index_mean+moca_executive_index_sd)) +
  geom_signif(comparisons = list(c('PCA', 'PPA_L')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

attention_index <- ggplot(t1_dx_means, aes(x=diagnostic_cohort, y=moca_attention_index_mean, fill = diagnostic_cohort)) + 
  geom_point(data = moca_t1, aes(x=diagnostic_cohort, y = moca_attention_index)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(data = t1_dx_means, aes(ymin = moca_attention_index_mean-moca_attention_index_sd, ymax = moca_attention_index_mean+moca_attention_index_sd)) +
  geom_signif(comparisons = list(c('PCA', 'PPA_L')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

language_index <- ggplot(t1_dx_means, aes(x=diagnostic_cohort, y=moca_language_index_mean, fill = diagnostic_cohort)) + 
  geom_point(data = moca_t1, aes(x=diagnostic_cohort, y = moca_language_index)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(data = t1_dx_means, aes(ymin = moca_language_index_mean-moca_language_index_sd, ymax = moca_language_index_mean+moca_language_index_sd)) +
  geom_signif(comparisons = list(c('PCA', 'PPA_L')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

visuospatial_index <- ggplot(t1_dx_means, aes(x=diagnostic_cohort, y=moca_visuospatial_index_mean, fill = diagnostic_cohort)) + 
  geom_point(data = moca_t1, aes(x=diagnostic_cohort, y = moca_visuospatial_index)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(data = t1_dx_means, aes(ymin = moca_visuospatial_index_mean-moca_visuospatial_index_sd, ymax = moca_visuospatial_index_mean+moca_visuospatial_index_sd)) +
  geom_signif(comparisons = list(c('PCA', 'PPA_L')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

orientation_index <- ggplot(t1_dx_means, aes(x=diagnostic_cohort, y=moca_orientation_index_mean, fill = diagnostic_cohort)) + 
  geom_point(data = moca_t1, aes(x=diagnostic_cohort, y = moca_orientation_index)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(data = t1_dx_means, aes(ymin = moca_orientation_index_mean-moca_orientation_index_sd, ymax = moca_orientation_index_mean+moca_orientation_index_sd)) +
  geom_signif(comparisons = list(c('PCA', 'PPA_L')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

grid.arrange(memory_index,executive_index,attention_index,language_index,visuospatial_index,orientation_index, nrow = 3)
g <- arrangeGrob(memory_index,executive_index,attention_index,language_index,visuospatial_index,orientation_index, nrow = 3)
ggsave("barchart_t1_moca_index_scores.pdf", g)

#index scores by cohort and timepoint (boxplot)

memory_index <- ggplot(timepoint_dx_means, aes(x=diagnostic_cohort, y=moca_memory_index_mean, fill = redcap_event_name)) + 
  geom_col(alpha = 0.5, position = 'dodge', width = 0.5) +
  geom_errorbar(aes(ymin = moca_memory_index_mean-moca_memory_index_sd, ymax = moca_memory_index_mean+moca_memory_index_sd), position = 'dodge', width = 0.5) +
  geom_signif(comparisons = list(c('PCA', 'PPA_L')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

executive_index <- ggplot(timepoint_dx_means, aes(x=diagnostic_cohort, y=moca_executive_index_mean, fill = redcap_event_name)) + 
  geom_col(alpha = 0.5, position = 'dodge', width = 0.5) +
  geom_errorbar(aes(ymin = moca_executive_index_mean-moca_executive_index_sd, ymax = moca_executive_index_mean+moca_executive_index_sd), position = 'dodge', width = 0.5) +
  geom_signif(comparisons = list(c('PCA', 'PPA_L')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

attention_index <- ggplot(timepoint_dx_means, aes(x=diagnostic_cohort, y=moca_attention_index_mean, fill = redcap_event_name)) + 
  geom_col(alpha = 0.5, position = 'dodge', width = 0.5) +
  geom_errorbar(aes(ymin = moca_attention_index_mean-moca_attention_index_sd, ymax = moca_attention_index_mean+moca_attention_index_sd), position = 'dodge', width = 0.5) +
  geom_signif(comparisons = list(c('PCA', 'PPA_L')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

language_index <- ggplot(timepoint_dx_means, aes(x=diagnostic_cohort, y=moca_language_index_mean, fill = redcap_event_name)) + 
  geom_col(alpha = 0.5, position = 'dodge', width = 0.5) +
  geom_errorbar(aes(ymin = moca_language_index_mean-moca_language_index_sd, ymax = moca_language_index_mean+moca_language_index_sd), position = 'dodge', width = 0.5) +
  geom_signif(comparisons = list(c('PCA', 'PPA_L')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

visuospatial_index <- ggplot(timepoint_dx_means, aes(x=diagnostic_cohort, y=moca_visuospatial_index_mean, fill = redcap_event_name)) + 
  geom_col(alpha = 0.5, position = 'dodge', width = 0.5) +
  geom_errorbar(aes(ymin = moca_visuospatial_index_mean-moca_visuospatial_index_sd, ymax = moca_visuospatial_index_mean+moca_visuospatial_index_sd), position = 'dodge', width = 0.5) +
  geom_signif(comparisons = list(c('PCA', 'PPA_L')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

orientation_index <- ggplot(timepoint_dx_means, aes(x=diagnostic_cohort, y=moca_orientation_index_mean, fill = redcap_event_name)) + 
  geom_col(alpha = 0.5, position = 'dodge', width = 0.5) +
  geom_errorbar(aes(ymin = moca_orientation_index_mean-moca_orientation_index_sd, ymax = moca_orientation_index_mean+moca_orientation_index_sd), position = 'dodge', width = 0.5) +
  geom_signif(comparisons = list(c('PCA', 'PPA_L')), map_signif_level = TRUE) +
  theme_bw()+ theme(legend.position = 'none')

grid.arrange(memory_index,executive_index,attention_index,language_index,visuospatial_index,orientation_index, nrow = 3)
