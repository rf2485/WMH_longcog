library(tidyverse)

cohort <- read_csv('RedCap Reports/Cohort_WMHLongCog_DATA_2020-06-03_1048.csv')
mri_bl_scans <- read_csv('RedCap Reports/BaselineMRIScans_WMHLongCog_DATA_2020-06-03_1050.csv')

master <- right_join(cohort, mri_bl_scans) %>% select(-redcap_event_name) %>%
  select(subject_id, scan_event_name, everything())

write_csv(master, 'RedCap Reports/Master_WMHLongCog_2020-06-03.csv')

#Scans_WMHLongCog RedCap report is filtered in Excel to only include baseline MRI scans done under the 
##Tau Study protocol in Bay 3. Redcap_event_name is changed to scan_event_name and MRI event column is
##deleted. Might write up script in R in the future for better reproducability and faster report updates.