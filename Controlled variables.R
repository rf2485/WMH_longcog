#run WMH and Index Scores first (before models)

#check to see if there is a difference in age or thickness. If so, regress out age/thickness.

#check for statisitcal significance
wmh_age <- lm(formula = corrected_wmh ~ age, data = wmh_and_moca)
wmh_thick <- lm(formula = corrected_wmh ~ global_MeanThickness, data = wmh_and_moca)
memory_age <- lm(formula = moca_memory_index ~ age, data = wmh_and_moca)
memory_thick <- lm(formula = moca_memory_index ~ global_MeanThickness, data = wmh_and_moca)
exec_age <- lm(formula = moca_executive_index ~ age, data = wmh_and_moca)
exec_thick <- lm(formula = moca_executive_index ~ global_MeanThickness, data = wmh_and_moca)
#attn_age <- lm(formula = moca_attention_index ~ age, data = wmh_and_moca)
attn_thick <- lm(formula = moca_attention_index ~ global_MeanThickness, data = wmh_and_moca)
lang_age <- lm(formula = moca_language_index ~ age, data = wmh_and_moca)
lang_thick <- lm(formula = moca_language_index ~ global_MeanThickness, data = wmh_and_moca)
vis_age <- lm(formula = moca_visuospatial_index ~ age, data = wmh_and_moca)
vis_thick <- lm(formula = moca_visuospatial_index ~ global_MeanThickness, data = wmh_and_moca)
ori_age <- lm(formula = moca_orientation_index ~ age, data = wmh_and_moca)
ori_thick <- lm(formula = moca_orientation_index ~ global_MeanThickness, data = wmh_and_moca)
z_thick <- lm(formula = moca_zscore ~ global_MeanThickness, data = wmh_and_moca)
tots_age <- lm(formula = mocatots ~ age, data = wmh_and_moca)
tots_thick <- lm(formula = mocatots ~ global_MeanThickness, data = wmh_and_moca)
age_thick <- lm(formula = age ~ global_MeanThickness, data = wmh_and_moca)

#regress out
wmh_and_moca$age_controlled_wmh <- wmh_and_moca$corrected_wmh - predict(wmh_age)
wmh_and_moca$age_controlled_exec <- wmh_and_moca$moca_executive_index - predict(exec_age)
wmh_and_moca$thick_controlled_wmh <- wmh_and_moca$corrected_wmh - predict(wmh_thick)
wmh_and_moca$thick_controlled_attn <- wmh_and_moca$moca_attention_index - predict(attn_thick)
wmh_and_moca$thick_controlled_vis <- wmh_and_moca$moca_visuospatial_index - predict(vis_thick)
wmh_and_moca$thick_controlled_ori <- wmh_and_moca$moca_orientation_index - predict(ori_thick)
wmh_and_moca$thick_controlled_age <- wmh_and_moca$age - predict(age_thick)

#test the models with the corrected variables
summary(lm(formula = age_controlled_wmh ~ moca_memory_index, data = wmh_and_moca))
summary(lm(formula = age_controlled_wmh ~ age_controlled_exec, data = wmh_and_moca))
thick_controlled_attn_wmh_age <- lm(formula = thick_controlled_wmh ~ thick_controlled_attn * thick_controlled_age, data = wmh_and_moca)
summary(lm(formula = age_controlled_wmh ~ thick_controlled_vis, data = wmh_and_moca))
summary(lm(formula = age_controlled_wmh ~ thick_controlled_ori, data = wmh_and_moca))

anova(wmh_attn_age_interaction, thick_controlled_attn_wmh_age)

library('sjPlot')
plot_model(thick_controlled_attn_wmh_age, type = 'int')

