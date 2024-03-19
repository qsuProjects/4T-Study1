############################## MICE ############################################

vars_imp <- c("age.at.onset", "cgm_sex", "Race.Ethnicity",     
              "cgm_dka_onset", "cgm_hba1c_onset", "cgm_insurance_type", 
              "cgm_language")

dat_gmi <- dat_cgm.4T1 %>% 
  filter(Week.of.Visit %in% wks & !is.na(GMI)) %>% 
  select(Record.ID, Week.of.Visit, GMI) %>%
  rename(gmi_week = Week.of.Visit) %>%
  mutate(record_id = paste0("4T.1-", sprintf("%03d", Record.ID))) %>%
  spread(key = gmi_week, value = GMI, sep = "") %>%
  select(-Record.ID)

df_wide <- dat_long.all %>% filter(cohort == "4T.1") %>% 
  select(record_id, studyday, cgm_hba1c, all_of(vars_imp)) %>%
  rename(a1c_day = studyday) %>%
  mutate(Public = ifelse(cgm_insurance_type == "Public", 1, 0),
         Hispanic = case_when(Race.Ethnicity == "Hispanic" ~ 1, 
                              Race.Ethnicity == "Unknown / Declined to state" ~ NA,
                              T ~ 0)) %>%   
  select(-Race.Ethnicity, -cgm_insurance_type) %>%
  spread(key = a1c_day, 
         value = cgm_hba1c,
         sep = "") %>%
  mutate_if(is.character, as.factor) %>%
  mutate_at(c("cgm_dka_onset", "Hispanic", "Public"), as.factor) %>%
  select(-a1c_day0, -c(a1c_day362:a1c_day391)) %>%
  left_join(dat_gmi, by = "record_id") 

imp_method <- c(rep("", 4), "pmm", rep("", 2), "logreg",        
                rep("pmm", (df_wide %>% select(a1c_day1:last_col()) %>% length()))) 

mi = mice(data = df_wide, m = n.imp, maxit = 5, print = T, seed = seed, method = imp_method)

save(list = c('mi'), file = "NMED-A128661-2023_mi.Rdata")