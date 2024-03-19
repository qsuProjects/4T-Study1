############################## Import data #####################################

mon.fu = 12

dat.hist0 <- readRDS(paste0("REDCap_pull_hist_", cutdate, ".rds"))  
dat.p4T0 <- readRDS(paste0("REDCap_pull_p4T_", cutdate, ".rds"))
dat.4T0 <- readRDS(paste0("REDCap_pull_4T1_", cutdate, ".rds"))
  
ID_exclude.p4T <- dat.p4T0 %>% filter(ltfu_oth %in% "not T1D") %>%
  pull(record_id) 
ID_exclude.4T1 <- c(53, 29, 38, 40, 74) 

dat_cgm.p4T <- readRDS("NMED-A128661-2023_cgm_p4T.rds")
dat_cgm.4T1 <- readRDS("NMED-A128661-2023_cgm_4T1.rds")

dat_tir.p4T <- readRDS("NMED-A128661-2023_tir_p4T.rds")
dat_tir.4T1 <- readRDS("NMED-A128661-2023_tir_4T1.rds")

dat_wear <- rbind(readRDS("NMED-A128661-2023_wear_p4T.rds"),
                  readRDS("NMED-A128661-2023_wear_4T1.rds")) 

############################## Historical ######################################

vars_decode.hist <- c("cgm_sex", "cgm_education_type", "cgm_insurance_type", 
                      "cgm_visit_type", "cgm_provider_type", "pump_type",
                      "a1c_collection_method")

for (i in 1:length(vars_decode.hist)){
  dat.hist0[, vars_decode.hist[i]] <- add_labels(dataName = dat.hist0, varName = vars_decode.hist[i])
}

dat.hist0$record_id <- as.numeric(dat.hist0$record_id)

dat_base <- dat.hist0 %>% 
  filter(is.na(redcap_repeat_instrument)) %>%   
  select(record_id, cgm_dob:demographics_complete) %>%
  mutate_at(c("cgm_dob", "cgm_diab_onsetdt", "cgm_start_date"), as.Date) %>%
  mutate(
    Primary.Language = case_when(cgm_language___1 == 1 ~ "English",
                                 cgm_language___2 == 1 | cgm_language___3 == 1 ~ "Non-English"),
    
    cgm_sex = recode(cgm_sex, `Trans-male` = "Male",
                              `Trans-female` = "Female"),
    
    Race.Ethnicity = case_when(cgm_ethnicity %in% 2 ~ "Hispanic",
                               
                               cgm_ethnicity %in% 1 & cgm_race_update___2 == 1 |
                                 cgm_race_update___5 == 1 ~ "Asian or Pacific Islander",
                               cgm_ethnicity %in% 1 & cgm_race_update___3 == 1 ~ "Non-Hispanic Black",
                               cgm_ethnicity %in% 1 & cgm_race_update___4 == 1 ~ "American Indian or Alaska Native",
                               cgm_ethnicity %in% 1 & cgm_race_update___1 == 1 ~ "Non-Hispanic White",
                               cgm_ethnicity %in% c(1, 99) & cgm_race_update___99 == 1 ~ "Other", 
 
                               cgm_race_update___98 == 1 &
                               cgm_race_update___1 == 0 & 
                               cgm_race_update___2 == 0 & 
                               cgm_race_update___3 == 0 & 
                               cgm_race_update___4 == 0 & 
                               cgm_race_update___5 == 0 & 
                               cgm_race_update___99 == 0 |

                               cgm_race_update___1 == 0 & 
                               cgm_race_update___2 == 0 & 
                               cgm_race_update___3 == 0 & 
                               cgm_race_update___4 == 0 & 
                               cgm_race_update___5 == 0 & 
                               cgm_race_update___99 == 0 &
                               cgm_race_update___98 == 0 |
                               
                                 cgm_ethnicity == 99 ~ "Unknown / Declined to state"),
    
    Race.Ethnicity = factor(Race.Ethnicity, levels = c("Non-Hispanic White", 
                                                       "Non-Hispanic Black",
                                                       "Asian or Pacific Islander", 
                                                       "American Indian or Alaska Native",
                                                       "Hispanic", "Other", 
                                                       "Unknown / Declined to state")),
    
    cgm_diab_onsetdt = replace(cgm_diab_onsetdt, record_id == 231, as.POSIXct("2016-09-16"))
    
  ) %>% left_join((baselineds %>% select(record_id, a1c, Insurancetype)), by = "record_id") %>%
  
  mutate(cgm_hba1c_onset = ifelse(is.na(cgm_hba1c_onset), a1c, cgm_hba1c_onset)) %>%
  rename(cgm_insurance_type.new = cgm_insurance_type,
         cgm_insurance_type = Insurancetype)

dat_pump <- dat.hist0 %>% 
  filter(redcap_repeat_instrument %in% "pump_start_form") %>% 
  select(record_id, redcap_repeat_instance, 
         pump_yesno:pump_start_form_complete) %>%
  mutate_at("pump_start_dt", as.Date) %>%
  mutate(pump_type = case_when(pump_type %in% c("Medtronic, 630G, 530G",
                                                "Medtronic 670G no HCL", 
                                                "Omnipod Eros/Classic", 
                                                "Omnipod Dash", "Omnipod 5", 
                                                "Tandem with no Sensor",
                                                "Other (older pumps, Animas, etc.)") ~ "Pump.OpenLoop",
                               pump_type %in% c("Tandem with Basal IQ") ~ "Pump.PLGS",
                               pump_type %in% c("Medtronic 670G HCL") ~ "Pump.HybridCL", 
                               pump_type %in% c("Tandem with Control IQ",
                                                "Loop/Open APS") ~ "Pump.AdvHybridCL")) %>%
  filter(!(record_id %in% 231 & redcap_repeat_instance %in% 1))

dat_pump.w <- dat_pump %>% filter(pump_yesno %in% 1) %>% 
  arrange(record_id, pump_type, pump_start_dt) %>%
  distinct(record_id, pump_type, .keep_all = T) %>% 
  select(record_id, pump_yesno, pump_type, pump_start_dt) %>%
  spread(key = 'pump_type',   
         value = 'pump_start_dt') %>% 
  arrange(record_id) %>%
  rename(Pump.Unknown = `<NA>`) %>% 
  rename_with(.fn = ~ paste0(.x, ".dt"), starts_with("Pump."))

dat_base.hist0 <- dat_base %>% left_join(dat_pump.w, by = "record_id") %>%
  
  mutate(Pump.OpenLoop = ifelse(!is.na(Pump.OpenLoop.dt), 1, 0), 
         Pump.PLGS = ifelse(!is.na(Pump.PLGS.dt), 1, 0), 
         Pump.HybridCL = ifelse(!is.na(Pump.HybridCL.dt), 1, 0), 
         Pump.AdvHybridCL = ifelse(!is.na(Pump.AdvHybridCL.dt), 1, 0),
         
         pump_start_dt = pmin(Pump.OpenLoop.dt, Pump.PLGS.dt, Pump.HybridCL.dt, Pump.AdvHybridCL.dt,
                              Pump.Unknown.dt, 
                              na.rm = T),
         cgm_height = NA,
         cgm_weight = NA,
         
         eos_date = NA,
         cohort = "hist")

dat_hist <- dat.hist0 %>% filter(redcap_repeat_instrument %in% "hba1c_tracking") %>%
  mutate_at(c("a1c_collection_date", "a1c_lab_result_date"), as.Date)

############################## Pilot 4T ########################################

vars_decode.p4T <- c("sex", "education_type", "insurance_type", "visit_type", 
  "hba1c_method", "provider_type", "pump_type", "cgm_pump_type", "pump_type_885620",     
  "a1c_collection_method") 

for (i in 1:length(vars_decode.p4T)){
  dat.p4T0[, vars_decode.p4T[i]] <- add_labels(dataName = dat.p4T0, varName = vars_decode.p4T[i])
}

dat_p4T <- subset(dat.p4T0, !(record_id %in% ID_exclude.p4T)) %>% 
  mutate(record_id = as.numeric(record_id))

dat_cgm.r <- dat_p4T %>% filter(redcap_event_name == "remote_monitoring_arm_1b") %>%
  select(record_id, date_contact:remote_monitoring_initiation_complete) 

ID_remote <- unique(dat_cgm.r$record_id) 

dat_base <- dat_p4T %>% filter(redcap_event_name == "baseline_data_arm_1" & 
                                 is.na(redcap_repeat_instrument)) %>%
  select(record_id:baseline_information_complete) %>%
  mutate_at(c("height_baseline", "weight_baseline", "bmi_baseline"), as.numeric) %>%
  mutate_at(c("dob", "diab_onset", "cgm_start_date"), as.Date) %>%
  mutate(Primary.Language = case_when(language___3 == 1 ~ "Other",
                                      language___2 == 1 ~ "Spanish",
                                      language___1 == 1 ~ "English")
  ) %>%
  
  left_join(dat_p4T %>% filter(redcap_event_name == "end_of_study_arm_1") %>%
              select(record_id, ltfu:no_longer_following_complete), by = "record_id")

dat_base$race_updated___9 <- with(dat_base, replace(race_updated___9, 
                                                    rowSums(data.frame(race_updated___0,
                                                                       race_updated___1,
                                                                       race_updated___2, 
                                                                       race_updated___3,
                                                                       race_updated___4,
                                                                       race_updated___98)) > 0, 0))
dat_base %<>%
  mutate(Hispanic = ifelse(ethnicity == 99, NA, ethnicity),
         Race.Ethnicity = case_when(Hispanic %in% 1 ~ "Hispanic",
                                    
                                    Hispanic %in% 0 & race_updated___1 == 1 |
                                      race_updated___4 == 1 ~ "Asian or Pacific Islander", 
                                    Hispanic %in% 0 & race_updated___2 == 1 ~ "Non-Hispanic Black",
                                    Hispanic %in% 0 & race_updated___3 == 1 ~ "American Indian or Alaska Native",
                                    Hispanic %in% 0 & race_updated___0 == 1 ~ "Non-Hispanic White",
                                    
                                    Hispanic %in% c(0, NA) & race_updated___98 == 1 ~ "Other", 
                                    
                                    race_updated___9 == 1 &
                                      race_updated___0 == 0 & 
                                      race_updated___1 == 0 & 
                                      race_updated___2 == 0 & 
                                      race_updated___3 == 0 & 
                                      race_updated___4 == 0 & 
                                      race_updated___98 == 0 |
                                      
                                      race_updated___0 == 0 & 
                                      race_updated___1 == 0 & 
                                      race_updated___2 == 0 & 
                                      race_updated___3 == 0 & 
                                      race_updated___4 == 0 & 
                                      race_updated___98 == 0 &
                                      race_updated___9 == 0 ~ "Unknown / Declined to state")) %>%
  
  mutate(Race.Ethnicity = factor(factor(Race.Ethnicity, 
                                        levels = c("Non-Hispanic White",
                                                   "Non-Hispanic Black",
                                                   "Asian or Pacific Islander",
                                                   "American Indian or Alaska Native",
                                                   "Hispanic", 
                                                   "Other",
                                                   "Unknown / Declined to state"))))

dat_pump <- dat_p4T %>% filter(redcap_event_name == "baseline_data_arm_1" &
                                 redcap_repeat_instrument %in% "pump_start_form") %>%
  select(record_id, pump_yesno:pump_start_form_complete) %>%
  rename(pump_type = pump_type_885620) %>%
  mutate(pump_start_dt = as.Date(pump_start_dt),
         pump_type = case_when(pump_type %in% c("Medtronic 670G no HCL", 
                                                "Omnipod Eros/Classic", 
                                                "Omnipod Dash", 
                                                "Tandem with no Sensor") ~ "Pump.OpenLoop",
                               pump_type %in% c("Medtronic, 630G", 
                                                "Tandem with Basal IQ") ~ "Pump.PLGS",
                               pump_type %in% c("Medtronic 670G HCL", 
                                                "Omnipod 5") ~ "Pump.HybridCL", 
                               pump_type %in% c("Omnipod Horizon",
                                                "Tandem with Control IQ",
                                                "Loop/Open APS") ~ "Pump.AdvHybridCL"))

dat_pump.w <- dat_pump %>% filter(pump_yesno %in% 1) %>% 
  arrange(record_id, pump_type, pump_start_dt) %>%
  distinct(record_id, pump_type, .keep_all = T) %>% 
  select(record_id, pump_yesno, pump_type, pump_start_dt) %>%
  spread(key = 'pump_type',  
         value = 'pump_start_dt') %>% 
  arrange(record_id) %>%
  rename_with(.fn = ~ paste0(.x, ".dt"), starts_with("Pump."))

dat_base.p4T0 <- dat_base %>% left_join(dat_pump.w, by = "record_id") %>%
  
  mutate_at("end_study_date", as.Date) %>%
  
  mutate(Pump.OpenLoop = ifelse(!is.na(Pump.OpenLoop.dt), 1, 0), 
         Pump.PLGS = ifelse(!is.na(Pump.PLGS.dt), 1, 0), 
         Pump.HybridCL = ifelse(!is.na(Pump.HybridCL.dt), 1, 0),  
         Pump.AdvHybridCL = ifelse(!is.na(Pump.AdvHybridCL.dt), 1, 0),
         pump_start_dt = pmin(Pump.OpenLoop.dt, Pump.PLGS.dt, Pump.AdvHybridCL.dt,
                              na.rm = T),
         cohort = "p.4T")

############################## 4T Study 1 ######################################

vars_decode.4T1 <- c("sex", "cgm_sex", "cgm_education_type", "cgm_insurance_type", 
                     "cgm_visit_type", "cgm_provider_type", "cgm_pump_type",
                     "pump_type", "dsdecod", "a1c_collection_method")

for (i in 1:length(vars_decode.4T1)){
  dat.4T0[, vars_decode.4T1[i]] <- add_labels(dataName = dat.4T0, varName = vars_decode.4T1[i])
}

dat_4T1 <- subset(dat.4T0, !(record_id %in% ID_exclude.4T1)) %>%
  mutate(record_id = as.numeric(record_id)) %>%
  filter(record_id %in% dat.4T0$record_id[dat.4T0$consent %in% 1])

dat_screen <- dat_4T1 %>% filter(redcap_event_name == "screening_arm_1" & is.na(redcap_repeat_instrument)) %>%
  select(record_id, center:screening_complete)

dat_cgm <- dat_4T1 %>% filter(redcap_event_name == "baseline_arm_1") %>%
  select(record_id, cgm_visitdt:cgm_visit_form_complete) 

dat_eos <- dat_4T1 %>% filter(redcap_event_name == "end_of_study_arm_1" & is.na(redcap_repeat_instrument)) %>%
  select(record_id, eos_demos_resend:end_of_study_complete)

dat_pump <- dat_4T1 %>% 
  filter(redcap_event_name == "end_of_study_arm_1" &
           redcap_repeat_instrument %in% "pump_start_form") %>% 
  select(record_id, pump_yesno:pump_start_form_complete) %>% 
  mutate(pump_start_dt = as.Date(pump_start_dt),
         pump_type = as.character(pump_type), 
         pump_type = case_when(pump_type %in% c("Medtronic 670G no HCL", 
                                                "Omnipod Eros/Classic", 
                                                "Omnipod Dash", 
                                                "Tandem with no Sensor") ~ "Pump.OpenLoop",
                               pump_type %in% c("Medtronic, 630G", 
                                                "Tandem with Basal IQ") ~ "Pump.PLGS",
                               pump_type %in% c("Medtronic 670G HCL", 
                                                "Omnipod 5 (HCL)") ~ "Pump.HybridCL", 
                               pump_type %in% c("Omnipod Horizon",
                                                "Tandem with Control IQ",
                                                "Loop/Open APS") ~ "Pump.AdvHybridCL"))

dat_pump.w <- dat_pump %>% filter(pump_yesno %in% 1) %>%
  arrange(record_id, pump_type, pump_start_dt) %>%
  distinct(record_id, pump_type, .keep_all = T) %>% 
  select(record_id, pump_yesno, pump_type, pump_start_dt) %>%
  spread(key = 'pump_type',  
         value = 'pump_start_dt') %>% 
  arrange(record_id) %>%
  rename_with(.fn = ~ paste0(.x, ".dt"), starts_with("Pump."))

dat_base <- list(dat_screen, dat_cgm %>% select(record_id, cgm_dob:cgm_provider_type), 
                 dat_pump.w, dat_eos) %>% reduce(full_join, by = 'record_id')

dat_base %<>%
    mutate_at(c("cgm_dob", "cgm_diab_onsetdt", "cgm_start_date",
                "dscomdt", "dsdiscdt"), as.Date) %>%
  mutate(
    Pump.OpenLoop = ifelse(!is.na(Pump.OpenLoop.dt), 1, 0), 
    Pump.PLGS = ifelse(!is.na(Pump.PLGS.dt), 1, 0), 
    Pump.HybridCL = ifelse(!is.na(Pump.HybridCL.dt), 1, 0), 
    Pump.AdvHybridCL = ifelse(!is.na(Pump.AdvHybridCL.dt), 1, 0),

    pump_start_dt = pmin(Pump.OpenLoop.dt, Pump.PLGS.dt, Pump.HybridCL.dt, Pump.AdvHybridCL.dt,
                         na.rm = T),
    
    cgm_language = case_when(cgm_language___3 == 1 ~ "Other",
                             cgm_language___2 == 1 ~ "Spanish",
                             cgm_language___1 == 1 ~ "English"),
    
    cgm_sex = ifelse(is.na(cgm_sex), as.character(sex), as.character(cgm_sex)),
    cgm_sex = recode(cgm_sex, `Trans-male` = "Male",
                     `Trans-female` = "Female"),
    
    disc.within.12mon = as.integer(dsdiscdt - cgm_diab_onsetdt)/30 <= mon.fu,
    eos_date = as.Date(pmin(dscomdt, dsdiscdt, na.rm = T)),  
    
    cohort = "4T.1" 
  )

dat_base$cgm_race_update___98 <- with(dat_base, replace(cgm_race_update___98, 
                                                        rowSums(data.frame(cgm_race_update___1, 
                                                                           cgm_race_update___2,
                                                                           cgm_race_update___3, 
                                                                           cgm_race_update___4,
                                                                           cgm_race_update___5, 
                                                                           cgm_race_update___99)) > 0, 0))

dat_base$cgm_race_update___98 <- with(dat_base, replace(cgm_race_update___98, 
                                                        rowSums(data.frame(cgm_race_update___1, 
                                                                           cgm_race_update___2,
                                                                           cgm_race_update___3, 
                                                                           cgm_race_update___4,
                                                                           cgm_race_update___5, 
                                                                           cgm_race_update___98, 
                                                                           cgm_race_update___99)) == 0, 1))

dat_base$race_multi.other <- with(dat_base, ifelse(rowSums(data.frame(cgm_race_update___1, cgm_race_update___2,
                                                                      cgm_race_update___3, cgm_race_update___4,
                                                                      cgm_race_update___5)) > 1 | cgm_race_update___99 == 1, 1, 0))

dat_base$cgm_race <- with(dat_base, ifelse(dat_base$race_multi.other == 1, 99,
                                           1*cgm_race_update___1 + 2*cgm_race_update___2 + 3*cgm_race_update___3 + 
                                           4*cgm_race_update___4 + 5*cgm_race_update___5 + 
                                           98*cgm_race_update___98))  

dat_base %<>% mutate(cgm_race = recode(cgm_race, `1` = "White or Caucasian",
                                                 `2` = "Asian or Pacific Islander",
                                                 `3` = "Black or African American",
                                                 `4` = "American Indian or Alaska Native",
                                                 `5` = "Asian or Pacific Islander", 
                                                 `99` = "Other",
                                                 `98` = "Unknown / Declined to state"),
                     Race.Ethnicity = case_when(cgm_ethnicity %in% 2 ~ "Hispanic",
                                                T ~ cgm_race),
                     Race.Ethnicity = recode(Race.Ethnicity, `White or Caucasian` = "Non-Hispanic White",
                                                             `Black or African American` = "Non-Hispanic Black"))

get_pros_re <- function(event, raceVar){
  
  dat_pros <- with(dat_4T1, dat_4T1[redcap_event_name == event & 
                                      is.na(redcap_repeat_instrument),  
                                    c(1, which(names(dat_4T1) == "race_child_eng___1"):
                                         which(names(dat_4T1) == "race_child_eng___99") ,
                                         which(names(dat_4T1) == "race_child_esp___1"):
                                         which(names(dat_4T1) == "race_child_esp___99"))])  
  
  dat_pros$race_multi.other <- with(dat_pros, ifelse(rowSums(data.frame(race_child_eng___1, 
                                                                        race_child_eng___3, race_child_eng___4,
                                                                        race_child_eng___5)) > 1 | 
                                                       race_child_eng___99 == 1, 1, 0))
  
  dat_pros$race <- with(dat_pros, ifelse(race_multi.other == 1, 99,
                                         1*race_child_eng___1 + 2*race_child_eng___2 + 3*race_child_eng___3 + 
                                         4*race_child_eng___4 + 5*race_child_eng___5 + 
                                        99*race_child_eng___99))  
  
  dat_pros$race <- with(dat_pros, ifelse(race_child_eng___2 == 1 | race_child_esp___2 == 1, 2, race))

  dat_pros$race <- with(dat_pros, recode(race, `1` = "Non-Hispanic Black",
                                               `2` = "Hispanic",
                                               `3` = "American Indian or Alaska Native",
                                               `4` = "Asian or Pacific Islander",
                                               `5` = "Non-Hispanic White",
                                              `99` = "Other",
                                               `0` = "Unknown / Declined to state"))
        
  out <- dat_pros[,c("record_id", "race")]
  names(out)[2] <- raceVar  
  
  return(out)
  
}

dat_pros.race <- merge(get_pros_re(event = "baseline_arm_1", raceVar = "race.pros.base"),
                       get_pros_re(event = "end_of_study_arm_1", raceVar = "race.pros.eos"), 
                       by = "record_id", all = T)
dat_pros.race$Race.Ethnicity.pros <- with(dat_pros.race, ifelse(race.pros.base == "Unknown / Declined to state", 
                                                                race.pros.eos, 
                                                                race.pros.base))

dat_base <- merge(dat_base, dat_pros.race[,c("record_id", "Race.Ethnicity.pros")], by = "record_id", all.x = T)

dat_base$Race.Ethnicity <- with(dat_base, ifelse(Race.Ethnicity.pros == "Unknown / Declined to state",
                                                 Race.Ethnicity, Race.Ethnicity.pros)) 

ID_update.dob <- c(49, 77, 128, 155)

dat_base %<>%
  
  left_join((dat_4T1 %>% select(record_id, redcap_event_name, dob_child_eng) %>% 
               filter(redcap_event_name == "baseline_arm_1" & record_id %in% ID_update.dob)), by = "record_id") %>%
  
  mutate(
    cgm_dob = case_when(record_id %in% ID_update.dob ~ dob_child_eng,
                        T ~ cgm_dob))

dat_base.4T0 <- dat_base

############################## Combine cohorts #################################

vars_4T1 <- c("record_id", "cgm_sex", "Race.Ethnicity", "cgm_hba1c_onset",
              "cgm_insurance_type", "cgm_education_type",
              
              "cgm_dob",
              "cgm_diab_onsetdt",  
              "cgm_start_date",
              "pump_start_dt",
              "eos_date",  
              
              "cgm_dka_onset", "lowest_ph", "low_bicarb",  
              "cgm_height", "cgm_weight",                  
              "pump_yesno", "Pump.OpenLoop", "Pump.PLGS",
              "Pump.HybridCL", "Pump.AdvHybridCL", 
              
              "Pump.OpenLoop.dt", "Pump.PLGS.dt",
              "Pump.HybridCL.dt", "Pump.AdvHybridCL.dt", 
              
              "cgm_language", "cohort")

dat_base.hist <- dat_base.hist0[,c("record_id", "cgm_sex", "Race.Ethnicity", 
                                   "cgm_hba1c_onset", "cgm_insurance_type", 
                                   "cgm_education_type",
                                   
                                   "cgm_dob", 
                                   "cgm_diab_onsetdt",  
                                   "cgm_start_date",
                                   "pump_start_dt",
                                   "eos_date",  
                                   
                                   "cgm_dka_onset", "lowest_ph", "low_bicarb", 
                                   "cgm_height", "cgm_weight",                 
                                   "pump_yesno", "Pump.OpenLoop", "Pump.PLGS",
                                   "Pump.HybridCL", "Pump.AdvHybridCL",
                                   
                                   "Pump.OpenLoop.dt", "Pump.PLGS.dt",
                                   "Pump.HybridCL.dt", "Pump.AdvHybridCL.dt", 
                                   
                                   "Primary.Language", "cohort")]

dat_base.p4T <- dat_base.p4T0[,c("record_id", "sex", "Race.Ethnicity", "hba1c_onset",
                                 "insurance_type", "education_type",
                                 
                                 "dob", 
                                 "diab_onset",  
                                 "cgm_start_date",
                                 "pump_start_dt",
                                 "end_study_date",  
                                 
                                 "dka_onsest", "lowest_ph", "low_bicarb",  
                                 "height_baseline", "weight_baseline",     
                                 "pump_yesno", "Pump.OpenLoop", "Pump.PLGS",
                                 "Pump.HybridCL", "Pump.AdvHybridCL",
                                 
                                 "Pump.OpenLoop.dt", "Pump.PLGS.dt",
                                 "Pump.HybridCL.dt", "Pump.AdvHybridCL.dt", 
                                 
                                 "Primary.Language", "cohort")]

dat_base.4T1 <- dat_base.4T0[, vars_4T1]

names(dat_base.hist) <- names(dat_base.p4T) <- names(dat_base.4T1) <- vars_4T1

dat_base.all <- rbind(dat_base.hist, dat_base.p4T, dat_base.4T1) %>%
  
  mutate(
    record_id = paste0(cohort, "-", sprintf("%03d", record_id)),
    
    age.at.onset = as.numeric(cgm_diab_onsetdt - cgm_dob)/365.25,
    age.at.onset.cat = cut(age.at.onset, breaks=c(min(age.at.onset, na.rm = T), 6, 12.991, 
                                                  max(age.at.onset, na.rm = T)), 
                           include.lowest = T),
    
    CGM.initiated = ifelse(!is.na(cgm_start_date), 1, 0),
    daystoCGM = as.integer(cgm_start_date - cgm_diab_onsetdt),
    CGM30days = ifelse(daystoCGM <= 30, 1, 0),
    daystoPump = as.integer(pump_start_dt - cgm_diab_onsetdt),
    daystoCGM = replace(daystoCGM, daystoCGM > mon.fu*30, NA),
    CGM.initiated = replace(CGM.initiated, is.na(daystoCGM), 0),
    Pump.OpenLoop = ifelse(Pump.OpenLoop & 
                             as.integer(Pump.OpenLoop.dt - cgm_diab_onsetdt) <= mon.fu*30, 1, 0), 
    
    Pump.PLGS = ifelse(Pump.PLGS & as.integer(Pump.PLGS.dt - cgm_diab_onsetdt) <= mon.fu*30, 1, 0), 
    Pump.HybridCL = ifelse(Pump.HybridCL & as.integer(Pump.HybridCL.dt - cgm_diab_onsetdt) <= mon.fu*30, 1, 0), 
    Pump.AdvHybridCL = ifelse(Pump.AdvHybridCL & as.integer(Pump.AdvHybridCL.dt - cgm_diab_onsetdt) <= mon.fu*30, 1, 0), 
    daystoPump = replace(daystoPump, daystoPump > mon.fu*30, NA),
    pump_yesno = replace(pump_yesno, is.na(daystoPump), 0), 
    
    remote.monitor = ifelse(cohort == "p.4T" & as.numeric(substr(record_id, 6, nchar(record_id))) %in% ID_remote |
                            cohort == "4T.1", 1, 0),
    
    Hispanic = case_when(Race.Ethnicity == "Hispanic" ~ 1,
                         Race.Ethnicity != "Hispanic" & Race.Ethnicity != "Unknown / Declined to state" ~ 0),
    
    Public = ifelse(cgm_insurance_type == "Public", 1, 0), 
    
    cgm_language = case_when(cgm_language %in% c("Spanish", "Other") ~ "Non-English",
                             T ~ cgm_language),
    cgm_sex = factor(cgm_sex),
    
    Race.Ethnicity = factor(Race.Ethnicity, levels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic",
                                                       "Asian or Pacific Islander", "American Indian or Alaska Native",
                                                       "Other", "Unknown / Declined to state")),
    
    cgm_insurance_type = factor(cgm_insurance_type, levels = c("Private", "Public", "Both", "No Insurance")),
    
    cohort = factor(cohort, levels = c("hist", "p.4T", "4T.1")),
    eos_date = as.Date(ifelse(!is.na(eos_date), eos_date, as.Date(cutdate)), origin = "1970-01-01"),
    CGM30days = ifelse(!is.na(CGM30days), CGM30days, 0),
    Pump.OpenLoop = ifelse(!is.na(Pump.OpenLoop), Pump.OpenLoop, 0),
    Pump.PLGS = ifelse(!is.na(Pump.PLGS), Pump.PLGS, 0),
    Pump.HybridCL = ifelse(!is.na(Pump.HybridCL), Pump.HybridCL, 0),
    Pump.AdvHybridCL = ifelse(!is.na(Pump.AdvHybridCL), Pump.AdvHybridCL, 0),
    AdvHybridCL.date1 = as.integer(Pump.AdvHybridCL.dt - cgm_diab_onsetdt),
    
  ) %>% 
  arrange(cohort, record_id)

############################## Combine outcome #################################

dat_a1c.hist <- process_a1c(cohort = "hist", event = "none")[[2]] %>% 
  select(record_id, a1c_result, a1c_collection_date, a1c_collection_method) %>%
  arrange(record_id, a1c_collection_date) %>%
  mutate(record_id = paste0("hist-", sprintf("%03d", record_id))) %>% 
  rename(cgm_hba1c = a1c_result) 

dat_a1c.p4T <- process_a1c(cohort = "p.4T", event = "baseline_data_arm_1")[[2]] %>% 
  select(record_id, a1c_result, a1c_collection_date, a1c_collection_method) %>%
  arrange(record_id, a1c_collection_date) %>%
  mutate(record_id = paste0("p.4T-", sprintf("%03d", record_id))) %>% 
  rename(cgm_hba1c = a1c_result) 

dat_a1c.4T1 <- process_a1c(cohort = "4T.1", event = "extension_phase_arm_1")[[2]] %>% 
  select(record_id, a1c_result, a1c_collection_date, a1c_collection_method) %>%
  arrange(record_id, a1c_collection_date) %>%
  mutate(record_id = paste0("4T.1-", sprintf("%03d", record_id))) %>% 
  rename(cgm_hba1c = a1c_result) 

dat_a1c.epic <- rbind(dat_a1c.hist, dat_a1c.p4T, dat_a1c.4T1) 

ID_addA1c <- setdiff(dat_base.all$record_id, 
                     (dat_a1c.epic %>% left_join(dat_base.all[,c("record_id", "cgm_diab_onsetdt")]) %>%             
                        filter(cgm_diab_onsetdt == a1c_collection_date) %>%
                        pull(record_id))) 

dat_a1c.cgm <- dat_base.all %>% filter(record_id %in% ID_addA1c) %>% 
  select(record_id, cgm_hba1c_onset, cgm_diab_onsetdt) %>%
  mutate(a1c_collection_method = "Clinical Point-of-care (POC)") %>% 
  rename(cgm_hba1c = cgm_hba1c_onset,
         a1c_collection_date = cgm_diab_onsetdt)

dat_a1c.all <- rbind(dat_a1c.epic, dat_a1c.cgm) %>% 
  arrange(record_id, a1c_collection_date) 

dat_base.all %<>% left_join(dat_a1c.all[, c("record_id", "cgm_hba1c", "a1c_collection_date")], by = "record_id") %>%
  filter(cgm_diab_onsetdt == a1c_collection_date | 
           is.na(cgm_diab_onsetdt)) %>%  
  mutate(diff_a1c0 = cgm_hba1c - cgm_hba1c_onset,
         cgm_hba1c_onset = case_when(!is.na(diff_a1c0) & diff_a1c0 != 0 ~ cgm_hba1c,
                                     T ~ cgm_hba1c_onset)) %>%
  select(-cgm_hba1c, -a1c_collection_date) %>%
  left_join(dat_wear, by = "record_id")


dat_long.all <- merge(dat_base.all, dat_a1c.all, by = "record_id", all = T) %>% 
  
  filter(a1c_collection_date <= pmin(as.Date(cutdate), 
                                     as.Date(eos_date, origin = "1970-01-01"), na.rm = T)) %>%
  mutate(
    a1c_collection_method = factor(a1c_collection_method, 
                                   levels = c("Clinical Point-of-care (POC)",
                                              "At Home HbA1c Kit", "Lab")),
    
    studyday = as.numeric(a1c_collection_date - cgm_diab_onsetdt) + 1,
    studymon = studyday/30,
    
    OpenLoopEver = ifelse(record_id %in% (dat_base.all %>% 
                                            filter(Pump.OpenLoop == 1) %>% pull(record_id)), 1, 0)) %>%
  
  arrange(record_id, studyday, desc(cgm_hba1c)) %>%
  distinct(record_id, studyday, .keep_all = T) %>%
  
  subset(dat_long.all, studyday >= 0) %>%
  mutate(month = cut(studymon, breaks = c(0, 1, 4.5, 7.5, 10.5, 13.5), 
                     labels = c(1, 3, 6, 9, 12), include.lowest = T),
         target.lt7 = ifelse(cgm_hba1c < 7, 1, 0),
         target.lt7.5 = ifelse(cgm_hba1c < 7.5, 1, 0))
