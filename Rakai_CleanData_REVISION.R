#################################################################################################
# Rakai data analysis for EMOD inputs - Data cleaning
# Nov 29, 2021
# Adam Akullian
#################################################################################################
library("ggplot2")
library("readstata13")
library("data.table")
library("dplyr")
library("haven")
library("tidyr")
library("tidyverse")
library("mgcv")
library("zoo")
library("metR")
library("tableone")
library("mgcv")

options(scipen = 999)

################################################################################################
#Bring in data
################################################################################################

setwd("RakaiData_23Oct2023")

verif_1 <- read_dta("verif_4.dta")
verif_1$research_id <- verif_1$study_id
length(unique(verif_1$research_id)) #41392 unique individuals in the cohort

subset(verif_1) %>%
  group_by(round) %>%
  summarise(count=length(unique(research_id)))

table(verif_1$resident)

# #restrict to only residents
# verif_1 <- subset(verif_1, resident==1)
# length(unique(verif_1$research_id)) #41388 unique individuals in the cohort (residents) hence this dataset is only residents
# table(verif_1$resident)
# class(verif_1$round)
# names(verif_1)

hivstatus_vlcopies_1 <- read_dta("hivstatus_vlcopies_4.dta")

nrow(hivstatus_vlcopies_1)
hivstatus_vlcopies_1$research_id <- hivstatus_vlcopies_1$study_id
hivstatus_vlcopies_1$hivdate <- hivstatus_vlcopies_1$intdate
length(unique(hivstatus_vlcopies_1$research_id)) #40,876 unique indiv before round 19

subset(hivstatus_vlcopies_1) %>%
  group_by(round) %>%
  summarise(count=length(unique(research_id)))

hivstatus_vlcopies_1_prev <- hivstatus_vlcopies_1 %>%
  group_by(round, finalhiv) %>%
  summarise(count=length(unique(research_id)))

spread(hivstatus_vlcopies_1_prev, finalhiv, count)

#length(unique(hivstatus_vlcopies_1$research_id[hivstatus_vlcopies_1$round=="R019"])) #37,396 unique indiv before round 19

ids <- unique(verif_1$research_id)
hivstatus_vlcopies_1_verif <- subset(hivstatus_vlcopies_1, research_id %in% ids) 

hivincidence_1 <- read_dta("hivincidence_4.dta")
quest_analysisvars_1 <- read_dta("quest_analysisvars_4.dta")
quest_analysisvars_1$research_id <- quest_analysisvars_1$study_id
quest_analysisvars_1 <- quest_analysisvars_1[order(quest_analysisvars_1$study_id, quest_analysisvars_1$round),] #sort by round
quest_1 <- read_dta("quest_4.dta")
quest_1$research_id <- quest_1$study_id

names(quest_1)

names(quest_1)
head(quest_1$rltnage1)
head(quest_1$rltn1)

tab.indiv.round <- hivstatus_vlcopies_1 %>%
  group_by(round) %>%
  summarise(count=length(unique(research_id)))

tab.indiv.round.hiv <- hivstatus_vlcopies_1 %>%
  group_by(round, finalhiv) %>%
  summarise(count=length(unique(research_id)))
spread(tab.indiv.round.hiv, finalhiv, count)
################################################################################################
#clean data
###############################################################################################
#code HIV status
table(hivstatus_vlcopies_1$finalhiv)
table(hivstatus_vlcopies_1$finalhiv, hivstatus_vlcopies_1$round)
hivstatus_vlcopies_1$hivstatus=NA
hivstatus_vlcopies_1$hivstatus[hivstatus_vlcopies_1$finalhiv=="N"]=0
hivstatus_vlcopies_1$hivstatus[hivstatus_vlcopies_1$finalhiv=="P" | hivstatus_vlcopies_1$finalhiv=="S"]=1
table(hivstatus_vlcopies_1$hivstatus)
table(hivstatus_vlcopies_1$hivstatus, hivstatus_vlcopies_1$round)

#code visit dates
table(hivstatus_vlcopies_1$study_id)
table(hivstatus_vlcopies_1$hivdate)
hivstatus_vlcopies_1$hivdate <- as.Date(hivstatus_vlcopies_1$hivdate, format='%Y-%m-%d')
summary(hivstatus_vlcopies_1$hivdate)
hivstatus_vlcopies_1$round <- as.numeric(gsub(".*?([0-9]+).*", "\\1", hivstatus_vlcopies_1$round)) #get numeric round     
hivstatus_vlcopies_1 <- hivstatus_vlcopies_1[order(hivstatus_vlcopies_1$study_id, hivstatus_vlcopies_1$round),] #sort by round
hivstatus_vlcopies_1$year <- format(hivstatus_vlcopies_1$hivdate, format="%Y")
table(hivstatus_vlcopies_1$year)
hivstatus_vlcopies_1$hivdate[hivstatus_vlcopies_1$year < 1999] <- NA
hivstatus_vlcopies_1$year[hivstatus_vlcopies_1$year < 1999] <- NA
table(hivstatus_vlcopies_1$year)
table(hivstatus_vlcopies_1$round,hivstatus_vlcopies_1$year)
hivstatus_vlcopies_1$research_id <- hivstatus_vlcopies_1$study_id

#Get mean hivdate by round and add in missing rounds between present rounds to be able to estimate incidence at missing rounds
hivdate_round <- hivstatus_vlcopies_1 %>% group_by(round) %>% summarise(hivdate=mean(hivdate, na.rm=T))

#fill in missing - gaps - in HIV status (those who serially tested negative or positive but had missing test status in the middle) 
hivstatus_vlcopies_1$hivstatus_impute_up=hivstatus_vlcopies_1$hivstatus
hivstatus_vlcopies_1$hivstatus_impute_down=hivstatus_vlcopies_1$hivstatus
hiv <- hivstatus_vlcopies_1 %>%
  group_by(research_id) %>%
  complete(round = min(round):max(round)) %>% #add in missing rounds between early and late rounds
  fill(hivstatus_impute_down, .direction = "down") %>% #impute missing negative/positive values between early late present values
  fill(hivstatus_impute_up, .direction = "up") %>%
  mutate(hivstatus_imputed=ifelse(hivstatus_impute_down==hivstatus_impute_up, hivstatus_impute_down, 
                                  ifelse(hivstatus_impute_down==0 & hivstatus_impute_up==1, NA, NA))) %>% #do not impute values between a negative and positive
  dplyr::select(-c(hivstatus_impute_down, hivstatus_impute_up)) %>%
  left_join(y=hivdate_round, by="round") %>% #assign mean value of hivdate for each round for missing rounds
  mutate(hivdate=ifelse(is.na(hivdate.x)==F, hivdate.x, hivdate.y))  %>%
  mutate(hivdate=as.Date(hivdate), year=ifelse(is.na(year)==F, year, format(hivdate, format='%Y'))) %>%
  dplyr::select(-c(hivdate.x, hivdate.y))

cohort_gaps <- hiv %>%
  group_by(research_id) %>%
  mutate(round_diff = round - dplyr::lag(round)) %>%
  ungroup()

#any jump larger than 1 means a missing round
gaps <- cohort_gaps %>%
  filter(round_diff > 1) %>%
  select(research_id, round, round_diff)

# how many individuals have gaps?
length(unique(gaps$research_id))

################################################################################################
#Create the seroconverter cohort for incidence analysis
#Also fill in missing co-variates (gaps) in variables (e.g., HIV status, ever sex, sex, birthdate)
    #where an individual's values is 0 at one time step and 0 at another with a gap
    #or is 1 at one time step and 1 at a subsequent with a gap in between
################################################################################################

# #Create dataset to indicate first neg, last neg, and first positive to generate the cohort - this dataset only includes those with at least one test
firstpos <- hiv %>%
  group_by(research_id) %>%
  arrange(research_id, round) %>%
  slice(match(1, hivstatus_imputed)) %>%
  select(research_id, round, hivdate) %>%
  mutate(first_pos=round) %>% select(-hivdate,-round)
lastneg <- hiv %>%
  group_by(research_id) %>%
  arrange(research_id, desc(round)) %>%
  slice(match(0, hivstatus_imputed)) %>%
  select(research_id, round, hivdate) %>%
  mutate(last_neg=round) %>% select(-hivdate,-round) %>%
  full_join(firstpos, by="research_id")
join_df <- hiv %>%
  group_by(research_id) %>%
  arrange(research_id, round) %>%
  slice(match(0, hivstatus_imputed)) %>%
  select(research_id, round, hivdate) %>%
  mutate(first_neg=round) %>% select(-hivdate,-round) %>%
  full_join(lastneg, by="research_id")
nrow(join_df) #40,876 unique id's
join_df <- join_df %>% mutate(cohortclass=ifelse(first_neg!=last_neg & is.na(first_neg)==F & is.na(first_pos)==T, "Serially negative",
                                                 ifelse(first_neg==last_neg & is.na(first_neg)==F & is.na(first_pos)==T, "One test negative",
                                                        ifelse(is.na(first_neg)==T & is.na(first_pos)==F, "serially positive",
                                                               ifelse(is.na(first_neg)==F & is.na(first_pos)==F, "Seroconversion","No Testing")))))
table(join_df$cohortclass)

#Impute a random round for seroconversion that is between the last negative and first positive, including the first positive.
set.seed(1984)

bad_bounds <- join_df %>%
  filter(cohortclass == "Seroconversion") %>%
  mutate(last_neg = as.integer(last_neg),
         first_pos = as.integer(first_pos)) %>%
  filter(!is.na(last_neg), !is.na(first_pos),
         first_pos <= last_neg) %>%
  select(research_id, last_neg, first_pos)

bad_bounds
n_distinct(bad_bounds$research_id)  # how many unique IDs

seroconverters <- join_df %>%
  filter(!research_id %in% bad_bounds$research_id) %>% #remove bad bounds from join_df
  filter(cohortclass == "Seroconversion") %>%
  select(research_id, last_neg, first_pos, cohortclass) 

seroconverters <- seroconverters %>%
  rowwise() %>%
  mutate(
    seroconv = if (is.na(last_neg) || is.na(first_pos) || first_pos <= last_neg) {
      NA_integer_
    } else if (first_pos == last_neg + 1) {
      first_pos
    } else {
      sample(seq.int(last_neg + 1L, first_pos), 1L)
    }
  ) %>%
  ungroup()

with(seroconverters, all(seroconv > last_neg | is.na(seroconv)))
with(seroconverters, all(seroconv <= first_pos | is.na(seroconv)))

join_df <- join_df %>%
  filter(!research_id %in% bad_bounds$research_id) %>% #remove bad bounds from join_df
  select(-matches("^seroconv$")) %>%                # drop any old seroconv col if present
  left_join(seroconverters %>% select(research_id, seroconv) %>%
              rename(seroconv_imp = seroconv),
            by = "research_id")

hiv_m_seroconv <- hiv %>%
  # bring in cohort labels + imputed seroconv round
  left_join(join_df, by = "research_id") %>%
  arrange(research_id, hivdate, round) %>%          # ensure within-ID order
  group_by(research_id) %>%
  mutate(
    # event assignment that matches your original intent:
    # - Serially negative: use imputed status (0/1 per row)
    # - Seroconversion: 1 at the imputed event round; 0 before; NA after
    hivinc = ifelse(
      cohortclass == "Serially negative", as.integer(hivstatus_imputed),
      ifelse(cohortclass == "Seroconversion" & round == seroconv_imp, 1L,
             ifelse(cohortclass == "Seroconversion" & round <  seroconv_imp, 0L, NA_integer_))
    ),
    # person-time assigned to the later visit in each interval
    py = ifelse(hivinc %in% c(0L,1L),
                as.numeric(hivdate - dplyr::lag(hivdate)) / 365.25,
                NA_real_)
  ) %>%
  ungroup()

n_events <- dplyr::n_distinct(hiv_m_seroconv$research_id[hiv_m_seroconv$hivinc == 1])
n_events

# incident count (should match 1,353 cohort-class seroconverters)
n_events <- dplyr::n_distinct(hiv_m_seroconv$research_id[hiv_m_seroconv$hivinc == 1])
n_events

# py should be present on event rows (unless previous date truly missing)
sum(is.na(hiv_m_seroconv$py[hiv_m_seroconv$hivinc == 1]))

# optional: overall incidence estimate
events <- sum(hiv_m_seroconv$hivinc == 1, na.rm = TRUE)
total_py <- sum(hiv_m_seroconv$py, na.rm = TRUE)
100 * events / total_py

################################################################################################
#merge in age, gender (and other demographics) to HIV data
################################################################################################
quest_1_demog <- quest_1[, c('research_id', 'round', 'comm_num','arvmed','cuarvmed','circum', 'evermarr','currmarr','eversex','occup1','occup2','ag1stsex','sexp1yr','educate', 'educyrs')]
verif_1_demog <- verif_1[, c('research_id', 'round', 'sex', 'ageyrs','birthdate','resident','mobility')]
verif_1_demog$round <- as.numeric(gsub(".*?([0-9]+).*", "\\1", verif_1_demog$round)) #get numeric round    
quest_1_demog$round <- as.numeric(gsub(".*?([0-9]+).*", "\\1", quest_1_demog$round)) #get numeric round    

ids <- unique(verif_1_demog$research_id)
subset(hiv_m_seroconv, research_id %in% ids)

hiv_m_d <- merge(hiv_m_seroconv, verif_1_demog, by = c("research_id", "round"), all.x=T)
nrow(hiv_m_d)
hiv_m_d <- merge(hiv_m_d, quest_1_demog, by = c("research_id", "round"), all.x=T)
#hiv_m_d <- merge(hiv_m_d, quest_analysisvars_1[,c(1,2,5)], by = c("research_id", "round"), all.x=T)
hiv_m_d <- hiv_m_d[order(hiv_m_d$research_id, hiv_m_d$round),] #sort by round
hiv_m_d$resident_imputed <- hiv_m_d$resident
table(hiv_m_d$resident) #data set only includes those who are residents

hiv_m_d <- hiv_m_d[order(hiv_m_d$research_id, hiv_m_d$round),] #sort by round
table(hiv_m_d$hivinc) #1,352 seroconversions
table(hiv_m_d$hivinc, hiv_m_d$resident)

################################################################################################
#Impute missing circumcision values, ever sex, sex, birthdate values (circumcision and eversex are only imputed between present values)
#Clean and recode other variables
################################################################################################

#circumcision status - 1 if circum==1 and 0 otherwise if circum==2
hiv_m_d$circ_imputed=NA
hiv_m_d$circ_imputed[hiv_m_d$circum==2]=0
hiv_m_d$circ_imputed[hiv_m_d$circum==1]=1
table(hiv_m_d$circ_imputed[hiv_m_d$sex=="M"])
table(hiv_m_d$circ_imputed[hiv_m_d$sex=="F"])

#fill in missing - gaps - in circumcision status (those who have a gap between two 0's or two 1's) 
hiv_m_d$circ_impute_up=hiv_m_d$circ_imputed
hiv_m_d$circ_impute_down=hiv_m_d$circ_imputed

hiv_m_d <- hiv_m_d %>%
  group_by(research_id) %>%
  fill(circ_impute_down, .direction = "down") %>% #impute missing negative/positive values between early late present values
  fill(circ_impute_up, .direction = "up") %>%
  mutate(circ_imputed=ifelse(circ_impute_down==circ_impute_up, circ_impute_down, 
                                  ifelse(circ_impute_down==0 & circ_impute_up==1, NA, NA))) %>% #do not impute values between a negative and positive
  dplyr::select(-c(circ_impute_down, circ_impute_up))

hiv_m_d <- hiv_m_d %>%
  group_by(research_id) %>%
  mutate(circ_imputed_lag = lag(circ_imputed))
ids <- hiv_m_d$research_id[hiv_m_d$circ_imputed_lag==1 & hiv_m_d$circ_imputed==0]
length(unique(ids)) #249 individuals have circumcision values that are out of order, we will mark all circ values here as NA
hiv_m_d$circ_imputed_fix = hiv_m_d$circ_imputed
hiv_m_d$circ_imputed_fix[hiv_m_d$research_id %in% ids] <- NA 
table(hiv_m_d$circ_imputed_fix[hiv_m_d$sex=="M"])
table(hiv_m_d$circ_imputed_fix[hiv_m_d$sex=="F"])

#fill in missing - gaps - in ever had sex (those who have a gap between two 0's or two 1's) 
table(hiv_m_d$eversex)
hiv_m_d$eversex_clean <- NA
hiv_m_d$eversex_clean[hiv_m_d$eversex==1 | hiv_m_d$eversex==8] <- 1
hiv_m_d$eversex_clean[hiv_m_d$eversex==2] <- 0
table(hiv_m_d$eversex_clean)
prop.table(table(hiv_m_d$eversex_clean, hiv_m_d$sex),2)
names(hiv_m_d)

hiv_m_d$eversex_clean_up=hiv_m_d$eversex_clean
hiv_m_d$eversex_clean_down=hiv_m_d$eversex_clean

hiv_m_d <- hiv_m_d %>%
  group_by(research_id) %>%
  fill(eversex_clean_down, .direction = "down") %>% #impute missing negative/positive values between early late present values
  fill(eversex_clean_up, .direction = "up") %>%
  mutate(eversex_clean_imputed=ifelse(eversex_clean_down==eversex_clean_up, eversex_clean_down, 
                             ifelse(eversex_clean_down==0 & eversex_clean_up==1, NA, NA))) %>% #do not impute values between a negative and positive
  dplyr::select(-c(eversex_clean_down, eversex_clean_up))

#fill in missing - gaps - in sex, birthdate and resident status

hiv_m_d$birthdate <- as.Date(hiv_m_d$birthdate, "%d/%m/%Y")

hiv_m_d <- hiv_m_d %>% group_by(research_id) %>% 
  fill(sex, .direction = "downup") %>% 
  fill(birthdate, .direction="downup") %>% 
  fill(resident_imputed, .direction="downup") #impute missing resident values ***Check with Kate*** change this to not impute

table(hiv_m_d$hivinc) #1352 incident infections in all, 1255 among residents

hiv_m_d$ageyrs2 <- floor(as.numeric(hiv_m_d$hivdate-hiv_m_d$birthdate)/365.25) #age based on birthday rounded down
hiv_m_d$ageyrs2 <- ifelse(is.na(hiv_m_d$ageyrs2), hiv_m_d$ageyrs, hiv_m_d$ageyrs2)
hiv_m_d<- subset(hiv_m_d, ageyrs2 >= 15 & ageyrs2 <= 50)

table(hiv_m_d$hivinc) #1317 incident infections in all, 1255 among residents

sum(is.na(hiv_m_d$ageyrs2)) #no NA values

#clean other variables
table(hiv_m_d$hivinc)
names(hiv_m_d)
table(hiv_m_d$mobility)

length(unique(hiv_m_d$research_id)) #40,800 unique individuals
length(unique(hiv_m_d$research_id[hiv_m_d$sex=="F"])) #23,000 unique individuals
length(unique(hiv_m_d$research_id[hiv_m_d$sex=="M"])) #17,800 unique individuals

################################################################################################
#Decide wither to restrict to residents or allow for non-resident person-time to be included
################################################################################################

table(hiv_m_d$round, hiv_m_d$resident_imputed) #data set only includes those who are residents
hiv_m_d_r <- hiv_m_d #no resident restriction
#hiv_m_d_r <- hiv_m_d[hiv_m_d$resident_imputed==1,] #restrict to residents==1

length(unique(hiv_m_d_r$research_id)) #40,800 unique individuals

#descriptive stats on cohort
table(hiv_m_d_r$hivinc) #1317 incident infections in all, 
table(hiv_m_d_r$hivinc[hiv_m_d_r$sex=="F"]) #795 incident infections in all,
table(hiv_m_d_r$hivinc[hiv_m_d_r$sex=="M"]) #522 incident infections in all, 
table(hiv_m_d_r$hivinc, hiv_m_d_r$year) #    

table(hiv_m_d_r$hivinc[hiv_m_d_r$sex=="F"]) #795 incident infections in all,
table(hiv_m_d_r$hivinc[hiv_m_d_r$sex=="M"]) #522 incident infections in all, 

#1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020
#0 3865 5442 4233 5667 6396 4720 3989 5102 4067 5188 5068 4850 6879 4961 3894 6596 5056 6991 5834 4714 2854 3035
#1   28   54   71   72   82   71   81   91   78   79   91   71   80   69   49   69   31   64   39   24    8   12

sum(hiv_m_d_r$py, na.rm=T) #141537.3 person-years of follow-up
sum(hiv_m_d_r$py[hiv_m_d_r$sex=="F"], na.rm=T) #76055.88 person-years of follow-up
sum(hiv_m_d_r$py[hiv_m_d_r$sex=="M"], na.rm=T) #65481.37 person-years of follow-up
length(unique(hiv_m_d_r$research_id[hiv_m_d_r$py>0])) #22017

#crude incidence rates
100*sum(hiv_m_d_r$hivinc, na.rm=T)/sum(hiv_m_d_r$py, na.rm=T) #0.93 per 100 py
100*sum(hiv_m_d_r$hivinc[hiv_m_d_r$sex=="F"], na.rm=T)/sum(hiv_m_d_r$py[hiv_m_d_r$sex=="F"], na.rm=T) #1.04 per 100 py
100*sum(hiv_m_d_r$hivinc[hiv_m_d_r$sex=="M"], na.rm=T)/sum(hiv_m_d_r$py[hiv_m_d_r$sex=="M"], na.rm=T) #0.80 per 100 py

#crude incidence rates by age and gender
100*sum(hiv_m_d_r$hivinc[hiv_m_d_r$sex=="F" & hiv_m_d_r$ageyrs2 <25], na.rm=T)/sum(hiv_m_d_r$py[hiv_m_d_r$sex=="F" & hiv_m_d_r$ageyrs2 <25], na.rm=T) #1.23982
100*sum(hiv_m_d_r$hivinc[hiv_m_d_r$sex=="M"], na.rm=T)/sum(hiv_m_d_r$py[hiv_m_d_r$sex=="M"], na.rm=T) #0.80 per 100 py

hist(hiv_m_d_r$ageyrs2)
hist(hiv_m_d_r$py)
summary(hiv_m_d_r$py) #check low and high values: Answer, for those who seroconvert, the timing of serconversion is randomly assigned and so could occur right after the last negative visit hence giving low values

names(hiv_m_d_r)
View(hiv_m_d_r[,c(1,2,7,14,17,32,37)])

#Create table for reporting in manuscript

library(dplyr)

# --- SETTINGS (edit these if your column names differ) ---
df <- hiv_m_d_r
sex_var <- "sex"     # e.g. "sex" or "gender"
age_var <- "ageyrs2"     # e.g. "age" or "ageyrs2"

# Helper to refer to columns by string names
sx  <- rlang::sym(sex_var)
age <- rlang::sym(age_var)

# Keep rows that contribute to incidence (have PY and hivinc coded)
dat <- df %>%
  ungroup() %>%                                  # <-- important
  filter(hivinc %in% c(0, 1), !is.na(py)) %>%
  mutate(
    age_band = cut(!!age,
                   breaks = c(15, 25, 35, 50),
                   labels = c("15-24","25-34","35-49"),
                   right = FALSE, include.lowest = TRUE)
  )

# Define age bands: (15–24, 25–34, 35–49)
dat <- dat %>%
  mutate(age_band = cut(!!age,
                        breaks = c(15, 25, 35, 50),
                        labels = c("15-24","25-34","35-49"),
                        right = FALSE, include.lowest = TRUE))

summarise_rates <- function(.data) {
  .data %>%
    summarise(
      person_years  = sum(py, na.rm = TRUE),
      cases         = sum(hivinc == 1, na.rm = TRUE),
      inc_per_100py = 100 * cases / person_years,
      lower_95 = 100 * (ifelse(cases == 0, 0,
                               qchisq(0.025, 2 * cases) / 2) / person_years),
      upper_95 = 100 * (qchisq(0.975, 2 * (cases + 1)) / 2) / person_years,
      .groups = "drop"
    ) %>%
    mutate(
      inc_95ci = sprintf("%.2f (%.2f–%.2f)", inc_per_100py, lower_95, upper_95)
    ) %>%
    ungroup()
}

# Apply summaries (assuming dat already has hivinc, py, sex, age_band)
overall_tbl <- dat %>% summarise_rates() %>%
  mutate(group = "Overall", subgroup = NA_character_) %>% relocate(group, subgroup)

by_sex_tbl <- dat %>%
  group_by(sex) %>%
  summarise_rates() %>%
  mutate(group = "Sex", subgroup = as.character(sex)) %>%
  select(group, subgroup, everything())

by_sex_age_tbl <- dat %>%
  filter(!is.na(age_band)) %>%
  group_by(sex, age_band) %>%
  summarise_rates() %>%
  mutate(group = "Sex × Age", subgroup = paste0(sex, " / ", age_band)) %>%
  select(group, subgroup, everything())

incidence_table <- bind_rows(overall_tbl, by_sex_tbl, by_sex_age_tbl) %>%
  mutate(person_years = round(person_years)) %>%            # round to nearest whole
  select(group, subgroup, person_years, cases, inc_95ci) %>%# keep only desired cols
  arrange(match(group, c("Overall","Sex","Sex × Age")), subgroup)

incidence_table

################################################################################################
#Rename hiv_m_d_r as hivstatus_vlcopies_1_demog and save Rds file
################################################################################################

hivstatus_vlcopies_1_demog = hiv_m_d_r
hivstatus_vlcopies_1_demog$agecat <- cut(hivstatus_vlcopies_1_demog$ageyrs2, c(15,20,25,30,35,40,45,50), right=F)
length(unique(hivstatus_vlcopies_1_demog$research_id))

setwd("C:\\Users\\adamak\\OneDrive - Gates Foundation\\Dropbox\\EMOD_rakai\\Data\\RakaiData_23Oct2023")
saveRDS(hivstatus_vlcopies_1_demog, file="hivstatus_vlcopies_1_demog_v4.rds")

