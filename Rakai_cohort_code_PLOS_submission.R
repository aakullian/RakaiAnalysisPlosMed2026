#################################################################################################
# Rakai cohort data analysis
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

#sets the directory to where this script is stored#
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

################################################################################################
#Read in cleaned dataset
################################################################################################

#brings in cleaned cohort data
hivstatus_vlcopies_1_demog <- readRDS("hivstatus_vlcopies_1_demog_v4.rds")

#create a year for each round 
hivdate_round <- hivstatus_vlcopies_1_demog %>% 
  dplyr::group_by(round) %>% 
  dplyr::summarise(hivdate=mean(hivdate, na.rm=T))
hivdate_round$Year <- as.numeric(substr(hivdate_round$hivdate, 1,4))
roundyear <- hivdate_round[,c(1,3)]
length(unique(hivstatus_vlcopies_1_demog$research_id))
length(unique(hivstatus_vlcopies_1_demog$research_id[is.na(hivstatus_vlcopies_1_demog$py)==F]))

#decide whether to restrict to residents
hivstatus_vlcopies_1_demog <- hivstatus_vlcopies_1_demog[hivstatus_vlcopies_1_demog$resident_imputed==1,] #restrict to residents==1
table(hivstatus_vlcopies_1_demog$hivinc) #1317 incident infections in all
table(hivstatus_vlcopies_1_demog$round)

hivstatus_vlcopies_1_demog %>%
  filter(hivinc==1) %>%
  group_by(sex) %>%
  summarise(n_participants = n_distinct(research_id))

hivstatus_vlcopies_1_demog %>%
  filter(hivinc==1, is.na(py)==F) %>%
  ungroup() %>%  # ensure no active grouping
  summarise(n_participants = n_distinct(research_id))

hivstatus_vlcopies_1_demog %>%
  filter(hivinc==1, ageyrs2 < 25) %>%
  group_by(sex) %>%
  summarise(n_participants = n_distinct(research_id))

# Step 1: find IDs with hivinc == 1 and py == NA
ids_with_na_py <- hivstatus_vlcopies_1_demog %>%
  filter(hivinc == 1, is.na(py)) %>%
  distinct(research_id) %>%
  pull(research_id)

# Step 2: show all records for those IDs
ids_with_na_py_all <- hivstatus_vlcopies_1_demog %>%
  filter(research_id %in% ids_with_na_py) %>%
  arrange(research_id, round)

#How many rounds did each participant complete?

# Count the number of rounds completed by each research_id
rounds_per_id <- hivstatus_vlcopies_1_demog %>%
  dplyr::group_by(research_id) %>%
  dplyr::summarise(n_rounds = n())

# Compute the average number of rounds per participant
average_rounds <- mean(rounds_per_id$n_rounds)

average_rounds

################################################################################################
#Incidence models
################################################################################################

table(hivstatus_vlcopies_1_demog$round)

hivstatus_vlcopies_1_demog %>%
  filter(!is.na(py)) %>%
  ungroup() %>%  # ensure no active grouping
  summarise(n_participants = n_distinct(research_id))

hivstatus_vlcopies_1_demog %>%
  filter(!is.na(py), ageyrs2<25) %>%
  ungroup() %>%  # ensure no active grouping
  summarise(n_participants = n_distinct(research_id))

#crude cases and py
table(hivstatus_vlcopies_1_demog$round)
crude_cases_py <- subset(hivstatus_vlcopies_1_demog, is.na(py)==F) %>%
  group_by(sex,round, ageyrs2) %>%
  summarise(pysum=sum(py), seroconvsum=sum(hivinc))

crude_cases_py_agecat <- subset(hivstatus_vlcopies_1_demog, is.na(py)==F) %>%
  group_by(sex,round, agecat) %>%
  summarise(pysum=sum(py), seroconvsum=sum(hivinc))  %>%
  mutate(strata=sex, strata=recode(strata, "F"="Female", "M"="Male")) %>%
  na.omit() 

table(hivstatus_vlcopies_1_demog$agecat)

table(hivstatus_vlcopies_1_demog$ageyrs2)
crude_cases_py_agecat3 <- subset(hivstatus_vlcopies_1_demog, is.na(py)==F) %>%
  mutate(agecat3 = cut(ageyrs2, breaks=c(0,24,34,49), labels=c("15-24","25-34","35-49"))) %>%
  group_by(sex,round, agecat3) %>%
  summarise(pysum=round(sum(py),0), seroconvsum=sum(hivinc))  %>%
  mutate(sex=recode(sex, "F"="Women", "M"="Men")) %>%
  rename(age=agecat3) %>%
  na.omit() 

subset(hivstatus_vlcopies_1_demog, is.na(py)==F) %>%
  mutate(agecat3 = cut(ageyrs2, breaks=c(0,24,34,49), labels=c("15-24","25-34","35-49"))) %>%
  group_by(sex, agecat3) %>%
  summarise(pysum=round(sum(py),0), seroconvsum=sum(hivinc))  %>%
  mutate(sex=recode(sex, "F"="Women", "M"="Men"), inc=seroconvsum/pysum * 100) %>%
  rename(age=agecat3) %>%
  na.omit() 

subset(hivstatus_vlcopies_1_demog, is.na(py)==F) %>%
  group_by(sex) %>%
  summarise(pysum=round(sum(py),0), seroconvsum=sum(hivinc))  %>%
  mutate(sex=recode(sex, "F"="Women", "M"="Men"), inc=seroconvsum/pysum * 100) %>%
  na.omit() 

subset(hivstatus_vlcopies_1_demog, is.na(py)==F) %>%
  ungroup() %>%
  summarise(pysum=round(sum(py),0), seroconvsum=sum(hivinc))  %>%
  mutate(inc=seroconvsum/pysum * 100) %>%
  na.omit() 

#run Poisson models to estimate incidence rates by year
glmfit <- glm(hivinc ~ as.factor(round) + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog,is.na(py)==F))
gamfit <- gam(hivinc ~ s(round, bs="tp") + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog,is.na(py)==F))
glmfit.m <- glm(hivinc ~ as.factor(round) + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog, sex=="M" & is.na(py)==F))
gamfit.m <- gam(hivinc ~ s(round, bs="tp") + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog, sex=="M" & is.na(py)==F))
glmfit.f <- glm(hivinc ~ as.factor(round) + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog, sex=="F" & is.na(py)==F))
gamfit.f <- gam(hivinc ~ s(round, bs="tp") + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog,sex=="F" & is.na(py)==F))

gamfit.m.age <- gam(hivinc ~ s(round) + s(ageyrs2) + ti(round, ageyrs2) + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog,sex=="M" & is.na(py)==F))
gamfit.f.age <- gam(hivinc ~ s(round) + s(ageyrs2) + ti(round, ageyrs2) + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog,sex=="F" & is.na(py)==F))

gamfit.m.1524 <- gam(hivinc ~ s(round, bs="tp") + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog,sex=="M" & is.na(py)==F & ageyrs2 <25))
gamfit.f.1524 <- gam(hivinc ~ s(round, bs="tp") + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog,sex=="F" & is.na(py)==F & ageyrs2 <25))
gamfit.m.2534 <- gam(hivinc ~ s(round, bs="tp") + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog,sex=="M" & is.na(py)==F & ageyrs2 >24 & ageyrs2 <35))
gamfit.f.2534 <- gam(hivinc ~ s(round, bs="tp") + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog,sex=="F" & is.na(py)==F & ageyrs2 >24 & ageyrs2 <35))
gamfit.m.3549 <- gam(hivinc ~ s(round, bs="tp") + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog,sex=="M" & is.na(py)==F & ageyrs2 >34 & ageyrs2 <50))
gamfit.f.3549 <- gam(hivinc ~ s(round, bs="tp") + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog,sex=="F" & is.na(py)==F & ageyrs2 >34 & ageyrs2 <50))

glmfit.m.1524 <- glm(hivinc ~ as.factor(round) + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog, sex=="M" & is.na(py)==F & ageyrs2 <25))
glmfit.f.1524 <- glm(hivinc ~ as.factor(round) + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog, sex=="F" & is.na(py)==F & ageyrs2 <25))
glmfit.m.2534 <- glm(hivinc ~ as.factor(round) + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog, sex=="M" & is.na(py)==F & ageyrs2 >24 & ageyrs2 <35))
glmfit.f.2534 <- glm(hivinc ~ as.factor(round) + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog, sex=="F" & is.na(py)==F & ageyrs2 >24 & ageyrs2 <35))
glmfit.m.3549 <- glm(hivinc ~ as.factor(round) + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog, sex=="M" & is.na(py)==F & ageyrs2 >34 & ageyrs2 <50))
glmfit.f.3549 <- glm(hivinc ~ as.factor(round) + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog, sex=="F" & is.na(py)==F & ageyrs2 >34 & ageyrs2 <50))

gamfit.m.agecat <- gam(hivinc ~ s(round, bs="tp") + agecat + s(round,by=agecat) + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog, sex=="M" & is.na(py)==F))
gamfit.f.agecat <- gam(hivinc ~ s(round, bs="tp") + agecat + s(round,by=agecat) + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog, sex=="F" & is.na(py)==F))
glmfit.m.agecat <- glm(hivinc ~ as.factor(round) + agecat + as.factor(round)*agecat + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog, sex=="M" & is.na(py)==F))
glmfit.f.agecat <- glm(hivinc ~ as.factor(round) + agecat + as.factor(round)*agecat + offset(log(py)), family = poisson, data = subset(hivstatus_vlcopies_1_demog, sex=="F" & is.na(py)==F))

summary(glmfit.m.agecat)
summary(glmfit.f.agecat)

gg <- expand.grid(round=seq(7,19), py=1) #create prediction grid
gg.age <- expand.grid(round=seq(7,19), py=1, ageyrs2=seq(15,49,1))

agecats <- unique(hivstatus_vlcopies_1_demog$agecat)
agecats <- agecats[!is.na(agecats)]
gg.agecat <- expand.grid(round=seq(7,19), py=1, agecat=agecats)

modelpreds_agecat3 <- rbind(
  cbind(sex="Men", age="15-24", model="glmfit",gg, as.data.frame(predict(glmfit.m.1524, gg, type="link", se.fit=TRUE))[,c(1:2)]),
  cbind(sex="Men", age="15-24", model="gamfit",gg, as.data.frame(predict(gamfit.m.1524, gg, type="link", se.fit=TRUE))),
  cbind(sex="Men", age="25-34", model="glmfit",gg, as.data.frame(predict(glmfit.m.2534, gg, type="link", se.fit=TRUE))[,c(1:2)]),
  cbind(sex="Men", age="25-34", model="gamfit",gg, as.data.frame(predict(gamfit.m.2534, gg, type="link", se.fit=TRUE))),
  cbind(sex="Men", age="35-49", model="glmfit",gg, as.data.frame(predict(glmfit.m.3549, gg, type="link", se.fit=TRUE))[,c(1:2)]),
  cbind(sex="Men", age="35-49", model="gamfit",gg, as.data.frame(predict(gamfit.m.3549, gg, type="link", se.fit=TRUE))),
  cbind(sex="Women", age="15-24", model="glmfit",gg, as.data.frame(predict(glmfit.f.1524, gg, type="link", se.fit=TRUE))[,c(1:2)]),
  cbind(sex="Women", age="15-24", model="gamfit",gg, as.data.frame(predict(gamfit.f.1524, gg, type="link", se.fit=TRUE))),
  cbind(sex="Women", age="25-34", model="glmfit",gg, as.data.frame(predict(glmfit.f.2534, gg, type="link", se.fit=TRUE))[,c(1:2)]),
  cbind(sex="Women", age="25-34", model="gamfit",gg, as.data.frame(predict(gamfit.f.2534, gg, type="link", se.fit=TRUE))),
  cbind(sex="Women", age="35-49", model="glmfit",gg, as.data.frame(predict(glmfit.f.3549, gg, type="link", se.fit=TRUE))[,c(1:2)]),
  cbind(sex="Women", age="35-49", model="gamfit",gg, as.data.frame(predict(gamfit.f.3549, gg, type="link", se.fit=TRUE)))
)
modelpreds_agecat3 <- modelpreds_agecat3 %>% mutate(inc=round(exp(fit)*100,2), se.fit=se.fit, lb=round(exp(fit - 1.96*(se.fit))*100,2), ub=round(exp(fit + 1.96*(se.fit))*100,2))
modelpreds_agecat3 <- merge(modelpreds_agecat3, roundyear, by="round")
modelpreds_agecat3 <- merge(modelpreds_agecat3, crude_cases_py_agecat3, by=c("sex","round","age"))

modelpreds.age <- rbind(cbind(strata="Male", model="tensorproductgam",gg.age, as.data.frame(predict(gamfit.m.age, gg.age, type="link", se.fit=TRUE))),
                        cbind(strata="Female", model="tensorproductgam",gg.age, as.data.frame(predict(gamfit.f.age, gg.age, type="link", se.fit=TRUE)))
)
modelpreds.age <- modelpreds.age %>% mutate(inc=exp(fit), se.fit=se.fit, lb=exp(fit - 1.96*(se.fit)), ub=exp(fit + 1.96*se.fit)) %>%
  group_by(ageyrs2, strata) %>%
  mutate(irr=inc/inc[round==7])

modelpreds.age <- merge(modelpreds.age, roundyear, by="round")

modelpreds.agecat <- rbind(
  cbind(strata="Male", model="tensorproductgam",gg.agecat, as.data.frame(predict(gamfit.m.agecat, gg.agecat, type="link", se.fit=TRUE))),
  cbind(strata="Female", model="tensorproductgam",gg.agecat, as.data.frame(predict(gamfit.f.agecat, gg.agecat, type="link", se.fit=TRUE)))
)

modelpreds.agecat <- modelpreds.agecat %>% mutate(inc=exp(fit), se.fit=se.fit, lb=exp(fit - 1.96*(se.fit)), ub=exp(fit + 1.96*se.fit)) %>%
  group_by(agecat, strata) %>%
  mutate(irr=inc/inc[round==7])

modelpreds.agecat <- merge(modelpreds.agecat, roundyear, by="round")
modelpreds.agecat <- merge(modelpreds.agecat, crude_cases_py_agecat, by=c("strata","round","agecat"))


#Add person-years for estimating number of incident infections in each group
py <- subset(hivstatus_vlcopies_1_demog,is.na(py)==F) %>%
  dplyr::group_by(round, ageyrs, sex)%>%
  dplyr::summarise(py=sum(py)) %>%
  mutate(ageyrs2=ageyrs, strata=sex, strata=recode(strata, "F"="Female", "M"="Male")) 

names(py)
names(modelpreds.age)
modelpreds.age2 <- modelpreds.age %>% left_join(py, by=c("round","ageyrs2","strata")) %>%
  mutate(newinf=inc*py.y, newinf_lb=lb*py.y, newinf_ub=ub*py.y)

################################################################################################
#Incidence plots
################################################################################################

modelpreds_agecat3_year <- modelpreds_agecat3 %>%
  left_join(hivdate_round, by="round")

ggplot() +
  geom_point(data=subset(modelpreds_agecat3, model=="glmfit" & sex!="All" & age!="All (15-49)"), aes(y=inc, x=Year, color=age),size=2)+ 
  geom_errorbar(data=subset(modelpreds_agecat3, model=="glmfit" & sex!="All" & age!="All (15-49)"), aes(x=Year, ymin=lb, ymax=ub, color=age), width=.3) +
  geom_line(data=subset(modelpreds_agecat3, model=="gamfit" & age!="All (15-49)" & sex!="All"), aes(y=inc, x=Year, color=age),size=1)+ 
  geom_ribbon(data=subset(modelpreds_agecat3, model=="gamfit" & age!="All (15-49)" & sex!="All"),aes(ymin=lb, ymax=ub, x=Year, fill=age), alpha = 0.3)+
  facet_grid(sex~age) +
  xlab("Year")+ ylab("Incidence rate per 100 py")+
  scale_x_continuous(breaks=seq(2000,2020,1)) + 
  scale_y_continuous(expand=c(0,0), trans="log2", breaks=c(0.03125, 0.0625,0.125,0.25,0.5,1,2),labels = scales::label_number(accuracy = 0.01)) +
  theme_bw(base_size=10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(legend.title=element_blank(), legend.position='none') + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

ggsave(
  filename = "incidence_by_gender_round_agecat_REVISION.tiff",
  plot = last_plot(),
  device = "tiff",
  width = 7,        # inches (≈2250 px at 300 dpi)
  height = 8.75,         # inches (≈2700 px at 300 dpi)
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

ggplot() +
  geom_point(data=subset(modelpreds_agecat3, model=="glmfit" & sex=="Women" & age!="All (15-49)"), aes(y=inc, x=Year, color=age),size=2)+ 
  geom_errorbar(data=subset(modelpreds_agecat3, model=="glmfit" & sex=="Women" & age!="All (15-49)"), aes(x=Year, ymin=lb, ymax=ub, color=age), width=.3) +
  geom_line(data=subset(modelpreds_agecat3, model=="gamfit" & age!="All (15-49)" & sex=="Women"), aes(y=inc, x=Year, color=age),size=1)+ 
  geom_ribbon(data=subset(modelpreds_agecat3, model=="gamfit" & age!="All (15-49)" & sex=="Women"),aes(ymin=lb, ymax=ub, x=Year, fill=age), alpha = 0.3)+
  facet_grid(sex~age) +
  xlab("Year")+ ylab("Incidence rate per 100 py")+
  scale_x_continuous(breaks=seq(2000,2020,1)) + 
  scale_y_continuous(expand=c(0,0), trans="log2", breaks=c(0.03125, 0.0625,0.125,0.25,0.5,1,2),labels = scales::label_number(accuracy = 0.01)) +
  theme_bw(base_size=10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(legend.title=element_blank(), legend.position='none') + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

ggsave(
  filename = "incidence_women_round_agecat_REVISION.tiff",
  plot = last_plot(),
  device = "tiff",
  width = 8.75,        # inches (≈2250 px at 300 dpi)
  height = 5.5,         # inches (≈2700 px at 300 dpi)
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

library(scales)

table(modelpreds.age$round)
ggplot(subset(modelpreds.age, Year %in% c(2000,2019) & strata == "Female")) +
  geom_line(aes(x=ageyrs2, y=inc*100, color=factor(Year), group=factor(Year)),size=1.5)+ 
  geom_ribbon(aes(x=ageyrs2, ymin=lb*100, ymax=ub*100, fill=factor(Year), group=factor(Year)),alpha=0.3,size=1.5)+ 
  #scale_colour_gradientn(colours=rainbow(3), guide="legend") +
  guides(colour = guide_legend(byrow=T, keywidth = 2, keyheight = 1, nrow=1)) +
  facet_grid(~strata) +
  scale_color_manual(values = c("2000" = "red", "2019" = "blue")) +
  scale_fill_manual(values  = c("2000" = "red", "2019" = "blue")) +
  theme(legend.position="bottom") +
  xlab("Age")+ylab("Inc/100py")+
  scale_x_continuous(expand=c(0,0), breaks=seq(15,50,5)) + 
  scale_y_continuous(
    expand = c(0, 0),
    trans = "log2",
    breaks = c(0.0156, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8),
    limits = c(0.008, 6),
    labels = number_format(accuracy = 0.01, decimal.mark = ".", trim = TRUE)
  )+
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme_bw(base_size=24) +
  theme(legend.position="bottom",legend.title=element_blank()) +
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black")) 

ggsave(
  filename = "incidence_by_gender_age_round7v19_women_REVISION_2.tiff",
  plot = last_plot(),
  device = "tiff",
  width = 7,        # inches (≈2250 px at 300 dpi)
  height = 7.5,         # inches (≈2700 px at 300 dpi)
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)


################################################################################################
#HIV prevalence estimates
################################################################################################

#Prevalence by gender over time

table(hivstatus_vlcopies_1_demog$sex)
table(hivstatus_vlcopies_1_demog$agecat)

#get population size of those with an HIV status and prevalence
tableprev <- subset(hivstatus_vlcopies_1_demog, hivstatus_imputed==1 | hivstatus_imputed==0) %>%
  group_by(year, sex) %>%
  summarize(prevcount = length(hivstatus_imputed[hivstatus_imputed==1 | hivstatus_imputed==0]), hivpos = length(hivstatus_imputed[hivstatus_imputed==1]),hivprev = hivpos/prevcount) %>%
  mutate(se=sqrt(hivprev*(1-hivprev)/prevcount), lb=hivprev-1.96*se, ub=hivprev+1.96*se)
#tableprev <- tableprev2[order(tableprev2$sex,tableprev2$year) , ]

tableprevagecat <- subset(hivstatus_vlcopies_1_demog, hivstatus_imputed==1 | hivstatus_imputed==0) %>%
  group_by(year, sex, agecat) %>%
  summarize(prevcount = length(hivstatus_imputed[hivstatus_imputed==1 | hivstatus_imputed==0]), hivpos = length(hivstatus_imputed[hivstatus_imputed==1]),hivprev = hivpos/prevcount) %>%
  mutate(se=sqrt(hivprev*(1-hivprev)/prevcount), lb=hivprev-1.96*se, ub=hivprev+1.96*se)

verif_1 <- verif_1[order(verif_1$research_id, verif_1$round),] #Eligibility data
hivstatus_vlcopies_1_demog$agecat <- cut(hivstatus_vlcopies_1_demog$ageyrs2, c(15,20,25,30,35,40,45,50), right=F)

#get population size of those with HIV status not missing
tablehiv_f <- subset(hivstatus_vlcopies_1_demog, hivstatus_imputed==1 | hivstatus_imputed==0 & sex=="F")
mytable.f <-table(tablehiv_f$agecat, tablehiv_f$year)
m.f <- matrix(mytable.f, dimnames=list(t(outer(colnames(mytable.f), rownames(mytable.f), FUN=paste)), NULL))
m.f <- as.data.frame(m.f)
library(tibble)
m.f <- tibble::rownames_to_column(m.f, c("VALUE"))
m.f <- data.frame(cbind(m.f$V1, do.call("rbind", strsplit(as.character(m.f$VALUE), " ", fixed = TRUE))))
colnames(m.f) <- c("effective_count", "Year", "AgeBin")
m.f$Gender <- "Female"
head(m.f)

tablehiv_m <- subset(hivstatus_vlcopies_1_demog, hivstatus_imputed==1 | hivstatus_imputed==0 & sex=="M")
mytable.m <-table(tablehiv_m$agecat, tablehiv_m$year)
m.m <- matrix(mytable.m, dimnames=list(t(outer(colnames(mytable.m), rownames(mytable.m), FUN=paste)), NULL))
m.m <- as.data.frame(m.m)
m.m <- tibble::rownames_to_column(m.m, c("VALUE"))
m.m <- data.frame(cbind(m.m$V1, do.call("rbind", strsplit(as.character(m.m$VALUE), " ", fixed = TRUE))))
colnames(m.m) <- c("effective_count", "Year", "AgeBin")
m.m$Gender <- "Male"
head(m.m)
m.merge <- rbind(m.f, m.m)

table(hivstatus_vlcopies_1_demog$hivstatus_imputed, useNA = "always")
#17679 prevalent infections in the longitudinal dataset
#118841 negative
#3968 without an HIV status - these could be imputed

gam.prev.m.age <- gam(hivstatus_imputed ~ s(round) + s(ageyrs2) + ti(round, ageyrs2), family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog, sex=="M"))
gam.prev.f.age <- gam(hivstatus_imputed ~ s(round) + s(ageyrs2) + ti(round, ageyrs2), family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog, sex=="F"))
glm.prev.m.age <- gam(hivstatus_imputed ~ round + agecat + round*agecat, family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog, sex=="M"))
glm.prev.f.age <- gam(hivstatus_imputed ~ round + agecat + round*agecat, family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog, sex=="F"))

summary(gam.prev.m.age)
summary(gam.prev.f.age)
summary(glm.prev.m.age)
summary(glm.prev.f.age)

gg.age <- expand.grid(round=seq(6,19), ageyrs2=seq(15,49,1))
gg.agecat <- expand.grid(round=seq(6,19), agecat=c("[15,20)", "[20,25)", "[25,30)", "[30,35)", "[35,40)", "[40,45)", "[45,50)"))

modelpreds.prev.age <- rbind(cbind(strata="Male", model="tensorproductgam",gg.age, as.data.frame(predict(gam.prev.m.age, gg.age, type="link", se.fit=TRUE))),
                             cbind(strata="Female", model="tensorproductgam",gg.age, as.data.frame(predict(gam.prev.f.age, gg.age, type="link", se.fit=TRUE)))
)

modelpreds.prev.agecat <- rbind(cbind(strata="Male", model="glm",gg.agecat, as.data.frame(predict(glm.prev.m.age, gg.agecat, type="link", se.fit=TRUE))),
                                cbind(strata="Female", model="glm",gg.agecat, as.data.frame(predict(glm.prev.f.age, gg.agecat, type="link", se.fit=TRUE)))
)

modelpreds.prev.age <- modelpreds.prev.age %>% mutate(odds=exp(fit), 
                                                      prob=odds/(1+odds), 
                                                      se.fit=se.fit, 
                                                      lb.odds=exp(fit - 1.96*se.fit), 
                                                      ub.odds=exp(fit + 1.96*se.fit),
                                                      lb.prob=lb.odds/(1+lb.odds), ub.prob=ub.odds/(1+ub.odds)) %>%
  group_by(ageyrs2, strata) %>%
  mutate(prevrr=prob/prob[round==6])

modelpreds.prev.agecat <- modelpreds.prev.agecat %>% mutate(odds=exp(fit), 
                                                            prob=odds/(1+odds), 
                                                            se.fit=se.fit, 
                                                            lb.odds=exp(fit - 1.96*se.fit), 
                                                            ub.odds=exp(fit + 1.96*se.fit),
                                                            lb.prob=lb.odds/(1+lb.odds), ub.prob=ub.odds/(1+ub.odds)) %>%
  group_by(agecat, strata) %>%
  mutate(prevrr=prob/prob[round==6])

summary(modelpreds.prev.age$prob)
summary(modelpreds.prev.agecat$prob)

#create prevalence by age dataset
hivdate_round$Year <- as.numeric(substr(hivdate_round$hivdate, 1,4))
roundyear <- hivdate_round[,c(1,3)]
names(modelpreds.prev.agecat)
modelpreds.prev.agecat <- merge(modelpreds.prev.agecat, roundyear, by="round")

#create prevalence by agecat dataset
modelpreds.prev.agecat.ingestform <- modelpreds.prev.agecat
modelpreds.prev.agecat.ingestform$Province = "All"
modelpreds.prev.agecat.ingestform$weight = ""
modelpreds.prev.agecat.ingestform$Prevalence = modelpreds.prev.agecat.ingestform$prob
names(modelpreds.prev.agecat.ingestform)
modelpreds.prev.agecat.ingestform$agecat <- gsub(",", ":", modelpreds.prev.agecat.ingestform$agecat)
m.merge$AgeBin <- gsub(",", ":", m.merge$AgeBin)
modelpreds.prev.agecat.ingestform$AgeBin <- modelpreds.prev.agecat.ingestform$agecat
modelpreds.prev.agecat.ingestform$Gender <- modelpreds.prev.agecat.ingestform$strata
names(modelpreds.prev.agecat.ingestform)
modelpreds.prev.agecat.ingestform <- merge(modelpreds.prev.agecat.ingestform, m.merge, by=c('Gender','Year', 'AgeBin'))
names(modelpreds.prev.agecat.ingestform)
modelpreds.prev.agecat.ingestform<- modelpreds.prev.agecat.ingestform[,c(2,17,1,3,19,18,20)]
modelpreds.prev.agecat.ingestform <- modelpreds.prev.agecat.ingestform[order(modelpreds.prev.agecat.ingestform$Gender),]

write.table(tableprev, "Rakai_prevalence_gender_R19_REVISION.csv",sep=",", row.names = F)
write.table(modelpreds.prev.agecat.ingestform, "Rakai_prevalence_preds_agecat_R19_REVISION.csv",sep=",", row.names = F)
write.table(modelpreds.prev.age[,c(1,2,3,4,8,11,12,13)], "Rakai_prevalence_preds_age_R19_REVISION.csv", col.names = c("Sex", "model", "round","age", "prevalence","lb","ub","prevRR"),sep=",", row.names = F)
names(modelpreds.prev.agecat)
write.table(modelpreds.prev.agecat[,c(3,2,1,14,4,8,11,12,13)], "Rakai_prevalence_preds_agecat_R19_REVISION.csv", col.names = c("model", "sex","round","year", "agecat", "prevalence","lb","ub","prevRR"),sep=",", row.names = F)

################################################################################################
#ART coverage and viral load
################################################################################################
table(hivstatus_vlcopies_1_demog$cuarvmed) #currently on ART
table(hivstatus_vlcopies_1_demog$arvmed) #ever on ART
table(hivstatus_vlcopies_1_demog$cuarvmed, hivstatus_vlcopies_1_demog$round)
table(hivstatus_vlcopies_1_demog$arvmed, hivstatus_vlcopies_1_demog$round) #ever on ART

ids=unique(hivstatus_vlcopies_1_demog$research_id[hivstatus_vlcopies_1_demog$hivstatus_imputed==1])

hivstatus_vlcopies_1_demog$onART <- 0
hivstatus_vlcopies_1_demog$onART[hivstatus_vlcopies_1_demog$arvmed==1] <- 1

tableartprevagecat <- subset(hivstatus_vlcopies_1_demog, hivstatus_imputed==1) %>%
  dplyr::group_by(year, sex, agecat) %>%
  dplyr::summarize(prevcount = length(hivstatus_imputed[onART==1 | onART==0]), artpos = length(hivstatus_imputed[onART==1]),artprev = artpos/prevcount) %>%
  dplyr::mutate(se=sqrt(artprev*(1-artprev)/prevcount), lb=artprev-1.96*se, ub=artprev+1.96*se)

hivstatus_vlcopies_1_demog$agecat_num <- unclass(hivstatus_vlcopies_1_demog$agecat)   
table(hivstatus_vlcopies_1_demog$agecat, hivstatus_vlcopies_1_demog$agecat_num)

ids_hiv_pos <- unique(hivstatus_vlcopies_1_demog$research_id[hivstatus_vlcopies_1_demog$hivstatus_imputed == 1])
#View(subset(hivstatus_vlcopies_1_demog, research_id %in% ids_hiv_pos))

table(hivstatus_vlcopies_1_demog$round)
#fit ART prev models
gam.artfit.m.age <- gam(onART ~ s(round) + s(ageyrs2) + ti(round, ageyrs2), family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog,sex=="M" & hivstatus_imputed==1))
gam.artfit.f.age <- gam(onART ~ s(round) + s(ageyrs2) + ti(round, ageyrs2), family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog,sex=="F" & hivstatus_imputed==1))
glm.artfit.m.agecat <- gam(onART ~ as.factor(round) + agecat + as.factor(round)*agecat, family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog, sex=="M" & hivstatus_imputed==1))
glm.artfit.f.agecat <- gam(onART ~ as.factor(round) + agecat + as.factor(round)*agecat, family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog, sex=="F" & hivstatus_imputed==1))

summary(gam.artfit.m.age)
summary(gam.artfit.f.age)
summary(glm.artfit.m.agecat)
summary(glm.artfit.f.agecat)

gg.age <- expand.grid(round=seq(6,19), ageyrs2=seq(15,49,1))
gg.agecat <- expand.grid(round=seq(6,19), agecat=c("[15,20)", "[20,25)", "[25,30)", "[30,35)", "[35,40)", "[40,45)", "[45,50)"))

preds.art.age <- rbind(cbind(strata="Male", model="tensorproductgam",gg.age, as.data.frame(predict(gam.artfit.m.age, gg.age, type="link", se.fit=TRUE))),
                       cbind(strata="Female", model="tensorproductgam",gg.age, as.data.frame(predict(gam.artfit.f.age, gg.age, type="link", se.fit=TRUE)))
)

preds.art.agecat <- rbind(cbind(strata="Male", model="glm",gg.agecat, as.data.frame(predict(glm.artfit.m.agecat, gg.agecat, type="link", se.fit=TRUE))),
                          cbind(strata="Female", model="glm",gg.agecat, as.data.frame(predict(glm.artfit.f.agecat, gg.agecat, type="link", se.fit=TRUE)))
)

preds.art.age <- preds.art.age %>% mutate(odds=exp(fit), 
                                          prob=odds/(1+odds), 
                                          se.fit=se.fit, 
                                          lb.odds=exp(fit - 1.96*se.fit), 
                                          ub.odds=exp(fit + 1.96*se.fit),
                                          lb.prob=lb.odds/(1+lb.odds), ub.prob=ub.odds/(1+ub.odds)) %>%
  group_by(ageyrs2, strata) %>%
  mutate(artrr=prob/prob[round==6])

preds.art.agecat <- preds.art.agecat %>% mutate(odds=exp(fit),
                                                prob=odds/(1+odds),
                                                se.fit=se.fit,
                                                lb.odds=exp(fit - 1.96*se.fit),
                                                ub.odds=exp(fit + 1.96*se.fit),
                                                lb.prob=lb.odds/(1+lb.odds), ub.prob=ub.odds/(1+ub.odds)) %>%
  group_by(agecat, strata) %>%
  mutate(artrr=prob/prob[round==6])

summary(preds.art.age$prob)
summary(preds.art.agecat$prob)

#fix out of range values for agecat models
preds.art.agecat$ub.prob[is.na(preds.art.agecat$ub.prob)==T] <- 0
preds.art.agecat$ub.prob[preds.art.agecat$prob<0.001] <- 0
summary(preds.art.agecat$prob)
summary(preds.art.agecat$lb.prob)
summary(preds.art.agecat$ub.prob)
preds.art.agecat$prob <- round(preds.art.agecat$prob, 3)
preds.art.agecat$lb.prob <- round(preds.art.agecat$lb.prob, 3)
preds.art.agecat$ub.prob <- round(preds.art.agecat$ub.prob, 3)

#fix out of range values for age continuous models
preds.art.age$ub.prob[is.na(preds.art.age$ub.prob)==T] <- 0
preds.art.age$ub.prob[preds.art.age$prob<0.01] <- 0
summary(preds.art.age$prob)
summary(preds.art.age$lb.prob)
summary(preds.art.age$ub.prob)
preds.art.age$prob <- round(preds.art.age$prob, 3)
preds.art.age$lb.prob <- round(preds.art.age$lb.prob, 3)
preds.art.age$ub.prob <- round(preds.art.age$ub.prob, 3)

#merge in year to round
preds.art.agecat <- merge(preds.art.agecat, roundyear, by="round")
preds.art.age <- merge(preds.art.age, roundyear, by="round")

################################################################################################
#VMMC
################################################################################################
table(hivstatus_vlcopies_1_demog$circ_imputed_fix, hivstatus_vlcopies_1_demog$round, hivstatus_vlcopies_1_demog$agecat)

#get population size and ART coverage of those with HIV status not missing
tablevmmc <- subset(hivstatus_vlcopies_1_demog, sex=="M" & is.na(circ_imputed_fix)==F) %>%
  group_by(agecat, round, .drop = FALSE) %>%
  summarize(VMMCcount = length(circ_imputed_fix),VMMCpos = length(circ_imputed_fix[circ_imputed_fix==1]),VMMCprev = mean(circ_imputed_fix)) %>%
  mutate(se=sqrt(VMMCprev*(1-VMMCprev)/VMMCcount), lb=VMMCprev-1.96*se, ub=VMMCprev+1.96*se, sex="Male")

write.table(tablevmmc, "Rakai_VMMC_proportions_R19.csv", row.names = F,sep=",")

glm.vmmc.m.agecat <- gam(circ_imputed_fix ~ s(round) + agecat + s(round, by=agecat), family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog, sex=="M"))
glm.vmmc.m.age <- gam(circ_imputed_fix ~ s(round) + s(ageyrs2) + ti(round, ageyrs2), family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog,sex=="M"))

summary(glm.vmmc.m.agecat)
summary(glm.vmmc.m.age)

gg.age <- expand.grid(round=seq(6,19), ageyrs2=seq(17,47,5))
gg.agecat <- expand.grid(round=seq(6,19), agecat=c("[15,20)", "[20,25)", "[25,30)", "[30,35)", "[35,40)", "[40,45)", "[45,50)"))

modelpreds.vmmc.agecat <- cbind(strata="Male", model="glm",gg.agecat, as.data.frame(predict(glm.vmmc.m.agecat, gg.agecat, type="link", se.fit=TRUE)))

modelpreds.vmmc.age <- rbind(cbind(strata="Male", model="tensorproductgam",gg.age, as.data.frame(predict(glm.vmmc.m.age, gg.age, type="link", se.fit=TRUE))),
                             cbind(strata="Female", model="tensorproductgam",gg.age, as.data.frame(predict(glm.vmmc.m.age, gg.age, type="link", se.fit=TRUE)))
)

modelpreds.vmmc.agecat <- modelpreds.vmmc.agecat %>% mutate(odds=exp(fit),
                                                            prob=odds/(1+odds),
                                                            se.fit=se.fit,
                                                            lb.odds=exp(fit - 1.96*se.fit),
                                                            ub.odds=exp(fit + 1.96*se.fit),
                                                            lb.prob=lb.odds/(1+lb.odds), ub.prob=ub.odds/(1+ub.odds)) %>%
  group_by(agecat, strata) %>%
  mutate(artrr=prob/prob[round==8])

modelpreds.vmmc.age <- modelpreds.vmmc.age %>% mutate(odds=exp(fit),
                                                      prob=odds/(1+odds),
                                                      se.fit=se.fit,
                                                      lb.odds=exp(fit - 1.96*se.fit),
                                                      ub.odds=exp(fit + 1.96*se.fit),
                                                      lb.prob=lb.odds/(1+lb.odds), ub.prob=ub.odds/(1+ub.odds)) %>%
  group_by(ageyrs2, strata) %>%
  mutate(artrr=prob/prob[round==8])

summary(modelpreds.vmmc.agecat$prob)
summary(modelpreds.vmmc.age$prob)

#merge in year to round
modelpreds.vmmc.agecat <- merge(modelpreds.vmmc.agecat, roundyear, by="round")
modelpreds.vmmc.age <- merge(modelpreds.vmmc.age, roundyear, by="round")

################################################################################################
#Age at sexual debut and ever had sex
################################################################################################

#NOTE: lots of missingness for eversex.  Individuals usually report age first sex at first round (usually round 6).  

names(hivstatus_vlcopies_1_demog)

hist(hivstatus_vlcopies_1_demog$ag1stsex)
table(hivstatus_vlcopies_1_demog$ag1stsex)
table(hivstatus_vlcopies_1_demog$eversex_clean_imputed)

eversex_table <- subset(hivstatus_vlcopies_1_demog) %>%
  dplyr::group_by(agecat, round, sex) %>%
  dplyr::summarize(eversexcount = length(eversex_clean_imputed), eversex = length(eversex_clean_imputed[eversex_clean_imputed==1]),eversexprop = eversex/eversexcount) %>%
  dplyr::mutate(se=sqrt((eversexprop*(1-eversexprop))/eversexcount), lb=eversexprop-1.96*se, ub=eversexprop+1.96*se)

#merge in Year
hivdate_round$Year <- as.numeric(substr(hivdate_round$hivdate, 1,4))
roundyear <- hivdate_round[,c(1,3)]
eversex_table <- merge(eversex_table, roundyear, by="round")
eversex_table

hivstatus_vlcopies_1_demog$ag1stsex_clean <- hivstatus_vlcopies_1_demog$ag1stsex
hivstatus_vlcopies_1_demog$ag1stsex_clean[hivstatus_vlcopies_1_demog$ag1stsex<10 | hivstatus_vlcopies_1_demog$ag1stsex > 96] <- NA
table(hivstatus_vlcopies_1_demog$ag1stsex_clean)
hist(hivstatus_vlcopies_1_demog$ag1stsex_clean)

age1stsex_round <- subset(hivstatus_vlcopies_1_demog, ageyrs2>20 & ageyrs2 <25) %>%
  group_by(round, sex) %>%
  summarise(median_age1stsex = median(ag1stsex_clean, na.rm=T))

#models
gam.eversexprev.m.age <- gam(eversex_clean_imputed ~ s(round) + s(ageyrs2) + ti(round, ageyrs2), family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog, sex=="M"))
gam.eversexprev.f.age <- gam(eversex_clean_imputed ~ s(round) + s(ageyrs2) + ti(round, ageyrs2), family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog, sex=="F"))

gg.age <- expand.grid(round=seq(6,19), ageyrs2=seq(15,49,1))

modelpreds.prev.age <- rbind(cbind(strata="Male", model="tensorproductgam",gg.age, as.data.frame(predict(gam.eversexprev.m.age, gg.age, type="link", se.fit=TRUE))),
                             cbind(strata="Female", model="tensorproductgam",gg.age, as.data.frame(predict(gam.eversexprev.f.age, gg.age, type="link", se.fit=TRUE)))
)


modelpreds.prev.age <- modelpreds.prev.age %>% mutate(odds=exp(fit), 
                                                      prob=odds/(1+odds), 
                                                      se.fit=se.fit, 
                                                      lb.odds=exp(fit - 1.96*se.fit), 
                                                      ub.odds=exp(fit + 1.96*se.fit),
                                                      lb.prob=lb.odds/(1+lb.odds), ub.prob=ub.odds/(1+ub.odds)) %>%
  group_by(ageyrs2, strata) %>%
  mutate(prevrr=prob/prob[round==6])

summary(modelpreds.prev.age$prob)

#merge in year
hivdate_round$Year <- as.numeric(substr(hivdate_round$hivdate, 1,4))
roundyear <- hivdate_round[,c(1,3)]
names(modelpreds.prev.age)
modelpreds.prev.age <- merge(modelpreds.prev.age, roundyear, by="round")

#ever sex modelled agecat
table(hivstatus_vlcopies_1_demog$agecat)
gam.eversexprev.m.agecat.1519 <- gam(eversex_clean_imputed ~ s(round), family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog, sex=="M" & agecat=="[15,20)"))
gam.eversexprev.f.agecat.1519 <- gam(eversex_clean_imputed ~ s(round), family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog, sex=="F" & agecat=="[15,20)"))
gam.eversexprev.m.agecat.2024 <- gam(eversex_clean_imputed ~ s(round), family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog, sex=="M" & agecat=="[20,25)"))
gam.eversexprev.f.agecat.2024 <- gam(eversex_clean_imputed ~ s(round), family = binomial, method = "REML", data = subset(hivstatus_vlcopies_1_demog, sex=="F" & agecat=="[20,25)"))

gg.agecat <- expand.grid(round=seq(6,19))

modelpreds.prev.agecat <- rbind(cbind(sex="Male", agecat="[15,20)", model="gam",gg.agecat, as.data.frame(predict(gam.eversexprev.m.agecat.1519, gg.agecat, type="link", se.fit=TRUE))),
                                cbind(sex="Female", agecat="[15,20)", model="gam",gg.agecat, as.data.frame(predict(gam.eversexprev.f.agecat.1519, gg.agecat, type="link", se.fit=TRUE))),
                                cbind(sex="Male", agecat="[20,25)", model="gam",gg.agecat, as.data.frame(predict(gam.eversexprev.m.agecat.2024, gg.agecat, type="link", se.fit=TRUE))),
                                cbind(sex="Female", agecat="[20,25)", model="gam",gg.agecat, as.data.frame(predict(gam.eversexprev.f.agecat.2024, gg.agecat, type="link", se.fit=TRUE)))
)


modelpreds.prev.agecat <- modelpreds.prev.agecat %>% mutate(odds=exp(fit), 
                                                            prob=odds/(1+odds), 
                                                            se.fit=se.fit, 
                                                            lb.odds=exp(fit - 1.96*se.fit), 
                                                            ub.odds=exp(fit + 1.96*se.fit),
                                                            lb.prob=lb.odds/(1+lb.odds), ub.prob=ub.odds/(1+ub.odds)) %>%
  group_by(agecat, sex) %>%
  mutate(prevrr=prob/prob[round==6])

summary(modelpreds.prev.agecat$prob)

#merge in year
hivdate_round$Year <- as.numeric(substr(hivdate_round$hivdate, 1,4))
roundyear <- hivdate_round[,c(1,3)]
names(modelpreds.prev.agecat)
modelpreds.prev.agecat <- merge(modelpreds.prev.agecat, roundyear, by="round")

head(modelpreds.prev.agecat)

eversex_modelpreds_prev_agecat <- modelpreds.prev.agecat

setwd("plots\\Figs_for_manuscript")

# --- 1) VMMC preds ---
vmmc_std <- modelpreds.vmmc.agecat %>%
  transmute(
    year   = as.integer(Year),
    strata = strata,
    agecat = as.character(agecat),
    prob   = prob,        # already on 0–1 scale per your head()
    lb     = lb.prob,
    ub     = ub.prob
  )

# --- 2) Ever-sex preds ---
eversex_std <- eversex_modelpreds_prev_agecat %>%
  transmute(
    year   = as.integer(Year),
    strata = sex,
    agecat = as.character(agecat),
    prob   = prob,
    lb     = lb.prob,
    ub     = ub.prob
  )

# --- 3) ART prevalence table ---
art_std <- tableartprevagecat %>%
  transmute(
    year   = as.integer(year),
    strata = sex,
    agecat = as.character(agecat),
    prob   = artprev,   # use artprev here
    lb     = lb,
    ub     = ub
  )

# --- bind, clean, optional recodes ---
out <- bind_rows(
  vmmc_std %>% mutate(source = "vmmc"),
  eversex_std %>% mutate(source = "eversex"),
  art_std %>% rename(strata = sex) %>% mutate(source = "art")   # rename here
) %>%
  mutate(
    strata = case_when(
      strata %in% c("F", "Female") ~ "Female",
      strata %in% c("M", "Male")   ~ "Male",
      TRUE ~ strata
    )
  ) %>%
  select(year, strata, agecat, prob, lb, ub, source) %>%  # drop extra sex column
  filter(!is.na(prob)) %>%
  arrange(year, strata, agecat)

out

library(ggplot2)
library(dplyr)
library(forcats)

# ----- optional: make sure age groups order correctly
age_order <- c("15-19","20-24","25-29","30-34","35-39","40-44","45-49")

plot_dat <- figdat %>%                       # <-- your data frame here
  mutate(
    age_cat   = factor(age_cat, levels = age_order),
    indicator = factor(indicator,
                       levels = c("Female ART","Male ART","VMMC",
                                  "Female initiated sex","Male initiated sex"))
  )

# Color palette (Okabe–Ito: colorblind-friendly; resembles your figure)
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Map to your 7 age groups in the order above
age_cols <- setNames(okabe_ito, age_order)

p <- ggplot(
  plot_dat,
  aes(x = survey_round, y = value, color = age_cat, group = age_cat)
) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  facet_wrap(~ indicator, nrow = 1, scales = "fixed") +
  scale_color_manual(values = age_cols, name = NULL) +
  scale_x_continuous(
    breaks = sort(unique(plot_dat$survey_round)),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(x = "Survey round", y = "%") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_minimal(base_size = 11, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 4, r = 0, b = 0, l = 0),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(face = "bold")
  )

p

library(magick)

# --- Settings ---
img_dir    <- "Fig5_elements"
out_file   <- file.path(img_dir, "Fig5_combined.tiff")
ncol       <- 2
target_w   <- 1200

label_size <- 40
label_color<- "black"
label_box  <- "white"

# original offsets (in pixels)
label_x_px <- 60
label_y_px <- 60

# Move labels UP and LEFT by these amounts (in inches)
move_up_in   <- 0.5   # up 0.5 inch  = 150 px
move_left_in <- 0.25  # left 0.25 in = 75 px

dpi <- 300
label_x_px <- max(0, label_x_px - round(move_left_in * dpi))  # move left
label_y_px <- max(0, label_y_px - round(move_up_in * dpi))    # move up

label_pad <- sprintf("+%d+%d", label_x_px, label_y_px)
# ---------------

dir.create(img_dir, showWarnings = FALSE)

paths <- list.files(img_dir, pattern = "\\.jpe?g$", full.names = TRUE)
stopifnot(length(paths) > 0)

order_vec <- c(4, 3, 2, 1)
if (all(order_vec %in% seq_along(paths))) paths <- paths[order_vec]

imgs <- lapply(paths, image_read)

if (!is.null(target_w)) {
  imgs <- lapply(imgs, image_scale, geometry = paste0(target_w))
}

annotate_panel <- function(im, lab) {
  image_annotate(
    im,
    text = lab,
    size = label_size,
    color = label_color,
    boxcolor = label_box,
    gravity = "northwest",
    location = label_pad,
    weight = 700
  )
}

labs        <- LETTERS[seq_along(imgs)]
img_labeled <- Map(annotate_panel, imgs, labs)

split_rows <- split(img_labeled, ceiling(seq_along(img_labeled) / ncol))
row_imgs   <- lapply(split_rows, function(row_list) image_append(image_join(row_list), stack = FALSE))
combined   <- image_append(image_join(row_imgs), stack = TRUE)

combined_white <- combined |> image_background("white", flatten = TRUE)

image_write(
  combined_white,
  path = out_file,
  format = "tiff",
  compression = "lzw",
  density = "300x300"
)

####
library(magick)

# --- Settings (same as yours) ---
img_dir    <- "Fig5_elements"
out_file   <- file.path(img_dir, "Fig5_combined.tiff")
ncol       <- 2
target_w   <- 1200

label_size <- 60
label_color<- "black"
label_box  <- "white"

# Base offsets (pixels from top-left)
base_x_px <- 60
base_y_px <- 60

# Global nudges (inches → pixels); use 0 if none
move_up_in   <- 0.0
move_left_in <- 0.0

# Extra nudge ONLY for panel D (inches → pixels)
extra_up_D_in <- 0.5

dpi <- 300
base_x_px <- base_x_px - round(move_left_in * dpi)   # allow negatives
base_y_px <- base_y_px - round(move_up_in   * dpi)   # allow negatives

# Helper to format "+x+y" with proper signs (e.g., "+60-90")
loc_str <- function(x, y) sprintf("%+d%+d", x, y)

# --------------------------------
dir.create(img_dir, showWarnings = FALSE)

paths <- list.files(img_dir, pattern = "\\.jpe?g$", full.names = TRUE)
stopifnot(length(paths) > 0)

# Your explicit order
order_vec <- c(4, 3, 2, 1)
if (all(order_vec %in% seq_along(paths))) paths <- paths[order_vec]

imgs <- lapply(paths, image_read)
if (!is.null(target_w)) imgs <- lapply(imgs, image_scale, geometry = paste0(target_w))

labs <- LETTERS[seq_along(imgs)]

# Annotate; panel D (i == 4) gets extra upward nudge
img_labeled <- Map(function(im, lab, i) {
  # start from base
  x <- base_x_px
  y <- base_y_px
  if (i == 4) {  # Panel D
    y <- y - round(extra_up_D_in * dpi)  # move further up (more negative)
  }
  image_annotate(
    im,
    text     = lab,
    size     = label_size,
    color    = label_color,
    boxcolor = label_box,
    gravity  = "northwest",
    location = loc_str(x, y),
    weight   = 700
  )
}, imgs, labs, seq_along(imgs))

# Build grid and write
split_rows <- split(img_labeled, ceiling(seq_along(img_labeled) / ncol))
row_imgs   <- lapply(split_rows, function(row_list) image_append(image_join(row_list), stack = FALSE))
combined   <- image_append(image_join(row_imgs), stack = TRUE)

combined_white <- combined |> image_background("white", flatten = TRUE)
image_write(
  combined_white,
  path = out_file,
  format = "tiff",
  compression = "lzw",
  density = "300x300"
)

cat("Saved:", out_file, "\n")
