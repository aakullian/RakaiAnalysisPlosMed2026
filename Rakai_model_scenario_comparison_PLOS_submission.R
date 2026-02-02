################################################################################
## Adam Akullian
## Compare simulation results across scenarios
## Uses summary outputs computed for each scenario across simulation runs
## Saves a dataset (scenario_comparison.RData) and plots (scenario_comparison_plots.pdf)
################################################################################

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

##############################################################################################
#Parameter distributions
##############################################################################################

#Bring in 250 parameter sets
infile.parms <- "rakai_samples_n100_ahripfa"
infile.parms <- read.csv(paste0(infile.parms, '.csv'))
names(infile.parms)
names(infile.parms)
infile.parms <- infile.parms[,c(1:26)]

##############################################################################################
#Reference data
##############################################################################################

#prevalence reference data
infileprev <- "Rakai_prevalence_gender_R19_REVISION"
data_prevalence_gender1549 <- read.csv(paste0('./', infileprev, '.csv'))
data_prevalence_gender1549 <- data_prevalence_gender1549 %>% 
  mutate(gender=sex, gender=recode(gender, 'F'='Women', 'M'='Men'))

infileprev <- "Rakai_prevalence_preds_agecat_R19_REVISION"
prevalence_age_gender <- read.csv(paste0('./', infileprev, '.csv'))

prevalence_age_gender <- prevalence_age_gender %>% mutate(age_group = cut(as.numeric(substr(agecat, 2, 3)), 
                                                                          breaks = c(14, 19, 24, 29, 34, 39, 44, 49), 
                                                                          labels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")),
                                                          gender=recode(sex, 'Female'='Women', 'Male'='Men'))

infileincage <- "Rakai_incpredictions_agecat3_R19_REVISION"
incidence_age_gender <- read.csv(paste0('./', infileincage, '.csv'))
head(incidence_age_gender)

##############################################################################################
## Plots set up
##############################################################################################

load("scenario_comparison.RData")
names(outputs)

#Figure Arguments
colors=c(scales::viridis_pal()(10))
example_colorts = c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF","#B4DE2CFF", "#FDE725FF")

##############################################################################################
#Calibration plots
##############################################################################################
setwd("..\\plots")

#reshape to create ggplots of histograms for paramter set values
infile.parms$ART.Efficacy = 1-infile.parms$ART.Efficacy
ggplot(gather(infile.parms), aes(value)) + 
  geom_histogram(bins = 7, fill="white", col="black") + 
  facet_wrap(~key, scales = 'free_x',ncol=5) +
  xlab("Parameter value") +
  ylab("Count") +
  theme(strip.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

options(scipen = 999)
quants <- c(0.25,0.50,0.75)
parm.quantiles <- as.data.frame(apply( infile.parms[1:length(infile.parms)],2  , quantile , probs = quants , na.rm = TRUE ))
parm.quantiles.t <- as.data.frame(t(parm.quantiles))

print(
        ggplot(data=subset(outputs$prev_adults_gender %>%
                             group_by(scenario, gender) %>%
                             arrange(scenario, gender, year) %>%
                             mutate(loess_prev_mean = predict(loess(prev_mean ~ year, span=0.25)),
                                    loess_prev_lb = predict(loess(prev_lower ~ year, span=0.25)),
                                    loess_prev_ub = predict(loess(prev_upper ~ year, span=0.25))), year>1980 & scenario=="Baseline-campaign_Rakai_Baseline-FullRun2050" )
        ) +
        geom_line(aes(x=year, y=loess_prev_mean, color = factor(scenario)), linewidth=1) +
        geom_ribbon(aes(x=year,ymin = loess_prev_lb, ymax = loess_prev_ub, fill = factor(scenario)), alpha = 0.3, show.legend = FALSE) +
        geom_point(data=data_prevalence_gender1549, aes(x=year, y=hivprev*100)) +
        geom_errorbar(data=data_prevalence_gender1549, aes(x=year, ymin=lb*100, ymax=ub*100)) +
        scale_y_continuous(expand = c(0,0), limits=c(0,20)) +
        scale_x_continuous(breaks = seq(1980,2050,10), limits=c(1985, 2050)) +
        labs(x = "Year", y = "Prevalence (%)") +
        facet_wrap(~gender) +
        ggtitle("Prevalence among adults age 15-49") +
          theme_bw(base_size=24) +
          theme(legend.position="none", legend.title=element_blank()) +
          theme(strip.background = element_rect(colour="black", fill="white"))+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_rect(colour = "black"))
)

print(
  ggplot(data=subset(outputs$prev_age_gender %>%
                       group_by(scenario, age_group, gender) %>%
                       arrange(scenario, age_group, gender, year) %>%
                       mutate(loess_prev_mean = predict(loess(prev_mean ~ year, span=0.25)),
                              loess_prev_lb = predict(loess(prev_lower ~ year, span=0.25)),
                              loess_prev_ub = predict(loess(prev_upper ~ year, span=0.25))), year %in% seq(1999,2020,1) & scenario=="Baseline-campaign_Rakai_Baseline-FullRun2050" )
  ) +
    geom_line(aes(x=age_group, y=loess_prev_mean, color = gender, group=gender), linewidth=1) +
    geom_ribbon(aes(x=age_group,ymin = loess_prev_lb, ymax = loess_prev_ub, fill = gender, group=gender), alpha = 0.3, show.legend = FALSE) +
    geom_point(data=prevalence_age_gender, aes(x=age_group, y=prevalence*100, color=gender)) +
    geom_errorbar(data=prevalence_age_gender, aes(x=age_group, ymin=lb*100, ymax=ub*100, color=gender)) +
    scale_y_continuous(expand = c(0,0)) +
    #scale_x_continuous(breaks = seq(1980,2050,10), limits=c(1980, 2049),expand = c(0,0)) +
    labs(x = "Year", y = "Prevalence (%)") +
    facet_wrap(~year) +
    ggtitle("Prevalence by age, gender, and year") +
    theme_bw(base_size=24) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position="none", legend.title=element_blank()) +
    theme(strip.background = element_rect(colour="black", fill="white"))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black"))
)

#Incidence 15-19 and 20-24
outputs$inc_age_3_gender <- outputs$inc_age_3_gender %>% 
  mutate(Age=age_group_3, Sex=gender, Year=year)

modelled_inc_1524 = subset(outputs$inc_age_3_gender, Age=="15-24" & Sex=="Women")
data_inc_1524 = subset(incidence_age_gender, Age=="15-24" & Sex=="Women")

data <- subset(modelled_inc_1524 %>%
         group_by(scenario, gender, age_group_3) %>%
         arrange(scenario, gender, age_group_3, year) %>%
         mutate(loess_inc_mean = predict(loess(inc_mean ~ year, span=0.25)),
                loess_inc_lb = predict(loess(inc_lower ~ year, span=0.25)),
                loess_inc_ub = predict(loess(inc_upper ~ year, span=0.25))), year>1980 & (scenario=="Baseline-campaign_Rakai_Baseline-FullRun2050") )


print(subset(modelled_inc_1524 %>%
               group_by(scenario, gender, age_group_3) %>%
               arrange(scenario, gender, age_group_3, year) %>%
               mutate(loess_inc_mean = predict(loess(inc_mean ~ year, span=0.25)),
                      loess_inc_lb = predict(loess(inc_lower ~ year, span=0.25)),
                      loess_inc_ub = predict(loess(inc_upper ~ year, span=0.25))), year>1980 & (scenario=="Baseline-campaign_Rakai_Baseline-FullRun2050") ) %>%
        ggplot() +
        geom_point(data=data_inc_1524, aes(x=Year, y=Incidence)) +
        geom_errorbar(data=data_inc_1524, aes(x=Year, ymin=lb,ymax=ub)) +
        geom_line(aes(x=Year, y=loess_inc_mean, color = factor(scenario)), linewidth=1) +
        geom_ribbon(aes(x=Year, ymin = loess_inc_lb, ymax = loess_inc_ub, fill = factor(scenario)), alpha = 0.3, show.legend = FALSE) +
        scale_color_manual(guide = guide_legend(title = ""),values=colors[c(1,8,10,9,2)],labels = c("Baseline","No ART","No change sexual debut","No VMMC")) +
        scale_fill_manual(values=colors[c(1,8,10,9,2)]) +
        #scale_y_continuous(expand = c(0,0), breaks=seq(0,2,0.5), limits=c(0,2)) +
        #scale_x_continuous(breaks = seq(2000,2050,10), limits=c(1980, 2049),expand = c(0,0)) +
        labs(x = "Year", y = "Incidence (per 100 py)") +
        #facet_grid(Age~Sex) +
        ggtitle("Incidence among women age 15-24") +
        theme_bw(base_size=24) +
        theme(legend.position="bottom") +
        theme(strip.background = element_rect(colour="black", fill="white"))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black"))
)

#ART coverage
outputs$art_age_gender$artcov_mean[is.nan(outputs$art_age_gender$artcov_mean)]<-0
table(outputs$art_age_gender$scenario)
print(
  ggplot(data=subset(outputs$art_age_gender %>%
                       group_by(scenario, gender, age_group) %>%
                       arrange(scenario, gender, age_group, year) %>%
                       mutate(loess_prev_mean = predict(loess(artcov_mean ~ year, span=0.25)),
                              loess_prev_lb = predict(loess(artcov_lower ~ year, span=0.25)),
                              loess_prev_ub = predict(loess(artcov_upper ~ year, span=0.25))), year>1980 & scenario=="Baseline-campaign_Rakai_Baseline-FullRun2050")
  ) +
    geom_line(aes(x=year, y=loess_prev_mean, color = gender), linewidth=1) +
    geom_ribbon(aes(x=year,ymin = loess_prev_lb, ymax = loess_prev_ub, fill = gender), alpha = 0.3, show.legend = FALSE) +
    #geom_point(data=data_prevalence_gender1549, aes(x=year, y=hivprev*100)) +
    #geom_errorbar(data=data_prevalence_gender1549, aes(x=year, ymin=lb*100, ymax=ub*100)) +
    scale_y_continuous(expand = c(0,0), limits=c(0,100)) +
    scale_x_continuous(breaks = seq(1980,2050,10), limits=c(1985, 2050)) +
    labs(x = "Year", y = "% on ART") +
    facet_grid(scenario~age_group) +
    ggtitle("ART coverage") +
    theme_bw(base_size=24) +
    theme(legend.position="none", legend.title=element_blank()) +
    theme(strip.background = element_rect(colour="black", fill="white"))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black"))
)

##############################################################################################
#Scenario plots
##############################################################################################
library(dplyr)
library(ggplot2)
library(tidyr)

#Incidence

#Incidence by scenario 15-19 and 20-24
model_inc_1524 <- print(subset(outputs$inc_age_gender %>%
               group_by(scenario, gender, age_group) %>%
               arrange(scenario, gender, age_group, year) %>%
               mutate(loess_inc_mean = predict(loess(inc_mean ~ year, span=0.25)),
                      loess_inc_lb = predict(loess(inc_lower ~ year, span=0.25)),
                      loess_inc_ub = predict(loess(inc_upper ~ year, span=0.25))), (age_group=="15-19" | age_group=="20-24") & gender=="Women" & year>2000 & (scenario=="Baseline-campaign_Rakai_Baseline-FullRun2050" |
                                                                                                                          scenario=="Baseline-campaign_Rakai_Nochangesexualdebut-FullRun2050" | 
                                                                                                                          scenario=="Baseline-campaign_Rakai_NoVMMC-FullRun2050"|
                                                                                                                          scenario=="Baseline-campaign_Rakai_NoART-FullRun2050") ) %>%
        ggplot(aes(x = year, y = loess_inc_mean)) +
        geom_line(aes(color = factor(scenario)), linewidth=1) +
        scale_color_manual(guide = guide_legend(title = ""),values=colors[c(1,8,10,9)],labels = c("Baseline","No ART","No change sexual debut","No VMMC")) +
        labs(x = "Year", y = "Incidence (per 100 py)") +
        scale_x_continuous(breaks=seq(2000,2050,10), limits=c(2000,2050)) +
        facet_grid(~age_group) +
        ggtitle("Incidence among adults age 15-19") +
        theme_bw(base_size=24) +
        theme(legend.position="bottom") +
        theme(strip.background = element_rect(colour="black", fill="white"))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black"))
)


#Cumulative infections averted 15-19 and 20-24
table(outputs$inf_averted_age_gender$scenario)
scenarios <- subset(outputs$inf_averted_age_gender)
baseline <- subset(outputs$inf_averted_age_gender, scenario=="Baseline-campaign_Rakai_Baseline-FullRun2050")
scenarios_comb <- merge(baseline, scenarios, by=c("gender", "age_group", "year"))
scenarios_comb$delta_inc_mean <- scenarios_comb$cuminc_mean.y - scenarios_comb$cuminc_mean.x
scenarios_comb$delta_inc_mean_pct <- (scenarios_comb$cuminc_mean.y - scenarios_comb$cuminc_mean.x)/ scenarios_comb$cuminc_mean.y
data <- subset(scenarios_comb, (age_group=="15-19" | age_group== "20-24") & gender=="Women" &  (scenario.y=="Baseline-campaign_Rakai_Baseline-FullRun2050" |
                                                                                                  scenario.y=="Baseline-campaign_Rakai_Nochangesexualdebut-FullRun2050" | 
                                                                                                  scenario.y=="Baseline-campaign_Rakai_NoVMMC-FullRun2050"|
                                                                                                  scenario.y=="Baseline-campaign_Rakai_NoART-FullRun2050"))

print(subset(data, scenario.y!="Baseline-campaign_Rakai_Nointerventions-FullRun2050") %>%
        ggplot(aes(x = year, y = delta_inc_mean_pct*100)) +
        geom_line(aes(color = factor(scenario.y)), linewidth=1) +
        scale_color_manual(guide = guide_legend(title = ""),values=colors[c(1,8,10,9)],labels = c("Baseline","ART","Increased AFS", "VMMC")) +
        labs(x = "Year", y = "Cumulative infections averted (%)") +
        scale_x_continuous(breaks=seq(2000,2050,10), limits=c(2000,2050)) +
        scale_y_continuous(limits=c(0,50)) +
        facet_grid(~age_group) +
        ggtitle("Cumulative infections averted") +
        theme_bw(base_size=24) +
        theme(legend.position="bottom") +
        theme(strip.background = element_rect(colour="black", fill="white"))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black"))
)

inc_data <- subset(outputs$inc_age_gender %>%
         group_by(scenario, gender, age_group) %>%
         arrange(scenario, gender, age_group, year) %>%
         mutate(loess_inc_mean = predict(loess(inc_mean ~ year, span=0.25)),
                loess_inc_lb = predict(loess(inc_lower ~ year, span=0.25)),
                loess_inc_ub = predict(loess(inc_upper ~ year, span=0.25))), (age_group=="15-19" | age_group=="20-24") & gender=="Women" & year>2000 & (scenario=="Baseline-campaign_Rakai_Baseline-FullRun2050" |
                                                                                                                                                          scenario=="Baseline-campaign_Rakai_Nochangesexualdebut-FullRun2050" | 
                                                                                                                                                          scenario=="Baseline-campaign_Rakai_NoVMMC-FullRun2050"| scenario=="Baseline-campaign_Rakai_NoART-FullRun2050") ) %>%
  select(age_group, gender, year, loess_inc_mean, scenario)

cum_inc_data <- subset(data, scenario.y!="Baseline-campaign_Rakai_Nointerventions-FullRun2050") %>%
  select(gender, age_group, year, delta_inc_mean_pct, scenario.y)

# 1) Harmonize column names first
inc_aligned <- inc_data %>%
  rename(scenario = scenario) %>%       # already 'scenario', this is a no-op for clarity
  mutate(delta_inc_mean_pct = NA_real_) # add missing col so rows bind cleanly

cum_aligned <- cum_inc_data %>%
  rename(scenario = scenario.y) %>%     # make 'scenario' consistent
  mutate(loess_inc_mean = NA_real_) %>% # add missing col so rows bind cleanly
  mutate(delta_inc_mean_pct = 100*delta_inc_mean_pct)

# 2) Row-bind, then pivot to long and rename measures
combined_long <- bind_rows(inc_aligned, cum_aligned) %>%
  pivot_longer(
    cols = c(loess_inc_mean, delta_inc_mean_pct),
    names_to = "measure",
    values_to = "value",
    values_drop_na = TRUE
  ) %>%
  mutate(
    measure = recode(measure,
                     loess_inc_mean = "incidence",
                     delta_inc_mean_pct = "cum_incidence")
  ) %>%
  select(age_group, gender, year, scenario, measure, value)

combined_long

print(subset(combined_long) %>%
        ggplot(aes(x = year, y = value)) +
        geom_line(aes(color = factor(scenario)), linewidth=1) +
        scale_color_manual(guide = guide_legend(title = ""),values=colors[c(1,8,10,9)],labels = c("Baseline","ART","Increased AFS", "VMMC")) +
        labs(x = "Year", y = "metric (%)") +
        scale_x_continuous(breaks=seq(2000,2050,10), limits=c(2000,2050)) +
        #scale_y_continuous(limits=c(0,50)) +
        facet_grid(measure~age_group, scales = "free_y") +
        ggtitle("") +
        theme_bw(base_size=24) +
        theme(legend.position="bottom") +
        theme(strip.background = element_rect(colour="black", fill="white"))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black"))
)

# Make sure facets ordered: incidence top, cumulative below
combined_long <- combined_long %>%
  mutate(measure = factor(measure, levels = c("incidence", "cum_incidence")))

p <-ggplot(combined_long, aes(x = year, y = value)) +
  geom_line(aes(color = factor(scenario)), linewidth = 1) +
  scale_color_manual(
    guide = guide_legend(title = ""),
    values = colors[c(1,8,10,9)],
    labels = c("Baseline", "ART", "Increased AFS", "VMMC")
  ) +
  labs(x = "Year", y = NULL) +
  scale_x_continuous(
    breaks = seq(2000, 2050, 10),
    limits = c(2000, 2050)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, NA)
  ) +
  facet_grid(
    rows = vars(measure),
    cols = vars(age_group),
    scales = "free_y",
    switch = "y",
    labeller = labeller(
      measure = c(
        cum_incidence = "Cumulative\nIncidence (%)",
        incidence = "Incidence (per 100 py)"
      )
    )
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "bottom",
    # remove strip box entirely (no border, no fill)
    strip.background = element_blank(),
    # keep vertical labels
    strip.text.y.left = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    strip.placement = "outside",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.background = element_blank(),
    axis.line = element_blank()
  )


