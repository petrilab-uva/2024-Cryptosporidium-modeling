#Title: Mixed Effects Modeling of Antibody Titers
#Name: Hannah So
#Date: 3/25/24

####Introduction###############################################################
#Calculating the mixed effects models for the fecal and antibody titers. 
#File involves formatting data to be fit for analysis. Mixed effects models were
#run and a data visualization was made from that data. 

####Load Packages and Data######################################################
#install.packages("tidyverse")
packageVersion("tidyverse") #'2.0.0'
library(tidyverse)

#install.packages("lmerTest")
packageVersion("lmerTest") #‘3.1.3’
library(lmerTest)

#install.packages("Matrix")
packageVersion("Matrix") # 1.6.5


Plasma.Antibody.Metadata <-read.csv("../Results/Plasma_Antibody_Metadata.csv")
Fecal.Antibody.Metadata <-read.csv("../Results/Fecal_Antibody_Metadata.csv")


####Formatting#########################################
#Check structure of data frames
#Fecal
str(Fecal.Antibody.Metadata)

Fecal.Antibody.Factored <- Fecal.Antibody.Metadata %>%
  dplyr::mutate(AgeInYears = AgeInDays/365.25) %>%
  dplyr::mutate(SID = factor(SID)) %>%
  dplyr::mutate(Season = factor(Season, levels = c("Dry", "Rainy"))) %>%
  dplyr::mutate(sex = factor(sex, levels = c("MALE", "FEMALE"))) %>%
  dplyr::mutate(DRAIN = factor(DRAIN, levels = c("2", "1"))) %>% #Open drain beside house, yes = 1, no = 2
  dplyr::select(-AgeInDays)

str(Fecal.Antibody.Factored) # Factor variables now correct.


#Plasma
str(Plasma.Antibody.Metadata)

Plasma.Antibody.Factored <- Plasma.Antibody.Metadata %>%
  dplyr::mutate(AgeInYears = AgeInDays/365.25) %>%
  dplyr::mutate(SID = factor(SID)) %>%
  dplyr::mutate(Season = factor(Season, levels = c("Dry", "Rainy"))) %>%
  dplyr::mutate(sex = factor(sex, levels = c("MALE", "FEMALE"))) %>%
  dplyr::mutate(DRAIN = factor(DRAIN, levels = c("2", "1"))) %>% #Open drain beside house, yes = 1, no = 2
  dplyr::select(-AgeInDays)

str(Plasma.Antibody.Factored) # Factor variables now correct.



####Mixed Effects Model - Fecal################################################

#1. Cp23 

lmer23 <- lmerTest::lmer(Cp23.RU ~ AgeInYears + Season + INCOME + sex + HAZ + DRAIN + 
                 EpisodeNumber + (1|SID), 
               data = Fecal.Antibody.Factored, REML = FALSE)
summary(lmer23)
confint(lmer23)

#Convert output into data frame
Cp23_fecal_summary <- as.data.frame(summary(lmer23)$coefficients[, c(1,2,5)]) 
Cp23_fecal_confint <- as.data.frame(confint(lmer23)[c(3:10), ])

#Combine Summary and Confidence Intervals
Cp23_fecal_summary <- tibble::rownames_to_column(Cp23_fecal_summary, "Features")
Cp23_fecal_confint <- tibble::rownames_to_column(Cp23_fecal_confint, "Features")
Cp23_fecal_output <-  dplyr::full_join(Cp23_fecal_summary, Cp23_fecal_confint, by = join_by(Features))

#Organizing Dataframe
Cp23_fecal_formatted <- Cp23_fecal_output %>%
  dplyr::filter(Features != "(Intercept)")%>%
  dplyr::rename(pvalue = "Pr(>|t|)") %>%
  dplyr::mutate(Significant = ifelse(pvalue < 0.05, "Yes", "No")) %>%
  dplyr::mutate(Significant = factor(Significant, levels = c("No", "Yes")))


#Change Features Names
#For dichotomous variables, the change is for the factor indicated in parentheses 
Cp23_fecal_formatted$Features[Cp23_fecal_formatted$Features == "AgeInYears"] <- "Age in Years"
Cp23_fecal_formatted$Features[Cp23_fecal_formatted$Features == "SeasonRainy"] <- "Season (Monsoon)"
Cp23_fecal_formatted$Features[Cp23_fecal_formatted$Features == "INCOME"] <- "Income"
Cp23_fecal_formatted$Features[Cp23_fecal_formatted$Features == "sexFEMALE"] <- "Sex (Female)"
Cp23_fecal_formatted$Features[Cp23_fecal_formatted$Features == "HAZ"] <- "HAZ Score"
Cp23_fecal_formatted$Features[Cp23_fecal_formatted$Features == "DRAIN1"] <- "Drain (Present)"
Cp23_fecal_formatted$Features[Cp23_fecal_formatted$Features == "EpisodeNumber"] <- "Episode Number"

Cp23_fecal_formatted$Features <- factor(Cp23_fecal_formatted$Features, 
                                        levels = c("Sex (Female)", "Income", "Drain (Present)", "HAZ Score", "Episode Number",  "Season (Monsoon)", "Age in Years"))

ggplot(Cp23_fecal_formatted, aes(x= Estimate, y = Features, color = Significant)) +
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_linerange(aes(xmax = `97.5 %`, xmin = `2.5 %`)) +
  geom_point(size = 4) +
  labs(y = NULL) +
  theme_bw() +
  theme(text = element_text(size = 25))

#ggsave("../Results/Figures/FecalCp23IgA_MEM.png", dpi = 300, width = 7.4)



#2. Cp17

lmer17 <- lmerTest::lmer(Cp17.RU ~ AgeInYears + Season + INCOME + sex + HAZ + DRAIN + 
                           EpisodeNumber + (1|SID), 
                         data = Fecal.Antibody.Factored, REML = FALSE)
summary(lmer17)
confint(lmer17)

#Convert output into data frame
Cp17_fecal_summary <- as.data.frame(summary(lmer17)$coefficients[, c(1,2,5)]) 
Cp17_fecal_confint <- as.data.frame(confint(lmer17)[c(3:10), ])

#Combine Summary and Confidence Intervals
Cp17_fecal_summary <- tibble::rownames_to_column(Cp17_fecal_summary, "Features")
Cp17_fecal_confint <- tibble::rownames_to_column(Cp17_fecal_confint, "Features")
Cp17_fecal_output <-  dplyr::full_join(Cp17_fecal_summary, Cp17_fecal_confint, by = join_by(Features))

#Organizing Dataframe
Cp17_fecal_formatted <- Cp17_fecal_output %>%
  dplyr::filter(Features != "(Intercept)")%>%
  dplyr::rename(pvalue = "Pr(>|t|)") %>%
  dplyr::mutate(Significant = ifelse(pvalue < 0.05, "Yes", "No")) %>%
  dplyr::mutate(Significant = factor(Significant, levels = c("No", "Yes")))


#Change Features Names
#For dichotomous variables, the change is for the factor indicated in parentheses 
Cp17_fecal_formatted$Features[Cp17_fecal_formatted$Features == "AgeInYears"] <- "Age in Years"
Cp17_fecal_formatted$Features[Cp17_fecal_formatted$Features == "SeasonRainy"] <- "Season (Monsoon)"
Cp17_fecal_formatted$Features[Cp17_fecal_formatted$Features == "INCOME"] <- "Income"
Cp17_fecal_formatted$Features[Cp17_fecal_formatted$Features == "sexFEMALE"] <- "Sex (Female)"
Cp17_fecal_formatted$Features[Cp17_fecal_formatted$Features == "HAZ"] <- "HAZ Score"
Cp17_fecal_formatted$Features[Cp17_fecal_formatted$Features == "DRAIN1"] <- "Drain (Present)"
Cp17_fecal_formatted$Features[Cp17_fecal_formatted$Features == "EpisodeNumber"] <- "Episode Number"

Cp17_fecal_formatted$Features <- factor(Cp17_fecal_formatted$Features, 
                                        levels = c("Sex (Female)", "Income", "Drain (Present)", "HAZ Score", "Episode Number",  "Season (Monsoon)", "Age in Years"))

ggplot(Cp17_fecal_formatted, aes(x= Estimate, y = Features, color = Significant)) +
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_linerange(aes(xmax = `97.5 %`, xmin = `2.5 %`)) +
  geom_point(size = 4) +
  labs(y = NULL) +
  theme_bw() +
  theme(text = element_text(size = 25))

#ggsave("../Results/Figures/FecalCp17IgA_MEM.png", dpi = 300, width = 7.4)



####Mixed Effects Model - Plasma#######################

#1. Cp23 IgG

lmerP23IgG <- lmerTest::lmer(Cp23.IgG.Plasma ~ AgeInYears + Season + INCOME + sex + HAZ + DRAIN + 
                               EpisodeNumber + (1|SID), 
                             data = Plasma.Antibody.Factored, REML = FALSE)
summary(lmerP23IgG)
confint(lmerP23IgG)

#Convert output into data frame
Cp23IgG_plasma_summary <- as.data.frame(summary(lmerP23IgG)$coefficients[, c(1,2,5)]) 
Cp23IgG_plasma_confint <- as.data.frame(confint(lmerP23IgG)[c(3:10), ])

#Combine Summary and Confidence Intervals
Cp23IgG_plasma_summary <- tibble::rownames_to_column(Cp23IgG_plasma_summary, "Features")
Cp23IgG_plasma_confint <- tibble::rownames_to_column(Cp23IgG_plasma_confint, "Features")
Cp23IgG_plasma_output <-  dplyr::full_join(Cp23IgG_plasma_summary, Cp23IgG_plasma_confint, by = join_by(Features))

#Organizing Dataframe
Cp23IgG_plasma_formatted <- Cp23IgG_plasma_output %>%
  dplyr::filter(Features != "(Intercept)")%>%
  dplyr::rename(pvalue = "Pr(>|t|)") %>%
  dplyr::mutate(Significant = ifelse(pvalue < 0.05, "Yes", "No")) %>%
  dplyr::mutate(Significant = factor(Significant, levels = c("No", "Yes")))


#Change Features Names
#For dichotomous variables, the change is for the factor indicated in parentheses 
Cp23IgG_plasma_formatted$Features[Cp23IgG_plasma_formatted$Features == "AgeInYears"] <- "Age in Years"
Cp23IgG_plasma_formatted$Features[Cp23IgG_plasma_formatted$Features == "SeasonRainy"] <- "Season (Monsoon)"
Cp23IgG_plasma_formatted$Features[Cp23IgG_plasma_formatted$Features == "INCOME"] <- "Income"
Cp23IgG_plasma_formatted$Features[Cp23IgG_plasma_formatted$Features == "sexFEMALE"] <- "Sex (Female)"
Cp23IgG_plasma_formatted$Features[Cp23IgG_plasma_formatted$Features == "HAZ"] <- "HAZ Score"
Cp23IgG_plasma_formatted$Features[Cp23IgG_plasma_formatted$Features == "DRAIN1"] <- "Drain (Present)"
Cp23IgG_plasma_formatted$Features[Cp23IgG_plasma_formatted$Features == "EpisodeNumber"] <- "Episode Number"

Cp23IgG_plasma_formatted$Features <- factor(Cp23IgG_plasma_formatted$Features, 
                                            levels = c("Sex (Female)", "Income", "Drain (Present)", "HAZ Score", "Episode Number",  "Season (Monsoon)", "Age in Years"))

ggplot(Cp23IgG_plasma_formatted, aes(x= Estimate, y = Features, color = Significant)) +
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_linerange(aes(xmax = `97.5 %`, xmin = `2.5 %`)) +
  geom_point(size = 4) +
  labs(y = NULL) +
  theme_bw() +
  theme(text = element_text(size = 25))

#ggsave("../Results/Figures/PlasmaCp23IgG_MEM.png", dpi = 300, width = 7.4)




#2. Cp23 IgA

lmerP23IgA <- lmerTest::lmer(Cp23.IgA.Plasma ~ AgeInYears + Season + INCOME + sex + HAZ + DRAIN + 
                               EpisodeNumber + (1|SID), 
                             data = Plasma.Antibody.Factored, REML = FALSE)
summary(lmerP23IgA)
confint(lmerP23IgA)

#Convert output into data frame
Cp23IgA_plasma_summary <- as.data.frame(summary(lmerP23IgA)$coefficients[, c(1,2,5)]) 
Cp23IgA_plasma_confint <- as.data.frame(confint(lmerP23IgA)[c(3:10), ])

#Combine Summary and Confidence Intervals
Cp23IgA_plasma_summary <- tibble::rownames_to_column(Cp23IgA_plasma_summary, "Features")
Cp23IgA_plasma_confint <- tibble::rownames_to_column(Cp23IgA_plasma_confint, "Features")
Cp23IgA_plasma_output <-  dplyr::full_join(Cp23IgA_plasma_summary, Cp23IgA_plasma_confint, by = join_by(Features))

#Organizing Dataframe
Cp23IgA_plasma_formatted <- Cp23IgA_plasma_output %>%
  dplyr::filter(Features != "(Intercept)")%>%
  dplyr::rename(pvalue = "Pr(>|t|)") %>%
  dplyr::mutate(Significant = ifelse(pvalue < 0.05, "Yes", "No")) %>%
  dplyr::mutate(Significant = factor(Significant, levels = c("No", "Yes")))


#Change Features Names
#For dichotomous variables, the change is for the factor indicated in parentheses 
Cp23IgA_plasma_formatted$Features[Cp23IgA_plasma_formatted$Features == "AgeInYears"] <- "Age in Years"
Cp23IgA_plasma_formatted$Features[Cp23IgA_plasma_formatted$Features == "SeasonRainy"] <- "Season (Monsoon)"
Cp23IgA_plasma_formatted$Features[Cp23IgA_plasma_formatted$Features == "INCOME"] <- "Income"
Cp23IgA_plasma_formatted$Features[Cp23IgA_plasma_formatted$Features == "sexFEMALE"] <- "Sex (Female)"
Cp23IgA_plasma_formatted$Features[Cp23IgA_plasma_formatted$Features == "HAZ"] <- "HAZ Score"
Cp23IgA_plasma_formatted$Features[Cp23IgA_plasma_formatted$Features == "DRAIN1"] <- "Drain (Present)"
Cp23IgA_plasma_formatted$Features[Cp23IgA_plasma_formatted$Features == "EpisodeNumber"] <- "Episode Number"

Cp23IgA_plasma_formatted$Features <- factor(Cp23IgA_plasma_formatted$Features, 
                                            levels = c("Sex (Female)", "Income", "Drain (Present)", "HAZ Score", "Episode Number",  "Season (Monsoon)", "Age in Years"))

ggplot(Cp23IgA_plasma_formatted, aes(x= Estimate, y = Features, color = Significant)) +
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_linerange(aes(xmax = `97.5 %`, xmin = `2.5 %`)) +
  geom_point(size = 4) +
  labs(y = NULL) +
  theme_bw() +
  theme(text = element_text(size = 25))

#ggsave("../Results/Figures/PlasmaCp23IgA_MEM.png", dpi = 300, width = 7.4)



#3. Cp17 IgG 

lmerP17IgG <- lmerTest::lmer(Cp17.IgG.Plasma ~ AgeInYears + Season + INCOME + sex + HAZ + DRAIN + 
                               EpisodeNumber + (1|SID), 
                             data = Plasma.Antibody.Factored, REML = FALSE)
summary(lmerP17IgG)
confint(lmerP17IgG)

#Convert output into data frame
Cp17IgG_plasma_summary <- as.data.frame(summary(lmerP17IgG)$coefficients[, c(1,2,5)]) 
Cp17IgG_plasma_confint <- as.data.frame(confint(lmerP17IgG)[c(3:10), ])

#Combine Summary and Confidence Intervals
Cp17IgG_plasma_summary <- tibble::rownames_to_column(Cp17IgG_plasma_summary, "Features")
Cp17IgG_plasma_confint <- tibble::rownames_to_column(Cp17IgG_plasma_confint, "Features")
Cp17IgG_plasma_output <-  dplyr::full_join(Cp17IgG_plasma_summary, Cp17IgG_plasma_confint, by = join_by(Features))

#Organizing Dataframe
Cp17IgG_plasma_formatted <- Cp17IgG_plasma_output %>%
  dplyr::filter(Features != "(Intercept)")%>%
  dplyr::rename(pvalue = "Pr(>|t|)") %>%
  dplyr::mutate(Significant = ifelse(pvalue < 0.05, "Yes", "No")) %>%
  dplyr::mutate(Significant = factor(Significant, levels = c("No", "Yes")))


#Change Features Names
#For dichotomous variables, the change is for the factor indicated in parentheses 
Cp17IgG_plasma_formatted$Features[Cp17IgG_plasma_formatted$Features == "AgeInYears"] <- "Age in Years"
Cp17IgG_plasma_formatted$Features[Cp17IgG_plasma_formatted$Features == "SeasonRainy"] <- "Season (Monsoon)"
Cp17IgG_plasma_formatted$Features[Cp17IgG_plasma_formatted$Features == "INCOME"] <- "Income"
Cp17IgG_plasma_formatted$Features[Cp17IgG_plasma_formatted$Features == "sexFEMALE"] <- "Sex (Female)"
Cp17IgG_plasma_formatted$Features[Cp17IgG_plasma_formatted$Features == "HAZ"] <- "HAZ Score"
Cp17IgG_plasma_formatted$Features[Cp17IgG_plasma_formatted$Features == "DRAIN1"] <- "Drain (Present)"
Cp17IgG_plasma_formatted$Features[Cp17IgG_plasma_formatted$Features == "EpisodeNumber"] <- "Episode Number"

Cp17IgG_plasma_formatted$Features <- factor(Cp17IgG_plasma_formatted$Features, 
                                            levels = c("Sex (Female)", "Income", "Drain (Present)", "HAZ Score", "Episode Number",  "Season (Monsoon)", "Age in Years"))

ggplot(Cp17IgG_plasma_formatted, aes(x= Estimate, y = Features, color = Significant)) +
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_linerange(aes(xmax = `97.5 %`, xmin = `2.5 %`)) +
  geom_point(size = 4) +
  labs(y = NULL) +
  theme_bw() +
  theme(text = element_text(size = 25))

#ggsave("../Results/Figures/PlasmaCp17IgG_MEM.png", dpi = 300, width = 7.4)



#4. Cp17 IgA 

lmerP17IgA <- lmerTest::lmer(Cp17.IgA.Plasma ~ AgeInYears + Season + INCOME + sex + HAZ + DRAIN + 
                               EpisodeNumber + (1|SID), 
                             data = Plasma.Antibody.Factored, REML = FALSE)
summary(lmerP17IgA)
confint(lmerP17IgA)

#Convert output into data frame
Cp17IgA_plasma_summary <- as.data.frame(summary(lmerP17IgA)$coefficients[, c(1,2,5)]) 
Cp17IgA_plasma_confint <- as.data.frame(confint(lmerP17IgA)[c(3:10), ])

#Combine Summary and Confidence Intervals
Cp17IgA_plasma_summary <- tibble::rownames_to_column(Cp17IgA_plasma_summary, "Features")
Cp17IgA_plasma_confint <- tibble::rownames_to_column(Cp17IgA_plasma_confint, "Features")
Cp17IgA_plasma_output <-  dplyr::full_join(Cp17IgA_plasma_summary, Cp17IgA_plasma_confint, by = join_by(Features))

#Organizing Dataframe
Cp17IgA_plasma_formatted <- Cp17IgA_plasma_output %>%
  dplyr::filter(Features != "(Intercept)")%>%
  dplyr::rename(pvalue = "Pr(>|t|)") %>%
  dplyr::mutate(Significant = ifelse(pvalue < 0.05, "Yes", "No")) %>%
  dplyr::mutate(Significant = factor(Significant, levels = c("No", "Yes")))


#Change Features Names
#For dichotomous variables, the change is for the factor indicated in parentheses 
Cp17IgA_plasma_formatted$Features[Cp17IgA_plasma_formatted$Features == "AgeInYears"] <- "Age in Years"
Cp17IgA_plasma_formatted$Features[Cp17IgA_plasma_formatted$Features == "SeasonRainy"] <- "Season (Monsoon)"
Cp17IgA_plasma_formatted$Features[Cp17IgA_plasma_formatted$Features == "INCOME"] <- "Income"
Cp17IgA_plasma_formatted$Features[Cp17IgA_plasma_formatted$Features == "sexFEMALE"] <- "Sex (Female)"
Cp17IgA_plasma_formatted$Features[Cp17IgA_plasma_formatted$Features == "HAZ"] <- "HAZ Score"
Cp17IgA_plasma_formatted$Features[Cp17IgA_plasma_formatted$Features == "DRAIN1"] <- "Drain (Present)"
Cp17IgA_plasma_formatted$Features[Cp17IgA_plasma_formatted$Features == "EpisodeNumber"] <- "Episode Number"

Cp17IgA_plasma_formatted$Features <- factor(Cp17IgA_plasma_formatted$Features, 
                                            levels = c("Sex (Female)", "Income", "Drain (Present)", "HAZ Score", "Episode Number",  "Season (Monsoon)", "Age in Years"))

ggplot(Cp17IgA_plasma_formatted, aes(x= Estimate, y = Features, color = Significant)) +
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_linerange(aes(xmax = `97.5 %`, xmin = `2.5 %`)) +
  geom_point(size = 4) +
  labs(y = NULL) +
  theme_bw() +
  theme(text = element_text(size = 25))

#ggsave("../Results/Figures/PlasmaCp17IgA_MEM.png", dpi = 300, width = 7.4)



####Correlations of Antibody and Age###########################################

#Check correlations between fecal antibody levels and age, as these were 
#significantly associated in the Mixed Effects Model.

#Fecal Data
#Cp17
ggplot(Fecal.Antibody.Factored, aes(AgeInYears, Cp17.RU)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Fecal anti-Cp17 IgA (RU)", x = "Age (Years)") +
  scale_x_continuous(breaks = seq(365, 1460, by = 365)) + 
  theme_classic() + 
  theme(text = element_text(size = 25)) +
  geom_hline(yintercept=10.06, linetype="dashed", linewidth = 1, color = "red")

cor.test(Fecal.Antibody.Factored$AgeInYears, Fecal.Antibody.Factored$Cp17.RU)

#ggsave("../Results/Figures/FecalCp17CorrelationPlot.png", dpi = 300, width = 6)


#Cp23
ggplot(Fecal.Antibody.Factored, aes(AgeInYears, Cp23.RU)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Fecal anti-Cp23 IgA (RU)", x = "Age (Years)") +
  scale_x_continuous(breaks = seq(365, 1460, by = 365)) + 
  theme_classic() + 
  theme(text = element_text(size = 25)) +
  geom_hline(yintercept=10.06, linetype="dashed", linewidth = 1, color = "red")

cor.test(Fecal.Antibody.Factored$AgeInYears, Fecal.Antibody.Factored$Cp23.RU)
#ggsave("../Results/Figures/FecalCp23CorrelationPlot.png", dpi = 300, width = 6)


