#Title: [4] Calculating Proportion of 1 vs 2+ PCR Positive Episodes
#Name: G. Brett Moreau
#Date: 5/13/2024

####Introduction################################################################
# One remaining question regarding the antibody data is whether repeated 
# Cryptosporidium-positive episodes are associated with shorter episodes. I will 
# be using 1 PCR positive vs 2+ PCR positive samples within an episode as my 
# comparison. Because surveillance sample collection is arbitrarily located 
# within an episode, episode length as measured from first positive to first
# negative can be arbitrary and misleading. Instead, I will compare episode 
# length by looking at episodes that had only one PCR positive sample versus 
# those with multiple PCR positive episodes.


####Load Data and Packages#####################################################

#install.packages("tidyverse")
library(tidyverse)
packageVersion("tidyverse") #'2.0.0'

#install.packages("readxl") 
library(readxl)
packageVersion("readxl") #‘1.4.2’

#install.packages("lmerTest")
library(lmerTest)
packageVersion("lmerTest") #‘3.1.3’


Crypto_episodes <- read.csv("../Results/FecalSample_CryptoEpisodeNumber.csv")




####Recategorize According to 1 or 2+ Positive PCRs per Episode#################


# Filter episode table to include only Group 1 individuals during the first 4yrs.
Crypto_episodes_grp1 <- Crypto_episodes %>%
  dplyr::filter(SID >=2059 & SID <=2250) %>%
  dplyr::filter(AgeInDays <= 1461)

length(unique(Crypto_episodes_grp1$SID)) # 192 SIDs, as expected.




####Evaluate Subclinical Episode Duration Over Repeated Episodes################

# We're focusing on subclinical episodes only, as diarrheal episodes are 
# collected more frequently and more likely to be of longer duration regardless.
# However, the total number of episodes (including both diarrheal and 
# subclinical) may influence episode duration. Therefore I'll use the absolute
# episode number, but filter to look only at subclinical episode duration.


# Group episodes into 1 or 2+ PCR positive episodes.
positives_per_episode_subclinical_total <- Crypto_episodes_grp1 %>%
  dplyr::ungroup() %>%
  dplyr::mutate(SID = factor(SID)) %>%
  dplyr::mutate(EpisodeNumber = factor(EpisodeNumber)) %>%
  dplyr::mutate(CryptoStatus = factor(CryptoStatus)) %>%
  dplyr::filter(EpisodeType == "Subclinical") %>%
  dplyr::group_by(SID, EpisodeNumber, .drop = FALSE) %>%
  dplyr::count(CryptoStatus) %>%
  dplyr::mutate(Positive_per_Episode_Category = ifelse(n == 0, "0", 
                                              ifelse(n == 1, "1", "2+"))) %>%
  plyr::rename(replace = c("n" = "Positive_per_Episode")) %>%
  dplyr::filter(CryptoStatus == "POSITIVE") %>%
  dplyr::select(-CryptoStatus) %>%
  dplyr::filter(Positive_per_Episode != "0") %>%
  dplyr::mutate(Positive_per_Episode_Category = factor(Positive_per_Episode_Category, levels = c("1", "2+"))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(EpisodeNumber, .drop = FALSE) %>%
  dplyr::count(Positive_per_Episode_Category) %>%
  dplyr::filter(EpisodeNumber != 0)

str(positives_per_episode_subclinical_total)

# Organize data set in wide format
positives_per_episode_subclinical_total_wide <- tidyr::pivot_wider(positives_per_episode_subclinical_total, 
                                                                   names_from = Positive_per_Episode_Category, 
                                                                   values_from = n)


# Convert absolute counts to proportions
positives_per_episode_subclinical_prop <- positives_per_episode_subclinical_total %>%
  dplyr::group_by(EpisodeNumber) %>%
  dplyr::mutate(Total_Episodes = sum(n)) %>%
  dplyr::mutate(n_prop = n / Total_Episodes * 100) %>%
  dplyr::mutate(Positive_per_Episode_Category = factor(Positive_per_Episode_Category, 
                                              levels = c("1", "2+")))


# Visualize change in episode duration over time.
ggplot(positives_per_episode_subclinical_prop, aes(x = EpisodeNumber, 
                                              y = n_prop, 
                                              fill = Positive_per_Episode_Category)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Episode Number", 
       y = "Percent Subclinical Episodes", 
       fill = "Duration (# of PCR Positive Samples)") +
  scale_fill_manual(values = c("blue", "white")) +
  theme_classic() +
  theme(legend.position = "bottom", 
        axis.title.x = element_text(face = "bold"),  
        axis.title.y = element_text(face = "bold"),  
        text = element_text(size = 18))

#ggsave("../Results/Figures/PCR_positive_subclinical.png", width = 7, height = 6, dpi = 300)

#write.csv(positives_per_episode_subclinical_total_wide, 
#          file = "../Results/PCR positve samples per episode_subclinical_total.csv", 
#          row.names = FALSE)
  

  

####Evaluate Episode Duration on a Per-SID Basis################################

# The previous code looked at the proportion of longer (2+ qPCR-positive samples) 
# episodes. However, I want to analyze this on a per-subject basis, which would
# allow for individual to be incorporated into a mixed effect model as a random
# effect. This will allow us to test whether there is a statistical association
# between episode duration (number of qPCR-positive samples per episode) and
# episode number.


# Group episodes into 1 or 2+ PCR positive episodes.
positives_per_episode <- Crypto_episodes_grp1 %>%
  dplyr::ungroup() %>%
  dplyr::mutate(SID = factor(SID)) %>%
  dplyr::mutate(EpisodeNumber = factor(EpisodeNumber)) %>%
  dplyr::mutate(CryptoStatus = factor(CryptoStatus)) %>%
  dplyr::filter(EpisodeType == "Subclinical") %>%
  dplyr::group_by(SID, EpisodeNumber, .drop = FALSE) %>%
  dplyr::count(CryptoStatus) %>%
  dplyr::mutate(Positive_per_Episode_Category = ifelse(n == 0, "0", 
                                                       ifelse(n == 1, "1", "2+"))) %>%
  plyr::rename(replace = c("n" = "Positive_per_Episode")) %>%
  dplyr::filter(CryptoStatus == "POSITIVE") %>%
  dplyr::filter(Positive_per_Episode != "0") %>%
  dplyr::select(-CryptoStatus) %>%
  dplyr::mutate(EpisodeNumber = as.character(EpisodeNumber)) %>%
  dplyr::mutate(EpisodeNumber = as.numeric(EpisodeNumber)) %>%
  ungroup() %>%
  as.data.frame()

#write.csv(positives_per_episode, 
#          file = "../Results/qPCR positives per episode_Subclinical.csv", 
#          row.names = FALSE)

str(positives_per_episode)


# Perform mixed effects model
lmer <- lmerTest::lmer(Positive_per_Episode ~ EpisodeNumber + (1|SID), 
                       data = positives_per_episode, REML = FALSE)

summary(lmer) # Significant association (p=0.00017)


# Add HAZ and WHZ to episode number model.
anthro <- read.csv("../Raw_data/Metadata_ForIndividuals.csv")
anthro <- dplyr::mutate(anthro, SID = factor(SID))

positives_per_episode_anthro <- dplyr::left_join(positives_per_episode, anthro, by = "SID")

# Perform mixed effects model
lmer_anthro <- lmerTest::lmer(Positive_per_Episode ~ EpisodeNumber + HAZ + WHZ + (1|SID), 
                       data = positives_per_episode_anthro, REML = FALSE)

summary(lmer_anthro) # Remains significant when controlling for HAZ and WHZ.

