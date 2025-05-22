#Title: [5] Association Between Cryptosporidium Burden and Episode Number
#Name: G. Brett Moreau
#Date: 6/17/2024

####Introduction################################################################

# This code evaluates the association between Cryptosporidium burden (as defined
# by Crypt Cq values) and episode number.


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


#### Organize Crypto Episode and Burden Data ###################################

# Import episode number information.Crypto burden values.
Crypto_episodes_all <- read.csv("../Results/FecalSample_CryptoEpisodeNumber.csv")


# Filter episode table to include only Group 1 individuals during the first 4yrs.
Crypto_episodes_grp1 <- Crypto_episodes_all %>%
  dplyr::filter(SID >=2059 & SID <=2250) %>%
  dplyr::filter(AgeInDays <= 1461)

length(unique(Crypto_episodes_grp1$SID)) # 192 SIDs, as expected.


# Filter to remove all sample except for the first positive of each episode.
Crypto_first_pos <- Crypto_episodes_grp1 %>%
  dplyr::filter(CryptoStatus == "POSITIVE") %>%
  dplyr::group_by(SID, EpisodeNumber) %>%
  dplyr::slice_min(AgeInDays) %>%
  dplyr::mutate(SID = factor(SID))


# Perform mixed effects model for association.

Episode_CT <- lmer(EpisodeNumber ~ CRYCTRT + (1|SID), data = Crypto_first_pos)

summary(Episode_CT) # p = 2.73x10^-9

ggplot(Crypto_first_pos, aes(x = EpisodeNumber, y = -CRYCTRT)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()


# Add WHZ information to control for malnutrition at enrollment.

anthro <- read.csv("../Raw_data/Metadata_ForIndividuals.csv")
anthro <- dplyr::mutate(anthro, SID = factor(SID))

Crypto_first_pos_anthro <- dplyr::left_join(Crypto_first_pos, anthro, by = join_by(SID))

str(Crypto_first_pos_anthro)

Episode_CT_anthro <- lmer(EpisodeNumber ~ CRYCTRT + WHZ + HAZ + (1|SID), data = Crypto_first_pos_anthro)

summary(Episode_CT_anthro) # Association between Crypto burden and episode number
# remains significant when controlling for HAZ and WHZ at baseline.
               