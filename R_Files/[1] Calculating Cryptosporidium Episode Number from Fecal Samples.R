#Title: [1] Calculating Cryptosporidium Episode Number from Fecal Samples
#Name: Hannah So
#Date: 3/14/24

####Introduction###############################################################

#Main goal of the project is to perform statistical analysis on associations 
#with Cryptosporidium-positive episodes. To do so, it was necessary to organize 
#the data files to obtain the variables of interest. This file focuses on 
#calculating Cryptosporidium episode number for diarrheal and surveillance stool 
#samples on a per child basis. 


####Load Data and Packages#####################################################
#install.packages("tidyverse")
packageVersion("tidyverse") #'2.0.0'
library(tidyverse)

#install.packages("readxl") 
packageVersion("readxl") #‘1.4.2’
library(readxl)


ENR_Fecal_Data <- read.csv("../Raw_Data/StoolSample_CryptoStatus.csv")


####Determining Infection Number###############################################

#Make data frame of all positive samples and of all first infection samples
Pos.Crypto <- dplyr::filter(ENR_Fecal_Data, CryptoStatus == "POSITIVE")

#Separating first crypto episode for each individual
Pos.FirstCrypto <- Pos.Crypto %>%
  dplyr::group_by(SID) %>%
  dplyr::slice_min(n=1, order_by = AgeInDays) %>%
  dplyr::mutate(EpisodeType = SampleType)


#Determine difference in sample ages on a per child basis
Inf.Num.Crypto <- Pos.Crypto %>%
  dplyr::group_by(SID) %>%
  dplyr::arrange(AgeInDays, .by_group = TRUE) %>%
  dplyr::mutate(PreviousPositive.AgeInDays = lag(AgeInDays, n=1)) %>%
  dplyr::mutate(SampleAgeDiff = AgeInDays - PreviousPositive.AgeInDays)

#Filter for infection episode criteria
Episode.Num.Crypto <- Inf.Num.Crypto %>%
  dplyr:: filter(SampleAgeDiff >= 65) %>% # Removes first episode for each SID
  dplyr:: group_by(SID) %>%
  dplyr::arrange(AgeInDays, .by_group = TRUE) %>%
  dplyr::mutate(EpisodeNumber = row_number()) %>%
  dplyr::mutate(EpisodeType = SampleType)

#Cut-off of 65 days because that is at least 2 collected stool samples
#that were negative 


#Change to reflect actual crypto infection episode number (adjust for removed 1st episode)
Corrected.Episode.Num.Crypto <- Episode.Num.Crypto 
Corrected.Episode.Num.Crypto$EpisodeNumber[Corrected.Episode.Num.Crypto$EpisodeNumber == 7] <- 8
Corrected.Episode.Num.Crypto$EpisodeNumber[Corrected.Episode.Num.Crypto$EpisodeNumber == 6] <- 7
Corrected.Episode.Num.Crypto$EpisodeNumber[Corrected.Episode.Num.Crypto$EpisodeNumber == 5] <- 6
Corrected.Episode.Num.Crypto$EpisodeNumber[Corrected.Episode.Num.Crypto$EpisodeNumber == 4] <- 5
Corrected.Episode.Num.Crypto$EpisodeNumber[Corrected.Episode.Num.Crypto$EpisodeNumber == 3] <- 4
Corrected.Episode.Num.Crypto$EpisodeNumber[Corrected.Episode.Num.Crypto$EpisodeNumber == 2] <- 3
Corrected.Episode.Num.Crypto$EpisodeNumber[Corrected.Episode.Num.Crypto$EpisodeNumber == 1] <- 2

#This table is data that excludes the very first infection, so what it
#considers the first infection is actually the second infection. 
#We are adjusting it here so it accurately reflects this



#Check that labeling is accurate 
table(Episode.Num.Crypto$EpisodeNumber, useNA = "ifany")
table(Corrected.Episode.Num.Crypto$EpisodeNumber, useNA = "ifany")


#Reformat data frames for rbind
Corrected.Episode.Num.Crypto <- dplyr::select(Corrected.Episode.Num.Crypto, 
                                              -PreviousPositive.AgeInDays,
                                              -SampleAgeDiff)


#Combine first episodes with subsequent episodes
Pos.FirstCrypto$EpisodeNumber <- 1

CryptoEpisodes.OnlyPositive <- rbind(Corrected.Episode.Num.Crypto, Pos.FirstCrypto)


#Separating Diarrheal and Subclinical Episode Numbers
DiarrhealEpisodes.OnlyPositive <- CryptoEpisodes.OnlyPositive %>%
  dplyr::filter(EpisodeType == "Diarrheal") %>%
  dplyr::group_by(SID) %>%
  dplyr::arrange(EpisodeNumber, .by_group = TRUE) %>%
  dplyr::mutate(DiarrhealEpisodeNumber = row_number()) %>%
  dplyr::mutate(SubclinicalEpisodeNumber = NA) 

SubclinicalEpisodes.OnlyPositive <- CryptoEpisodes.OnlyPositive %>%
  dplyr::filter(EpisodeType == "Subclinical") %>%
  dplyr::group_by(SID) %>%
  dplyr::arrange(EpisodeNumber, .by_group = TRUE) %>%
  dplyr::mutate(SubclinicalEpisodeNumber = row_number()) %>%
  dplyr::mutate(DiarrhealEpisodeNumber = NA)

CombinedEpisodes.OnlyPositive <- rbind(DiarrhealEpisodes.OnlyPositive, SubclinicalEpisodes.OnlyPositive)


#Include all crypto positive samples that are not distinct episodes
Lessthan65 <- Inf.Num.Crypto %>%
  dplyr::filter(SampleAgeDiff < 65) %>%
  dplyr::mutate(EpisodeNumber = NA) %>%
  dplyr::mutate(EpisodeType = NA) %>%
  dplyr::mutate(SubclinicalEpisodeNumber = NA) %>%
  dplyr::mutate(DiarrhealEpisodeNumber = NA) %>%
  dplyr::select(-PreviousPositive.AgeInDays, -SampleAgeDiff)

#Combine all positive samples
CryptoEpisodes.AllPositives <- rbind(CombinedEpisodes.OnlyPositive, Lessthan65)




####Incorporate Negative Samples################################################

#Merge Positives with Negatives
Neg.Crypto <- ENR_Fecal_Data %>%
  dplyr::filter(CryptoStatus == "NEGATIVE") %>%
  dplyr::mutate(EpisodeType = NA) %>%
  dplyr::mutate(SubclinicalEpisodeNumber = NA) %>%
  dplyr::mutate(DiarrhealEpisodeNumber = NA) %>%
  dplyr::mutate(EpisodeNumber = NA)

AllSamples <- rbind(CryptoEpisodes.AllPositives, Neg.Crypto)


#Label number of infection for NAs 
AllSamples_CryptoEpisode <- AllSamples %>%
  dplyr::group_by(SID) %>%
  dplyr::arrange(AgeInDays, .by_group = TRUE) %>%
  tidyr::fill(EpisodeNumber, .direction = "down") %>%
  tidyr::fill(EpisodeType, .direction= "down") %>%
  tidyr::fill(DiarrhealEpisodeNumber, .direction= "down") %>%
  tidyr::fill(SubclinicalEpisodeNumber, .direction= "down")
  
#Confirm all positive samples have been imputed
AllSamples %>%
  dplyr::filter(CryptoStatus == "POSITIVE") %>%
  dplyr::group_by(EpisodeNumber) %>%
  dplyr::tally()

AllSamples_CryptoEpisode %>%
  dplyr::filter(CryptoStatus == "POSITIVE") %>%
  dplyr::group_by(EpisodeNumber) %>%
  dplyr::tally() # No NAs after filling

#Add episode number for negative samples before first infection
AllSamples_CryptoEpisode$EpisodeNumber[is.na(AllSamples_CryptoEpisode$EpisodeNumber)] <- 0
AllSamples_CryptoEpisode$DiarrhealEpisodeNumber[is.na(AllSamples_CryptoEpisode$DiarrhealEpisodeNumber)] <- 0
AllSamples_CryptoEpisode$SubclinicalEpisodeNumber[is.na(AllSamples_CryptoEpisode$SubclinicalEpisodeNumber)] <- 0

AllSamples_CryptoEpisode <- dplyr::select(AllSamples_CryptoEpisode,
                                         -SITE, -SampleType)

#write.csv(AllSamples_CryptoEpisode,"../Results/FecalSample_CryptoEpisodeNumber.csv",row.names = FALSE)






