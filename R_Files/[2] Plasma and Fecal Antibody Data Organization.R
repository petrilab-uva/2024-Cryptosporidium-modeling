#Title: Plasma and Fecal Antibody Data Organization 
#Name: Hannah So
#Date: 3/22/24

####Introduction################################################################
#This file organizes antibody data frames so that they 
#can be merged with metadata and crypto episode number data. 


####Load Data and Packages######################################################
##Packages##
#install.packages("tidyverse")
packageVersion("tidyverse") # '2.0.0'
library(tidyverse)

Fecal.Antibody <- read.csv("../Raw_Data/Fecal_Antibody_RU.csv")
Fecal.EpisodeNumber <- read.csv("../Results/FecalSample_CryptoEpisodeNumber.csv")

Plasma.Antibody <- read.csv("../Raw_Data/Plasma_Antibody_RU.csv")

Metadata <- read.csv("../Raw_Data/Metadata_ForIndividuals.csv")



####Data Organization#########################################################

#Add Metadata to Fecal Antibody Data 
Fecal.Antibody.Metadata <- Fecal.Antibody %>%
  dplyr::left_join(., Fecal.EpisodeNumber, by = join_by(SPECID, SID)) %>%
  dplyr::left_join(., Metadata, by = join_by(SID))

length(unique(Fecal.Antibody.Metadata$SID)) # 54 unique SIDs

#write.csv(Fecal.Antibody.Metadata, "../Results/Fecal_Antibody_Metadata.csv", 
#          row.names = FALSE)


#Add Columns to Allow Binding of Metadata and Plasma Antibody Data
Plasma.Antibody.Rbind <- Plasma.Antibody %>%
  dplyr::mutate(CryptoStatus = NA) %>%
  dplyr::mutate(EpisodeNumber = NA) %>%
  dplyr::mutate(EpisodeType = NA) %>%
  dplyr::mutate(DiarrhealEpisodeNumber = NA) %>%
  dplyr::mutate(SubclinicalEpisodeNumber = NA) 

EpisodeNumber.Rbind <- Fecal.EpisodeNumber %>%
  dplyr::mutate(Season = NA) %>%
  dplyr::mutate(Cp17.IgG.Plasma = NA) %>%
  dplyr::mutate(Cp17.IgA.Plasma = NA) %>%
  dplyr::mutate(Cp23.IgG.Plasma = NA) %>%
  dplyr::mutate(Cp23.IgA.Plasma = NA) %>%
  dplyr::select(-CRYCTRT)
  

#Combine plasma antibody and metadata
Plasma.Antibody.Metadata <- rbind(Plasma.Antibody.Rbind, EpisodeNumber.Rbind)
Plasma.Antibody.Metadata <- dplyr::left_join(Plasma.Antibody.Metadata, Metadata, by = join_by(SID))


#Fill Episode Number and remove samples without Ab values.
Plasma.Antibody.Metadata.Filled <- Plasma.Antibody.Metadata %>%
  dplyr::group_by(SID) %>%
  dplyr::arrange(AgeInDays, .by_group = TRUE) %>%
  tidyr::fill(EpisodeNumber, .direction = "down") %>%
  tidyr::fill(EpisodeType, .direction = "down") %>%
  tidyr::fill(DiarrhealEpisodeNumber, .direction = "down") %>%
  tidyr::fill(SubclinicalEpisodeNumber, .direction = "down") %>%
  dplyr::select(-CryptoStatus) %>%
  tidyr::drop_na(Cp17.IgA.Plasma)

table(Plasma.Antibody.Metadata.Filled$Cp17.IgG.Plasma, useNA = "ifany") # No NAs

length(unique(Plasma.Antibody.Metadata.Filled$SID)) # 54 unique SIDs

#write.csv(Plasma.Antibody.Metadata.Filled, "../Results/Plasma_Antibody_Metadata.csv", row.names = FALSE)

