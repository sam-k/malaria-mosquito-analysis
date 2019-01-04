# ----------------------------------------- #
#         Spat21 Data Set Merging           #
#              Mosquito Data                #
#             January 4, 2018               #
#            K. Sumner, S. Kim              #
# ----------------------------------------- #

#### --------- load packages ----------------- ####
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(magrittr)


#### --------- set up environment ----------------- ####
wd <- "~/Projects/Malaria collab/Spatial R21 projects/Spat21 cleaning, analysis/"
CLEANED_FP <- paste0(wd, "Data/Data Sets/cleaned_data.Rdata")


#### --------- read in mosquito data ----------------- ####
load(CLEANED_FP)  # allspecies_data, anopheles_data, qpcr_data


#### ------------- merge mosquito data sets ---------------- ####

# Extract sample IDs from anopheles_data.
anopheles_data$sample.id <- gsub("\\s?H\\s?", " ", anopheles_data$sample.id.head)

# Group qpcr_data by sample ID.
qpcr_data$Sample.ID <- gsub("\\s?[AH]\\s?", " ", qpcr_data$Sample.Name)
qpcr_data$Head.Abd  <- gsub("[^AH]", "", qpcr_data$Sample.Name)
qpcr_data %<>%
  mutate(Head.Abd = factor(Head.Abd)) %>%
  arrange(Sample.ID, Head.Abd)
qpcr_groupeddata <- as.data.frame(matrix(nrow=nrow(qpcr_data), ncol=17), stringsAsFactors=FALSE)  # overshoot # of rows
names(qpcr_groupeddata) <- c("sample.id",
                             "H.HbtubCT1","H.HbtubCT2","H.pfr364CT1","H.pfr364CT2","H.pfr364Q1","H.pfr364Q2","H.Has.Hb","H.Has.Pf",
                             "A.HbtubCT1","A.HbtubCT2","A.pfr364CT1","A.pfr364CT2","A.pfr364Q1","A.pfr364Q2","A.Has.Hb","A.Has.Pf")
temp_count <- 0
for(i in 1:nrow(qpcr_data)) {
  if(qpcr_data[[i, "Sample.ID"]] != ifelse(i>1, qpcr_data[[i-1, "Sample.ID"]], "")) {
    temp_count <- temp_count + 1
  }
  qpcr_groupeddata[[temp_count, "sample.id"]] <- qpcr_data[[i, "Sample.ID"]]
  if(qpcr_data[[i, "Head.Abd"]] == "H") {
    qpcr_groupeddata[temp_count, 2:9]   <- c(qpcr_data[i, 4:9], qpcr_data[i, 25:26])
  } else if(qpcr_data[[i, "Head.Abd"]] == "A") {
    qpcr_groupeddata[temp_count, 10:17] <- c(qpcr_data[i, 4:9], qpcr_data[i, 25:26])
  }
}
qpcr_groupeddata %<>% filter(!is.na(sample.id))  # trim empty rows

# Merge anopheles_data with qpcr_data.
merged_data <- left_join(anopheles_data, qpcr_groupeddata, by="sample.id")


#### -------- clean up environment ----------------- ####
rm(i, list=ls(pattern="^temp_"))