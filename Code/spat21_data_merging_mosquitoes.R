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
VALIDATED_FP  <- paste0(wd, "Data/Data Sets/validated_data.Rdata")
MERGED_CSV_FP <- paste0(wd, "Data/Data Sets/merged_mosquito_data.csv")
MERGED_RDS_FP <- paste0(wd, "Data/Data Sets/merged_mosquito_data.rds")
LOG_FP        <- paste0(wd, "Code/spat21_data_merging_mosquitoes.log")
close(file(LOG_FP, open="w"))  # clear log file
write.log <- function(...) {
  for(temp_output in list(...)) {
    write(temp_output, file=LOG_FP, append=TRUE)
  }
  write("", file=LOG_FP, append=TRUE)
}


#### --------- read in mosquito data ----------------- ####
load(VALIDATED_FP)  # allspecies_data, anopheles_data, qpcr_data


#### ------------- merge mosquito data sets ---------------- ####

write.log("# ------ MERGE MOSQUITO DATA ------ #")

# Group qpcr_data by sample ID.
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
write.log("Converted qPCR data to wide format by sample ID")

# Merge anopheles_data with qpcr_data.
merged_data <- left_join(anopheles_data, qpcr_groupeddata, by="sample.id")
write.log("Merged anopheles descriptive data with wide qPCR data")


#### ------------- validate merged data ---------------- ####

# Check if any entries were not merged.
unmerged_anoph_data <- anopheles_data[-which(anopheles_data$sample.id %in% merged_data$sample.id), ]
write.log("All descriptive entries were present in qPCR data and merged correctly")
unmerged_qpcr_data  <- qpcr_groupeddata[-which(qpcr_groupeddata$sample.id %in% merged_data$sample.id), ]
write.log(paste("qPCR entries", paste(unmerged_qpcr_data$sample.id, collapse=", "), "were absent from descriptive data and did not merge"))

#### -------- export merged data ----------------- ####
write.csv(merged_data, file=MERGED_CSV_FP, row.names=FALSE)
saveRDS(merged_data, file=MERGED_RDS_FP)


#### -------- clean up environment ----------------- ####
rm(i, list=ls(pattern="^temp_"))