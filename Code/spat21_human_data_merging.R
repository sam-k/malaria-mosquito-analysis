# ----------------------------------------- #
#         Spat21 Data Set Merging           #
#                Human Data                 #
#              March 4, 2019                #
#                  S. Kim                   #
# ----------------------------------------- #

#### ------------------ load packages ------------------ ####
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(magrittr)


#### ---------------- set up environment --------------- ####
.wd <- "~/Projects/Malaria collab/Spatial R21 projects/Spat21 cleaning, analysis/"
CLEANED_FP    <- paste0(.wd, "Data/Data Sets/human_cleaned_data.Rdata")
MERGED_FP     <- paste0(.wd, "Data/Data Sets/human_merged_data.Rdata")
LOG_FP        <- paste0(.wd, "Code/spat21_human_data_merging.log")
close(file(LOG_FP, open="w"))  # clear log file
write.log <- function(...) {
  for(output in list(...)) {
    write(output, file=LOG_FP, append=TRUE)
  }
  write("", file=LOG_FP, append=TRUE)
}


#### --------------- read in human data --------------- ####

load(CLEANED_FP)  # monthly_data, qpcr_data
# monthly_data %<>% arrange(collection.date, household.id)
# qpcr_data    %<>% arrange(Sample.ID, Head.Abd)


#### -------------- merge human datasets -------------- ####

write.log("# ------ MERGE HUMAN DATA ------ #")

merged_data <- left_join(monthly_data, qpcr_data, by="Sample Name")


#### -------------- validate merged data --------------- ####

write.log("# ------ VALIDATE MERGING ------ #")

unmerged_monthly <- setdiff(monthly_data$`Sample Name`, qpcr_data$`Sample Name`)
unmerged_qpcr    <- setdiff(qpcr_data$`Sample Name`, monthly_data$`Sample Name`)

unmerged_monthly_prefixes <- sort(unique(str_extract(unmerged_monthly, "^[KMS]\\d{2}-\\d{6}")))
unmerged_qpcr_prefixes    <- sort(unique(str_extract(unmerged_qpcr,    "^[KMS]\\d{2}-\\d{6}")))
length(unmerged_monthly_prefixes) <- max(length(unmerged_monthly_prefixes), length(unmerged_qpcr_prefixes))
length(unmerged_qpcr_prefixes)    <- max(length(unmerged_monthly_prefixes), length(unmerged_qpcr_prefixes))
unmerged_prefixes <- cbind(unmerged_monthly_prefixes, unmerged_qpcr_prefixes)


#### ---------------- export merged data --------------- ####

# Export data.
save(merged_data, file=MERGED_FP)