# ----------------------------------------- #
#        Spat21 Data Set Importing          #
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
MONTHLY_FP  <- paste0(.wd, "Data/Human Data Sets/hum_monthly_merged_with_table_and_sick_with_exclusion_data_27FEB2019.RDS")
QPCR_FP     <- paste0(.wd, "Data/Human Data Sets/spat21_qpcr_data_clean_human_dbs_16JAN2019.RDS")
IMPORTED_FP <- paste0(.wd, "Data/Human Data Sets/human_imported_data.Rdata")
LOG_FP      <- paste0(.wd, "Code/spat21_human_data_import.log")
close(file(LOG_FP, open="w"))  # clear log file
write.log <- function(...) {
  for(output in list(...)) {
    write(output, file=LOG_FP, append=TRUE)
  }
  write("", file=LOG_FP, append=TRUE)
}


#### --------------- read in human data ---------------- ####

write.log("# ------ IMPORT RAW DATA ------ #")

# Read in the monthly visit data set.
monthly_data <- read_rds(MONTHLY_FP)
# Read in the human qPCR data set.
qpcr_data    <- read_rds(QPCR_FP)

# Look at summaries of all data sets.
summary(monthly_data)
summary(qpcr_data)
str(monthly_data, list.len=ncol(monthly_data))  # avoid truncating output
str(qpcr_data)
write.log("monthly_data dims:", paste(ncol(monthly_data), "vars"), paste(nrow(monthly_data), "obs"))
write.log("qpcr_data dims:",    paste(ncol(qpcr_data), "vars"),    paste(nrow(qpcr_data), "obs"))



#### ------------------ reformat data ------------------ ####

write.log("# ------ REFORMAT MONTHLY DATA ------ #")

monthly_data %<>%
  mutate_at(c("memID"), as.integer) %>%
  mutate_at(c("to_where","med_oth_hum_monthly_data"), as.character) %>%
  mutate_at(c("village_all_data"), factor)

write.log("# ------ REFORMAT QPCR DATA ------ #")

qpcr_data[qpcr_data=="Undetermined"] <- NA
qpcr_data %<>%
  mutate_at(c("pfr364Std8a","pfr364Std8b","pfr364Std10b"), as.numeric)


#### ------------ export reformatted data -------------- ####
save(monthly_data, qpcr_data, file=IMPORTED_FP)