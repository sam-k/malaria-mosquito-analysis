# ----------------------------------------- #
#         Spat21 Data Set Cleaning          #
#                Human Data                 #
#              March 4, 2019                #
#                  S. Kim                   #
# ----------------------------------------- #

#### ------------------ load packages ------------------ ####
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(magrittr)

#### ---------------- set up environment --------------- ####
.wd <- "~/Projects/Malaria collab/Spatial R21 projects/Spat21 cleaning, analysis/"
IMPORTED_FP <- paste0(.wd, "Data/Data Sets/human_imported_data.Rdata")
CLEANED_FP  <- paste0(.wd, "Data/Data Sets/human_cleaned_data.Rdata")
LOG_FP      <- paste0(.wd, "Code/spat21_human_data_cleaning.log")
close(file(LOG_FP, open="w"))  # clear log file
write.log <- function(...) {
  for(output in list(...)) {
    write(output, file=LOG_FP, append=TRUE)
  }
  write("", file=LOG_FP, append=TRUE)
}
.village_dict <- list(Kinesamo="K", Maruti="M", Sitabicha="S")  # codes for villages


#### ---------------- read in human data --------------- ####
load(IMPORTED_FP)  # monthly_data, qpcr_data


#### ---------------- clean monthly data --------------- ####

write.log("# ------ CLEAN MONTHLY DATA ------ #")

discr_m_monthly_ID <- monthly_data %>%
  filter(!is.na(monthly_unq_memID)) %>%
  select(monthly_unq_memID, unq_memID, today_hum_monthly_data,
         month_hum_monthly_data, year_hum_monthly_data, month_year_combo_monthly_data,
         HH_ID, memID, village_name_hum_monthly_data, village_all_data, `Sample Name`) %>%
  filter(!grepl("^[KMS]\\d{2}-\\d{6}-[1-9]?\\d$", monthly_unq_memID) |
           !grepl("^[KMS]\\d{2}_[1-9]?\\d$", unq_memID) |
           !grepl("^[1-9]?\\d-\\d{4}$", month_year_combo_monthly_data) |
           !grepl("^[KMS]\\d{2}-\\d{6}-[1-9]?\\d$", `Sample Name`) |
           !startsWith(monthly_unq_memID, HH_ID) |
           !endsWith(monthly_unq_memID, as.character(memID)) |
           month_hum_monthly_data != strtoi(str_extract(monthly_unq_memID, "(?<=-\\d{2})\\d{2}")) |
           substr(year_hum_monthly_data, 3, 4) != str_extract(monthly_unq_memID, "(?<=-\\d{4})\\d{2}") |
           paste0(month_hum_monthly_data, "-", year_hum_monthly_data) != month_year_combo_monthly_data |
           !startsWith(monthly_unq_memID, HH_ID) |
           !endsWith(monthly_unq_memID, as.character(memID)) |
           !startsWith(monthly_unq_memID, as.character(village_all_data)) |
           format(today_hum_monthly_data, "%d%m%y") != str_extract(monthly_unq_memID, "(?<=-)\\d{6}") |
           village_all_data != .village_dict[village_name_hum_monthly_data]) %>%
  arrange(monthly_unq_memID)
write.log("All healthy monthly entries are formatted correctly")

discr_m_sick_ID <- monthly_data %>%
  filter(is.na(monthly_unq_memID)) %>%
  select(unq_memID, HH_ID, memID, village_name_hum_sick_data, today_hum_sick_data, sick_unq_memID, `Sample Name`) %>%
  filter(!grepl("^[KMS]\\d{2}_[1-9]?\\d$", unq_memID) |
           !grepl("^[KMS]\\d{2}-\\d{6}-[1-9]?\\d-R$", sick_unq_memID) |
           !grepl("^[KMS]\\d{2}-\\d{6}-[1-9]?\\d-R$", `Sample Name`) |
           !startsWith(unq_memID, HH_ID) |
           !endsWith(unq_memID, as.character(memID)) |
           substr(HH_ID, 1, 1) != .village_dict[village_name_hum_sick_data] |
           !startsWith(sick_unq_memID, HH_ID) |
           !endsWith(sick_unq_memID, paste0(memID, "-R")) |
           format(today_hum_sick_data, "%d%m%y") != str_extract(sick_unq_memID, "(?<=-)\\d{6}") |
           sick_unq_memID != `Sample Name`) %>%
  arrange(unq_memID)
write.log("All sick monthly entries are formatted correctly")


#### ----------------- clean qPCR data ----------------- ####

write.log("# ------ CLEAN QPCR DATA ------ #")

# Check if sample IDs follow the correct format (X## X#####).
discr_q_ID <- qpcr_data %>%
  select(`Sample Name`) %>%
  filter(!grepl("^[KMS]\\d{2}-\\d{6}-[1-9]?\\d(-R)?$", `Sample Name`)) %>%
  arrange(`Sample Name`)
write.table(discr_q_ID, row.names=FALSE, col.names=FALSE, file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
.incorrect_ids <- c("K07-030817-08","K07-030817-09","M15-311017-P-R","M16-270618-P-R","M03-26018-2")
.correct_ids   <- c("K07-030817-8", "K07-030817-9", "M15-311017-6-R","M16-270618-4-R","M03-26018-2-R")
for(.i in 1:length(.incorrect_ids)) {
  qpcr_data$`Sample Name`[qpcr_data$`Sample Name`==.incorrect_ids[[.i]]] <- .correct_ids[[.i]]
}
write.log("K07-030817-8, K07-030817-9 were miscoded and were corrected",
          "M15-311017-P-R, M16-270618-P-R were corrected to M15-311017-6-R, M16-270618-4-R",
          "M03-26018-2 was corrected to M03-26018-2-R")


#### --------------- export cleaned data --------------- ####
save(monthly_data, qpcr_data, file=CLEANED_FP)