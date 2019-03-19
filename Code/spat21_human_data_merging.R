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
CLEANED_FP    <- paste0(.wd, "Data/Human Data Sets/human_cleaned_data.Rdata")
MERGED_FP     <- paste0(.wd, "Data/Human Data Sets/human_merged_data.Rdata")
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


#### ------------ standardize IDs by date ------------- ####

# Remove duplicate sick/non-sick entries, within each and across both datasets.

monthly_truedata <- monthly_data %>%
  mutate(sickID = str_extract(`Sample Name`, "R$")) %>%
  mutate(date   = coalesce(today_hum_monthly_data, today_hum_sick_data)) %>%
  arrange(date, HH_ID, memID)
qpcr_truedata <- qpcr_data %>%
  mutate(date   = dmy(str_extract(`Sample Name`, "(?<=-)\\d{6}(?=-)"))) %>%
  mutate(HH_ID  = str_extract(`Sample Name`, "^[KMS]\\d{2}")) %>%
  mutate(memID  = str_extract(`Sample Name`, "(?<=-)[1-9]?\\d(?=[AB-]|$)")) %>%
  mutate(altID  = str_extract(`Sample Name`, "(?<=\\d)[AB]")) %>%
  mutate(sickID = str_extract(`Sample Name`, "R$")) %>%
  arrange(date, HH_ID, memID, `Sample Name`)

x_monthly_data <- monthly_truedata %>%  # temporary: testing w/ just IDs
  select(`Sample Name`, HH_ID, memID, sickID, date) %>%
  mutate(sampleID = str_extract(`Sample Name`, ".*[^-R]"))
x_qpcr_data <- qpcr_truedata %>%
  select(`Sample Name`, HH_ID, memID, altID, sickID, date) %>%
  mutate(sampleID = str_extract(`Sample Name`, ".*[^-R]"))

.dup_monthly_entries <- !duplicated(x_monthly_data$sampleID, fromLast=TRUE)  # filter out monthly duplicates thru sorting
.dup_qpcr_entries    <- !duplicated(x_qpcr_data$sampleID, fromLast=TRUE)
x_monthly_data %<>% filter(.dup_monthly_entries)
x_qpcr_data    %<>% filter(.dup_qpcr_entries)
write.log(sprintf("%d monthly subjects had both sick and non-sick entries, latter were discarded",
                  sum(!.dup_monthly_entries)),
          sprintf("%d qPCR samples had both sick and non-sick entries, latter were discarded",
                  sum(!.dup_qpcr_entries)))

.sample_ids <- unique(c(x_monthly_data$`Sample Name`, x_qpcr_data$`Sample Name`)) %>%
  as.data.frame() %>%
  `colnames<-`(c("Sample Name")) %>%
  mutate(sampleID = str_extract(`Sample Name`, ".*[^-R]"))

.dup_names <- .sample_ids$sampleID[duplicated(.sample_ids$sampleID)]
.dup_monthly_entries <- x_monthly_data$sampleID %in% .dup_names & is.na(x_monthly_data$sickID)
.dup_qpcr_entries    <- x_qpcr_data$sampleID %in% .dup_names & is.na(x_qpcr_data$sickID)
x_monthly_data$`Sample Name`[.dup_monthly_entries] %<>% paste0("-R")
x_qpcr_data$`Sample Name`[.dup_qpcr_entries]       %<>% paste0("-R")
write.log(sprintf("%d monthly subjects were marked sick to match qPCR data",
                  sum(.dup_monthly_entries)),
          sprintf("%d qPCR samples were marked sick to match monthly data",
                  sum(.dup_qpcr_entries)))


# #### -------------- merge human datasets -------------- ####
# 
# write.log("# ------ MERGE HUMAN DATA ------ #")
# 
# merged_data <- left_join(monthly_data, qpcr_data, by=c("True Sample Name","Sample Name"))
# 
# 
# #### -------------- validate merged data --------------- ####
# 
# write.log("# ------ VALIDATE MERGING ------ #")
# 
# unmerged_monthly <- setdiff(monthly_data$`Sample Name`, qpcr_data$`Sample Name`)
# unmerged_qpcr    <- setdiff(qpcr_data$`Sample Name`, monthly_data$`Sample Name`)
# 
# unmerged_monthly_prefixes <- sort(unique(str_extract(unmerged_monthly, "^[KMS]\\d{2}-\\d{6}")))
# unmerged_qpcr_prefixes    <- sort(unique(str_extract(unmerged_qpcr,    "^[KMS]\\d{2}-\\d{6}")))
# length(unmerged_monthly_prefixes) <- max(length(unmerged_monthly_prefixes), length(unmerged_qpcr_prefixes))
# length(unmerged_qpcr_prefixes)    <- max(length(unmerged_monthly_prefixes), length(unmerged_qpcr_prefixes))
# unmerged_prefixes <- cbind(unmerged_monthly_prefixes, unmerged_qpcr_prefixes)
# 
# 
# #### ---------------- export merged data --------------- ####
# 
# # Export data.
# save(merged_data, file=MERGED_FP)