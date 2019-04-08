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
library(stringr)
library(lubridate)
library(magrittr)
library(reshape2)


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


#### ------------ standardize IDs by date ------------- ####

# Remove duplicate sick/non-sick entries within each dataset.

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
  select(`Sample Name`, HH_ID, memID, altID, sickID, date, pf_pcr_infection_status) %>%
  mutate(sampleID = str_extract(`Sample Name`, ".*[^-R]"))

# c_monthly_data <- x_monthly_data  # temporary: copy original datasets, just in case
# c_qpcr_data    <- x_qpcr_data
# 
# .dup_monthly <- duplicated(x_monthly_data$sampleID, fromLast=TRUE)  # filter out monthly duplicates thru sorting
# .dup_qpcr    <- duplicated(x_qpcr_data$sampleID, fromLast=TRUE)
# dup_monthly <- x_monthly_data %>%
#   filter(.dup_monthly | duplicated(.$sampleID)) %>%
#   arrange(sampleID)
# dup_qpcr <- x_qpcr_data %>%
#   filter(.dup_qpcr | duplicated(.$sampleID)) %>%
#   arrange(sampleID)
# x_monthly_data %<>% filter(!.dup_monthly)
# x_qpcr_data    %<>% filter(!.dup_qpcr)
# write.log(sprintf("%d monthly subjects had both sick and non-sick entries, latter were discarded", sum(.dup_monthly)),
#           sprintf("%d qPCR samples had both sick and non-sick entries, latter were discarded", sum(.dup_qpcr)))

# Remove duplicate sick/non-sick entries across both datasets.

# .sample_ids <- unique(c(x_monthly_data$`Sample Name`, x_qpcr_data$`Sample Name`)) %>%
#   as.data.frame() %>%
#   `colnames<-`(c("Sample Name")) %>%
#   mutate(sampleID = str_extract(`Sample Name`, ".*[^-R]"))
# 
# .dup_names <- .sample_ids$sampleID[duplicated(.sample_ids$sampleID)]
# .dup_monthly_entries <- x_monthly_data$sampleID %in% .dup_names & is.na(x_monthly_data$sickID)
# .dup_qpcr_entries    <- x_qpcr_data$sampleID %in% .dup_names & is.na(x_qpcr_data$sickID)
# x_monthly_data$`Sample Name`[.dup_monthly_entries] %<>% paste0("-R")
# x_qpcr_data$`Sample Name`[.dup_qpcr_entries]       %<>% paste0("-R")
# write.log(sprintf("%d monthly subjects were marked sick to match qPCR data", sum(.dup_monthly_entries)),
#           sprintf("%d qPCR samples were marked sick to match monthly data", sum(.dup_qpcr_entries)))
# 
# x_monthly_data %<>%
#   mutate(c.date = date) %>%
#   mutate(`c.Sample Name` = `Sample Name`) %>%
#   mutate(sickID = str_extract(`Sample Name`, "R$")) %>%
#   arrange(HH_ID, memID, date, `Sample Name`)
# x_qpcr_data %<>%
#   mutate_at(c("memID"), as.integer) %>%
#   mutate(c.date = date) %>%
#   mutate(`c.Sample Name` = `Sample Name`) %>%
#   mutate(sickID = str_extract(`Sample Name`, "R$")) %>%
#   arrange(HH_ID, memID, altID, date, `Sample Name`)

# Match monthly visit dates with qPCR DBS dates by crawling both datasets.

.i <- 1  # monthly dataset index
.j <- 1  # qPCR dataset index
.dup <- FALSE  # were there multiple visits for a single qPCR date?
.misdated_monthly      <- vector("list", nrow(x_monthly_data))  # misdated monthly entries, first occurrence
.misdated_corr_qpcr    <- vector("list", nrow(x_qpcr_data))     # correspondng qPCR entries
.misdated_dup_monthly  <- vector("list", nrow(x_monthly_data))  # misdated monthly entries, subseq. occurrences
.misdated_dup_corr_pcr <- vector("list", nrow(x_qpcr_data))     # corresponding qPCR entries
while(.i <= nrow(x_monthly_data) & .j <= nrow(x_qpcr_data)) {
  if(x_monthly_data$HH_ID[[.i]] == x_qpcr_data$HH_ID[[.j]]) {
    if(x_monthly_data$memID[[.i]] == x_qpcr_data$memID[[.j]]) {
      if(x_monthly_data$date[[.i]] > x_qpcr_data$date[[.j]])    {
        .monthly_date_i <- x_monthly_data$date[[.i]]
        .qpcr_date_j    <- x_qpcr_data$date[[.j]]
        for(.k in 1:6) {  # monthly visit date can be 0-6 days later than qPCR DBS date
          if(.monthly_date_i == .qpcr_date_j + .k) {
            if(.dup) {
              .misdated_dup_monthly[[.i]]  <- x_monthly_data$`Sample Name`[[.i]]
              .misdated_dup_corr_pcr[[.j]] <- x_qpcr_data$`Sample Name`[[.j]]
            } else {
              .misdated_monthly[[.i]]   <- x_monthly_data$`Sample Name`[[.i]]
              .misdated_corr_qpcr[[.j]] <- x_qpcr_data$`Sample Name`[[.j]]
              x_monthly_data$date[[.i]]          <- .qpcr_date_j
              x_monthly_data$sampleID[[.i]]      <- x_qpcr_data$sampleID[[.j]]
              x_monthly_data$`Sample Name`[[.i]] <- x_qpcr_data$`Sample Name`[[.j]]
              .dup <- TRUE
            }
            break
          }
        }
        .j <- .j+1
      } else {
        if(x_monthly_data$date[[.i]] < x_qpcr_data$date[[.j]]) { .i <- .i+1 } else { .j <- .j+1 }
        .dup <- FALSE
      }
    } else {
      if(x_monthly_data$memID[[.i]] < x_qpcr_data$memID[[.j]]) { .i <- .i+1 } else { .j <- .j+1 }
      .dup <- FALSE
    }
  } else {
    if(x_monthly_data$HH_ID[[.i]] < x_qpcr_data$HH_ID[[.j]]) { .i <- .i+1 } else { .j <- .j+1 }
    .dup <- FALSE
  }
}
.misdated_monthly %<>% .[lengths(.) != 0] %>%
  { do.call(rbind, Map(cbind, ., str_extract(., "(?<=-)\\d{6}(?=-)"))) } %>%
  as.data.frame() %>%
  `colnames<-`(c("Sample Name.m","date.m")) %>%
  mutate_at(c("Sample Name.m"), as.character) %>%
  mutate_at(c("date.m"), dmy) %>%
  mutate(HH_ID.m = str_extract(`Sample Name.m`, "^[KMS]\\d{2}")) %>%
  mutate(memID.m = str_extract(`Sample Name.m`, "(?<=-)[1-9]?\\d(?=[AB-]|$)")) %>%
  arrange(HH_ID.m, memID.m, date.m)
.misdated_corr_qpcr %<>% .[lengths(.) != 0] %>%
  { do.call(rbind, Map(cbind, ., str_extract(., "(?<=-)\\d{6}(?=-)"))) } %>%
  as.data.frame() %>%
  `colnames<-`(c("Sample Name.q","date.q")) %>%
  mutate_at(c("Sample Name.q"), as.character) %>%
  mutate_at(c("date.q"), dmy) %>%
  mutate(HH_ID.q = str_extract(`Sample Name.q`, "^[KMS]\\d{2}")) %>%
  mutate(memID.q = str_extract(`Sample Name.q`, "(?<=-)[1-9]?\\d(?=[AB-]|$)")) %>%
  arrange(HH_ID.q, memID.q, date.q)
misdated_merged <- cbind(.misdated_monthly, .misdated_corr_qpcr) %>%
  select(`Sample Name.m`, `date.m`, `Sample Name.q`, `date.q`) %>%
  mutate(diff = as.integer(date.m - date.q))


#### -------------- merge human datasets -------------- ####

x_merged_data <- full_join(x_monthly_data, x_qpcr_data, by="sampleID", suffix=c(".m",".q"))

unmerged_monthly <- x_merged_data[is.na(x_merged_data$`Sample Name.q`), ] %>%
  select(`Sample Name.m`, HH_ID.m, memID.m, sickID.m, date.m, sampleID)
unmerged_qpcr    <- x_merged_data[is.na(x_merged_data$`Sample Name.m`), ] %>%
  select(`Sample Name.q`, HH_ID.q, memID.q, sickID.q, date.q, sampleID)

x_merged_data %<>% filter(!is.na(`Sample Name.m`) & !is.na(`Sample Name.q`))


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