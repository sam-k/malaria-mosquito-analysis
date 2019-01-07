# ----------------------------------------- #
#        Spat21 Data Set Validation         #
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
LOG_FP     <- paste0(wd, "Code/spat21_data_checks_mosquitoes.log")
close(file(LOG_FP, open="w"))  # clear log file
write.log <- function(...) {
  for(temp_output in list(...)) {
    write(temp_output, file=LOG_FP, append=TRUE)
  }
  write("", file=LOG_FP, append=TRUE)
}


#### --------- read in mosquito data ----------------- ####
load(CLEANED_FP)  # allspecies_data, anopheles_data, qpcr_data


#### -------- validate data ----------------- ####

# Perform data checks for anopheles_data.
write.log("# ------ VALIDATE ANOPH. DESCRIPTIVE DATA ------ #")
# Check if village names are consistent.
anoph_village_discrepancies <- anopheles_data %>%
  select(household.id, village, sample.id.head, sample.id.abdomen) %>%
  filter((substr(household.id,1,1)!=substr(sample.id.head,1,1)) & (substr(household.id,1,1)!=substr(sample.id.abdomen,1,1)) &
         (substr(household.id,1,1)!=substr(village,1,1)))
write.log("All village names match")
# Check if heads/abdomens are missing.
anoph_counts <- anopheles_data %>%
  select(village, sample.id.head, sample.id.abdomen) %>%
  group_by(village) %>%
  summarize(H=sum(!is.na(sample.id.head)), A=sum(!is.na(sample.id.abdomen))) %>%
  as.data.frame()
rownames(anoph_counts) <- substr(anoph_counts$village, 1, 1)
anoph_counts$village <- NULL
write.table(anoph_counts, col.names=NA, file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
write.log("S02 A00001 is missing (as noted in comment)",
          "K05 A00034, M01 A00058, S07 A00008 are commented as missing, but present in data")
# Check if abdominal statuses are correct.
anoph_abd_discrepancies <- anopheles_data %>%
  select()

# Perform data checks for qpcr_data.
write.log("# ------ VALIDATE QPCR DATA ------ #")


#### -------- export validated data ----------------- ####




#### -------- clean up environment ----------------- ####
rm(list=ls(pattern="^temp_"))