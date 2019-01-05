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

# Perform data checks for qpcr_data.
write.log("# ------ VALIDATE QPCR DATA ------ #")