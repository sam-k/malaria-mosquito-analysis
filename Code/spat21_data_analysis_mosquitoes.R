# ----------------------------------------- #
#         Spat21 Data Set Analysis          #
#              Mosquito Data                #
#             January 29, 2018              #
#                  S. Kim                   #
# ----------------------------------------- #

#### ------------------ load packages ------------------ ####
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(magrittr)
library(tableone)


#### ---------------- set up environment --------------- ####
.wd <- "~/Projects/Malaria collab/Spatial R21 projects/Spat21 cleaning, analysis/"
MERGED_FP <- paste0(.wd, "Data/Data Sets/merged_mosquito_data.Rdata")
LOG_FP    <- paste0(.wd, "Code/spat21_data_analysis_mosquitoes.log")
close(file(LOG_FP, open="w"))  # clear log file
write.log <- function(...) {
  for(output in list(...)) {
    write(output, file=LOG_FP, append=TRUE)
  }
  write("", file=LOG_FP, append=TRUE)
}


#### -------------- read in mosquito data -------------- ####
load(MERGED_FP)  # allspecies_data, merged_data


#### ----------- analyze all species dataset ----------- ####

tab_allsp_am <- allspecies_data %>%
  select(collection.date, village,
         anoph.unfed, anoph.bloodfed, anoph.halfgravid, anoph.gravid, anoph.undetermined, anoph.total) %>%
  group_by(village) %>%
  summarize(n=n())
