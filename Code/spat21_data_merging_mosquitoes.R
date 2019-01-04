# ----------------------------------------- #
#        Spat21 Data Set Cleaning           #
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
CLEANED_FP = paste0(wd, "Data/Data Sets/cleaned_data.Rdata")
# zero = 1e-6  # threshold for zero CT value


#### --------- read in mosquito data ----------------- ####

load(CLEANED_FP)


#### ------------- clean each variable in mosquito data sets ---------------- ####


## -------- allspecies_data ----------------- ####
