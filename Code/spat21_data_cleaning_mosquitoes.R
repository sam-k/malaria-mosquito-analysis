# ----------------------------------------- #
#        Spat21 Data Set Cleaning           #
#              Mosquito Data                #
#            December 18, 2018              #
#            K. Sumner, S. Kim              #
# ----------------------------------------- #

#### --------- load packages ----------------- ####
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(haven)


#### --------- set up environment ----------------- ####
wd <- "~/Projects/Malaria collab/Spatial R21 projects/Spat21 cleaning, analysis/"
allspecies_fp = paste0(wd, "Data/Data Sets/MOZZIECollectionSummary_June2017_July2018.csv")
anopheles_fp  = paste0(wd, "Data/Data Sets/MOZZIEFemaleAnophele_June2017_July2018.csv")
anopheles2_fp = paste0(wd, "Data/Data Sets/Individual_female_anoph_long.dta")
qpcr_fp       = paste0(wd, "Data/Data Sets/Mozzie mosquito compiled detection results 18Dec2018.csv")
zero = 1e-6  # threshold for zero CT value


#### --------- read in mosquito data ----------------- ####

# Read in the mosquito descriptive data sets.
# Read in the data set with all mosquito species.
allspecies_data <- read.csv(allspecies_fp, stringsAsFactors=FALSE)
# Read in the wide data set with only anopheles mosquitoes.
anopheles_widedata <- read.csv(anopheles_fp, stringsAsFactors=FALSE)
# anopheles2_data <- read_dta(anopheles2_fp)  # STATA
# Read in the mosquito qPCR data sets.
qpcr_data <- read.csv(qpcr_fp, stringsAsFactors=FALSE)

# Clean column names.
names(anopheles_widedata) <- tolower(gsub("..", ".", names(anopheles_widedata), fixed=TRUE))
names(anopheles_widedata) <- tolower(gsub("\\.$", "", names(anopheles_widedata)))  # remove trailing periods

# Look at summaries of all the data sets.
summary(allspecies_data)
summary(anopheles_widedata)
summary(qpcr_data)
str(allspecies_data)
str(anopheles_widedata)
str(qpcr_data)

# Output a CSV file of all the variable names.
allnames <- c(names(allspecies_data), names(anopheles_widedata), names(qpcr_data))
allnames <- data.frame(allnames)
write_csv(allnames, paste0(wd, "Data/Data Dictionary/spat21_data_mosquito_dictionary.csv"))


#### ------------- clean each variable in mosquito data sets ---------------- ####

# Positive, negative, missing counts for each data set.
counts <- matrix(rep(0,15), nrow=3, ncol=5, dimnames=list(c("allspecies","anopheles","qpcr"),
                                                          c("hb_positive","hb_negative","pf_positive","pf_negative","missing")))

# Rename and reformat all species data columns.
names(allspecies_data) <- c("household.id","repeat.instrument","repeat.instance","collection.date","collection.time","village","collection.done.by",
                            "anoph.unfed","anoph.bloodfed","anoph.halfgravid","anoph.gravid","anoph.undetermined","anoph.total","num.male.anoph",
                            "culex.unfed","culex.bloodfed","culex.halfgravid","culex.gravid","culex.undetermined","culex.total","num.male.culex",
                            "form.checked.by","form.checked.date","form.entered.by","form.entered.date","complete")
temp_cols <- c("household.id","repeat.instrument","village","collection.done.by","form.checked.by","form.entered.by","complete")
allspecies_data[temp_cols]        <- lapply(allspecies_data[temp_cols], factor)
temp_cols <- c("repeat.instance",
               "anoph.unfed","anoph.bloodfed","anoph.halfgravid","anoph.gravid","anoph.undetermined","anoph.total","num.male.anoph",
               "culex.unfed","culex.bloodfed","culex.halfgravid","culex.gravid","culex.undetermined","culex.total","num.male.culex")
allspecies_data[temp_cols]        <- lapply(allspecies_data[temp_cols], as.integer)
temp_cols <- c("collection.date","form.checked.date","form.entered.date")
allspecies_data$collection.date   <- mdy(allspecies_data$collection.date)
allspecies_data$form.checked.date <- mdy(allspecies_data$form.checked.date)
allspecies_data$form.entered.date <- mdy(allspecies_data$form.entered.date)
allspecies_data$collection.time   <- as.logical(allspecies_data$collection.time)

# Reformat anopheles data columns.
temp_cols <- c("household.id","repeat.instrument","village","collection.done.by","samples.prepared.by","species.id.done.by",
               "form.checked.by","form.entered.by","complete")
for(i in 1:16) {
  temp_cols <- unlist(list(temp_cols, paste0(
    c("sample.id.head.","sample.id.abdomen.","abdominal.status.","species.type.","specify.species.","comment."), i)))
}
anopheles_widedata[temp_cols]      <- lapply(anopheles_widedata[temp_cols], factor)
temp_cols <- c("repeat.instance","total.number.of.mosquitos.in.the.household")
anopheles_widedata[temp_cols]      <- lapply(anopheles_widedata[temp_cols], as.integer)
temp_cols <- c("collection.date","form.checked.date","form.entered.on")
anopheles_widedata[temp_cols]      <- lapply(anopheles_widedata[temp_cols], mdy)
anopheles_widedata$collection.time <- as.logical(anopheles_widedata$collection.time)
# STATA operations, same as above.
# temp_cols <- c("household_id","redcap_repeat_instrument","sampleid","specify","comment")
# anopheles2_data[temp_cols]      <- lapply(anopheles2_data[temp_cols], factor)
# temp_cols <- c("redcap_repeat_instance","index","village2","staffname2","sampleby","speciesid","totalnumber","abdominal","speciestype")
# anopheles2_data[temp_cols]      <- lapply(anopheles2_data[temp_cols], as.integer)
# anopheles2_data$collectiondate2 <- mdy(anopheles2_data$collectiondate2)
# anopheles2_data$collecttime     <- as.logical(anopheles2_data$collecttime)
# convert wide format to long format
anopheles_data <- data.frame(matrix(nrow=16*nrow(anopheles_widedata), ncol=21))  # overshoot # of rows
names(anopheles_data) <- c("household.id","repeat.instrument","repeat.instance","collection.date","collection.time","village",
                           "collection.done.by","samples.prepared.by","species.id.done.by","total.number.of.mosquitos.in.the.household",
                           "sample.id.head","sample.id.abdomen","abdominal.status","species.type","specify.species","comment",
                           "form.checked.by","form.checked.date","form.entered.by","form.entered.on","complete")
# temp_count <- 1
# for(i in 1:nrow(anopheles_widedata)) {
#   for(j in 1:16) {
#     if()
#   }
# }


# Reformat qPCR data columns.
qpcr_data[qpcr_data == "Undetermined"] <- NA
temp_cols <- c("Sample.Name","Experiment.Name")
qpcr_data[temp_cols] <- lapply(qpcr_data[temp_cols], factor)
temp_cols <- c("HbtubCT1","HbtubCT2","pfr364CT1","pfr364CT2","pfr364Std5a","pfr364Std5b","pfr364Std6a","pfr364Std6b")
qpcr_data[temp_cols] <- lapply(qpcr_data[temp_cols], as.numeric)
qpcr_data[is.na(qpcr_data)] <- NA  # correct NaNs to NAs
# Process qPCR CT values.
qpcr_data$Has.Hb <- FALSE
qpcr_data$Has.Pf <- FALSE
qpcr_data$Has.Hb[which(qpcr_data$HbtubCT1<=zero & qpcr_data$HbtubCT2<=zero & qpcr_data$pfr364CT1<=zero & qpcr_data$pfr364CT2<=zero)] <- NA
qpcr_data$Has.Pf[which(is.na(qpcr_data$Has.Hb))] <- NA
temp_cols <- c("HbtubCT1","HbtubCT2","pfr364CT1","pfr364CT2")
qpcr_data[temp_cols] <- lapply(qpcr_data[temp_cols], function(x) { ifelse(x<=zero, NA, x) })
qpcr_data$Has.Hb[which(qpcr_data$HbtubCT1>0  | qpcr_data$HbtubCT2>0)]  <- TRUE
qpcr_data$Has.Pf[which(qpcr_data$pfr364CT1>0 | qpcr_data$pfr364CT2>0)] <- TRUE
counts["qpcr","missing"]     <- sum(is.na(qpcr_data$Has.Hb))
counts["qpcr","hb_positive"] <- sum(qpcr_data$Has.Hb, na.rm=TRUE)
counts["qpcr","pf_positive"] <- sum(qpcr_data$Has.Pf, na.rm=TRUE) - 1  # no parasitemia for M06 A0026
counts["qpcr","hb_negative"] <- nrow(qpcr_data) - counts["qpcr","hb_positive"] - counts["qpcr","missing"]
counts["qpcr","pf_negative"] <- nrow(qpcr_data) - counts["qpcr","pf_positive"] - counts["qpcr","missing"] + 1  # no parasitemia for M06 A0026


## -------- allspecies_data ----------------- ####

# Household ID
table(anopheles_widedata$household.id, useNA="always")
str(anopheles_widedata$household.id)
attr(anopheles_widedata$household.id, "labels")