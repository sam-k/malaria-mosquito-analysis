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
ALLSPECIES_FP         <- paste0(wd, "Data/Data Sets/MOZZIECollectionSummary_June2017_July2018.csv")
ANOPHELES_FP          <- paste0(wd, "Data/Data Sets/MOZZIEFemaleAnophele_June2017_July2018.csv")
QPCR_FP               <- paste0(wd, "Data/Data Sets/Mozzie mosquito compiled detection results 18Dec2018.csv")
DATA_DICT_FP          <- paste0(wd, "Data/Data Dictionary/spat21_data_mosquito_dictionary.csv")
CLEANED_ALLSPECIES_FP <- paste0(wd, "Data/Data Sets/cleaned_allspecies_data.csv")
CLEANED_ANOPHELES_FP  <- paste0(wd, "Data/Data Sets/cleaned_anopheles_data.csv")
CLEANED_QPCR_FP       <- paste0(wd, "Data/Data Sets/cleaned_qpcr_data.csv")
CLEANED_FP            <- paste0(wd, "Data/Data Sets/cleaned_data.Rdata")
zero <- 1e-6  # threshold for zero CT value


#### --------- read in mosquito data ----------------- ####

# Read in the mosquito descriptive data sets.
# Read in the data set with all mosquito species.
allspecies_data    <- read.csv(ALLSPECIES_FP, stringsAsFactors=FALSE)
# Read in the wide data set with only anopheles mosquitoes.
anopheles_widedata <- read.csv(ANOPHELES_FP, stringsAsFactors=FALSE)
# Read in the mosquito qPCR data sets.
qpcr_data          <- read.csv(QPCR_FP, stringsAsFactors=FALSE)

# Clean column names.
names(anopheles_widedata) <- tolower(gsub("..", ".", names(anopheles_widedata), fixed=TRUE))
names(anopheles_widedata) <- tolower(gsub("\\.$", "", names(anopheles_widedata)))  # remove trailing periods
anopheles_widedata %<>% rename(form.entered.date = form.entered.on)  # consistent name

# Look at summaries of all the data sets.
summary(allspecies_data)
summary(anopheles_widedata)
summary(qpcr_data)
str(allspecies_data)
str(anopheles_widedata)
str(qpcr_data)

# Output a CSV file of all the variable names.
allnames <- data.frame(c(names(allspecies_data), names(anopheles_widedata), names(qpcr_data)))
write_csv(allnames, DATA_DICT_FP)


#### ------------- clean each variable in mosquito data sets ---------------- ####

# Positive, negative, missing counts for each data set.
counts <- matrix(rep(0,15), nrow=3, ncol=5, dimnames=list(c("allspecies","anopheles","qpcr"),
                                                          c("hb_positive","hb_negative","pf_positive","pf_negative","missing")))

# Rename and reformat all_species_data columns.
names(allspecies_data) <- c("household.id","repeat.instrument","repeat.instance","collection.date","collection.time","village","collection.done.by",
                            "anoph.unfed","anoph.bloodfed","anoph.halfgravid","anoph.gravid","anoph.undetermined","anoph.total","num.male.anoph",
                            "culex.unfed","culex.bloodfed","culex.halfgravid","culex.gravid","culex.undetermined","culex.total","num.male.culex",
                            "form.checked.by","form.checked.date","form.entered.by","form.entered.date","complete")
allspecies_data %<>%
  mutate_at(c("household.id","repeat.instrument","village","collection.done.by","form.checked.by","form.entered.by","complete"), factor) %>%
  mutate_at(c("repeat.instance",
              "anoph.unfed","anoph.bloodfed","anoph.halfgravid","anoph.gravid","anoph.undetermined","anoph.total","num.male.anoph",
              "culex.unfed","culex.bloodfed","culex.halfgravid","culex.gravid","culex.undetermined","culex.total","num.male.culex"), as.integer) %>%
  mutate_at(c("collection.date","form.checked.date","form.entered.date"), mdy) %>%
  mutate(collection.time = as.logical(collection.time))

# Reformat anopheles_data columns from wide to long.
anopheles_data <- as.data.frame(matrix(nrow=16*nrow(anopheles_widedata), ncol=21), stringsAsFactors=FALSE)  # long data, overshooting # of rows
names(anopheles_data) <- c("household.id","repeat.instrument","repeat.instance","collection.date","collection.time","village",
                           "collection.done.by","samples.prepared.by","species.id.done.by","total.number.of.mosquitos.in.the.household",
                           "sample.id.head","sample.id.abdomen","abdominal.status","species.type","specify.species","comment",
                           "form.checked.by","form.checked.date","form.entered.by","form.entered.date","complete")
temp_count <- 1
for(i in 1:nrow(anopheles_widedata)) {
  header <- anopheles_widedata[i, 1:10]
  footer <- anopheles_widedata[i, 107:111]
  for(j in 1:16) {
    if(anopheles_widedata[[i, 5+6*j]] != "") {  # first column of j-th "block"
      anopheles_data[temp_count, ] <- c(header, anopheles_widedata[i, (5+6*j):(10+6*j)], footer)
      temp_count <- temp_count + 1
    }
  }
}
anopheles_data %<>% filter(!is.na(household.id))  # trim empty rows
anopheles_data[anopheles_data==""] <- NA
# Rename and reformat anopheles_data columns.
anopheles_data %<>%
  mutate_at(c("household.id","repeat.instrument","village","collection.done.by","samples.prepared.by","species.id.done.by",
              "abdominal.status","species.type","specify.species","comment",
              "form.checked.by","form.entered.by","complete"), factor) %>%
  mutate_at(c("repeat.instance","total.number.of.mosquitos.in.the.household"), as.integer) %>%
  mutate_at(c("collection.date","form.checked.date","form.entered.date"), mdy) %>%
  mutate(collection.time = as.logical(collection.time))


# Reformat qpcr_data columns.
qpcr_data[qpcr_data == "Undetermined"] <- NA
qpcr_data %<>%
  mutate_at(c("Sample.Name","Experiment.Name"), factor) %>%
  mutate_at(c("HbtubCT1","HbtubCT2","pfr364CT1","pfr364CT2","pfr364Std5a","pfr364Std5b","pfr364Std6a","pfr364Std6b"), as.numeric)
qpcr_data[is.na(qpcr_data)] <- NA  # correct NaNs to NAs
# Process qPCR CT values.
qpcr_data$Has.Hb <- FALSE
qpcr_data$Has.Pf <- FALSE
qpcr_data$Has.Hb[which(qpcr_data$HbtubCT1<=zero & qpcr_data$HbtubCT2<=zero & qpcr_data$pfr364CT1<=zero & qpcr_data$pfr364CT2<=zero)] <- NA
qpcr_data$Has.Pf[which(is.na(qpcr_data$Has.Hb))] <- NA
qpcr_data %<>%
  mutate_at(c("HbtubCT1","HbtubCT2","pfr364CT1","pfr364CT2"), function(x) { ifelse(x<=zero, NA, x) })
qpcr_data$Has.Hb[which(qpcr_data$HbtubCT1>0  | qpcr_data$HbtubCT2>0)]  <- TRUE
qpcr_data$Has.Pf[which(qpcr_data$pfr364CT1>0 | qpcr_data$pfr364CT2>0)] <- TRUE
counts["qpcr","missing"]     <- sum(is.na(qpcr_data$Has.Hb))
counts["qpcr","hb_positive"] <- sum(qpcr_data$Has.Hb, na.rm=TRUE)
counts["qpcr","pf_positive"] <- sum(qpcr_data$Has.Pf, na.rm=TRUE) - 1  # no parasitemia for M06 A0026
counts["qpcr","hb_negative"] <- nrow(qpcr_data) - counts["qpcr","hb_positive"] - counts["qpcr","missing"]
counts["qpcr","pf_negative"] <- nrow(qpcr_data) - counts["qpcr","pf_positive"] - counts["qpcr","missing"] + 1  # no parasitemia for M06 A0026


## -------- export cleaned data ----------------- ####

# Export cleaned data as CSV files.
write_csv(allspecies_data, CLEANED_ALLSPECIES_FP)
write_csv(anopheles_data, CLEANED_ANOPHELES_FP)
write_csv(qpcr_data, CLEANED_QPCR_FP)
# Export cleaned data as a single Rdata file.
save(allspecies_data, anopheles_data, qpcr_data, file=CLEANED_FP)