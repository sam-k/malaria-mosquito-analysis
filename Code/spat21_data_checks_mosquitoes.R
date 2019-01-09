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
CLEANED_FP   <- paste0(wd, "Data/Data Sets/cleaned_data.Rdata")
VALIDATED_FP <- paste0(wd, "Data/Data Sets/validated_data.Rdata")
LOG_FP       <- paste0(wd, "Code/spat21_data_checks_mosquitoes.log")
close(file(LOG_FP, open="w"))  # clear log file
write.log <- function(...) {
  for(temp_output in list(...)) {
    write(temp_output, file=LOG_FP, append=TRUE)
  }
  write("", file=LOG_FP, append=TRUE)
}


#### --------- read in mosquito data ----------------- ####
load(CLEANED_FP)  # allspecies_data, anopheles_data, qpcr_data


#### -------- validate descriptive data ----------------- ####

write.log("# ------ VALIDATE ANOPH. DESCRIPTIVE DATA ------ #")

# Check if village names are consistent.
discr_an_village <- anopheles_data %>%
  select(household.id, village, sample.id.head, sample.id.abdomen) %>%
  filter((substr(household.id,1,1)!=substr(sample.id.head,1,1)) & (substr(household.id,1,1)!=substr(sample.id.abdomen,1,1)) &
         (substr(household.id,1,1)!=substr(village,1,1))) %>%
  arrange(household.id, sample.id.head)
write.log("All village names match")

# Check if household names are consistent.
discr_an_hh <- anopheles_data %>%
  select(household.id, sample.id.head, sample.id.abdomen) %>%
  filter((substr(household.id,2,3)!=substr(sample.id.head,2,3)) & (substr(household.id,2,3)!=substr(sample.id.abdomen,2,3))) %>%
  arrange(household.id, sample.id.head)
write.table(discr_an_hh, row.names=FALSE, col.names=c("HH ID","Sample ID H","Sample ID A"),
            file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
temp_h_ids <- sort(discr_an_hh$sample.id.head)
anopheles_data$household.id <- substr(anopheles_data$sample.id.head, 1, 3)
write.log(paste("Household IDs for", paste(gsub("H","",temp_h_ids), collapse=", "), "did not match sample IDs and were overridden"))

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
          "K05 A00034, M01 A00058, S07 A00008 are present in data, but commented as missing")

# Check if last five numbers of sample IDs are consistent.
discr_an_id <- anopheles_data %>%
  select(sample.id.head, sample.id.abdomen) %>%
  filter(substr(sample.id.head,6,10)!=substr(sample.id.abdomen,6,10)) %>%
  arrange(sample.id.head)
write.table(discr_an_id, row.names=FALSE, col.names=c("Sample ID H","Sample ID A"),
            file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
anopheles_data$sample.id.abdomen[anopheles_data$sample.id.abdomen=="M09 A0097"] <- "M09 A00097"
write.log("K05 H00008 and K05 A00005 are grouped together in the same entry",
          "M09 A0097 appeared to be formatted incorrectly and was corrected")

# Extracted sample IDs.
anopheles_data$sample.id <- gsub("\\s?[AH]\\s?", " ", anopheles_data$sample.id.head)
anopheles_data$sample.id[is.na(anopheles_data$sample.id.head)] <-
  gsub("\\s?[AH]\\s?", " ", anopheles_data$sample.id.abdomen[is.na(anopheles_data$sample.id.head)])  # in case sample H ID is NA
write.log("Extracted sample IDs")

# Check if abdominal statuses and species types are correct.
temp_discr_an_status <- anopheles_data %>%
  select(sample.id, abdominal.status) %>%
  filter(!(abdominal.status %in% c("Blood Fed","Gravid","Half Gravid","Undetermined","Unfed"))) %>%
  as.data.frame() %>%
  mutate_at(c("abdominal.status"), as.character) %>%
  replace(is.na(.), "NA")
temp_discr_an_species <- anopheles_data %>%
  select(sample.id, species.type) %>%
  filter(!grepl("An. ", species.type) & !(species.type %in% c("Other, Specify","Un-identified"))) %>%
  as.data.frame() %>%
  mutate_at(c("species.type"), as.character) %>%
  replace(is.na(.), "NA")
discr_an_abdsp <- merge(temp_discr_an_status, temp_discr_an_species, by="sample.id", all=TRUE)
for(temp_col in c("abdominal.status","species.type")) {  # fill in corresponding info for discrepancies
  temp_na_col <- which(is.na(discr_an_abdsp[, temp_col]))
  discr_an_abdsp[temp_na_col, temp_col] <- as.character(anopheles_data[temp_na_col, temp_col])
}
discr_an_abdsp[discr_an_abdsp=="NA"] <- NA
write.table(discr_an_abdsp, row.names=FALSE, col.names=c("Sample ID","Abd Status","Species Type"),
            file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
temp_ids <- c("K01 00026","M07 00011","M07 00031","M07 00047","M07 00062","M07 00092","M07 00109",
              "M07 00128","M07 00185","M09 00018","M09 00038","M09 00103","M14 00018","M14 00031")
temp_status  <- anopheles_data$abdominal.status[which(anopheles_data$sample.id %in% temp_ids)]
temp_species <- anopheles_data$species.type[which(anopheles_data$sample.id %in% temp_ids)]
anopheles_data$species.type[anopheles_data$sample.id %in% temp_ids]     <- temp_status
anopheles_data$abdominal.status[anopheles_data$sample.id %in% temp_ids] <- temp_species
anopheles_data$abdominal.status[is.na(anopheles_data$abdominal.status)] <- "Undetermined"
write.log(paste("Abd statuses and species for", paste(temp_ids, collapse=", "), "appeared swapped and were corrected"),
          "Abd statuses for K05 00038, K14 00041 were missing and were corrected to Undetermined")


#### -------- validate qPCR data ----------------- ####

write.log("# ------ VALIDATE QPCR DATA ------ #")

# Perform data checks for qpcr_data.
qpcr_data$Sample.ID <- gsub("\\s?[AH]\\s?", " ", qpcr_data$Sample.Name)
qpcr_data$Head.Abd  <- gsub("[^AH]", "", qpcr_data$Sample.Name)
write.log("Extracted sample IDs and heads/abdomens")
qpcr_data %<>%
  mutate(Head.Abd = factor(Head.Abd)) %>%
  arrange(Sample.ID, Head.Abd)
write.log("No parasitemia for M06 A00026, considered missing")


#### -------- export validated data ----------------- ####
save(allspecies_data, anopheles_data, qpcr_data, file=VALIDATED_FP)


#### -------- clean up environment ----------------- ####
rm(list=ls(pattern="^temp_"))