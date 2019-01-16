# ----------------------------------------- #
#         Spat21 Data Set Cleaning          #
#              Mosquito Data                #
#             January 4, 2018               #
#            K. Sumner, S. Kim              #
# ----------------------------------------- #

#### ------------------ load packages ------------------ ####
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(magrittr)

#### ---------------- set up environment --------------- ####
.wd <- "~/Projects/Malaria collab/Spatial R21 projects/Spat21 cleaning, analysis/"
IMPORTED_FP <- paste0(.wd, "Data/Data Sets/imported_mosquito_data.Rdata")
CLEANED_FP  <- paste0(.wd, "Data/Data Sets/cleaned_mosquito_data.Rdata")
LOG_FP      <- paste0(.wd, "Code/spat21_data_cleaning_mosquitoes.log")
close(file(LOG_FP, open="w"))  # clear log file
write.log <- function(...) {
  for(.output in list(...)) {
    write(.output, file=LOG_FP, append=TRUE)
  }
  write("", file=LOG_FP, append=TRUE)
}
.village_dict <- list(K="Kinesamo", M="Maruti", S="Sitabicha")  # codes for villages


#### -------------- read in mosquito data -------------- ####
load(IMPORTED_FP)  # allspecies_data, anopheles_data, qpcr_data


#### -------------- clean all species data ------------- ####

write.log("# ------ CLEAN ALL SPP. DESCRIPTIVE DATA ------ #")

# Check if household IDs follow the correct format (X##).
discr_sp_hhformat <- allspecies_data %>%
  select(household.id) %>%
  filter(!grepl("^[KMS]\\d{2}$", household.id)) %>%
  arrange(household.id)
write.log("All household IDs are formatted correctly")

# Check if village names are consistent with household IDs.
discr_sp_villageid <- allspecies_data %>%
  select(household.id, village) %>%
  as.data.frame() %>%
  mutate_at(c("village"), as.character) %>%
  filter(village != .village_dict[substr(household.id,1,1)]) %>%
  arrange(household.id, village)
write.table(discr_sp_villageid, row.names=FALSE, col.names=c("HH","Village"),
            file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
allspecies_data$village <- .village_dict[substr(allspecies_data$household.id, 1, 1)]
write.log(paste("Village names for", nrow(discr_sp_villageid), "samples did not match household/sample IDs and were corrected"))

# Sort dataset.
allspecies_data %<>% arrange(collection.date, household.id)


#### -------------- clean descriptive data ------------- ####

write.log("# ------ CLEAN ANOPH. DESCRIPTIVE DATA ------ #")

# Check if household/sample IDs follow the correct format (X## X#####).
discr_an_hhformat <- anopheles_data %>%
  select(household.id) %>%
  filter(!grepl("^[KMS]\\d{2}$", household.id)) %>%
  arrange(household.id)
write.log("All household IDs are formatted correctly")
discr_an_idformat <- anopheles_data %>%
  select(sample.id.head, sample.id.abdomen) %>%
  filter(!(grepl("^[KMS]\\d{2}\\sH\\d{5}$", sample.id.head) & grepl("^[KMS]\\d{2}\\sA\\d{5}$", sample.id.abdomen))) %>%
  arrange(sample.id.head, sample.id.abdomen)
write.table(discr_an_idformat, row.names=FALSE, col.names=c("Sample ID H","Sample ID A"),
            file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
.h_incorrect_ids <- c( "K1 H00027","M05 A00002","M06 A00016","S06 A00012")
.h_correct_ids   <- c("K14 H00027","M05 H00002","M06 H00016","S06 H00012")
for(.i in 1:length(.h_incorrect_ids)) {
  anopheles_data$sample.id.head[anopheles_data$sample.id.head==.h_incorrect_ids[[.i]]] <- .h_correct_ids[[.i]]
}
.a_incorrect_ids <- c( "M09 A0097","M14 H00067","S06 H00012")
.a_correct_ids   <- c("M09 A00097","M14 A00067","S06 A00012")
for(.i in 1:length(.a_incorrect_ids)) {
  anopheles_data$sample.id.abdomen[anopheles_data$sample.id.abdomen==.a_incorrect_ids[[.i]]] <- .a_correct_ids[[.i]]
}
write.log("M05 00002, M06 00016, M09 00097, M14 00067, S06 00012 had incorrect head/abd designations and were corrected",
          "M09 A00097 was written as M09 A0097 and was corrected",
          "K14 H00027 was written as K1 H00027 and was corrected",
          "S02 A00001 is missing (as noted in comment)",
          "K05 A00034, M01 A00058, S07 A00008 are present in data, but commented as missing")

# Check if sample IDs are unique and consistent with each other.
discr_an_idid <- anopheles_data %>%
  select(sample.id.head, sample.id.abdomen) %>%
  filter((gsub("[AH]","",sample.id.head) != gsub("[AH]","",sample.id.abdomen))
       | duplicated(sample.id.head)    | duplicated(sample.id.head, fromLast=TRUE)
       | duplicated(sample.id.abdomen) | duplicated(sample.id.abdomen, fromLast=TRUE)) %>%
  arrange(sample.id.head, sample.id.abdomen)
write.table(discr_an_idid, row.names=FALSE, col.names=c("Sample ID H","Sample ID A"),
            file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
anopheles_data$sample.id.abdomen[anopheles_data$sample.id.head=="K05 H00008"] <- "K05 A00008"
anopheles_data %<>% .[-which(duplicated(.$sample.id.head, fromLast=TRUE)
                           | duplicated(.$sample.id.abdomen, fromLast=TRUE)), ]
write.log("K05 A00008 was written as K05 A00005 and was corrected",
          "Duplicate entries of M13 00035, M03 00021 were removed")

# Extract sample IDs.
anopheles_data$sample.id <- gsub("\\s*[AH]\\s*", " ", anopheles_data$sample.id.head)
anopheles_data$sample.id[is.na(anopheles_data$sample.id.head)] <-
  gsub("\\s*[AH]\\s*", " ", anopheles_data$sample.id.abdomen[is.na(anopheles_data$sample.id.head)])  # in case sample ID H is NA
write.log("Extracted sample IDs")

# Check if village names are consistent with household/sample IDs.
discr_an_villageid <- anopheles_data %>%
  select(household.id, village, sample.id) %>%
  as.data.frame() %>%
  mutate_at(c("village"), as.character) %>%
  filter((village!=.village_dict[substr(household.id,1,1)]) | (village!=.village_dict[substr(sample.id,1,1)])) %>%
  arrange(household.id, sample.id, village)
write.table(discr_an_villageid, row.names=FALSE, col.names=c("HH","Village","Sample ID"),
            file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
for(.i in 1:nrow(anopheles_data)) {
  anopheles_data$village[[.i]] <- .village_dict[substr(anopheles_data$household.id[[.i]], 1, 1)]
}
write.log(paste("Village names for", nrow(discr_an_villageid), "samples did not match household/sample IDs and were corrected"))

# Check if household names are consistent with sample IDs.
discr_an_hhid <- anopheles_data %>%
  select(household.id, sample.id.head, sample.id.abdomen) %>%
  filter((substr(household.id,2,3)!=substr(sample.id.head,2,3)) & (substr(household.id,2,3)!=substr(sample.id.abdomen,2,3))) %>%
  arrange(household.id, sample.id.head)
write.table(discr_an_hhid, row.names=FALSE, col.names=c("HH","Sample ID H","Sample ID A"),
            file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
anopheles_data$household.id <- substr(anopheles_data$sample.id.head, 1, 3)
write.log(paste("Household IDs for", nrow(discr_an_hhid), "samples did not match sample IDs and were overridden"))

# Check if abdominal statuses and species types are correct.
.discr_an_status <- anopheles_data %>%
  select(sample.id, abdominal.status) %>%
  filter(!(abdominal.status %in% c("Blood Fed","Gravid","Half Gravid","Undetermined","Unfed"))) %>%
  as.data.frame() %>%
  mutate_at(c("abdominal.status"), as.character) %>%
  replace(is.na(.), "NA")
.discr_an_species <- anopheles_data %>%
  select(sample.id, species.type) %>%
  filter(!grepl("^An\\. ", species.type) & !(species.type %in% c("Other, Specify","Un-identified"))) %>%
  as.data.frame() %>%
  mutate_at(c("species.type"), as.character) %>%
  replace(is.na(.), "NA")
discr_an_abdsp <- merge(.discr_an_status, .discr_an_species, by="sample.id", all=TRUE)
for(.col in c("abdominal.status","species.type")) {  # fill in corresponding info for discrepancies
  .na_col <- which(is.na(discr_an_abdsp[, .col]))
  discr_an_abdsp[.na_col, .col] <- as.character(anopheles_data[.na_col, .col])
}
discr_an_abdsp[discr_an_abdsp=="NA"] <- NA
write.table(discr_an_abdsp, row.names=FALSE, col.names=c("Sample ID","Abd Status","Species Type"),
            file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
.ids     <- c("K01 00026","M07 00011","M07 00031","M07 00047","M07 00062","M07 00092","M07 00109",
              "M07 00128","M07 00185","M09 00018","M09 00038","M09 00103","M14 00018","M14 00031")
.status  <- anopheles_data$abdominal.status[which(anopheles_data$sample.id %in% .ids)]
.species <- anopheles_data$species.type[which(anopheles_data$sample.id %in% .ids)]
anopheles_data$species.type[anopheles_data$sample.id %in% .ids]         <- .status
anopheles_data$abdominal.status[anopheles_data$sample.id %in% .ids]     <- .species
anopheles_data$abdominal.status[is.na(anopheles_data$abdominal.status)] <- "Undetermined"
write.log(paste("Abd statuses and species for", length(.ids), "samples appeared swapped and were corrected"),
          "Abd statuses for K05 00038, K14 00041 were missing and were corrected to Undetermined")

# Check if specified species and comments are correct.
discr_an_sppcomment <- anopheles_data %>%
  select(sample.id, species.type, specify.species, comment) %>%
  filter((!is.na(species.type)    & !grepl("^An\\. ", species.type) & species.type!="Un-identified")
       | (!is.na(specify.species) & !grepl("^An\\. ", specify.species))
       | (!is.na(comment)         &  grepl("^An\\. ", comment))) %>%
  as.data.frame() %>%
  arrange(species.type, sample.id, specify.species, comment)
write.table(discr_an_sppcomment, row.names=FALSE, col.names=c("Sample ID","Species Type","Specify Sp.","Comment"),
            file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
.which_ids <- which(anopheles_data$sample.id %in% c("M07 00013","M07 00094","M14 00020"))
anopheles_data$comment[.which_ids] <- as.character(anopheles_data$specify.species[.which_ids])
anopheles_data$specify.species[.which_ids] <- NA
.which_ids <- which(anopheles_data$sample.id %in% c("K02 00028","M03 00025","M07 00135","M07 00136","M09 00086","M14 00038","M14 00059",
                                                    "M15 00027","M15 00033","M15 00061","M15 00067","M16 00007","S01 00016"))
anopheles_data %<>% mutate_at(c("species.type"), as.character)
anopheles_data$species.type[.which_ids] <- as.character(anopheles_data$specify.species[.which_ids])
anopheles_data$specify.species[.which_ids] <- NA
anopheles_data %<>% mutate_at(c("species.type"), factor)
write.log("Specified species and comments for M07 00013, M07 00094, M14 00020 appeared swapped and were corrected",
          paste("Duplicate species rows for", nrow(discr_an_sppcomment)-3, "samples were merged (An. pretoriensis)"))

# Sort dataset and reorder columns.
.h_id_col <- which(names(anopheles_data)=="sample.id.head")
anopheles_data %<>% .[c(names(.)[1:(.h_id_col-1)], "sample.id", names(.)[.h_id_col:(ncol(.)-1)])] %>%
  droplevels() %>%  # remove empty levels
  arrange(sample.id)


#### ----------------- clean qPCR data ----------------- ####

write.log("# ------ CLEAN QPCR DATA ------ #")

# Check if sample IDs follow the correct format (X## X#####).
discr_qp_hhformat <- qpcr_data %>%
  select(Sample.Name) %>%
  filter(!grepl("^[KMS]\\d{2}\\s[AH]\\d{5}$", Sample.Name)) %>%
  arrange(Sample.Name)
write.log("All sample IDs are formatted correctly")

# Extract sample IDs and head/abd statuses.
qpcr_data$Sample.ID <- gsub("\\s?[AH]\\s?", " ", qpcr_data$Sample.Name)
qpcr_data$Head.Abd  <- gsub("[^AH]", "", qpcr_data$Sample.Name)
qpcr_data %<>% mutate(Head.Abd=factor(Head.Abd))
write.log("Extracted sample IDs and heads/abdomens")

# Sort dataset and reorder columns.
qpcr_data %<>% .[c("Sample.ID", "Head.Abd", names(qpcr_data)[1:26])] %>% arrange(Sample.ID, Head.Abd)


#### --------------- export cleaned data --------------- ####
save(allspecies_data, anopheles_data, qpcr_data, file=CLEANED_FP)