# ----------------------------------------- #
#         Spat21 Data Set Merging           #
#              Mosquito Data                #
#             January 4, 2018               #
#                  S. Kim                   #
# ----------------------------------------- #

#### ------------------ load packages ------------------ ####
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(magrittr)
library(tibble)
library(openxlsx)


#### ---------------- set up environment --------------- ####
.wd <- "~/Projects/Malaria collab/Spatial R21 projects/Spat21 cleaning, analysis/"
CLEANED_FP    <- paste0(.wd, "Data/Data Sets/cleaned_mosquito_data.Rdata")
MERGED_CSV_FP <- paste0(.wd, "Data/Data Sets/merged_mosquito_data.csv")
MERGED_RDS_FP <- paste0(.wd, "Data/Data Sets/merged_mosquito_data.rds")
MERGED_TAB_FP <- paste0(.wd, "Data/Merged mosquito tabulations.xlsx")
LOG_FP        <- paste0(.wd, "Code/spat21_data_merging_mosquitoes.log")
close(file(LOG_FP, open="w"))  # clear log file
write.log <- function(...) {
  for(output in list(...)) {
    write(output, file=LOG_FP, append=TRUE)
  }
  write("", file=LOG_FP, append=TRUE)
}


#### -------------- read in mosquito data -------------- ####

load(CLEANED_FP)  # allspecies_data, anopheles_data, qpcr_data
allspecies_data %<>% arrange(collection.date, household.id)
anopheles_data  %<>% arrange(sample.id)
qpcr_data       %<>% arrange(Sample.ID, Head.Abd)

#### ------------- merge mosquito data sets ------------ ####

write.log("# ------ MERGE MOSQUITO DATA ------ #")

# Group qpcr_data by sample ID.
qpcr_groupeddata <- as.data.frame(matrix(nrow=nrow(qpcr_data), ncol=17), stringsAsFactors=FALSE)  # overshoot # of rows
names(qpcr_groupeddata) <- c("sample.id",
                             "H.HbtubCT1","H.HbtubCT2","H.has.Hb","H.pfr364CT1","H.pfr364CT2","H.pfr364Q1","H.pfr364Q2","H.has.Pf",
                             "A.HbtubCT1","A.HbtubCT2","A.has.Hb","A.pfr364CT1","A.pfr364CT2","A.pfr364Q1","A.pfr364Q2","A.has.Pf")
.count <- 0
for(.i in 1:nrow(qpcr_data)) {
  if(qpcr_data[[.i, "Sample.ID"]] != ifelse(.i>1, qpcr_data[[.i-1, "Sample.ID"]], "")) {
    .count <- .count + 1
  }
  qpcr_groupeddata[[.count, "sample.id"]] <- qpcr_data[[.i, "Sample.ID"]]
  if(qpcr_data[[.i, "Head.Abd"]] == "H") {
    qpcr_groupeddata[.count, 2:9]   <- qpcr_data[.i, 6:13]
  } else if(qpcr_data[[.i, "Head.Abd"]] == "A") {
    qpcr_groupeddata[.count, 10:17] <- qpcr_data[.i, 6:13]
  }
}
qpcr_groupeddata %<>% filter(!is.na(sample.id))  # trim empty rows
write.log("Converted qPCR data to wide format by sample ID")

# Combine head/abdomen data for Hb/Pf statuses.
qpcr_groupeddata$any.has.Hb <- qpcr_groupeddata$H.has.Hb | qpcr_groupeddata$A.has.Hb
qpcr_groupeddata$any.has.Hb[is.na(qpcr_groupeddata$any.has.Hb)] <- FALSE  # NAs should be false
qpcr_groupeddata$any.has.Pf <- qpcr_groupeddata$H.has.Pf | qpcr_groupeddata$A.has.Pf
qpcr_groupeddata$any.has.Pf[is.na(qpcr_groupeddata$any.has.Pf)] <- FALSE  # NAs should be false

# Merge anopheles descriptive data with qPCR data.
merged_data <- left_join(anopheles_data, qpcr_groupeddata, by="sample.id") %>%
  mutate_at(c("village"), as.character) %>%
  select(-c(repeat.instrument, repeat.instance, collection.date, collection.time, total.number.of.mosquitos.in.the.household,
            collection.done.by, samples.prepared.by, species.id.done.by, sample.id.head, sample.id.abdomen,
            specify.species, comment, form.checked.by, form.checked.date, form.entered.by, form.entered.date, complete))
write.log("Merged anopheles descriptive data with wide qPCR data")


#### -------------- validate merged data --------------- ####

write.log("# ------ VALIDATE MERGING ------ #")

# Check if any anopheles descriptive entries were not merged.
unmerged_anoph <- merged_data %>%
  select(sample.id, any.has.Hb, H.has.Hb, A.has.Hb, any.has.Pf, H.has.Pf, A.has.Pf) %>%
  filter(is.na(any.has.Hb) | is.na(any.has.Pf)) %>%
  arrange(sample.id) %>%
  as.data.frame()
write.table(unmerged_anoph, row.names=FALSE, col.names=c("Sample ID","Any Hb","H Hb","A Hb","Any Pf","H Pf","A Pf"),
            file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
write.log("If any.has.XX is NA, that entry was not present in the anopheles descriptive data",
          paste("From the anopheles descriptive dataset,", nrow(unmerged_anoph), "entries did not merge and were discarded from the data"))

# Check if any qPCR entries were not merged.
unmerged_qpcr <- qpcr_groupeddata %>%
  select(sample.id, any.has.Hb, H.has.Hb, A.has.Hb, any.has.Pf, H.has.Pf, A.has.Pf) %>%
  filter(!(qpcr_groupeddata$sample.id %in% merged_data$sample.id)) %>%
  arrange(sample.id) %>%
  as.data.frame()
write.table(unmerged_qpcr, row.names=FALSE, col.names=c("Sample ID","Any Hb","H Hb","A Hb","Any Pf","H Pf","A Pf"),
            file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()
write.log(paste("From the qPCR dataset,", nrow(unmerged_qpcr), "entries were absent in the descriptive data and did not merge"),
          "Samples were absent from shipments and were discarded from the data")


#### -------------- tabulate merged data --------------- ####

write.log("# ------ TABULATE MERGED DATA ------ #")

# Tabulate village and abdominal status.
tab_village_abd <- rbind(table(merged_data$village),
                         table(factor(merged_data$abdominal.status, levels=c("Blood Fed","Half Gravid","Gravid","Unfed","Undetermined")),
                               merged_data$village)) %>%
  cbind(Total=rowSums(.))
row.names(tab_village_abd)[1] <- "Total female anoph collected"
tab_village_abd[["Total female anoph collected", "Total"]] <-
  sum(tab_village_abd["Total female anoph collected", 1:(which(colnames(tab_village_abd)=="Total")-1)])
write.table(tab_village_abd, col.names=NA, file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()

# Tabulate village and species.
tab_village_spp <- table(merged_data$species.type, merged_data$village) %>%
  cbind(Total=rowSums(.)) %>%
  as.data.frame() %>%
  rownames_to_column("Species.Type") %>%  # preserve rownames
  arrange(desc(Total)) %>%
  column_to_rownames("Species.Type")
.unid_row <- which(rownames(tab_village_spp)=="Un-identified")
tab_village_spp %<>%
  .[c(rownames(.)[1:(.unid_row-1)], rownames(.)[(.unid_row+1):nrow(.)], "Un-identified"), ] %>%  # reorder rows
  rbind(Other=colSums(.[5:(nrow(.)-1), ])) %>%
  .[c(rownames(.)[1:4], "Other", "Un-identified"), ]
write.table(tab_village_spp, col.names=NA, file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()

# Tabulate species counts per abdominal and infection status.
.tab_a_abd <- table(merged_data$species.type, merged_data$abdominal.status) %>%
  cbind(`Gravid or Half Gravid`=rowSums(.[, c("Gravid","Half Gravid")])) %>%
  as.data.frame() %>%
  select(`Blood Fed`, `Gravid or Half Gravid`)
make_a_inf_table <- function(row_name, col_name, crit_name, crit=TRUE) {
  table(merged_data[, row_name], merged_data[, col_name]) %>%
    as.data.frame() %>%
    `colnames<-`(c(row_name, col_name, crit_name)) %>%
    filter(.[, col_name]==crit) %>%
    column_to_rownames(row_name) %>%
    select(crit_name)
}
tab_infection_spp <- cbind(.tab_a_abd,
                           make_a_inf_table("species.type", "any.has.Hb", "Human Fed"),
                           make_a_inf_table("species.type", "H.has.Pf", "Infected (Head)"),
                           make_a_inf_table("species.type", "A.has.Pf", "Infected (Abdomen)"),
                           make_a_inf_table("species.type", "any.has.Pf", "Infected (Mosquito)")) %>%
  rownames_to_column("Species Type") %>%  # preserve rownames
  arrange(desc(`Blood Fed`)) %>%
  column_to_rownames("Species Type")
.unid_row <- which(rownames(tab_infection_spp)=="Un-identified")
tab_infection_spp %<>%
  .[c(rownames(.)[1:(.unid_row-1)], rownames(.)[(.unid_row+1):nrow(.)], "Un-identified"), ] %>%  # reorder rows
  rbind(Other=colSums(.[5:(nrow(.)-1), ])) %>%
  .[c(rownames(.)[1:4], "Other", "Un-identified"), ]
write.table(tab_infection_spp, col.names=NA, file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()


#### ---------------- export merged data --------------- ####

# Export data.
write.csv(merged_data, file=MERGED_CSV_FP, row.names=FALSE)
saveRDS(merged_data, file=MERGED_RDS_FP)

# Export tabulations.
.wb <- createWorkbook()
addWorksheet(.wb, "SK_tables")
writeData(.wb, sheet="SK_tables", x=tab_village_abd, rowNames=TRUE)
writeData(.wb, sheet="SK_tables", x=tab_village_spp, rowNames=TRUE, startRow=nrow(tab_village_abd)+3)
writeData(.wb, sheet="SK_tables", x=tab_infection_spp, rowNames=TRUE, startCol=ncol(tab_village_abd)+3)
saveWorkbook(.wb, MERGED_TAB_FP, overwrite=TRUE)