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
library(tibble)
library(openxlsx)


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


#### ---------- tabulate all species dataset ----------- ####

write.log("# ------ TABULATE ALL SPP. DATA ------ #")
allspecies_data["month"] <- month(allspecies_data$collection.date)

# Tabulate abdominal statuses by village.
tab_allsp_a_v <- allspecies_data %>%
  select(village, anoph.unfed, anoph.bloodfed, anoph.halfgravid, anoph.gravid, anoph.undetermined, anoph.total) %>%
  group_by(village) %>%
  summarize_all(sum)

# Tabulate abdominal statuses by collection month.
tab_allsp_a_m <- allspecies_data %>%
  select(month, anoph.unfed, anoph.bloodfed, anoph.halfgravid, anoph.gravid, anoph.undetermined, anoph.total) %>%
  group_by(month) %>%
  summarize_all(sum)

# Tabulate collection months by village.
tab_allsp_m_v <- allspecies_data %>%
  select(village, month, anoph.total) %>%
  group_by(village, month) %>%
  summarize_all(sum)


#### -------------- tabulate merged data --------------- ####

write.log("# ------ TABULATE MERGED DATA ------ #")

# Tabulate village and abdominal status.
tab_village_abd <- rbind(table(merged_data$village),
                         table(factor(merged_data$abdominal.status, levels=c("Blood Fed","Half Gravid","Gravid","Unfed","Undetermined")),
                               merged_data$village)) %>%
  cbind(Total=rowSums(.)) %>%
  as.data.frame()
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


#### ----------------- export analyses ----------------- ####

.wb <- createWorkbook()
addWorksheet(.wb, "SK_tables")
writeData(.wb, sheet="SK_tables", x=tab_village_abd, rowNames=TRUE)
writeData(.wb, sheet="SK_tables", x=tab_village_spp, rowNames=TRUE, startRow=nrow(tab_village_abd)+3)
writeData(.wb, sheet="SK_tables", x=tab_infection_spp, rowNames=TRUE, startCol=ncol(tab_village_abd)+3)
saveWorkbook(.wb, MERGED_TAB_FP, overwrite=TRUE)