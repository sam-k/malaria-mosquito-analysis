# ----------------------------------------- #
#         Spat21 Data Set Analysis          #
#              Mosquito Data                #
#             January 29, 2019              #
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
library(ggplot2)
library(reshape2)


#### ---------------- set up environment --------------- ####
.wd <- "~/Projects/Malaria collab/Spatial R21 projects/Spat21 cleaning, analysis/"
MERGED_FP <- paste0(.wd, "Data/Mosquito Data Sets/moz_merged_data.Rdata")
LOG_FP    <- paste0(.wd, "Code/spat21_moz_data_analysis.log")
TAB_FP <- paste0(.wd, "Data/mosquito_tabulations.xlsx")
close(file(LOG_FP, open="w"))  # clear log file
write.log <- function(...) {
  for(output in list(...)) {
    write(output, file=LOG_FP, append=TRUE)
  }
  write("", file=LOG_FP, append=TRUE)
}
MONTHS <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")


#### -------------- read in mosquito data -------------- ####
load(MERGED_FP)  # allspecies_data, merged_data


#### ---------- tabulate all species dataset ----------- ####

write.log("# ------ TABULATE ALL SPP. DATA ------ #")
allspecies_data["month"] <- month(allspecies_data$collection.date)

# Tabulate abdominal statuses by village.
tab_allsp_a_v <- allspecies_data %>%
  select(village, anoph.unfed, anoph.bloodfed, anoph.halfgravid, anoph.gravid, anoph.undetermined, anoph.total) %>%
  group_by(village) %>%
  summarize_all(sum) %>%
  as.data.frame() %>%
  column_to_rownames(var="village")
write.table(tab_allsp_a_v, col.names=NA, file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()

# Tabulate abdominal statuses by collection month.
tab_allsp_a_m <- allspecies_data %>%
  select(month, anoph.unfed, anoph.bloodfed, anoph.halfgravid, anoph.gravid, anoph.undetermined, anoph.total) %>%
  group_by(month) %>%
  summarize_all(sum) %>%
  as.data.frame() %>%
  column_to_rownames(var="month")
rownames(tab_allsp_a_m) <- MONTHS
write.table(tab_allsp_a_m, col.names=NA, file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()

# Tabulate collection months by village.
tab_allsp_m_v <- allspecies_data %>%
  select(village, month, anoph.total) %>%
  group_by(village, month) %>%
  summarize_all(sum) %>%
  spread(month, anoph.total) %>%
  as.data.frame() %>%
  column_to_rownames(var="village")
colnames(tab_allsp_m_v) <- MONTHS
write.table(tab_allsp_m_v, col.names=NA, file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()


#### -------------- tabulate merged data --------------- ####

write.log("# ------ TABULATE MERGED DATA ------ #")

# Tabulate abdominal statuses by village.
tab_merged_a_v <- rbind(table(merged_data$village),
                        table(factor(merged_data$abdominal.status, levels=c("Blood Fed","Half Gravid","Gravid","Unfed","Undetermined")),
                              merged_data$village)) %>%
  cbind(Total=rowSums(.)) %>%
  as.data.frame()
row.names(tab_merged_a_v)[1] <- "Total female anoph collected"
tab_merged_a_v[["Total female anoph collected", "Total"]] <-
  sum(tab_merged_a_v["Total female anoph collected", 1:(which(colnames(tab_merged_a_v)=="Total")-1)])
write.table(tab_merged_a_v, col.names=NA, file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()

# Tabulate species by village.
tab_merged_s_v <- table(merged_data$species.type, merged_data$village) %>%
  cbind(Total=rowSums(.)) %>%
  as.data.frame() %>%
  rownames_to_column("Species.Type") %>%  # preserve rownames
  arrange(desc(Total)) %>%
  column_to_rownames("Species.Type")
.unid_row <- which(rownames(tab_merged_s_v)=="Un-identified")
tab_merged_s_v %<>%
  .[c(rownames(.)[1:(.unid_row-1)], rownames(.)[(.unid_row+1):nrow(.)], "Un-identified"), ] %>%  # reorder rows
  rbind(Other=colSums(.[5:(nrow(.)-1), ])) %>%
  .[c(rownames(.)[1:4], "Other", "Un-identified"), ]
write.table(tab_merged_s_v, col.names=NA, file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()

# Tabulate species counts per abdominal and infection status.
.tab_a_abd <- table(merged_data$species.type, merged_data$abdominal.status) %>%
  cbind(`Gravid or Half Gravid`=rowSums(.[, c("Gravid","Half Gravid")])) %>%
  as.data.frame() %>%
  select(`Blood Fed`, `Gravid or Half Gravid`)
make_ai_table <- function(row_name, col_name, crit_name, crit=TRUE) {
  table(merged_data[, row_name], merged_data[, col_name]) %>%
    as.data.frame() %>%
    `colnames<-`(c(row_name, col_name, crit_name)) %>%
    filter(.[, col_name]==crit) %>%
    column_to_rownames(row_name) %>%
    select(crit_name)
}
tab_merged_s_ai <- cbind(.tab_a_abd,
                         make_ai_table("species.type", "any.has.Hb", "Human Fed"),
                         make_ai_table("species.type", "H.has.Pf", "Infected (Head)"),
                         make_ai_table("species.type", "A.has.Pf", "Infected (Abdomen)"),
                         make_ai_table("species.type", "any.has.Pf", "Infected (Mosquito)")) %>%
  rownames_to_column("Species Type") %>%  # preserve rownames
  arrange(desc(`Blood Fed`)) %>%
  column_to_rownames("Species Type")
.unid_row <- which(rownames(tab_merged_s_ai)=="Un-identified")
tab_merged_s_ai %<>%
  .[c(rownames(.)[1:(.unid_row-1)], rownames(.)[(.unid_row+1):nrow(.)], "Un-identified"), ] %>%  # reorder rows
  rbind(Other=colSums(.[5:(nrow(.)-1), ])) %>%
  .[c(rownames(.)[1:4], "Other", "Un-identified"), ]
write.table(tab_merged_s_ai, col.names=NA, file=LOG_FP, append=TRUE, quote=FALSE, sep="\t")
write.log()


#### ----------------- visualize data ------------------ ####

write.log("# ------ VISUALIZE DATA ------ #")

# Correlation between parasitemias in replicates.
.temp_dat <- merged_data %>%
  select(village, H.pfr364CT1, H.pfr364CT2, A.pfr364CT1, A.pfr364CT2)
plot_pfr_corr <- ggplot(na.omit(.temp_dat)) +
  geom_point(aes(x=H.pfr364CT1, y=H.pfr364CT2, color="head"), size=1) +
  geom_smooth(aes(x=H.pfr364CT1, y=H.pfr364CT2, color="head"), method=lm) +
  geom_point(aes(x=A.pfr364CT1, y=A.pfr364CT2, color="abd"), size=1) +
  geom_smooth(aes(x=A.pfr364CT1, y=A.pfr364CT2, color="abd"), method=lm) +
  labs(x="Parasitemia CT rep1", y="Parasitemia CT rep2") +
  scale_color_manual(name="", values=c("head"="blue","abd"="red"), labels=c("Head","Abd"))
plot(plot_pfr_corr)

# Distribution of parasitemia for each village.
.temp_dat <- merged_data %>%
  select(village, H.pfr364CT1, H.pfr364CT2, A.pfr364CT1, A.pfr364CT2) %>%
  mutate_at(c("village"), factor)
plot_pfr_village <- ggplot(na.omit(.temp_dat)) +
  geom_boxplot(aes(x=village, y=H.pfr364CT1)) +
  labs(x="Village", y="Parasitemia CT")
plot(plot_pfr_village)

# Log-risk regression model to predict malaria infection prevalence.
model_regr <- glm(any.has.Pf~abdominal.status+village+species.type, family=binomial("logit"), data=merged_data)
summary(model_regr)
plot(model_regr)


#### ----------------- export analyses ----------------- ####

.wb <- createWorkbook()
addWorksheet(.wb, "allsp_tables")
writeData(.wb, sheet="allsp_tables", x=tab_allsp_a_m, rowNames=TRUE)
writeData(.wb, sheet="allsp_tables", x=tab_allsp_a_v, rowNames=TRUE, startCol=ncol(tab_allsp_a_m)+3)
writeData(.wb, sheet="allsp_tables", x=tab_allsp_m_v, rowNames=TRUE, startCol=ncol(tab_allsp_a_m)+3, startRow=nrow(tab_allsp_a_v)+3)
addWorksheet(.wb, "anoph_tables")
writeData(.wb, sheet="anoph_tables", x=tab_merged_a_v,  rowNames=TRUE)
writeData(.wb, sheet="anoph_tables", x=tab_merged_s_v,  rowNames=TRUE, startRow=nrow(tab_merged_a_v)+3)
writeData(.wb, sheet="anoph_tables", x=tab_merged_s_ai, rowNames=TRUE, startCol=ncol(tab_merged_a_v)+3)
saveWorkbook(.wb, TAB_FP, overwrite=TRUE)