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

.month_names <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
.old_abd_stats <- c("bloodfed","gravid","halfgravid","unfed","undetermined")
.new_abd_stats <- c("Blood-fed","Gravid","Half Gravid","Unfed","Undetermined")
allspecies_data %<>%
  mutate(year.month = paste0("'", substr(year(collection.date), 3, 4), " ", .month_names[as.integer(month(collection.date))]))

# Visualize summarized Anopheles data.
tab_allsp_anoph <- allspecies_data %>%
  select(year.month, village, anoph.unfed, anoph.bloodfed, anoph.halfgravid, anoph.gravid, anoph.undetermined) %>%
  group_by(village, year.month) %>%
  summarize_all(sum) %>%
  as.data.frame() %>%
  mutate(village = substr(village, 1, 1))
melted_allsp_anoph <- tab_allsp_anoph %>%
  melt(c("village","year.month"),
       c("anoph.bloodfed","anoph.gravid","anoph.halfgravid","anoph.unfed","anoph.undetermined")) %>%
  `colnames<-`(c("Village","Month","variable","value")) %>%
  mutate_at(c("variable"), as.character)
for(.i in 1:length(.old_abd_stats)) {
  melted_allsp_anoph[melted_allsp_anoph==paste0("anoph.",.old_abd_stats[.i])] <- .new_abd_stats[.i]
}
melted_allsp_anoph$variable %<>% factor(levels=rev(c("Blood-fed","Gravid","Half Gravid","Unfed","Undetermined")))
melted_allsp_anoph$Month    %<>% factor(levels=c(paste("'17",.month_names), paste("'18",.month_names)))
plot_allsp_anoph <- ggplot(melted_allsp_anoph, aes(x=Village, y=value, fill=variable)) + 
  geom_bar(stat="identity", position="stack") +
  facet_grid(~Month) +
  labs(x="Month", y="Number of mosquitoes") +
  scale_fill_brewer("Abdominal status", palette="YlOrRd") +
  ggtitle("Anopheles") +
  theme(plot.title=element_text(hjust=0.5))
plot(plot_allsp_anoph)

# Visualize summarized Culex data.
tab_allsp_culex <- allspecies_data %>%
  select(year.month, village, culex.unfed, culex.bloodfed, culex.halfgravid, culex.gravid, culex.undetermined) %>%
  group_by(village, year.month) %>%
  summarize_all(sum) %>%
  as.data.frame() %>%
  mutate(village = substr(village, 1, 1))
melted_allsp_culex <- tab_allsp_culex %>%
  melt(c("village","year.month"),
       c("culex.bloodfed","culex.gravid","culex.halfgravid","culex.unfed","culex.undetermined")) %>%
  `colnames<-`(c("Village","Month","variable","value")) %>%
  mutate_at(c("variable"), as.character)
for(.i in 1:length(.old_abd_stats)) {
  melted_allsp_culex[melted_allsp_culex==paste0("culex.",.old_abd_stats[.i])] <- .new_abd_stats[.i]
}
melted_allsp_culex$variable %<>% factor(levels=rev(c("Blood-fed","Gravid","Half Gravid","Unfed","Undetermined")))
melted_allsp_culex$Month    %<>% factor(levels=c(paste("'17",.month_names), paste("'18",.month_names)))
plot_allsp_culex <- ggplot(melted_allsp_culex, aes(x=Village, y=value, fill=variable)) + 
  geom_bar(stat="identity", position="stack") +
  facet_grid(~Month) +
  labs(x="Month", y="Number of mosquitoes") +
  scale_fill_brewer("Abdominal status", palette="PuRd") +
  ggtitle("Culex") +
  theme(plot.title=element_text(hjust=0.5))
plot(plot_allsp_culex)

# Visualize correlations between parasite densities (Q's) in replicates.
.merged_h <- merged_data %>%
  select(H.pfr364Q1, H.pfr364Q2) %>%
  `colnames<-`(c("pfr364Q1","pfr364Q2")) %>%
  mutate(Head.Abd = "H")
.merged_h %<>% .[!is.na(.$pfr364Q1) & !is.na(.$pfr364Q2), ]
.merged_a <- merged_data %>%
  select(A.pfr364Q1, A.pfr364Q2) %>%
  `colnames<-`(c("pfr364Q1","pfr364Q2")) %>%
  mutate(Head.Abd = "A")
.merged_a %<>% .[!is.na(.$pfr364Q1) & !is.na(.$pfr364Q2), ]
melted_reps <- rbind(.merged_h, .merged_a) %>%
  as.data.frame()
melted_reps[melted_reps=="H"] <- "Head"
melted_reps[melted_reps=="A"] <- "Abdomen"
plot_reps <- ggplot(melted_reps, aes(x=pfr364Q1, y=pfr364Q2, color=Head.Abd)) + 
  geom_point(alpha=0.5) +
  geom_smooth(method=lm, size=0.75) +
  coord_fixed(ratio=1, xlim=c(0, 0.002), ylim = c(0, 0.003)) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_color_brewer("Mosquito Part", palette="Set2")
plot(plot_reps)

# Correlation between parasitemias (CTs) in replicates.
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
table(merged_data$abdominal.status, useNA = "always")
exposure <- ifelse(merged_data$abdominal.status %in% c("Blood Fed", "Gravid", "Half Gravid"), "Blood Fed", "Unfed")
table(exposure, useNA="always")
merged_data$exposure <- exposure
merged_data$exposure = as.factor(merged_data$exposure)
merged_data$village = as.factor(merged_data$village)
merged_data$species.type = as.factor(merged_data$species.type)
merged_data$exposure = relevel(merged_data$exposure,"Unfed")
model_regr <- glm(any.has.Pf~exposure+village+species.type, family=binomial("logit"), data=merged_data)
summary(model_regr)
# plot(model_regr)
confint(model_regr)


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