#common event  

##################
path <- "~/stage/results/"
##################

common_events <- read.delim(paste0(path,"common_events.txt"))
common_raw_events <- read.delim(paste0(path, "raw_common_events.txt"))
EX.vast <- read.delim(paste0(path, "vast-tools/EX_significant.txt"))
VastDB_EX <- read.delim(paste0(path, "vast-tools/VASTDB_EX.txt"))
all_SE_significant_rmats <- read.delim(paste0(path, "rmats/SE_significant.MATS.JC.txt"))


common_events <- merge(common_events, EX.vast[, c("EVENT", "coordinate","diff")], by.x = "Vast_ID", by.y = "EVENT", all.x = TRUE)
common_raw_events <- merge(common_raw_events, VastDB_EX[, c("EVENT", "coordinate")], by.x = "Vast_ID", by.y = "EVENT", all.x = TRUE)

colnames(common_events)[colnames(common_events) == "coordinate"] <- "Vast_coordinates"
colnames(common_events)[colnames(common_events) == "diff"] <- "Vast_diff"
colnames(common_raw_events)[colnames(common_raw_events) == "coordinate"] <- "Vast_coordinates"

common_events_filtered <- common_events[, c("rMats_event_ID", "Vast_coordinates", "Vast_diff","Vast_ID")]
common_raw_events_filtered <- common_raw_events[, c("rMats_event_ID", "Vast_coordinates", "Vast_ID")]

common_se_rmats <- merge(all_SE_significant_rmats, common_events_filtered, by.x = "ID", by.y = "rMats_event_ID")
common_raw_se_rmats <- merge(all_SE_significant_rmats, common_raw_events_filtered, by.x = "ID", by.y = "rMats_event_ID")

additional_insights <- c(paste0("Number of differential alternative splicing events detected by rMATS annotated in VASTDB: ", nrow(common_raw_se_rmats), " (before merging on VASTDB alternative exon coordinates) / ", nrow(unique(common_raw_se_rmats[,c(25,26)])), " (after merging)"),
                         paste0("Number of differential alternative splicing events detected in common: ", nrow(common_se_rmats), " with ", nrow(common_se_rmats[common_se_rmats$diff==common_se_rmats$Vast_diff,]), " common events sharing the same deltaPSI (before merging on VASTDB alternative exon coordinates) / ", nrow(unique(common_se_rmats[,c(3,25,26)])), " (after merging)" ))

insights <- file(paste0(path,"complementary_insights.txt"), open ="a") 
writeLines(additional_insights, insights)
close(insights)

#vérification
#length(intersect(unique(common_raw_se_rmats[,27]),unique(common_se_rmats[,28])))

