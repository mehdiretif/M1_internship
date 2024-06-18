#common event  

common_events <- read.delim("~/stage/common_events.txt")
common_raw_events <- read.delim("~/stage/common_events_based_non_significant_vast.txt")

common_events <- merge(common_events, EX.vast[, c("EVENT", "coordinate","diff")], by.x = "Vast_ID", by.y = "EVENT", all.x = TRUE)
common_raw_events <- merge(common_raw_events, EX_raw_inclusion_file.DIFF[, c("EVENT", "coordinate")], by.x = "Vast_ID", by.y = "EVENT", all.x = TRUE)

colnames(common_events)[colnames(common_events) == "coordinate"] <- "Vast_coordinates"
colnames(common_events)[colnames(common_events) == "diff"] <- "Vast_diff"
colnames(common_raw_events)[colnames(common_raw_events) == "coordinate"] <- "Vast_coordinates"

common_events_filtered <- common_events[, c("rMats_event_ID", "Vast_coordinates", "Vast_diff","Vast_ID")]
common_raw_events_filtered <- common_raw_events[, c("rMats_event_ID", "Vast_coordinates", "Vast_ID")]

common_se_rmats <- merge(all_significant_rmats$SE, common_events_filtered, by.x = "ID", by.y = "rMats_event_ID")
common_raw_se_rmats <- merge(all_significant_rmats$SE, common_raw_events_filtered, by.x = "ID", by.y = "rMats_event_ID")

print("nombre d'evenements identifiés par rMATS retouvés dans base de données VASTDB, avant fusion des evements")
nrow(common_raw_se_rmats) 
print("après")
nrow(unique(common_raw_se_rmats[,c(25,26)]))

print("nombre d'evenements en commun avant fusion des evenements identifié par rMATS sur la base des coordonnées de VAST-tools")
nrow(common_se_rmats) 
print("après")
nrow(unique(common_se_rmats[,c(3,25,26)]))
print("avec correspondance au niveau du sens du differentiel")
nrow(common_se_rmats[common_se_rmats$diff==common_se_rmats$Vast_diff,])

#View(unique(common_se_rmats[,c(3,25,26)]))



