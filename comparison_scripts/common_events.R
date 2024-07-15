#common event  

library(dplyr)
library(ggplot2)
library(tidyr)

##################Paths
path <- "~/stage/results_hg19/"
##################

event_rmats <- c("A3SS","A5SS","RI","SE","MXE")
rmats_results <- list()

for (event in event_rmats){
  rmats_results[[event]] <- read.delim(paste0(path,"rmats/", event, "_significant.MATS.JC.txt"))
}

event_vast <- c("ALTA", "ALTD", "INT", "EX")
vast_results <- list()

for (event in event_vast){
  vast_results[[event]] <- read.delim(paste0(path,"vast-tools/", event, "_significant.txt"))
}

##Number of DAS genes in common between the tools

a5ss_altd <- length(intersect(rmats_results[["A5SS"]]$geneSymbol, vast_results[["ALTD"]]$GENE))
a3ss_alta <- length(intersect(rmats_results[["A3SS"]]$geneSymbol, vast_results[["ALTA"]]$GENE))
ri_int <- length(intersect(rmats_results[["RI"]]$geneSymbol, vast_results[["INT"]]$GENE))
se_ex <- length(intersect(rmats_results[["SE"]]$geneSymbol, vast_results[["EX"]]$GENE))

list_genes_MXE_SE <- unique(c(rmats_results[["SE"]]$geneSymbol, rmats_results[["MXE"]]$geneSymbol))
se_mxe_ex <- length(intersect(list_genes_MXE_SE, vast_results[["EX"]]$GENE))

list_all_genes_rmats <- unique(c(rmats_results[["A5SS"]]$geneSymbol, rmats_results[["A3SS"]]$geneSymbol, rmats_results[["RI"]]$geneSymbol, rmats_results[["SE"]]$geneSymbol, rmats_results[["MXE"]]$geneSymbol))
list_all_genes_vast <- unique(c(vast_results[["ALTD"]]$GENE, vast_results[["ALTA"]]$GENE, vast_results[["RI"]]$GENE, vast_results[["EX"]]$GENE))

all <- length(intersect(list_all_genes_rmats, list_all_genes_vast))

common_summary_table <- data.frame(
  Common_DAS_genes = c(a3ss_alta,a5ss_altd,ri_int,se_ex,se_mxe_ex,all)
)
rownames(common_summary_table) <- c("A3SS_ALTA","A5SS_ALTD","RI_INT","SE_EX","SE+MXE_EX","ALL")

write.table(x = common_summary_table, file=paste0(path, "common_summary.tab"), sep = "\t", quote = FALSE, row.names=TRUE)

###################Comparison 

######plot

rmats_summary_table <- read.delim(paste0(path, "rmats_summary.tab"))
vast_summary_table <- read.delim(paste0(path, "vast_summary.tab"))

histo_table <- data.frame(
  event = c("A3SS","A5SS","RI","ES|EX","ES+MXE|EX"),
  rMATS = c(rmats_summary_table[-5,1],(rmats_summary_table[4,1]+rmats_summary_table[5,1])),
  VAST = c(vast_summary_table[,1],vast_summary_table[4,1])
)
histo_table <- histo_table %>% pivot_longer(!event, names_to="splicing_tool", values_to="number_of_events")

pdf(file = paste0(path,"events_histo.pdf"))
histo <-ggplot(data = histo_table) +
  #ggtitle("Number of differential alternative splicing events \n identified according to the tool used")+
  geom_bar(aes(x = event, y = number_of_events, fill = splicing_tool),stat = "identity", position="dodge")+
  scale_x_discrete(name = "Alternative splicing patterns")+
  scale_y_continuous(name = "Number of differential alternative splicing events")+
  scale_fill_discrete(name = "Splicing tool", labels = c("rMATS","Vast-tools"))
print(histo)
dev.off()

#################Common Differential Alternative Splicing events 

common_events <- read.delim(paste0(path,"common_events.txt"))
common_raw_events <- read.delim(paste0(path, "raw_common_events.txt"))
VastDB_EX <- read.delim(paste0(path, "vast-tools/VASTDB_EX.txt"))

common_events <- merge(common_events, vast_results[["EX"]][, c("EVENT", "coordinate","diff")], by.x = "Vast_ID", by.y = "EVENT", all.x = TRUE)
common_raw_events <- merge(common_raw_events, VastDB_EX[, c("EVENT", "coordinate")], by.x = "Vast_ID", by.y = "EVENT", all.x = TRUE)

colnames(common_events)[colnames(common_events) == "coordinate"] <- "Vast_coordinates"
colnames(common_events)[colnames(common_events) == "diff"] <- "Vast_diff"
colnames(common_raw_events)[colnames(common_raw_events) == "coordinate"] <- "Vast_coordinates"

common_events_filtered <- common_events[, c("rMats_event_ID", "Vast_coordinates", "Vast_diff","Vast_ID")]
common_raw_events_filtered <- common_raw_events[, c("rMats_event_ID", "Vast_coordinates", "Vast_ID")]

common_se_rmats <- merge(rmats_results[["SE"]], common_events_filtered, by.x = "ID", by.y = "rMats_event_ID")
common_raw_se_rmats <- merge(rmats_results[["SE"]], common_raw_events_filtered, by.x = "ID", by.y = "rMats_event_ID")

additional_insights <- c(paste0("Number of differential alternative splicing events detected by rMATS annotated in VASTDB: ", nrow(common_raw_se_rmats), " (before merging on VASTDB alternative exon coordinates) / ", nrow(unique(common_raw_se_rmats[,c(25,26)])), " (after merging)"),
                         paste0("Number of differential alternative splicing events detected in common: ", nrow(common_se_rmats), " with ", nrow(common_se_rmats[common_se_rmats$diff==common_se_rmats$Vast_diff,]), " common events sharing the same deltaPSI (before merging on VASTDB alternative exon coordinates) / ", nrow(unique(common_se_rmats[,c(3,25,26)])), " (after merging)" ))

insights <- file(paste0(path,"complementary_insights.txt"), open ="a") 
writeLines(additional_insights, insights)
close(insights)

#vérification
#length(intersect(unique(common_raw_se_rmats[,27]),unique(common_se_rmats[,28])))


