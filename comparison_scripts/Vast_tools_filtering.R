#############################Inputs
path_to_vast_out <- "~/stage/vast/vast_out_hg38/"
output_path <- "~/stage/results_hg38/" 
reference_genome <- "hg38"
number_of_sample <- "6" #total from all the conditions
conserve_events_without_coordinates <- "yes" #("yes"/"no") WARNING: if "yes", the identify_common_events.py script can encounter errors. 
#############################

#Vast filtering

dir.create(file.path(output_path, "vast-tools"), recursive = TRUE)

raw_vast_all <- read.delim(paste0(path_to_vast_out,"final_INCLUSION_FILE.DIFF.txt"))

vast_all <- raw_vast_all[!is.na(raw_vast_all$MV.dPsi._at_0.95),] #lines with NA suppression
vast_all <- vast_all[vast_all$MV.dPsi._at_0.95>=0.1,]
if(conserve_events_without_coordinates=="no"){
  vast_all <- vast_all[vast_all$coordinate !="",]
}

event_vast <- c("ALTA","ALTD","INT","EX")
event_type <- c()
significant_events <- c()
das_genes <- c()

for (event in event_vast) {
  subset_data <- vast_all[grep(event, vast_all$EVENT),]
  
  subset_data$diff <- ifelse(subset_data$E.dPsi.>0,'+','-')
  assign(paste0(event, ".vast"), subset_data)
  
  write.table(x = subset_data, file=paste0(output_path, "vast-tools/", event, "_significant.txt"), sep = "\t", quote = FALSE, row.names=FALSE)
  
  event_type <- c(event_type, event)
  significant_events <- c(significant_events, nrow(subset_data))
  das_genes <- c(das_genes, length(unique(subset_data$GENE)))
}

vast_summary_table <- data.frame(
  Significant_Events = significant_events,
  DAS_Genes = das_genes
)  

rownames(vast_summary_table) <- event_type

write.table(x = vast_summary_table, file=paste0(output_path, "vast_summary.tab"), sep = "\t", quote = FALSE, row.names=TRUE)

DAS_vast <- length(unique(sort(vast_all$GENE)))

#EX events from the inclusion table

combine_table <- read.delim(paste0(path_to_vast_out, "vast_out/INCLUSION_LEVELS_FULL-",reference_genome,"-", number_of_sample,".tab"))

EX_raw_inclusion_file.DIFF <- combine_table[grep("EX",combine_table$EVENT),]
all_EX_events <- raw_vast_all[raw_vast_all$EVENT %in% EX_raw_inclusion_file.DIFF$EVENT,]
all_EX_events <- all_EX_events[all_EX_events$coordinate !="",]

write.table(x = all_EX_events, file=paste0(output_path,"vast-tools/VASTDB_EX.txt"), sep = "\t", quote = FALSE, row.names=FALSE)

detectable_EX_events <- nrow(EX_raw_inclusion_file.DIFF) 

insights <- c(paste0(DAS_vast, " differentially alternatively spliced genes detected by Vast-tools"),
              paste0(detectable_EX_events, " alternative exon skipping events are annotated in VASTDB"))

file_path <- paste0(output_path,"complementary_insights.txt")

if(file.exists(file_path)){
  based_insights <- file(file_path, open ="a") 
  writeLines(insights, based_insights)
  close(based_insights)
} else{
  writeLines(insights, file_path)}



