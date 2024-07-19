# SPDX-FileCopyrightText: 2024 Marchand Mehdi <mehdi.retif@gmail.com>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

############################Inputs
path_to_rmats_raw_results_directory <- "~/stage/rmats/"
output_path <- "~/stage/results_hg19/" 
############################

dir.create(file.path(output_path, "rmats"), recursive = TRUE)

all_significant_rmats <- list()

event_rmats <- c("A3SS","A5SS","RI","SE","MXE")
event_type <- c()
significant_events <- c()
das_genes <- c()
list_all_genes_rmats <- c()

for (event in event_rmats)
{
  files_list <- list.files(path=path_to_rmats_raw_results_directory, pattern=paste0(event,".MATS.JC.txt"), recursive=TRUE, full.names=TRUE)
  for (pattern in files_list)
  {
    input <- read.delim(pattern, sep="\t", header=TRUE)
    
    input$basemean <- apply(input, 1, function(row) {
      samples <- c(as.numeric(unlist(strsplit(row['IJC_SAMPLE_1'], ','))),
                   as.numeric(unlist(strsplit(row['SJC_SAMPLE_1'], ','))),
                   as.numeric(unlist(strsplit(row['IJC_SAMPLE_2'], ','))),
                   as.numeric(unlist(strsplit(row['SJC_SAMPLE_2'], ','))))
      mean(samples, na.rm = TRUE)
    })
    
    input$diff <- ifelse(input$IncLevelDifference>0,'+','-')
    
    input <- input[input$basemean >=20,]
    input <- input[input$FDR <= 0.05,]
    input <- input[abs(input$IncLevelDifference)>=0.1,]
    
    all_significant_rmats[[event]] <- input
    write.table(x = all_significant_rmats[[event]], file=paste0(output_path,"rmats/", event, "_significant.MATS.JC.txt"), sep = "\t", quote = FALSE, row.names=FALSE)
  }
  
  event_type <- c(event_type, event)
  significant_events <- c(significant_events, nrow(all_significant_rmats[[event]]))
  das_genes <- c(das_genes, length(unique(all_significant_rmats[[event]]$geneSymbol)))
  #comptage das genes pour all rmats events
  list_all_genes_rmats <- c(list_all_genes_rmats, all_significant_rmats[[event]]$geneSymbol)
}

rmats_summary_table <- data.frame(
  Significant_Events = significant_events,
  DAS_Genes = das_genes
)

rownames(rmats_summary_table) <- event_type
write.table(x = rmats_summary_table, file=paste0(output_path, "rmats_summary.tab"), sep = "\t", quote = FALSE, row.names=TRUE)


##rMATS strict duplicates filtering 

rmats_SE_wthout_duplicats <- all_significant_rmats$SE[!duplicated(all_significant_rmats$SE[,-c(1,8:24)]),]
rmats_MXE_wthout_duplicats <- all_significant_rmats$MXE[!duplicated(all_significant_rmats$MXE[,-c(1,10:26)]),]

final_rmats_summary <- rmats_summary_table
final_rmats_summary["SE","Significant_Events"] <- nrow(rmats_SE_wthout_duplicats)
final_rmats_summary["MXE","Significant_Events"] <- nrow(rmats_MXE_wthout_duplicats)

final_rmats_summary["SE","DAS_Genes"] <- length(unique(rmats_SE_wthout_duplicats$geneSymbol))
final_rmats_summary["MXE","DAS_Genes"] <- length(unique(rmats_MXE_wthout_duplicats$geneSymbol))

write.table(x = final_rmats_summary, file=paste0(output_path, "final_rmats_summary.tab"), sep = "\t", quote = FALSE, row.names=TRUE)

rmats_DAS <- length(unique(sort(list_all_genes_rmats)))
insight <- paste0(rmats_DAS, " differentially alternatively spliced genes detected by rMATS")

file_path <- paste0(output_path,"complementary_insights.txt")

if(file.exists(file_path)){
  based_insights <- file(file_path, open ="a") 
  writeLines(insight, based_insights)
  close(based_insights)
} else{
  writeLines(insight, file_path)}
       
  
       

