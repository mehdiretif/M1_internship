library(dplyr)
library(ggplot2)
library(tidyr)

#rMATS output filtering 

path_rmats <- "~/stage/rmats/"
all_significant_rmats <- list()

event_rmats <- c("A3SS","A5SS","RI","SE","MXE")
event_type <- c()
significant_events <- c()
das_genes <- c()
list_all_genes_rmats <- c()

for (event in event_rmats)
{
  files_list <- list.files(path=path_rmats, pattern=paste0(event,".MATS.JC.txt"), recursive=TRUE, full.names=TRUE)
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
    write.table(x = all_significant_rmats[[event]], file=paste0(path_rmats, event, "_significant.MATS.JC.txt"), sep = "\t", quote = FALSE, row.names=FALSE)
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

print("das_all_genes:")
print(length(unique(sort(list_all_genes_rmats))))

#rMATS strict duplicates filtering 

rmats_SE_wthout_duplicats <- all_significant_rmats$SE[!duplicated(all_significant_rmats$SE[,-c(1,8:24)]),]
rmats_MXE_wthout_duplicats <- all_significant_rmats$MXE[!duplicated(all_significant_rmats$MXE[,-c(1,10:26)]),]

final_rmats_summary <- rmats_summary_table
final_rmats_summary["SE","Significant_Events"] <- nrow(rmats_SE_wthout_duplicats)
final_rmats_summary["MXE","Significant_Events"] <- nrow(rmats_MXE_wthout_duplicats)

final_rmats_summary["SE","DAS_Genes"] <- length(unique(rmats_SE_wthout_duplicats$geneSymbol))
final_rmats_summary["MXE","DAS_Genes"] <- length(unique(rmats_MXE_wthout_duplicats$geneSymbol))

#Vast filtering

path_vast <- "~/stage/vast/results/nf_results/"

vast_all <- read.delim(paste0(path_vast,"final_INCLUSION_FILE.DIFF.txt"))

vast_all <- vast_all[!is.na(vast_all$MV.dPsi._at_0.95),] #suppressed lines with na
vast_all <- vast_all[vast_all$MV.dPsi._at_0.95>=0.1,]

event_vast <- c("ALTA","ALTD","INT","EX")
event_type <- c()
significant_events <- c()
das_genes <- c()

for (event in event_vast) {
  subset_data <- vast_all[grep(event, vast_all$EVENT),]
  
  subset_data$diff <- ifelse(subset_data$E.dPsi.>0,'+','-')
  assign(paste0(event, ".vast"), subset_data)
  
  write.table(x = subset_data, file=paste0(path_vast, event, "_significant.txt"), sep = "\t", quote = FALSE, row.names=FALSE)
  
  event_type <- c(event_type, event)
  significant_events <- c(significant_events, nrow(subset_data))
  das_genes <- c(das_genes, length(unique(subset_data$GENE)))
}

vast_summary_table <- data.frame(
  Significant_Events = significant_events,
  DAS_Genes = das_genes
)  

rownames(vast_summary_table) <- event_type

print("all_das_genes_identified_by_vast")
print(length(unique(sort(vast_all$GENE))))

#number of DAS genes in common between the tools

a5ss_altd <- length(intersect(all_significant_rmats$A5SS$geneSymbol, ALTD.vast$GENE))
a3ss_alta <- length(intersect(all_significant_rmats$A3SS$geneSymbol, ALTA.vast$GENE))
ri_int <- length(intersect(all_significant_rmats$RI$geneSymbol, INT.vast$GENE))
se_ex <- length(intersect(rmats_SE_wthout_duplicats$geneSymbol, EX.vast$GENE))

list_genes_MXE_SE <- c(rmats_SE_wthout_duplicats$geneSymbol, rmats_MXE_wthout_duplicats$geneSymbol)
se_mxe_ex <- length(intersect(list_genes_MXE_SE, EX.vast$GENE))

all <- length(intersect(unique(list_all_genes_rmats), unique(sort(vast_all$GENE))))

common_summary_table <- data.frame(
  Common_DAS_genes = c(a3ss_alta,a5ss_altd,ri_int,se_ex,se_mxe_ex,all)
)
rownames(common_summary_table) <- c("a3ss_alta","a5ss_altd","ri_int","se_ex","se_mxe_ex","all")

setdiff(intersect(list_genes_MXE_SE, EX.vast$GENE),intersect(rmats_SE_wthout_duplicats$geneSymbol, EX.vast$GENE))

#View(rmats_summary_table)
#View(common_summary_table)
#View(final_rmats_summary)
#View(vast_summary_table)

#EX events from the inclusion table

final_INCLUSION_FILE.DIFF <- read.delim("~/stage/vast/results/INCLUSION_LEVELS_FULL-hg19-6.tab")

EX_raw_inclusion_file.DIFF <- final_INCLUSION_FILE.DIFF[grep("EX",final_INCLUSION_FILE.DIFF$EVENT),]
EX_raw_inclusion_file.DIFF <- EX_raw_inclusion_file.DIFF[EX_raw_inclusion_file.DIFF$COORD !="",] 
nrow(EX_raw_inclusion_file.DIFF) #nombre d'EX detectable
#write.table(x = EX_raw_inclusion_file.DIFF, file="~/stage/vast/results/nf_results/EX_raw_inclusion_file.DIFF.txt", sep = "\t", quote = FALSE, row.names=FALSE)

#plot

histo_table <- data.frame(
  event = c("A3SS","A5SS","RI","ES|EX","ES+MXE|EX"),
  rMATS = c(rmats_summary_table[-5,1],(rmats_summary_table[4,1]+rmats_summary_table[5,1])),
  VAST = c(vast_summary_table[,1],vast_summary_table[4,1])
)
histo_table <- histo_table %>% pivot_longer(!event, names_to="splicing_tool", values_to="number_of_events")

pdf(file = "~/stage/events_histo.pdf")
histo <-ggplot(data = histo_table) +
  #ggtitle("Number of differential alternative splicing events \n identified according to the tool used")+
  geom_bar(aes(x = event, y = number_of_events, fill = splicing_tool),stat = "identity", position="dodge")+
  scale_x_discrete(name = "Alternative splicing patterns")+
  scale_y_continuous(name = "Number of differential alternative splicing events")+
  scale_fill_discrete(name = "Splicing tool", labels = c("rMATS","Vast-tools"))
print(histo)
dev.off()



