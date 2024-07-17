## Vast-tools Nextflow pipeline

### Options

```
    --reads     '/path/to/read/*_files.fastq'
    --db        '/path/to/VASTDB' 
    --groups    '/path/to/file/containing/insights/about/samples' (refer to the groups.txt file) 
    --groupA    'conditionA' 
    --groupB    'conditionB' 
    --gen_ref   'genome_reference' 
    --cores     'number_of_cores' 
    --output    '/path/to/output/directory'
    --expr      'yes/no' (differential expression analysis)
``` 

The databases are available on the [vast-tools github](https://github.com/vastgroup/vast-tools).

### Example

```
nextflow run /Xnfs/abc/nf_scratch/mmarchand/vast_projet/vast_tools.nf \
-profile psmn \
-c /Xnfs/abc/nf_scratch/mmarchand/vast_projet/nextflow.config \
--reads '/Xnfs/abc/nf_scratch/mmarchand/vast_projet/trimmed_fastq_test/*_R{1,2}_cutadapt_match.fastq.gz' \
--db '/Xnfs/abc/nf_scratch/mmarchand/vast_projet/VASTDB/' \
--groups '/Xnfs/abc/nf_scratch/mmarchand/vast_projet/groups.txt' \
--groupA 'siDDX5_17' \
--groupB 'siGL2' \
--gen_ref 'hg19' \
--cores '16' \
--output '/Xnfs/abc/nf_scratch/mmarchand/vast_projet/vast_out/'
--expr 'yes'
```

Please use ```_``` after ```*``` for the paths of the fastq files. If the directory contains more fastq files than those of interest, you can create a new directory with symbolic link directing to the fastq files of interest. 

The groups file contains the names of the fastq files as first column (example: 5Y_siDDX5_17_B1 for 5Y_siDDX5_17_B1_R1_cutadapt_match.fastq.gz). The second column contains condition IDs (groupA and groupB).

The nextflow pipeline generate two final files, namely **final_INCLUSION_FILE.DIFF.txt** for the differential alternative splicing analysis, and the **DiffGE.tab** for the differential expression analysis (with the --expr option).

## Command lines used to perform Differential Alternative Splicing (DAS) analysis with Vast-tools

#### align

```
vast-tools align path/to/5Y_siDDX5_17_B1_R1_cutadapt_match.fastq.gz path/to/5Y_siDDX5_17_B1_R2_cutadapt_match.fastq.gz -sp genome_reference (hg19/hg38...) --dbDir path/to/VASTDB --cores 16 -o path/to/output/directory --expr (for DE analysis)
```

The module has to be performed for each pair of fastq files. 

#### combine

```
vast-tools combine -o /path/to/output/directory (same that for align) -sp genome_reference --cores 16 --dbDir path/to/VASTDB -C (for DE analysis)
```

#### diff

```
vast-tools diff -a 5Y_siDDX5_17_B1_R1_cutadapt_match,5Y_siDDX5_17_B2_R1_cutadapt_match,5Y_siDDX5_17_B3_R1_cutadapt_match -b 5Y_siGL2_B1_R1_cutadapt_match,5Y_siGL2_B2_R1_cutadapt_match,5Y_siGL2_B3_R1_cutadapt_match --sampleNameA=siDDX5_17 --sampleNameB=siGL2 -i path/to/output/directory/INCLUSION_LEVELS.tab 
```

In paired-end context, the paired replicates are merged under the name of the first file (R1) in the INCLUSION_LEVELS.tab (it was not tested in with single-end datasets).

Check the name of the INCLUSION_LEVELS.tab (it depends on the genome reference and the number of replicates used). 

#### compare_expr (for Differential Expression analysis)

```
vast-tools compare_expr path/to/output/directory/cRPKM_and_count_tab -a 5Y_siDDX5_17_B1_R1_cutadapt_match,5Y_siDDX5_17_B2_R1_cutadapt_match,5Y_siDDX5_17_B3_R1_cutadapt_match -b 5Y_siGL2_B1_R1_cutadapt_match,5Y_siGL2_B2_R1_cutadapt_match,5Y_siGL2_B3_R1_cutadapt_match -name_A siDDX5_17 -name_B siGL2
```

## Comparison between Vast-tools and rMATS outputs

The **rMATS_filtering.R**  and **Vast_tools_filtering.R** scripts filter significant events and compare the results between the tools. The first one takes in input the path to raw rMATS results (JC) and output directory. The second requires the path to vast_out directory (Vast-tools Nextflow pipeline results) and output directory (use the same that the one used for **rMATS_filtering.R**) in addition to the reference genome (hg19,hg38...), the number of samples and whether the events without coordinates have to be conserved or not. These events cause issues for executing the **identify_common_events.py** script. These options can be modified directly at the beginning of the R scripts.  

The **identify_common_events.py** script identifies the corresponding events between rMATS differential alternative splicing events and Vast-tools differential events or events annotated in VASTDB from the previous script outputs. It takes in input the output path used for the **rMATS_filtering.R** and  **Vast-tools_filtering.R** scripts.

```
python3 identify_common_events.py ~/stage/results_hg19/
```

Finally, the **common_events.R** use the python script outputs to quantify the number of common events and events detected by rMATS that could be detected by Vast-tools (annotated in VASTDB). It takes in input the output path of the **rMATS_filtering.R** and  **Vast-tools_filtering.R** scripts (input can be modified at the beginning of the R script).

### Results

rMATS detects exon skipping (ES), mutually exclusive exons (MXE), retained intron (RI), alternative 5' and 3' splice sites (A5SS and A3SS) alternative splicing events. In other hand, Vast-tools detects and groups the ES, microexon skipping (MIC in event_type column) and MXE patterns into alternative exon skipping (EX) events. The RI, A5SS and A3SS patterns are also detected by the toolset.  
The tools are described in more detail in the M1-MEMOIRE-Marchand_Mehdi_2023-2024.pdf.

Four tables can be found in the result directory:  
  - rmats_summary.tab: number of differential alternative splicing (DAS) events and differentially alternatively spliced (DAS) genes according to the type of alternative splicing event (generated with **rMATS_filtering.R**).  
  - final_rmats_summary.tab: rmats_summary.tab without strict duplicates (events with strictly the same alternative exon coordinates) (generated with **rMATS_filtering.R**).  
  - vast_summary.tab: number of DAS events and DAS genes according to the type of alternative splicing event (generated with **Vast-tools_filtering.R**). 
  - common_summary.tab: number of differentially alternatively spliced (DAS) genes identified by both tools (generated with **common_events.R**).

These results are illustrated in the events_histo.pdf (see below).

![DAS_plot](https://github.com/mehdiretif/M1_internship/blob/main/image/example_events_histo.png)


The common_events.txt file (generated with **common_events.R**) corresponds to the differential ES events identified by both tools. In other hand, the raw_common_events.txt file corresponds to the rMATS differential ES events annotated in VASTDB (according to the human genome version). **WARNING: some duplicates are present in these files due to rMATS results redundancy**. These files were obtained through the **identify_common_events.py** script that identify the Vast-tools events (differential ES events or ES annotated in VASTDB) that correspond to the rMATS differential ES events with a margin of alternative splicing exon coordinates of 5 nucleotides (the flanking exon coordinates are not taken into account).  
The number of differential ES events detected in common without duplicates are given in the complementary_insights.txt with some other insights. 

