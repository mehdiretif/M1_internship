## Vast-tools Nextflow pipeline

### Options

```
    --reads     '/path/to/read/files' 
    --db        '/path/to/VASTDB' 
    --groups    '/path/to/file/containing/insights/about/samples' (refer to the groups.txt file) 
    --groupA    'conditionA' 
    --groupB    'conditionB' 
    --gen_ref   'genome_reference' 
    --cores     'number_of_cores' 
    --paired    'yes'(paired-sequencing) or 'no' 
    --output    '/path/to/output/directory'
``` 

The databases are available on the [vast-tools github](https://github.com/vastgroup/vast-tools).

### Example

```
nextflow run /Xnfs/abc/nf_scratch/mmarchand/vast_projet/vast_tools.nf 
-profile psmn \
-c /Xnfs/abc/nf_scratch/mmarchand/vast_projet/nextflow.config \
--reads '/Xnfs/abc/nf_scratch/mmarchand/vast_projet/trimmed_fastq_test/*_R{1,2}_cutadapt_match.fastq.gz' \
--db '/Xnfs/abc/nf_scratch/mmarchand/vast_projet/VASTDB/' \
--groups '/Xnfs/abc/nf_scratch/mmarchand/vast_projet/groups.txt' \
--groupA 'siDDX5_17' \
--groupB 'siGL2' \
--gen_ref 'hg19' \
--cores '16' \
--paired 'yes' \
--output '/Xnfs/abc/nf_scratch/mmarchand/vast_projet/vast_out/' 
```

### Comparison between Vast-tools and rMATS outputs

The process_rMATS_Vast_outputs.R filters significant events and compare the results between the tools.

The paths can be modified in the script. The difference between the final_rmats_summary.tab corresponds to the rmats_summary.tab without strict duplicates (same alternative exon coordinates) for SE and MXE events. 

The identify_common_events.py script identifies the corresponding events between rMATS differential alternative splicing events and Vast-tools differential events or events annotated in VASTDB from the previous script outputs.

Finally, the process_rMATS_Vast_outputs.R use the python script outputs to quantify the number of common events and events detected by rMATS that could be detected by Vast-tools (annotated in VASTDB). The results are added to the complementary_insights.txt file. 