# Vast-tools Nextflow pipeline

## Options

--reads '/path/to/reads' \
--db '/path/to/VASTDB' \
--groups '/path/to/file/containing/insights/about/samples' (refer the groups.txt file) \
--groupA 'conditionA' \
--groupB 'conditionB' \
--gen_ref 'genome_reference' \
--cores 'number_of_cores' \
--paired 'yes'(paired-sequencing) or 'no' \
--output '/path/to/output/directory' \
 
The databases are available on the [vast-tools github](https://github.com/vastgroup/vast-tools)

## Example

nextflow run /Xnfs/abc/nf_scratch/mmarchand/vast_projet/vast_tools.nf /\\
-profile psmn /\\
-c /Xnfs/abc/nf_scratch/mmarchand/vast_projet/nextflow.config /\\
--reads '/Xnfs/abc/nf_scratch/mmarchand/vast_projet/trimmed_fastq_test/*_R{1,2}_cutadapt_match.fastq.gz' /\\
--db '/Xnfs/abc/nf_scratch/mmarchand/vast_projet/VASTDB/' /\\
--groups '/Xnfs/abc/nf_scratch/mmarchand/vast_projet/groups.txt' /\\
--groupA 'siDDX5_17' /\\
--groupB 'siGL2' /\\
--gen_ref 'hg19' /\\
--cores '16' /\\
--paired 'yes' /\\
--output '/Xnfs/abc/nf_scratch/mmarchand/vast_projet/vast_out/'

