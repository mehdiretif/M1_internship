#!/usr/bin/env nextflow

params.reads = "/Xnfs/abc/nf_scratch/mmarchand/vast_projet/trimmed_fastq/*_R{1,2}_cutadapt_match.fastq.gz"
params.db = "/Xnfs/abc/nf_scratch/mmarchand/vast_projet/VASTDB/"
params.groups = "/Xnfs/abc/nf_scratch/mmarchand/vast_projet/groups.txt" 
params.groupA = "siDDX5_17"
params.groupB = "siGL2"	
params.gen_ref = "hg19"
params.cores = "16"
params.paired = "TRUE"
params.output = "/Xnfs/abc/nf_scratch/mmarchand/vast_projet/nf_results/"

log.info "VAST-TOOLS NF"
log.info "============="
log.info "reads (FASTQ file)		: ${params.reads}"
log.info "groups (group  A : group B)	: ${params.groupA} : ${params.groupB}"
log.info "reference genome		: ${params.gen_ref}"
log.info "cores				: ${params.cores}"
log.info "paired (TRUE/FALSE)		: ${params.paired}"
log.info "output (path)			: ${params.output}"
log.info "\n"

if (params.paired == "TRUE") {
        diff_variable = "_R1_cutadapt_match"
        }
else {
        diff_variable = "_cutadapt_match"
        }

process align {
	tag "$pair_id"
	label "multi_thread"	
	container "marchandrm/vast-tools:2.5.1"
	publishDir params.output, mode: 'copy'
	
	input:
	tuple val(pair_id), path(read_path)
			
	output:
	path "vast_out/to_combine/${pair_id}${diff_variable}.eej2", emit: eej2
	path "vast_out/to_combine/${pair_id}${diff_variable}.exskX", emit: exskX
	path "vast_out/to_combine/${pair_id}${diff_variable}.info", emit: info
	path "vast_out/to_combine/${pair_id}${diff_variable}.IR2", emit: IR2
	path "vast_out/to_combine/${pair_id}${diff_variable}.IR.summary_v2.txt", emit: IR_summary
	path "vast_out/to_combine/${pair_id}${diff_variable}.micX", emit: micX
	path "vast_out/to_combine/${pair_id}${diff_variable}.MULTI3X", emit: MULTI3X
		
	script:	
	"""
	vast-tools align ${read_path[0]} ${read_path[1]} -sp ${params.gen_ref} --dbDir ${params.db} --cores ${params.cores} -o vast_out
	"""
}	

process combine_process {
	container "marchandrm/vast-tools:2.5.1"
	label "multi_thread"
	publishDir params.output, mode: 'copy'

	input:
	path eej2
	path exskX
	path info
	path IR2
	path IR_summary
	path micX
	path MULTI3X	

	output:
	path "vast_out/INCLUSION*.tab", emit: inclusion_tab  
	
	script:
	"""
	mkdir -p vast_out/to_combine
	mv -t vast_out/to_combine/. $eej2 $IR2 $IR_summary $micX $MULTI3X $exskX $info
	vast-tools combine -o ./vast_out -sp ${params.gen_ref} --cores ${params.cores} --dbDir ${params.db}
	"""
}

process diff {
	container "marchandrm/vast-tools:2.5.1"
	label "multi_thread"	
	publishDir params.output, mode: 'copy'

	input: 
	path inclusion_tab
	path eej2
	path exskX
	path info
	path IR2
	path IR_summary
	path micX
	path MULTI3X

	output:
	path "vast_out/*diff_output.tab", emit: diff_table
	path "vast_out/*diff_output.pdf"
	
	script:
	"""
	mkdir -p vast_out/to_combine
	mv -t vast_out/to_combine/. $eej2 $IR2 $IR_summary $micX $MULTI3X $exskX $info
	mv $inclusion_tab vast_out/
	a=\$(grep ${params.groupA} ${params.groups} | cut -f 1 |  sed "s/\$/${diff_variable}/" | paste -d, -s)
	b=\$(grep ${params.groupB} ${params.groups} | cut -f 1 |  sed "s/\$/${diff_variable}/" | paste -d, -s)
	
	vast-tools diff -a "\$a" -b "\$b" --sampleNameA=${params.groupA} --sampleNameB=${params.groupB} -i $inclusion_tab -c ${params.cores} -d diff_output
	"""
}

process new_table {
	label "multi_thread"
	publishDir params.output, mode: 'copy'

	input:
	path diff_table
	path inclusion_tab

	output:
	path "final*DIFF.txt"

	script:
	"""
	#Boucle pour obtenir le numero des colonnes d'interet 
	value=8 #initialisation
	list=("\$value")

	#sample_nb=\$(wc -l ${params.groups})

	i=1
	while [ "\$i" -lt 6 ]; do
		value=\$((value + 2))
		list+=("\$value")
		i=\$((i + 1))
	done

	###CREATION NOUVEAU FICHIER
	#Creation de la premiere ligne/colonne

	first_line=\$(head -n1 $diff_table)
	first_line+='\t'coordinate'\t'full_coordinates'\t'event_type

	#Ajout nom des colonnes
	for sample in "\${list[@]}"
		do
			col_sample=\$(head -n1 $inclusion_tab | cut -d '\t' -f"\$sample")
			first_line+='\t'"\$col_sample"
		done

	#Creation contenu du fichier
	new_file=\$(awk -v cols="\${list[*]}" '
		BEGIN {
			split(cols, colArray, " ") 
		}

		FNR == NR {
			if (NR > 1) {
				inclLevels[\$2] = \$0
			}
			next
		}
		{
			if (FNR == 1 && NR > 1) { next }
			id = \$2

			split(inclLevels[id], inclCols, "\t")
			coord = inclCols[3]
			full_coord = inclCols[5]
			event_type = inclCols[6]
			line = \$0 "\t" coord "\t" full_coord "\t" event_type
			for (i in colArray) {
				line = line "\t" inclCols[colArray[i]]
			}
			print line
		}
	' $inclusion_tab $diff_table)

	echo -e "\$first_line\n\$new_file" > final_INCLUSION_FILE.DIFF.txt

	"""
}

workflow {
	read_pairs = channel.fromFilePairs(params.reads, checkIfExists: true)
	align(read_pairs)
	ch_input_combine_eej2 = align.out.eej2
	ch_input_combine_exskX = align.out.exskX
	ch_input_combine_info = align.out.info
	ch_input_combine_IR2 = align.out.IR2
	ch_input_combine_IR_summary = align.out.IR_summary
	ch_input_combine_micX = align.out.micX
	ch_input_combine_MULTI3X = align.out.MULTI3X

	combine_process(ch_input_combine_eej2.collect(),ch_input_combine_exskX.collect(),ch_input_combine_info.collect(),ch_input_combine_IR2.collect(),ch_input_combine_IR_summary.collect(),ch_input_combine_micX.collect(),ch_input_combine_MULTI3X.collect())
	
	diff(combine_process.out,ch_input_combine_eej2.collect(),ch_input_combine_exskX.collect(),ch_input_combine_info.collect(),ch_input_combine_IR2.collect(),ch_input_combine_IR_summary.collect(),ch_input_combine_micX.collect(),ch_input_combine_MULTI3X.collect())

	new_table(diff.out.diff_table,combine_process.out)
}	

