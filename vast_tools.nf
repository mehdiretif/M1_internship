#!/usr/bin/env nextflow

log.info "VAST-TOOLS NF"
log.info "============="
log.info "reads (FASTQ file)		: ${params.reads}"
log.info "groups (group  A : group B)	: ${params.groupA} : ${params.groupB}"
log.info "reference genome		: ${params.gen_ref}"
log.info "cores				: ${params.cores}"
log.info "paired (yes/no)		: ${params.paired}"
log.info "output (path)			: ${params.output}"
log.info "\n"

if (params.paired == "yes") {
        params.diff_variable = "_R1_cutadapt_match"
        }
else {
        params.diff_variable = "_cutadapt_match"
        }


include { align } from './vast_nf_modules/vast_align'
include { combine_process } from './vast_nf_modules/vast_combine'
include { diff } from './vast_nf_modules/vast_diff'
include { new_table } from './vast_nf_modules/create_file'

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

