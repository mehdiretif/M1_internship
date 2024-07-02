process align {
	tag "$pair_id"
	label "multi_thread"	
	container "marchandrm/vast-tools:2.5.1"
	publishDir params.output, mode: 'copy'
	
	input:
	tuple val(pair_id), path(read_path)
			
	output:
	path "vast_out/to_combine/${pair_id}${params.diff_variable}.eej2", emit: eej2
	path "vast_out/to_combine/${pair_id}${params.diff_variable}.exskX", emit: exskX
	path "vast_out/to_combine/${pair_id}${params.diff_variable}.info", emit: info
	path "vast_out/to_combine/${pair_id}${params.diff_variable}.IR2", emit: IR2
	path "vast_out/to_combine/${pair_id}${params.diff_variable}.IR.summary_v2.txt", emit: IR_summary
	path "vast_out/to_combine/${pair_id}${params.diff_variable}.micX", emit: micX
	path "vast_out/to_combine/${pair_id}${params.diff_variable}.MULTI3X", emit: MULTI3X
		
	script:	
	"""
	vast-tools align ${read_path[0]} ${read_path[1]} -sp ${params.gen_ref} --dbDir ${params.db} --cores ${params.cores} -o vast_out
	"""
}	