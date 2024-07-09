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
	vast-tools combine -o ./vast_out -sp ${params.gen_ref} --cores ${params.cores} --dbDir ${params.db} -C
	"""
}