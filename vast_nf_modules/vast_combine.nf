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
	path resume
	path tmp
	path bias
	path cRPKM
	
	
	output:
	path "vast_out/INCLUSION*.tab", emit: inclusion_tab  
	path "vast_out/cRPKM_AND_COUNTS*.tab", emit: cRPKM_and_count_tab
	path "vast_out/COUNTS-*.tab", emit: count_tab
	path "vast_out/cRPKM-*.tab", emit: cRPKM_tab
	

	script:
	"""
	mkdir -p vast_out/to_combine
	mv -t vast_out/to_combine/. $eej2 $IR2 $IR_summary $micX $MULTI3X $exskX $info
	
	mkdir vast_out/expr_out
	mv -t vast_out/expr_out/. $bias $cRPKM

	mkdir vast_out/tmp
	mv -t	vast_out/tmp/. $tmp $resume

	vast-tools combine -o ./vast_out -sp ${params.gen_ref} --cores ${params.cores} --dbDir ${params.db} -C
	"""
}