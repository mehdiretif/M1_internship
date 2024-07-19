/*
 * SPDX-FileCopyrightText: 2024 Marchand Mehdi <mehdi.retif@gmail.com>
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 */

process diff {
	container "marchandrm/vast-tools_limma:2.5.1"
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
	a=\$(grep ${params.groupA} ${params.groups} | cut -f 1 |  sed "s/\$/${params.diff_variable}/" | paste -d, -s)
	b=\$(grep ${params.groupB} ${params.groups} | cut -f 1 |  sed "s/\$/${params.diff_variable}/" | paste -d, -s)
	
	vast-tools diff -a "\$a" -b "\$b" --sampleNameA=${params.groupA} --sampleNameB=${params.groupB} -i $inclusion_tab -c ${params.cores} -d diff_output
	"""
}