process expr {
    container "marchandrm/vast-tools:2.5.1"
	label "multi_thread"	
	publishDir params.output, mode: 'copy'

	input: 
	path cRPKM_and_count_tab
    path count_tab
    path cRPKM_tab
    path bias
    path cRPKM

    output:
	path "vast_out/DiffGE*.tab"

    script :
    """
    mkdir -p vast_out/expr_out
	mv -t vast_out/expr_out/. $bias $cRPKM
	mv -t vast_out/. $cRPKM_and_count_tab $count_tab $cRPKM_tab

    a=\$(grep ${params.groupA} ${params.groups} | cut -f 1 |  sed "s/\$/${params.diff_variable}/" | paste -d, -s)
	b=\$(grep ${params.groupB} ${params.groups} | cut -f 1 |  sed "s/\$/${params.diff_variable}/" | paste -d, -s)
	
    vast-tools compare_expr vast_out/$cRPKM_and_count_tab -a "\$a" -b "\$b" -name_A ${params.groupA} -name_B ${params.groupB}
    """

}

