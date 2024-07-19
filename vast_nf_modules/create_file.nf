/*
 * SPDX-FileCopyrightText: 2024 Marchand Mehdi <mehdi.retif@gmail.com>
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 */

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

	sample_nb=\$(wc -l < ${params.groups})

	i=1
	while [ "\$i" -lt "\$sample_nb" ]; do
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