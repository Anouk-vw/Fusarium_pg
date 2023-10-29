mkdir -p $1

ids=$(cat /home/anouk/anouk2/Cuban/BUSCO_tree/complete_busco_in_all.txt)

for dir in $(find Verticilliodes -type d -name "single_copy_busco_sequences"); do
	sppname=$(echo $dir | sed 's/\/run_hypocreales_odb10\/busco_sequences\/single_copy_busco_sequences//' | sed 's/\.\///');
	for id in ${ids}; do
		#echo ${id}
		file=$(echo ${dir}/${id}.faa)
		#f=$(echo $file | cut -d'/' -f 7)
		f=${id}
		echo ${file} cp sed cut mv to ${1}/${sppname}_${f}.1 ${1}/${sppname}_${f};
		cp $file ${1}/${sppname}_${f}
		sed -i 's/^>/>'${sppname}'|/g' ${1}/${sppname}_${f}
		cut -f 1 -d ":" ${1}/${sppname}_${f} | tr '[:lower:]' '[:upper:]' > ${1}/${sppname}_${f}.1
		mv ${1}/${sppname}_${f}.1 ${1}/${sppname}_$f;
		done
	done
