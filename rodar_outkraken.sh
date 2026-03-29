#!/bin/bash

set -euo pipefail

# checar argumento
if [[ $# -lt 1 ]]; then
    echo "Uso: $0 <nivel_taxonomico>"
    echo "Exemplo: $0 family"
    exit 1
fi

tax_level=$1

echo "nível taxonômico: $tax_level"

#rodando para loop para cada par de outputs de kraken2
for class_file in *_classification.tsv; do

	echo "processando $class_file"

	#define nome base 
	basename=${class_file%_classification.tsv}

	#procura report correspondente
	report_file="${basename}_report.tsv"

	#checa se report existe
	if [[ ! -f "$report_file" ]]; then
		echo "Arquivo de report nao encontrado"
		continue
	fi

	echo "processando $class_file e $report_file"

	out_table="${basename}_table.tsv"
	no_nulls_table="${basename}_nonulls_${tax_level}.tsv"
	plot_name="${basename}_plot"
	abundance_name="${basename}_abudance_${tax_level}.tsv"

	#criando tabela
	./outkraken.R table "$class_file" "$report_file" "$out_table"

	echo "tabela em $out_table"

	#rodando no_nulls
	./outkraken.R no_nulls "$out_table" "$tax_level" "$no_nulls_table"

	echo "tabela sem nulos em $no_nulls_table para o nivel $tax_level"

        #rodando plot e abudance
	./outkraken.R plot "$no_nulls_table" "$tax_level" "$plot_name"

	echo "grafico salvo em $plot_name para o nivel $tax_level"

	./outkraken.R abundance "$no_nulls_table" "$tax_level" "$abundance_name"

	echo "abudancia salva em $abundance_name para o nivel $tax_level"

done	


