                               -----outkraken-----

Programa outkraken.R utiliza os outputs do sofrware kraken2 (classification e report) e cria uma tabela organizando todos os níveis taxonômicos de cada read classificado. A partir dessa tabela, outkraken.R pode gerar gráficos de abundância da diversidade da amostra, um resumo com números relativos e absolutos da diversidade.

./outkraken.R help

--------------------------------------------------------------------------------------------
Kraken_outputs uso:

table (classification_kraken, report_kraken, output)
  --> gera tabela com táxons para cada read

plot (output de table, tax, output_prefixo)
  --> gera um plot com o nivel taxonomico informado

abundance (output de table, tax, prefixo output)
  --> gera resumos de abundância por táxon

no_nulls (out_put de table, tax, output_filtrado)
  --> retira nulos de um nível taxonomico

Niveis taxonômicos: phylum, class, order, family, genus, species

--------------------------------------------------------------------------------------------

O programa é útil para organizar por read classificado os níveis taxonômicos encontrados pelo kraken2, além do controle de nulos em análises de diversidade taxonômica a partir de dados NGS. Alternativamente recomendamos o uso do software krakentools para tratamento dos outputs do kraken2.


                             -----rodar_outkraken-----

Um bashscript que roda o outkraken para todas as duplas de report e classification em uma pasta. Os arquivos de input precisam ter terminação *_classification.tsv e *_report.tsv. Além disso, o programa pede o nível taxonômico desejado para o no_nulls, plot e abundance.

uso: 

./rodar_outkraken.sh <nível_taxonômico> (kingdom, phylum, class, order, family, genus, species)

No repositório há amostras de teste de outputa do kraken2, além de exemplos de outputs do outkraken.

em: /out_outkraken 
