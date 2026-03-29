#!/usr/bin/env Rscript

#instalando e carregando pacotes necessários
packages <- c("dplyr", "stringr", "readr", "ggplot2", "scales")

for (i in packages){
	if (!requireNamespace(i, quietly = TRUE)){
		install.packages(i)
	}
	library(i, character.only = TRUE)
}

#-----------
#função principal 
#função tem como inputs os ambos outputs do kraken2, a partir do id taxonômico
#organiza os táxons aninhados em colunas em um novo dataframe
kraken_to_table_df <- function(classification_df,
                               report_df,
                               keep_unclassified = FALSE) {
  
  message("Preparando report...")
  
  #cria dataframe com as informações de report
  report <- report_df %>%
    mutate(
      name = str_trim(V6),
      taxid = as.character(V5),
      rank = V4
    )
  
  #define os níveis taxonômicos a partir do report
  taxonomy_levels <- c(
    D = "domain",
    P = "phylum",
    C = "class",
    O = "order",
    F = "family",
    G = "genus",
    S = "species"
  )
  
  #cria vetor vazio com o tamanho de taxonomy_levels
  current_lineage <- setNames(rep(NA, length(taxonomy_levels)),
                              taxonomy_levels)
  
  #cria lista vazia
  tax_list <- list()
  
  message("Reconstruindo linhagens...")
  
  #para cada linha de report pega o nível taxonomico, o nome do táxon e o id (este para 
  #conferir com o classification)
  for (i in seq_len(nrow(report))) {
    
    rank  <- report$rank[i]
    taxid <- report$taxid[i]
    name  <- report$name[i]
    
    #verifica se a informação em rank está em taxonomy_levels
    if (rank %in% names(taxonomy_levels)) {
      
      #pega o nome do táxon
      level <- taxonomy_levels[rank]
      current_lineage[level] <- name
      
      below <- which(names(current_lineage) == level)
      if (below < length(current_lineage)) {
        current_lineage[(below+1):length(current_lineage)] <- NA
      }
      
      tax_list[[taxid]] <- current_lineage
    }
  }
  
  message("Montando tabela taxonômica...")
  
  #montando tabela com os niveis taxonomicos e seus nomes
  taxonomy_df <- bind_rows(
    lapply(names(tax_list), function(x) {
      data.frame(
        taxid = x,
        as.data.frame(t(tax_list[[x]])),
        stringsAsFactors = FALSE
      )
    })
  )
  
  message("Preparando classification...")
  
  #repetindo processo para o classification
  classification <- classification_df %>%
    rename(
      status  = 1,
      read_id = 2,
      taxid   = 3
    ) %>%
    mutate(taxid = as.character(taxid)) %>%
    filter(status == "C")
  
  message("Juntando...")
  
  #fazendo um join entre classification e report já estruturado com os niveis
  #taxonomicos em colunas, por taxaid (informação comum)
  final <- classification %>%
    left_join(taxonomy_df, by = "taxid")
  
  message("Finalizado.")
  
  return(final)
}

#---------
#comando: gerar tabela
#a função entra com o classification e report, e gera um output, utilizando
#a função anterior (kraken_to_table_df)

cmd_table <- function(classification_file, report_file, output_file){
  
  message("Lendo arquivos...")
  
  #define classification
  classification_df <- read.delim(classification_file,
                                  header = FALSE,
                                  sep = "\t")
  #define report
  report_df <- read.delim(report_file,
                          header = FALSE,
                          sep = "\t")
  
  #executa função 
  result <- kraken_to_table_df(classification_df, report_df)
  
  #salva tabela
  write.table(result,
              output_file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  message("Tabela salva em: ", output_file)
}

#--------
#comando: gerar gráficos
#cria gráfico de barras com abundância por nível taxonômico definido como argumento

cmd_plot <- function(input_table, tax_level, prefix){

  #le tabela	
  df <- read.delim(input_table, sep="\t")
  
  #verifica se tax é válido
  if(!(tax_level %in% colnames(df))){
    stop(paste("Nivel taxonômico inválido:", tax_level))
  }
  
  #cria gráfico
  bar <- ggplot(df, aes(x = reorder(.data[[tax_level]], .data[[tax_level]], length))) +
    geom_bar() +
    coord_flip() +
    scale_y_continuous(labels = comma)
  
  #salva gráfico
  ggsave(paste0(prefix,"_", tax_level, ".pdf"), plot = bar)
  
  message("Gráfico salvo.")
}

#---------
#comando: abundancia por taxon

cmd_abundance <- function(input_table, tax_level, output_file){
  
  df <- read.delim(input_table, sep = "\t")
  
  #verifica se o nível existe na tabela
  if(!(tax_level %in% colnames(df))){
    stop(paste("Nível taxonômico inválido:", tax_level))
  }
  
  #pega numero total de reads
  total_reads <- nrow(df)
  
  #calcula abundancia
  abundance <- df %>%
    count(.data[[tax_level]], name = "reads") %>%
    mutate(
      relative_abundance = reads / total_reads,
      percent = relative_abundance * 100
    ) %>%
    arrange(desc(reads))
  
  #cria coluna com valores de abundancia (relativa)
  colnames(abundance)[1] <- tax_level
  
  #salva arquivo com as informações por táxon
  write.table(
    abundance,
    output_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  message("Abundância salva em: ", output_file)
}

#-------------
#comando: retira nulos

cmd_no_nulls <- function(input_table, tax_level, output_file){
  
  #le tabela
  df <- read.delim(input_table, sep = "\t")
  
  #verifica se o nível existe
  if(!(tax_level %in% colnames(df))){
    stop(paste("Nível taxonômico inválido:", tax_level))
  }
  
  #filtra nivel taxonomico definino no argumento e retira nulos
  df_filtered <- df %>%
    filter(!is.na(.data[[tax_level]]) & .data[[tax_level]] != "")
  
  #escreve tabela
  write.table(
    df_filtered,
    output_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  message("Tabela filtrada salva em: ", output_file)
}


#---------------
# CLI

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0){
  stop("
Uso:

kraken_tools table <classification> <report> <output.tsv>
kraken_tools plot <tabela.tsv> <prefix>
")
}

#organizando argumentos para cada comando
cmd <- args[1]

switch(cmd,
       
       table = cmd_table(args[2], args[3], args[4]),
       
       plot  = cmd_plot(args[2], args[3], args[4]),
       
       abundance = cmd_abundance(args[2], args[3], args[4]),
       
       no_nulls = cmd_no_nulls(args[2], args[3], args[4]),
       
       help = message(
         "
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
"),
       
       stop("Comando inválido.")
)
