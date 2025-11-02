#!/usr/bin/env bash

# Contrôle de l'exécution : affichage des commandes exécutées, arrêt si erreur ou usage de variables non déclarées
set -eux 

#Activation de l'environnement "rnaseq" contenant les outils pour l'analyse
eval "$(conda shell.bash hook)"
conda activate rnaseq



########################################################################################################################################

################################################# DECLARATION DES VARIABLES CONTENANT LES REPERTOIRES ##################################

########################################################################################################################################

# '''''''''''''''''''''' '''''''''''''''''''''''''''''''''''''Rappel de l'arborescence ''''''''''''''''''''''''''''''''''''''''''''''''

# ───────────────────────────────
# sequence_data/
# │
# │
# Mapping/
# │
# ├── reference_genome/
# │
# ├── reference_genome_annotation/
# │
# ├── STAR_alignment_output/
# │
# │── reference_genome_index/
# │
# count_table/
# │
# │
# trimmomatic_results/
# ───────────────────────────────


echo -e "\n Définition des répertoires\n"
# Répertoire contenant les fichiers .fastq.gz de reads
sequence_data_dir=sequence_data

# Répertoire contenant les sorties du trimming des séquences par trimmomatic
trimming_results_dir=trimmomatic_results

# Répertoire contenant les fichiers de table de comptage
count_table_dir=count_table

# Répertoire contenant les fichiers utilisés pour l'indexation et l'alignement sur le génome de référence
mapping_dir=Mapping

# Répertoire contenant le génome de référence 
reference_genome_dir="${mapping_dir}/reference_genome"

# Répertoire contenant l'annotation du génome de référence
reference_genome_annotation_dir="${mapping_dir}/reference_genome_annotation"
reference_genome_annotation_file="gencode.v24lift37.basic.annotation.gtf"

# Répertoire contenant les sorties d'indexation du génome de référence
reference_genome_index_dir="${mapping_dir}/reference_genome_index"

# Répertoire contenant les sorties d'alignement sur le génome de référence
STAR_alignment_output_dir="${mapping_dir}/STAR_alignment_output"

echo -e "\n Fin Définition des répertoires\n"


########################################################################################################################################

############################################################# TRIMMING DES READS #######################################################

########################################################################################################################################

#Trimming des séquences avec trimmomatic
echo -e "\n Trimming : retrait/troncature des séquences de mauvaise qualité \n"

for R1 in "${sequence_data_dir}"/*R1.fastq.gz; do
  R2="${R1/.R1/.R2}"    # fabrique le nom du fichier fastq R2 en remplaçant R1 par R2
  sample=$(basename "$R1")  
  sample="${sample%.R1.fastq.gz}" #récupère le préfixe avant .R1
  echo -e "Trimming de $sample\n  
    Fichier R1 : $(basename "$R1")\n  
    Fichier R2 : $(basename "$R2")\n"

  trimmomatic PE "$R1" "$R2" \
    -baseout "${trimming_results_dir}/${sample}.fastq" \
    LEADING:25 TRAILING:25 MINLEN:50
done

echo -e "\n Fin du Trimming \n"


########################################################################################################################################

######################################## MAPPING DES READS & INDEXATION DES FICHIERS D'ALIGNEMENT ######################################

########################################################################################################################################

echo -e "\n Début du Mapping sur le génome de référence (chr18)"

for R1 in "${trimming_results_dir}"/*1P.fastq; do
  R2="${R1/1P/2P}"    # fabrique le chemin relatif du fichier 2P.fastq en remplaçant 1P par 2P
  sample=$(basename "$R1")  
  sample=${sample%_1P.fastq} # récupère le préfixe avant _1P
  echo -e "Mapping de $sample\n  
    Fichier 1P : $(basename "$R1")\n  
    Fichier 2P : $(basename "$R2")\n"

  STAR --runThreadN 4  --outFilterMultimapNmax 1\
  --genomeDir ${reference_genome_index_dir}\
  --outSAMattributes All --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "${STAR_alignment_output_dir}"/"${sample}" \
  --readFilesIn "${R1}" "${R2}"

  # Indexation avec des fichiers .bam avec samtools
  bam_file="${STAR_alignment_output_dir}/${sample}Aligned.sortedByCoord.out.bam" #fichier d'alignement bam pour la paire de reads R1 R2 en cours
  
  #Contrôle de l'existence du fichier bam
  if [[ -f "${bam_file}" ]]; then
    samtools index "${bam_file}"
  else
    echo "BAM : ${bam_file} manquant"
  fi
done

echo -e "\n Fin du Mapping sur le génome de référence (chr18)"



########################################################################################################################################

#################################################### GENERATION DES TABLES DE COMPTAGES ################################################

########################################################################################################################################

#Table de comptage
echo -e "\n Détermination des comptages\n"

#Génération d'une table "counts_table_gene_id.txt" contenant pour chaque échantillon, le décompte par gène, du nombre de reads 
#s'étant alignés sur les exons dans le génome de référence
counts_table_gene_id_file="table_gene_id_counts.txt"

featureCounts -p -t exon -g gene_id -a "${reference_genome_annotation_dir}"/"${reference_genome_annotation_file}" \
   -o ${count_table_dir}/${counts_table_gene_id_file} ${STAR_alignment_output_dir}/*.bam

#Génération d'un fichier d'équivalence gene_id gene_name(HUGO)
gene_id_gene_name_equivalence_file="gene_id_gene_name_equivalence.txt" #Nom du fichier des équivalences

perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' \
"${reference_genome_annotation_dir}"/"${reference_genome_annotation_file}" | sort \
| uniq > ${count_table_dir}/${gene_id_gene_name_equivalence_file}

#Génération de la table de comptage finale dans un fichier "table_gene_name_counts.txt"

count_table_gene_name_file="table_gene_name_counts.txt" # Nom du fichier contenant la table de comptage 
#avec pour colonnes: gene_name | gene_length | gene_counts_sample1 | gene_counts_sample2 | ....

#tri de la table de comptage au format:  gene_id | counts
sort ${count_table_dir}/${counts_table_gene_id_file} > ${count_table_dir}/${counts_table_gene_id_file/.txt/_sorted.txt}

#tri de la table d'équivalence gene_id | gene_name
sort ${count_table_dir}/${gene_id_gene_name_equivalence_file} > ${count_table_dir}/${gene_id_gene_name_equivalence_file/.txt/_sorted.txt}

#Jointure des deux tables triées précédentes dans un fichier temporaire temp3
join ${count_table_dir}/${counts_table_gene_id_file/.txt/_sorted.txt} \
${count_table_dir}/${gene_id_gene_name_equivalence_file/.txt/_sorted.txt} | grep "chr18" > ${count_table_dir}/temp3


#Table de comptage finale
awk '{print $NF " " $6 " " $7 " " $8 " " $9 " " $10  " " $11 " " $12 }' ${count_table_dir}/temp3 > ${count_table_dir}/${count_table_gene_name_file}
 
conda deactivate
echo -e "\nFin de la génération de la table de comptage\n"

#Vérification des comptages pour le gène CDH2
echo -e "\nVérification des comptages pour le gène CDH2\n"
grep "CDH2" ${count_table_dir}/${count_table_gene_name_file}