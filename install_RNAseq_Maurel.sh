#!/usr/bin/env bash

# Contrôle de l'exécution : affichage des commandes exécutées (-x), arrêt si erreur (-e) ou usage de variables non déclarées (-u)
set -eux 

########################################################################################################################################

############################################### ARBORESCENCE DU PROJET #################################################################

########################################################################################################################################

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

echo -e "\n Création des répertoires \n"

# Répertoire contenant les fichiers .fastq.gz (reads)
mkdir -p sequence_data   # option -p permet de créer le dossier si et seulement s'il n'existe pas
sequence_data_dir=sequence_data

# Répertoire contenant les sorties du trimming des séquences par trimmomatic
mkdir -p trimmomatic_results

# Répertoire contenant les fichiers de table de comptage
mkdir -p count_table

# Répertoire contenant les fichiers utilisés pour l'indexation et l'alignement sur le génome de référence
mkdir -p Mapping
mapping_dir=Mapping

# Répertoire contenant le génome de référence 
mkdir -p "${mapping_dir}"/reference_genome
reference_genome_dir="${mapping_dir}"/reference_genome

# Répertoire contenant l'annotation du génome de référence
mkdir -p "${mapping_dir}"/reference_genome_annotation
reference_genome_annotation_dir="${mapping_dir}"/reference_genome_annotation

# Répertoire contenant les sorties d'indexation du génome de référence
mkdir -p "${mapping_dir}"/reference_genome_index
reference_genome_index_dir="${mapping_dir}"/reference_genome_index

# Répertoire contenant les sorties d'alignement sur le génome de référence
mkdir -p "${mapping_dir}"/STAR_alignment_output

echo -e "\n Fin de la création des répertoires \n"

########################################################################################################################################

############################################### INSTALLATION DES OUTILS ################################################################

########################################################################################################################################

# Création d'un environnement "rnaseq" dédié 
# contenant les outils nécessaires pour l'analyse
eval "$(conda shell.bash hook)" #Activation des fonctions de conda dans le shell


echo -e "\n Création de l'environnement 'rnaseq' qui
contient les outils fastqc, trimmomatic, samtools, star, subread et perl \n"

conda create -y -q -n rnaseq -c conda-forge -c bioconda fastqc trimmomatic star samtools subread perl openjdk

echo -e "\n Fin de création de l'environnement rnaseq contenant les outils \n"



########################################################################################################################################

############################################### TELECHARGEMENTS ET EXTRACTION DES ARCHIVES #############################################

########################################################################################################################################

# Téléchargement de l'archive contenant les reads
echo -e "\n Téléchargement de l'archive contenant les reads sous le nom sequence_data.zip \n"

sequence_archive_name="sequence_data.zip" #nom de l'archive
sequence_data_link="https://zenodo.org/api/records/7525740/files-archive" #lien à partir duquel télécharger les archives
wget "${sequence_data_link}" -O "${sequence_archive_name}"

# Décompression de l'archive 
echo -e "\n Décompression de l'archive 'sequence_data.zip' contenant les reads \n
dans le dossier sequence_data \n"

unzip "${sequence_archive_name}" -d "${sequence_data_dir}"

# Téléchargement du génome de référence
echo -e "\n Téléchargement du génome de référence dans un dossier Mapping/reference_genome\n"

ref_genome_archive_name="chr18.fa.gz" #nom de l'archive
ref_genome_link="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz" #lien de téléchargement du génome de référence(chr18)
wget "${ref_genome_link}" -O "${reference_genome_dir}/${ref_genome_archive_name}"

# Extraction de l'archive contenant le génome de référence au format fasta
echo -e "\n Décompression de l'archive 'chr18.fa.gz' contenant le fichier fasta du chr18\n"

reference_genome_fasta_file_name="chr18.fa" # nom du fichier fasta après décompression
gunzip -c "${reference_genome_dir}/${ref_genome_archive_name}" > "${reference_genome_dir}/${reference_genome_fasta_file_name}"

# Téléchargement de l'annotation du chr 18
echo -e "\n Téléchargement et de l'annotation du génome de référence (chr 18) dans un dossier Mapping/reference_genome_annotation\n"

genome_annotation_archive_name="gencode.v24lift37.basic.annotation.gtf.gz" #nom de l'archive
genome_annotation_archive_link="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz"
wget "${genome_annotation_archive_link}" -O "${reference_genome_annotation_dir}/${genome_annotation_archive_name}"

# Extraction de l'archive
genome_annotation_gtf_file_name="gencode.v24lift37.basic.annotation.gtf" # nom du fichier gtf après décompression
gunzip -c "${reference_genome_annotation_dir}/${genome_annotation_archive_name}" \
> "${reference_genome_annotation_dir}/${genome_annotation_gtf_file_name}"


########################################################################################################################################

################################################### INDEXATION DU GENOME DE REFERENCE ##################################################

########################################################################################################################################

#Activation de l'environnement contenant l'outil dédié (STRAR)
conda activate rnaseq
echo -e "\n Indexation du génome de référence \n"

#Création de l'index
STAR --runMode genomeGenerate --runThreadN 4 \
  --genomeDir ${reference_genome_index_dir} \
  --genomeFastaFiles ${reference_genome_dir}/${reference_genome_fasta_file_name} \
  --sjdbGTFfile ${reference_genome_annotation_dir}/${genome_annotation_gtf_file_name}
conda deactivate 
