# TP RNA-seq — EMT - MAUREL DYLANE TCHUISSE II

**Cours** : Atelier NGS — Daniel Gautheret (I2BC), 15/09/2025  
**Jeu de données** : Données RNAseq d'une EMT induite par expression ectopique du gène Zeb1 dans des cellules de cancer du poumon non à petites cellules (lignée H358) - Yang *et al.*, 2016  
**Objectif** : Identifier les gènes **différentiellement exprimés** entre Day0 (cellules non induites) et Day7 (cellules après 7 jours d'induction).

> Les données ont été réduites réduites par échantillonange de ~0,5% des reads, du chromosome chr18 uniquement 

---

## Protocole d'exécution du présent pipeline

```bash
#Désarchiver l'archive
tar -xvf Maurel_Tchuisse_rna_seq.tar
#Se placer dans le dossier désarchivé
cd Maurel_Tchuisse_rna_seq

## installation des outils + téléchargement des données (génome de référence, annotation du génome de référence, données RNAseq) + indexation du génome de référence
bash install_RNAseq_Maurel.sh 

# Exécution du pipeline : Trimming → mapping → indexation BAMfiles → Table de comptage
bash run_RNAseq_Maurel.sh
```

---

## Arborescence du projet

```
sequence_data/
trimmomatic_results/
count_table/
Mapping/
  reference_genome/
  reference_genome_annotation/
  reference_genome_index/
  STAR_alignment_output/
```

### Description des dossiers

- **sequence_data/** — contient les fichiers FASTQ bruts téléchargés.  
- **trimmomatic_results/** — contient les FASTQ nettoyés après trimming.  
- **count_table/** — contient les tables finales de comptage par gène.  
- **Mapping/reference_genome/** — contient le fichier FASTA du génome de référence.  
- **Mapping/reference_genome_annotation/** — contient l’annotation du génome de référence.  
- **Mapping/reference_genome_index/** — contient l’index généré pour le génome de référence.  
- **Mapping/STAR_alignment_output/** — contient les fichiers d'alignement des reads sur le génome de référence et leurs index.

---

## Environnement créé automatiquement

---

## 🔧 Étapes du pipeline

### `install_RNAseq_Maurel.sh`
- Crée l’environnement conda
    **conda env : `rnaseq`**
    Outils installés automatiquement :FastQC, Trimmomatic, STAR, Samtools,   Subread/featureCounts, Perl
- Télécharge les FASTQ (Zenodo)
- Télécharge la séquence FASTA du chr18 + annotation (fichier GTF)
- Décompresse le FASTA & GTF du génome de référence
- Construit l’index STAR du génome de référence

### `run_RNAseq_Maurel.sh`
- Nettoyage des reads (Trimmomatic)
- Alignement STAR (BAM trié)
- Indexation BAM
- Comptage par gène (featureCounts)
- Conversion gene_id → gene_name
- Filtration chr18
