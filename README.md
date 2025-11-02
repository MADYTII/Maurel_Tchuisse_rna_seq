# TP RNA-seq — EMT - MAUREL DYLANE TCHUISSE II

**Cours** : Atelier NGS — Daniel Gautheret (I2BC), 15/09/2025  
**Jeu de données** : Yang *et al.*, 2016 — EMT (H358, induction ZEB1)  
**Objectif** : Identifier les gènes **différentiellement exprimés** entre Day0 et Day7.

> Données réduites (~0,5% des reads, uniquement chr18) → temps de calcul très réduit.  
> Sur données complètes : ressources ×200.

---

## Protocole d'exécution du présent pipeline

```bash
tar -xvf Maurel_Tchuisse_rna_seq.tar
cd Maurel_Tchuisse_rna_seq

# Si conda n'est pas activé automatiquement
# source ~/miniconda3/etc/profile.d/conda.sh

bash install_RNAseq_Maurel.sh     # installation + données + références + indexation
bash run_RNAseq_Maurel.sh         # trimming → mapping → index BAM → comptage
```

Sorties principales :
- BAM + index : `Mapping/STAR_alignment_output/`
- Table brute de comptages : `count_table/table_gene_id_counts.txt`
- Table finale annotée chr18 : `count_table/table_gene_name_counts.txt`

---

## 📁 Arborescence du projet

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
- **Mapping/reference_genome/** — contient le fichier FASTA du chromosome 18.  
- **Mapping/reference_genome_annotation/** — contient l’annotation Gencode.  
- **Mapping/reference_genome_index/** — contient l’index STAR généré.  
- **Mapping/STAR_alignment_output/** — contient les alignements BAM et leurs index.

---

## 🧰 Environnement créé automatiquement

**conda env : `rnaseq`**

Outils installés automatiquement :
FastQC, Trimmomatic, STAR, Samtools, Subread/featureCounts, Perl

---

## 🔧 Étapes du pipeline

### `install_RNAseq_Maurel.sh`
- Crée l’environnement conda
- Télécharge les FASTQ (Zenodo)
- Télécharge chr18 + annotation Gencode v24lift37
- Décompresse FASTA & GTF
- Construit l’index STAR

### `run_RNAseq_Maurel.sh`
- Nettoyage des reads (Trimmomatic)
- Alignement STAR (BAM trié)
- Indexation BAM
- Comptage par gène (featureCounts)
- Conversion gene_id → gene_name
- Filtration chr18

---

## ✅ Vérification rapide

```bash
ls Mapping/STAR_alignment_output/*.bam
samtools idxstats Mapping/STAR_alignment_output/*.bam | head
head count_table/table_gene_name_counts.txt
```

---

## 🛠️ Dépannage

| Problème | Solution |
|---|---|
Conda non activé | `source ~/miniconda3/etc/profile.d/conda.sh` |
Pairs FASTQ non détectées | Format obligatoire `*.R1.fastq.gz` / `*.R2.fastq.gz` |
Erreur annotation | Vérifier présence du GTF téléchargé |
Manque mémoire index STAR | OK pour chr18 (4–8 Go), ~32 Go pour genome complet |

---

## 🧬 Contexte biologique

- Cellules **H358** — EMT induite par **ZEB1**
- PolyA+ RNA-seq, paired-end 2×100 nt
- **CDH2** doit être plus exprimé en Day7 (marqueur EMT)
- Exemples de gènes d’intérêt : C18orf21, SLC39A6

---

## 📚 Références

Yang Y *et al.*, *Mol Cell Biol*, 2016  
Dataset SRA : SRP066794  
Références : UCSC hg19 chr18, Gencode v24lift37  
Outils : FastQC, Trimmomatic, STAR, Samtools, Subread/featureCounts, IGV

---

## 🏷️ Crédits

Scripts : Maurel & Tchuisse  
Encadrement : Daniel Gautheret — I2BC
