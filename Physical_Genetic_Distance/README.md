# Physical_Genetic_Distance 

## Présentation générale

Ce dossier contient les scripts relatifs au calcul et à l'analyse des distances physiques entre les gènes et des similarités de séquences protéiques associées.


## Scripts et utilisation

### Split_gff_by_fam.py
Script pour extraire les gènes d'un fichier GFF et les exporter dans un autre fichier GFF en fonction de la sous-famille (NLR, RLK, RLP).
Ce script extrait les gènes du fichier GFF donné en input (contenant des exons, CDS, etc.) et le split en 3 fichiers GFF correspondant chacun à une sous-famille.

**Ligne de commande :**

```bash
python Split_gff_by_fam.py input_file.gff output_prefix
```

### Split_fasta_by_fam.py  
  Script pour séparer un fichier FASTA en 3 fichiers individuels en fonction de la sous-famille du gène associé.
  Ce script extrait les IDs des gènes à partir d'un fichier GFF, cherche des correspondances dans le fichier FASTA donné afin d'assigner chaque séquence protéique à la sous-famille correspondante, puis split ce fichier en 3 fichiers FASTA en fonction de la sous-famille.

**Ligne de commande :**

```bash
python Split_fasta_by_fam.py input.fasta gff_file output_prefix
```

### GFF_to_dictionary.py
  Script pour parser un fichier GFF et le charger en un dictionnaire python.
  Ce script crée une clé pour chaque attribut du fichier GFF, calcule la distance physique (en paires de bases) entre deux gènes adjacents et vérifie si chacun des gènes de cette paire est canonique (fonctionnel) ou non.
  Les résultats sont exportés en un fichier CSV.
  
```bash
python GFF_to_dictionary.py input_file.gff output_file.csv
```

### Similarity_BLOSUM62.py 
  Script pour calculer des scores de similarité entre deux séquences protéiques adjacentes à partir du programme d'alignement MAFFT (Katoh et al. 2002 ; Katoh et al. 2013) et de la matrice de substitution BLOSUM62 (Henikoff & Henikoff 1992).
  Ce script aligne avec MAFFT les séquences adjacentes du fichier FASTA donné, calcule des scores de similarité en utilisant la matrice BLOSUM62, puis normalise ces scores pour obtenir des pourcentages de similarité.  

```bash  
python Similarity_BLOSUM62.py input.fasta output_file.csv
```
  
### BLASTP_all-vs-all.sh 
  Script servant à lancer un job sur le cluster IO pour créer une database à partir d'un fichier FASTA grâce à la commande MAKEBLASTDB de BLAST+ (Camacho et al. 2009), puis lancer un BLAST "all-against-all" de ce fichier FASTA sur la database créée.
  Ce script BLAST toutes les protéines de la database contre toutes les autres afin d'obtenir les hits BLAST de toutes les paires possibles. 

```bash  
bash BLASTP_all-vs-all.sh input.fasta output_prefix
``` 



