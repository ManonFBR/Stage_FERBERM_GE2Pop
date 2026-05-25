# Physical_Genetic_Distance 

## Présentation générale

Ce dossier contient les scripts relatifs au calcul et à l'analyse des distances physiques entre les gènes et des similarités de séquences protéiques associées.

## Ordre d'utilisation
**Indispensable :**
- ☆ *Split_gff_by_fam.py* avec le fichier GFF complet
- ☆ *Split_fasta_by_fam.py* avec le fichier FASTA complet et l'un des fichiers GFF de sortie du script précédent
  
**Scripts dépendants de ceux-ci :**  
- *GFF_to_dictionary.py* / ☆ *GFF_to_dictionary_allpairs.py* avec l'un des fichiers GFF de sortie de *Split_gff_by_fam.py*
- ☆ *New_SIM_SCRIPT.py* avec l'un des fichiers FAA de sortie de *Split_fasta_by_fam.py*
  
**Autres scripts :**  
Scripts pour lancer des jobs SLURM qui effectuent les scripts précédents ou d'autres.

## Scripts et utilisation

### ☆ Split_gff_by_fam.py
Script pour extraire les gènes d'un fichier GFF et les exporter dans un autre fichier GFF en fonction de la sous-famille (NLR, RLK, RLP).  
Ce script extrait les gènes du fichier GFF donné en input (contenant des exons, CDS, etc.) et le split en 3 fichiers GFF correspondant chacun à une sous-famille.  

**Ligne de commande :**

```bash
python Split_gff_by_fam.py input_file.gff output_prefix
```
**Attention !**  
Ligne 126 : *if attr.startswith('Fam='):*  
--> Remplacer "*Fam=*" par le terme utilisé pour la famille dans le GFF (ex : "*Family =*", "*Famille :*" etc.)


### ☆ Split_fasta_by_fam.py  
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

Si aucun nom de sortie n'est spécifié, l'output est automatiquement intitulé "distance_entre_genes_$basename.tsv". 

### ☆ GFF_to_dictionary_allpairs.py
  Script pour parser un fichier GFF et le charger en un dictionnaire python.  
  Ce script crée une clé pour chaque attribut du fichier GFF, calcule la distance physique (en paires de bases) entre toutes les paires de gènes possibles sur un même chromosome et vérifie si chacun des gènes de cette paire est canonique (fonctionnel) ou non.  
  Les gènes situés sur différents chromosomes sont également pris en compte, le fichier de sortie contient une colonne "same_chromosome" indiquant si les deux gènes de la paire se situent sur le même chromosome ou non (yes/no). Dans le cas où ces gènes se trouvent sur différents chromosomes, le script renvoie "-" pour la distance physique entre les gènes. 
  Les résultats sont exportés en un fichier CSV.  
  Ce script fonctionne exactement comme le précédent, à la seule différence qu'il calcule les distances pairwise pour toutes les paires de gènes et non pas seulement entre gènes adjacents.  
  
```bash
python GFF_to_dictionary_allpairs.py input_file.gff output_file.csv
```
Si aucun nom de sortie n'est spécifié, l'output est automatiquement intitulé "distance_entre_genes_$basename.tsv". 

### ☆ New_SIM_SCRIPT.py 
  Script retenu pour calculer des scores de similarité entre deux séquences protéiques adjacentes à partir du programme d'alignement MAFFT (Katoh et al. 2002 - https://doi.org/10.1093/nar/gkf436 ; Katoh et al. 2013 - https://doi.org/10.1093/molbev/mst010) et de la matrice de substitution BLOSUM62 (Henikoff & Henikoff 1992 - https://doi.org/10.1073/pnas.89.22.10915).  
  Ce script aligne avec MAFFT les séquences adjacentes du fichier FASTA donné, calcule des scores de similarité en utilisant la matrice BLOSUM62, puis normalise ces scores pour obtenir des pourcentages de similarité. Version améliorée et finale du script *Similarity_BLOSUM62.py*.   

```bash  
python New_SIM_SCRIPT.py input.fasta output_file.csv
```

### Similarity_BLOSUM62.py 
  Script pour calculer des scores de similarité entre deux séquences protéiques adjacentes à partir du programme d'alignement MAFFT (Katoh et al. 2002 - https://doi.org/10.1093/nar/gkf436 ; Katoh et al. 2013 - https://doi.org/10.1093/molbev/mst010) et de la matrice de substitution BLOSUM62 (Henikoff & Henikoff 1992 - https://doi.org/10.1073/pnas.89.22.10915).  
  Ce script aligne avec MAFFT les séquences adjacentes du fichier FASTA donné, calcule des scores de similarité en utilisant la matrice BLOSUM62, puis normalise ces scores pour obtenir des pourcentages de similarité.    

```bash  
python Similarity_BLOSUM62.py input.fasta output_file.csv
```
**Paramètres :**  
*Debug = False* permet de supprimer les fichiers temporaires d'alignement, ce qui réduit les bugs liés à la taille des fichiers (pour N protéines, N^2 alignements sont générés).  
Si vous souhaitez conserver les alignements pour les visualiser, changer par *Debug = True*.  
  
**Résolution de problèmes / dépannage :**  
  Voici les erreurs fréquentes et leurs solutions lors de l'utilisation de ce script :
#### Modules / commandes introuvables  
1. Erreur : Lmod has detected the following error: These module(s) or extension(s) exist but cannot be loaded as requested: "python/3.9"  
   Cause : le module Python demandé n'existe pas ou dépend d'un autre module.  
   Solution :
```bash  
module spider python/3.9
```
=> "You will need to load all module(s) on any one of the lines below before the "python/3.9" module is available to load."  
    Pour mon cas, le module parent manquant était bioinfo-ifb.
```bash  
module load bioinfo-ifb
```
2. Erreur : PermissionError: [Errno 13] Permission denied: 'mafft'  
            RuntimeError: MAFFT n'est pas installé ou pas fonctionnel  
   Cause : MAFFT n'est pas dans le path/nécessite d'être chargé/nécessite un module parent  
   Solution :
```bash  
module spider mafft
```
=> "You will need to load all module(s) on any one of the lines below before the "mafft/7.515" module is available to load." (encore bioinfo-ifb).  
```bash  
module purge
module load bioinfo-ifb
module load mafft/7.515
which mafft
mafft --version
```

#### Autres bugs possibles



### run_array_Sim_BLOSUM_job.sh
  Script servant à lancer un job sur le cluster IO pour faire tourner le script Python "Similarity_BLOSUM62.py".  
    
#SBATCH --output=/path_to_result_files/Sim_BLOSUM_%A_%a.out  
#SBATCH --error=/path_to_error_files/Sim_BLOSUM_%A_%a.err  
#SBATCH --array=1-N : N = Nombre de fichiers FASTA à traiter parallèlement.  

L19-22 : Activation de l'environnement pandas facultative.  
L27-28 : Remplacer "[...]" par le path du dossier contenant les fichiers FASTA (L27) et du dossier créé pour stocker les résultats (L28).  
L53 : Remplacer par le path du script python Similarity_BLOSUM62.py.  

```bash  
sbatch run_array_Sim_BLOSUM_job.sh
```

**Résolution de problèmes / dépannage :**
  Voici les erreurs fréquentes et leurs solutions lors de l'utilisation de ce script :  
#### Jobs SLURM bloqués ou annulés (PD)
1. Erreur :  
```bash  
squeue -u ton_user
```
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
       JOB_NUM[1] cpu-dedic Sim_BLOS  ton_user PD       0:00      1 (AssocMaxWallDurationPerJobLimit)
Cause : Limite de temps maximale pour ton compte / partition atteinte.  
Solution : Vérifier les limites du compte puis adapter la durée des jobs dans le script SLURM.  
```bash  
sacctmgr show assoc user=<ton_user> format=User,Account,Partition,MaxWall
```
2. Erreur (dans le fichier .err) : slurmstepd-io-cpu-09: error: *** JOB <ID> CANCELLED DUE TO TIME LIMIT ***  
   Cause : Le temps du job dans le script SLURM était inférieur au temps nécessaire pour traiter toutes les séquences.  
   Solution : Augmenter le #SBATCH -t dans le script SLURM.  


### BLASTP_all-vs-all.sh 
  Script servant à lancer un job sur le cluster IO pour créer une database à partir d'un fichier FASTA grâce à la commande MAKEBLASTDB de BLAST+ (Camacho et al. 2009 - https://doi.org/10.1186/1471-2105-10-421), puis lancer un BLAST "all-against-all" de ce fichier FASTA sur la database créée.  
  Ce script BLAST toutes les protéines de la database contre toutes les autres afin d'obtenir les hits BLAST de toutes les paires possibles.  

```bash  
sbatch BLASTP_all-vs-all.sh input.fasta output_prefix
``` 

### mafft_job_copy.sh 
  Script servant à lancer un job sur le cluster IO pour aligner les séquences d'un fichier FASTA à l'aide de MAFFT (Katoh et al. 2002 - https://doi.org/10.1093/nar/gkf436 ; Katoh et al. 2013 - https://doi.org/10.1093/molbev/mst010).  
  Ce script utilise MAFFT pour aligner toutes les séquences protéiques d'un fichier FASTA donné et retourne l'alignement en fichier FASTA qui sera utilisé en input du script suivant celui-ci (distmat_job.sh).

```bash  
sbatch mafft_job_copy.sh
``` 

### distmat_job.sh 
  Script servant à lancer un job sur le cluster IO pour calculer des distances pairwise entre séquences protéiques alignées (précédemment avec mafft_job_copy.sh) grâce à l'outil distmat du package EMBOSS (https://galaxy-iuc.github.io/emboss-5.0-docs/distmat.html).  
  Ce script utilise distmat pour sortir une matrice de distances génétiques entre paires de séquences protéiques alignées.

```bash  
sbatch distmat_job.sh
``` 





















