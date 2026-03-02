#!/bin/bash
#SBATCH --job-name=mafft_align_parallel
#SBATCH --output=mafft_align_%j.log
#SBATCH --error=mafft_align_%j.err
#SBATCH -p cpu-dedicated
#SBATCH -A dedicated-cpu@cirad-normal
#SBATCH -t 04:00:00             # adapter selon la taille des fichiers
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2       # 2 cœurs par MAFFT
#SBATCH --mem=64G                # Mémoire totale disponible pour le job

module load mafft

INPUT_DIR="FAA/INPUT_DATA_MAFFT"
OUTPUT_DIR="FAA/Aligned"
mkdir -p $OUTPUT_DIR

# Boucle sur tous les fichiers FAA
for f in $INPUT_DIR/*.faa; do
    base=$(basename $f .faa)
    echo "Alignement de $base ..."
    # Aligner le fichier avec MAFFT
    mafft --auto --thread $SLURM_CPUS_PER_TASK "$f" > "$OUTPUT_DIR/${base}_aligned.fasta"
    echo "Terminé : $base"
done

echo "Tous les alignements sont terminés."
