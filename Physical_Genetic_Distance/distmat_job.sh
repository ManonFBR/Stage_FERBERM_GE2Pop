#!/bin/bash
#SBATCH --job-name=distmat_parallel
#SBATCH --output=distmat_%j.log
#SBATCH --error=distmat_%j.err
#SBATCH -p cpu-dedicated
#SBATCH -A dedicated-cpu@cirad-normal
#SBATCH -t 02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module load emboss/6.6.0

INPUT_DIR="/storage/simple/users/ferberm/replicated/DATA_STAGE/FAA/Aligned"
OUTPUT_DIR="/storage/simple/users/ferberm/replicated/DATA_STAGE/Distmat"
mkdir -p $OUTPUT_DIR

# Boucle sur tous les fichiers alignés
for f in $INPUT_DIR/*_aligned.fasta; do
    base=$(basename $f _aligned.fasta)
    # Lancer chaque distmat en arrière-plan
    distmat -sequence $f -protmethod 2 -outfile $OUTPUT_DIR/${base}.distmat -auto &
done

# Attendre que tous les jobs en arrière-plan soient terminés
wait

echo "Tous les fichiers .distmat sont générés."
