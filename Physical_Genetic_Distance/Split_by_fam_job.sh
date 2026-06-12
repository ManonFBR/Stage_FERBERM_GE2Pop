#!/bin/bash
#SBATCH --job-name=Split_by_fam_job
#SBATCH --output=Split_by_fam_job_%j.out
#SBATCH --error=Split_by_fam_job_%j.err
#SBATCH -p cpu-dedicated
#SBATCH -A dedicated-cpu@cirad
#SBATCH -t 00:30:00                        # Durée max 30 min (à changer si besoin mais pour moi ça suffisait)
#SBATCH --cpus-per-task=1                  # 1 cœur
#SBATCH --mem-per-cpu=2G

# -----------------------------
# Charger Python
# -----------------------------
module load python/3.11

# -----------------------------
# Dossiers et fichiers d'entrée
# -----------------------------
GFF_FILE="/storage/simple/users/ferberm/replicated/DATA_STAGE/input.gff"
FASTA_FILE="/storage/simple/users/ferberm/replicated/DATA_STAGE/proteins.fasta"
OUTPUT_DIR="/storage/simple/users/ferberm/replicated/DATA_STAGE/SPLIT_RESULTS"
GFF_OUTPUT_DIR="$OUTPUT_DIR/GFF_by_fam"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$GFF_OUTPUT_DIR"

# -----------------------------
# Étape 1 : Split du GFF par famille
# -----------------------------
echo "[$(date)] Étape 1 : Split GFF par famille..."

python Split_gff_by_fam.py "$GFF_FILE" "$GFF_OUTPUT_DIR"

if [ $? -ne 0 ]; then
    echo "[ERREUR] Split_gff_by_fam.py a échoué. Abandon."
    exit 1
fi

echo "[$(date)] Split GFF terminé."

# -----------------------------
# Étape 2 : Split du FASTA par famille
# (utilise les GFF générés à l'étape 1)
# -----------------------------
echo "[$(date)] Étape 2 : Split FASTA par famille..."

python Split_fasta_by_fam.py "$FASTA_FILE" "$GFF_OUTPUT_DIR"

if [ $? -ne 0 ]; then
    echo "[ERREUR] Split_fasta_by_fam.py a échoué."
    exit 1
fi

echo "[$(date)] Split FASTA terminé."
echo "[$(date)] Job terminé avec succès. Résultats dans : $OUTPUT_DIR"
