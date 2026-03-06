#!/bin/bash
#SBATCH --job-name=Sim_BLOSUM_job
#SBATCH --output=/storage/simple/users/ferberm/replicated/DATA_STAGE/SIM_RESULTS/Sim_BLOSUM_%A_%a.out
#SBATCH --error=/storage/simple/users/ferberm/replicated/DATA_STAGE/SIM_RESULTS/Sim_BLOSUM_%A_%a.err
#SBATCH -p cpu-dedicated
#SBATCH -A dedicated-cpu@cirad-normal
#SBATCH -t 16:00:00                # Durée suffisante pour les gros fichiers
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-30                # Nombre de fichiers FASTA à traiter

# -----------------------------
# Purge anciens modules et charger les bons
# -----------------------------
module purge
module load bioinfo-ifb         
module load mafft/7.515
module load python/3.9

# -----------------------------
# Activer ton environnement pandas
# -----------------------------
source ~/env_pandas/bin/activate

# -----------------------------
# Dossiers
# -----------------------------
FASTA_DIR="/storage/simple/users/ferberm/replicated/DATA_STAGE/FASTAS_FILES"
OUTPUT_DIR="/storage/simple/users/ferberm/replicated/DATA_STAGE/SIM_RESULTS"
mkdir -p "$OUTPUT_DIR"

# -----------------------------
# Lister tous les fichiers FASTA
# -----------------------------
FASTAS=($FASTA_DIR/*.faa)

# -----------------------------
# Sélectionner le fichier correspondant à cette tâche
# -----------------------------
FASTA_FILE=${FASTAS[$SLURM_ARRAY_TASK_ID-1]}
BASE=$(basename "$FASTA_FILE" .faa)
OUT_FILE="$OUTPUT_DIR/${BASE}_similarity.tsv"

# -----------------------------
# Vérifications optionnelles
# -----------------------------
echo "Job $SLURM_ARRAY_TASK_ID: $FASTA_FILE -> $OUT_FILE"
echo "Python: $(which python) -- $(python --version)"
echo "MAFFT: $(which mafft) -- $(mafft --version)"

# -----------------------------
# Lancer le script Python
# -----------------------------
python /storage/simple/users/ferberm/replicated/DATA_STAGE/Similarity_BLOSUM62.py "$FASTA_FILE" "$OUT_FILE" --nproc 4
