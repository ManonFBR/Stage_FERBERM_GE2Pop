#!/bin/bash
#SBATCH --job-name=blastp_job
#SBATCH --output=blastp_job.out
#SBATCH --error=blastp_job.err
#SBATCH -p cpu-dedicated                   # Partition à utiliser
#SBATCH -A dedicated-cpu@cirad             # Account
#SBATCH -t 01:00:00                         # Durée max 1h
#SBATCH --cpus-per-task=8                  # Nombre de cœurs à utiliser
#SBATCH --mem-per-cpu=2G                   # Mémoire par cœur (2Go par défaut)

# Usage: ./blast_job.sh <fasta_file> [db_name] [output_file] [num_threads]

usage() {
    echo "Usage: $0 <fasta_file> [db_name] [output_file] [num_threads]"
    echo ""
    echo "Examples:"
    echo "  $0 proteins.faa"
    echo "  $0 proteins.faa my_db results.txt 16"
    exit 1
}

# Vérifier les arguments
if [[ $# -eq 0 ]]; then
    usage
fi

FASTA_FILE="$1"
DB_NAME="${2:-$(basename "$FASTA_FILE" .faa)_db}"
OUTPUT_FILE="${3:-$(basename "$FASTA_FILE" .faa)_blast_results.out}"
NUM_THREADS="${4:-8}"

# Vérifier que le fichier existe
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "Error: File '$FASTA_FILE' not found"
    exit 1
fi

# Charger le module BLAST
if command -v module &> /dev/null; then
    module load blast
fi

# Fonction de recherche BLAST
blast_search() {
    local fasta_file=$1
    local db_name=$2
    local output_file=$3
    local num_threads=$4

    echo "Running BLAST search..."
    echo "Input file: $fasta_file"
    echo "Database: $db_name"
    echo "Output: $output_file"
    echo "Threads: $num_threads"
    echo ""

    # Créer la base
    echo "Creating database..."
    makeblastdb -in "$fasta_file" -dbtype prot -out "$db_name"

    if [[ $? -ne 0 ]]; then
        echo "Error creating database"
        exit 1
    fi

    # Lancer blastp
    echo "Running BLAST search..."
    blastp -query "$fasta_file" \
           -db "$db_name" \
           -out "$output_file" \
           -outfmt 6 \
           -num_threads "$num_threads"

    if [[ $? -eq 0 ]]; then
        echo "Done. Results saved to: $output_file"
    else
        echo "Error running BLAST"
        exit 1
    fi
}

# Exécuter
blast_search "$FASTA_FILE" "$DB_NAME" "$OUTPUT_FILE" "$NUM_THREADS"
