"""
Script to calculate similarity scores between adjacent protein sequences using BLOSUM62 matrix and MAFFT alignment.

Aligns adjacent protein sequences using MAFFT.

Calculates similarity scores between sequences using BLOSUM62 matrix, and normalize the scores to a percentage.

Args:
    input_fasta (str): Path to FASTA file
    output_tsv (str): Path to TSV file
"""


from Bio import SeqIO
from Bio import AlignIO
import subprocess
import tempfile
import os
import shutil
import pandas as pd
import sys

# ----------------------------
# CHECK : MAFFT installé ?
# ----------------------------
try:
    subprocess.run(
        ["mafft", "--version"],
        capture_output=True,
        check=True
    )
except Exception:
    raise RuntimeError("MAFFT n'est pas installé ou pas fonctionnel")


# -------------------
# Paramètres de gaps
# -------------------
GAP_OPEN = -7
GAP_EXTEND = -1

# ----------------------
# Matrice BLOSUM62
# ----------------------
BLOSUM62 = {
    ('A','A'): 4,  ('A','R'): -1, ('A','N'): -2, ('A','D'): -2, ('A','C'): 0,
    ('A','Q'): -1, ('A','E'): -1, ('A','G'): 0,  ('A','H'): -2, ('A','I'): -1,
    ('A','L'): -1, ('A','K'): -1, ('A','M'): -1, ('A','F'): -2, ('A','P'): -1,
    ('A','S'): 1,  ('A','T'): 0,  ('A','W'): -3, ('A','Y'): -2, ('A','V'): 0,
    ('R','R'): 5,  ('R','N'): 0,  ('R','D'): -2, ('R','C'): -3,
    ('R','Q'): 1,  ('R','E'): 0,  ('R','G'): -2, ('R','H'): 0,  ('R','I'): -3,
    ('R','L'): -2, ('R','K'): 2,  ('R','M'): -1, ('R','F'): -3, ('R','P'): -2,
    ('R','S'): -1, ('R','T'): -1, ('R','W'): -3, ('R','Y'): -2, ('R','V'): -3,
    ('N','N'): 6,  ('N','D'): 1,  ('N','C'): -3,
    ('N','Q'): 0,  ('N','E'): 0,  ('N','G'): 0,  ('N','H'): 1,  ('N','I'): -3,
    ('N','L'): -3, ('N','K'): 0,  ('N','M'): -2, ('N','F'): -3, ('N','P'): -2,
    ('N','S'): 1,  ('N','T'): 0,  ('N','W'): -4, ('N','Y'): -2, ('N','V'): -3,
    ('D','D'): 6,  ('D','C'): -3,
    ('D','Q'): 0,  ('D','E'): 2,  ('D','G'): -1, ('D','H'): -1, ('D','I'): -3,
    ('D','L'): -4, ('D','K'): -1, ('D','M'): -3, ('D','F'): -3, ('D','P'): -1,
    ('D','S'): 0,  ('D','T'): -1, ('D','W'): -4, ('D','Y'): -3, ('D','V'): -3,
    ('C','C'): 9,
    ('C','Q'): -3, ('C','E'): -4, ('C','G'): -3, ('C','H'): -3, ('C','I'): -1,
    ('C','L'): -1, ('C','K'): -3, ('C','M'): -1, ('C','F'): -2, ('C','P'): -3,
    ('C','S'): -1, ('C','T'): -1, ('C','W'): -2, ('C','Y'): -2, ('C','V'): -1,
    ('Q','Q'): 5,  ('Q','E'): 2,  ('Q','G'): -2, ('Q','H'): 0,  ('Q','I'): -3,
    ('Q','L'): -2, ('Q','K'): 1,  ('Q','M'): 0,  ('Q','F'): -3, ('Q','P'): -1,
    ('Q','S'): 0,  ('Q','T'): -1, ('Q','W'): -2, ('Q','Y'): -1, ('Q','V'): -2,
    ('E','E'): 5,  ('E','G'): -2, ('E','H'): 0,  ('E','I'): -3,
    ('E','L'): -3, ('E','K'): 1,  ('E','M'): -2, ('E','F'): -3, ('E','P'): -1,
    ('E','S'): 0,  ('E','T'): -1, ('E','W'): -3, ('E','Y'): -2, ('E','V'): -2,
    ('G','G'): 6,  ('G','H'): -2, ('G','I'): -4,
    ('G','L'): -4, ('G','K'): -2, ('G','M'): -3, ('G','F'): -3, ('G','P'): -2,
    ('G','S'): 0,  ('G','T'): -2, ('G','W'): -2, ('G','Y'): -3, ('G','V'): -3,
    ('H','H'): 8,  ('H','I'): -3,
    ('H','L'): -3, ('H','K'): -1, ('H','M'): -2, ('H','F'): -1, ('H','P'): -2,
    ('H','S'): -1, ('H','T'): -2, ('H','W'): -2, ('H','Y'): 2,  ('H','V'): -3,
    ('I','I'): 4,
    ('I','L'): 2,  ('I','K'): -3, ('I','M'): 1,  ('I','F'): 0,  ('I','P'): -3,
    ('I','S'): -2, ('I','T'): -1, ('I','W'): -3, ('I','Y'): -1, ('I','V'): 3,
    ('L','L'): 4,  ('L','K'): -2, ('L','M'): 2,  ('L','F'): 0,  ('L','P'): -3,
    ('L','S'): -2, ('L','T'): -1, ('L','W'): -2, ('L','Y'): -1, ('L','V'): 1,
    ('K','K'): 5,  ('K','M'): -1, ('K','F'): -3, ('K','P'): -1,
    ('K','S'): 0,  ('K','T'): -1, ('K','W'): -3, ('K','Y'): -2, ('K','V'): -2,
    ('M','M'): 5,  ('M','F'): 0,  ('M','P'): -2,
    ('M','S'): -1, ('M','T'): -1, ('M','W'): -1, ('M','Y'): -1, ('M','V'): 1,
    ('F','F'): 6,  ('F','P'): -4,
    ('F','S'): -2, ('F','T'): -2, ('F','W'): 1,  ('F','Y'): 3,  ('F','V'): -1,
    ('P','P'): 7,
    ('P','S'): -1, ('P','T'): -1, ('P','W'): -4, ('P','Y'): -3, ('P','V'): -2,
    ('S','S'): 4,  ('S','T'): 1,  ('S','W'): -3, ('S','Y'): -2, ('S','V'): -2,
    ('T','T'): 5,  ('T','W'): -2, ('T','Y'): -2, ('T','V'): 0,
    ('W','W'): 11, ('W','Y'): 2,  ('W','V'): -3,
    ('Y','Y'): 7,  ('Y','V'): -1,
    ('V','V'): 4,
}

def blosum62_score(a, b):
    return BLOSUM62.get((a,b), BLOSUM62.get((b,a), -4))

# ---------------------------------
# Fonction pour aligner avec MAFFT
# ---------------------------------
def align_with_mafft_biopython(seq1, seq2, debug=True):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(f">{seq1.id}\n{seq1.seq}\n>{seq2.id}\n{seq2.seq}\n")
        temp_file = f.name

    result = subprocess.run(
        [
            'mafft',
            '--globalpair',
            '--maxiterate', '0',
            '--thread', '1',
            '--op', str(abs(GAP_OPEN)),
            '--ep', str(abs(GAP_EXTEND)),
            temp_file
        ],
        capture_output=True,
        text=True,
        check=True
    )

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f_out:
        f_out.write(result.stdout)
        output_file = f_out.name

    alignment = AlignIO.read(output_file, "fasta")
    if len(alignment) != 2:
        raise ValueError("Alignement invalide : nombre de séquences ≠ 2")

    os.unlink(temp_file) # Suppression des fichiers temporaires

    if debug:
        shutil.move(output_file, os.path.basename(output_file))
        return alignment, os.path.basename(output_file)
    else:
        os.unlink(output_file)
        return alignment, None

# -------------------------------------------------------------
# Calcul du score BLOSUM62 + GAPS puis normalisation [0 ; 100]
# -------------------------------------------------------------
def calculate_blosum_score_norm(alignment):
    seq1 = str(alignment[0].seq)
    seq2 = str(alignment[1].seq)

    total_score = 0
    score_max = 0
    score_min = 0

    MIN_BLOSUM = min(BLOSUM62.values())  # -4

    in_gap1 = False
    in_gap2 = False

    for a, b in zip(seq1, seq2):

        if a == '-':
            if not in_gap1:
                total_score += GAP_OPEN
                score_min += GAP_OPEN
            else:
                total_score += GAP_EXTEND
                score_min += GAP_EXTEND
            in_gap1 = True
            in_gap2 = False
            continue

        if b == '-':
            if not in_gap2:
                total_score += GAP_OPEN
                score_min += GAP_OPEN
            else:
                total_score += GAP_EXTEND
                score_min += GAP_EXTEND
            in_gap2 = True
            in_gap1 = False
            continue


        in_gap1 = in_gap2 = False

        # SUBSTITUTION
        s = blosum62_score(a, b)
        total_score += s

        score_max += blosum62_score(a, a)
        score_min += MIN_BLOSUM

    if score_max == score_min:
        return 0.0

    normalized_score = (total_score - score_min) / (score_max - score_min) * 100

    return max(0.0, min(100.0, normalized_score))

# --------------------------------------------------
# Fonction principale faisant appel aux précédentes
# --------------------------------------------------
def main(input_fasta, output_tsv, debug=True):
    genes = list(SeqIO.parse(input_fasta, "fasta"))
    results = []

    for i in range(len(genes) - 1):
        alignment, aln_file = align_with_mafft_biopython(genes[i], genes[i+1], debug)
        similarity = calculate_blosum_score_norm(alignment)

        results.append({
            "gene_1": genes[i].id,
            "gene_2": genes[i+1].id,
            "similarity": similarity,
            "alignment": aln_file
        })

    df = pd.DataFrame(results)
    df.to_csv(output_tsv, sep="\t", index=False)
    print("Done.")

# ---------------------------------
# Appel de la fonction principale
# ---------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python Similarity_BLOSUM62.py <input_fasta> <output_tsv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
