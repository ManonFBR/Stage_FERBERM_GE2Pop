import sys
import os
import re
import csv
from collections import defaultdict

#!/usr/bin/env python3
"""
count_exons_per_gene.py
-----------------------
Extrait le nombre d'exons par gène depuis un fichier GFF.

Usage:
    python count_exons_per_gene.py <fichier.gff> [fichier_sortie.csv]

Si aucun fichier de sortie n'est spécifié, le CSV sera créé dans le même
répertoire que le fichier d'entrée avec le suffixe '_exon_counts.csv'.
"""

def parse_attributes(attr_string):
    """
    Parse la colonne 9 du GFF pour extraire les valeurs de ID, Parent et Fam.
    Retourne un dictionnaire {clé: valeur}.
    """
    attrs = {}

    # Extraire ID
    id_match = re.search(r'ID=([^;]+)', attr_string)
    if id_match:
        attrs['ID'] = id_match.group(1).strip()

    # Extraire Parent
    parent_match = re.search(r'Parent=([^;]+)', attr_string)
    if parent_match:
        attrs['Parent'] = parent_match.group(1).strip()

    # Extraire Fam (NBS-LRR, LRR-RLK ou LRR-RLP)
    fam_match = re.search(r'Fam=([^;]+)', attr_string)
    if fam_match:
        attrs['Fam'] = fam_match.group(1).strip()

    return attrs


def count_exons_per_gene(gff_path, output_path=None):
    """
    Lit un fichier GFF et produit un CSV avec le nombre d'exons par gène.

    Paramètres
    ----------
    gff_path  : str  – chemin vers le fichier GFF en entrée
    output_path : str – chemin du CSV de sortie (optionnel)

    Retourne
    --------
    output_path : str – chemin du fichier CSV créé
    """
    if not os.path.isfile(gff_path):
        raise FileNotFoundError(f"Fichier introuvable : {gff_path}")

    # Définir le chemin de sortie par défaut
    if output_path is None:
        base = os.path.splitext(gff_path)[0]
        output_path = base + "_exon_counts.csv"

    # ------------------------------------------------------------------ #
    # Passe 1 : collecter les informations des lignes "gene"              #
    #   → chr, ID du gène, Fam                                            #
    # ------------------------------------------------------------------ #
    gene_info = {}   # gene_id -> {'chr': ..., 'fam': ...}

    # ------------------------------------------------------------------ #
    # Passe 2 : compter les exons                                         #
    #   Les exons ont un Parent = <mrna_id>                               #
    #   Les mRNA ont un Parent = <gene_id>                                #
    #   → on construit un mapping mrna_id -> gene_id                      #
    # ------------------------------------------------------------------ #
    mrna_to_gene = {}   # mrna_id  -> gene_id
    exon_counts  = defaultdict(int)  # gene_id -> nombre d'exons

    with open(gff_path, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')

            # Ignorer les commentaires et lignes vides
            if not line or line.startswith('#'):
                continue

            cols = line.split('\t')
            if len(cols) < 9:
                continue

            chrom      = cols[0]
            feat_type  = cols[2]
            attr_str   = cols[8]

            attrs = parse_attributes(attr_str)
            feat_id = attrs.get('ID', '')
            parent  = attrs.get('Parent', '')
            fam     = attrs.get('Fam', 'NA')

            if feat_type == 'gene':
                if feat_id:
                    gene_info[feat_id] = {'chr': chrom, 'fam': fam}

            elif feat_type in ('mRNA', 'transcript'):
                # Parent d'un mRNA = le gène
                if feat_id and parent:
                    mrna_to_gene[feat_id] = parent

            elif feat_type == 'exon':
                # Parent d'un exon = le mRNA
                # On remonte jusqu'au gène via mrna_to_gene
                gene_id = mrna_to_gene.get(parent)
                if gene_id:
                    exon_counts[gene_id] += 1

    # ------------------------------------------------------------------ #
    # Écriture du CSV de sortie                                           #
    # ------------------------------------------------------------------ #
    rows_written = 0
    with open(output_path, 'w', newline='') as out_fh:
        writer = csv.writer(out_fh, delimiter='\t')
        writer.writerow(['Chr', 'ID', 'Nb_exons', 'Family'])

        for gene_id, info in gene_info.items():
            nb_exons = exon_counts.get(gene_id, 0)
            writer.writerow([info['chr'], gene_id, nb_exons, info['fam']])
            rows_written += 1

    print(f"✓ {rows_written} gènes traités.")
    print(f"✓ Fichier de sortie : {output_path}")
    return output_path


# ---------------------------------------------------------------------- #
# Point d'entrée                                                          #
# ---------------------------------------------------------------------- #
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    gff_file    = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) >= 3 else None

    count_exons_per_gene(gff_file, output_file)
