"""
Script to split a FASTA file into separate files for each family.

Extracts genes IDs from a GFF file.

Checks for corresponding gene IDs in the FASTA file to assign them to the appropriate family.

Exports the results to separate FASTA files.

Args:
    gff_file (str): Path to GFF file
    fasta_file (str): Path to FASTA file
"""


import sys
import os

def get_gene_ids_from_gff(gff_file):
    gene_ids = set()

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            attributes = fields[8]

            for attr in attributes.split(";"):
                if attr.startswith("ID="):
                    gene_ids.add(attr.replace("ID=", ""))

    return gene_ids


def split_fasta_by_ids(fasta_file, ids, output_fasta):
    with open(fasta_file) as fin, open(output_fasta, "w") as fout:
        write = False
        kept = 0

        for line in fin:
            if line.startswith(">"):
                seq_id = line[1:].split()[0]
                write = seq_id in ids
                if write:
                    kept += 1
            if write:
                fout.write(line)

    return kept


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python Code_split_fasta.py <proteins.fasta> <gff_directory>")
        sys.exit(1)

    fasta = sys.argv[1]
    gff_dir = sys.argv[2]

    species = os.path.basename(fasta).replace("_proteins.fasta", "")
    os.makedirs("FAA", exist_ok=True)

    gff_files = [f for f in os.listdir(gff_dir) if f.endswith(".gff")]

    if not gff_files:
        print("[ERROR] Aucun fichier .gff trouvé dans", gff_dir)
        sys.exit(1)

    for gff in gff_files:
        family = gff.replace(f"{species}_LRR_genes_", "").replace(".gff", "")
        gff_path = os.path.join(gff_dir, gff)

        print(f"\nTraitement {family}")

        ids = get_gene_ids_from_gff(gff_path)
        print(f"  {len(ids)} IDs trouvés")

        output_faa = os.path.join("FAA", f"{species}_{family}.faa")
        kept = split_fasta_by_ids(fasta, ids, output_faa)

        print(f"  {kept} séquences écrites → {output_faa}")
