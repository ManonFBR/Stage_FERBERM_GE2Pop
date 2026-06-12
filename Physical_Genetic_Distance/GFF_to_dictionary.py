import sys
import re
import pandas as pd

"""
Script to parse a GFF file and load it into a python dictionary.
Script pour parser un fichier GFF et le charger en dictionnaire python.

Creates a key for each feature in the GFF file.
Crée une "clé" pour chaque attribut dans le fichier GFF.

Calculates the physical distance (base pairs) between two adjacent genes and checks the class type (canonical or non-canonical) of each gene of the pair.
Calcule la distance physique (en paires de bases) entre deux gènes adjacents et check la catégorie de fonctionnalité (canonique ou non) pour chaque gène de la paire.

Exports the results to a CSV file.
Exporte les résultats en fichier CSV.

Args:
    gff_file (str): Path to GFF file
"""

# --------------------------------------
# Fonction pour créer le dictionnaire
# --------------------------------------
def parse_gff_to_dict(gff_file, feature_type="gene", preview=3):

    dico_gff = {}

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): #ignorer les commentaires et lignes vides
                continue

            cols = line.strip().split("\t")
            if len(cols) < 9: #sécurité si ligne à moins de 9 colonnes (GFF buggé donc saute la ligne)
                continue

            chromosome, source, feature, start, end, score, strand, phase, attributes = cols

            if feature_type and feature != feature_type: #sécurité s'il reste des lignes =/= gènes
                continue

            raw_attrs = attributes.replace(" / ", ";").split(";") #uniformisation des séparateurs dans les attributs

            attr_dict = {}
            gene_class = None

            for attr in raw_attrs:
                attr = attr.strip()

                #extraction de Gene-Class
                match_gc = re.search(r"Gene-Class:([A-Za-z\-]+)", attr)
                if match_gc:
                    gene_class = match_gc.group(1) #garde seulement le groupe ()
                    attr = re.sub(r"Gene-Class:[A-Za-z\-]+", "", attr).strip() #évite les doublons

                #nettoyage de Origin-Class
                if attr.startswith("Origin-Class:"):
                    attr = attr.split(" ")[0]  # garde seulement Canonical / Non-canonical

                if "=" in attr:
                    key, value = attr.split("=", 1)
                    attr_dict[key] = value
                elif ":" in attr:
                    key, value = attr.split(":", 1)
                    attr_dict[key] = value

            gene_id = attr_dict.get("ID")
            if gene_id is None:
                continue

            dico_gff[gene_id] = {
                "chromosome": chromosome,
                "source": source,
                "feature": feature,
                "start": int(start),
                "end": int(end),
                "score": score,
                "strand": strand,
                "phase": phase,
                "Gene-Class": gene_class,
                **attr_dict
            }

    # Aperçu
    print(f"\nAperçu des {preview} premiers gènes dans {gff_file}:")
    for i, (gid, info) in enumerate(dico_gff.items()):
        if i >= preview:
            break
        print(f"{gid}: {info}")

    return dico_gff

# ------------------------------------------------------------
# Fonction pour calculer les distances + qualifier les paires
# ------------------------------------------------------------

def calculate_phy_dist_from_dict(dico_gff):

    results = []

    # transformer en DataFrame pour trier facilement
    df = pd.DataFrame.from_dict(dico_gff, orient="index")
    df["gene_id"] = df.index
    df = df.sort_values(by=["chromosome", "start"])

    prev_row = None

    for _, row in df.iterrows():
        if prev_row is not None and row["chromosome"] == prev_row["chromosome"]:

            distance = row["start"] - prev_row["end"]

            class_1 = prev_row["Gene-Class"]
            class_2 = row["Gene-Class"]

            if class_1 == "Canonical" and class_2 == "Canonical":
                pair_type = "Canonical"
            elif class_1 == "Non-canonical" and class_2 == "Non-canonical":
                pair_type = "Non-canonical"
            else:
                pair_type = "Mixed"

            results.append({
                "gene_1": prev_row["gene_id"],
                "gene_2": row["gene_id"],
                "chromosome": row["chromosome"],
                "distance": distance,
                "gene_class_1": class_1,
                "gene_class_2": class_2,
                "pair_type": pair_type,
                "start_1": prev_row["start"],
                "end_1": prev_row["end"],
                "start_2": row["start"],
                "end_2": row["end"]

            })

        prev_row = row

    return results

# -------------------------------------
# Fonction pour exporter les résultats
# -------------------------------------
def export_phy_dist_to_tsv(results, gff_filename):
    distance_df = pd.DataFrame(results)
    base_name = gff_filename.rsplit('.', 1)[0]
    output_filename = f"distance_entre_genes_{base_name}.tsv"
    distance_df.to_csv(output_filename, sep="\t", index=False)
    return distance_df


# -------------------------------
# Appel des fonctions
# -------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage : python GFF_to_dictionary.py <fichier.gff>")
        sys.exit(1)

    gff_filename = sys.argv[1]

    #Parser le GFF en dictionnaire
    dico_gff = parse_gff_to_dict(gff_filename)

    print(f"\nNombre total de gènes parsés : {len(dico_gff)}")

    #Calculer les distances physiques et types de paires
    results = calculate_phy_dist_from_dict(dico_gff)

    #Exporter les résultats en TSV
    distance_df = export_phy_dist_to_tsv(results, gff_filename)

    #Aperçu des résultats
    print("\nAperçu du fichier final :")
    print(distance_df.head())
