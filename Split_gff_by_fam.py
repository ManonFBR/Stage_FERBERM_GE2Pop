"""
Script to extract genes from a GFF file and export them to another GFF file according to family.

Extracts genes from the GFF file containing CDS, exons, etc.

Split the GFF file into separate files for each family.

Exports the results to a GFF file.

Args:
    gff_file (str): Path to GFF file
"""


def filter_gff_genes(gff_file):
    """
    Extrait uniquement les gènes du fichier GFF (colonne 3 = "gene")

    Args:
        gff_file (str): Chemin vers le fichier GFF à traiter

    Returns:
        str: Chemin vers le fichier GFF filtré contenant uniquement les gènes
    """
    import os

    # Créer le nom du fichier de sortie
    input_basename = os.path.splitext(os.path.basename(gff_file))[0]
    output_dir = os.path.dirname(gff_file) or "."
    filtered_file = os.path.join(output_dir, f"{input_basename}_genes.gff")

    gene_count = 0

    try:
        with open(gff_file, 'r') as f_in, open(filtered_file, 'w') as f_out:
            for line_num, line in enumerate(f_in, 1):
                # Ignorer les lignes vides et les commentaires
                if line.strip() == '' or line.startswith('#'):
                    f_out.write(line)
                    continue

                # Split la ligne GFF
                columns = line.strip().split('\t')

                # Vérifier que la ligne a 9 colonnes
                if len(columns) < 9:
                    print(f"Attention: Ligne {line_num} n'a pas 9 colonnes, ignorée")
                    continue

                # Extraire le type (colonne 3, index 2)
                feature_type = columns[2]

                # Garder uniquement les gènes
                if feature_type == "gene":
                    f_out.write(line)
                    gene_count += 1

        print(f"Fichier GFF filtré créé avec succès!")
        print(f"Nombre de gènes trouvés: {gene_count}")
        print(f"Fichier de sortie: {filtered_file}")

        return filtered_file

    except FileNotFoundError:
        print(f"Erreur: Le fichier {gff_file} n'existe pas")
        return None
    except Exception as e:
        print(f"Erreur lors du traitement: {e}")
        return None


def split_gff_by_family(gff_file, output_dir=None):
    """
    Split un fichier GFF en fonction de la famille (colonne 9, attributes 'Fam=')
    Ne traite que les gènes (filtrés avec filter_gff_genes)

    Args:
        gff_file (str): Chemin vers le fichier GFF à traiter (fichier filtré des gènes)
        output_dir (str): Répertoire de sortie (par défaut: même répertoire que le fichier)

    Returns:
        dict: Dictionnaire avec les familles en clés et les chemins de fichiers en valeurs
    """
    import os
    from pathlib import Path

    # Déterminer le répertoire de sortie
    if output_dir is None:
        output_dir = os.path.dirname(gff_file) or "."
    else:
        os.makedirs(output_dir, exist_ok=True)

    # Extraire le nom du fichier d'entrée sans extension
    input_basename = os.path.splitext(os.path.basename(gff_file))[0]

    # Dictionnaire pour stocker les fichiers ouverts par famille
    family_files = {}
    family_counts = {}

    try:
        with open(gff_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                # Ignorer les lignes vides et les commentaires
                if line.strip() == '' or line.startswith('#'):
                    continue

                # Split la ligne GFF
                columns = line.strip().split('\t')

                # Vérifier que la ligne a 9 colonnes
                if len(columns) < 9:
                    print(f"Attention: Ligne {line_num} n'a pas 9 colonnes, ignorée")
                    continue

                # Vérifier que c'est un gène (colonne 3, index 2)
                feature_type = columns[2]
                if feature_type != "gene":
                    continue

                # Extraire les attributes (colonne 9, index 8)
                attributes = columns[8]

                # Extraire la famille
                family = None
                for attr in attributes.split(';'):
                    if attr.startswith('Fam='):
                        family = attr.split('=', 1)[1].strip()
                        break

                # Si pas de famille trouvée, utiliser "Unknown"
                if family is None:
                    family = "Unknown"

                # Créer un fichier de sortie pour cette famille s'il n'existe pas
                if family not in family_files:
                    output_file = os.path.join(output_dir, f"{input_basename}_{family}.gff")
                    family_files[family] = open(output_file, 'w')
                    family_counts[family] = 0

                # Écrire la ligne dans le fichier de la famille
                family_files[family].write(line)
                family_counts[family] += 1

        # Fermer tous les fichiers
        for f in family_files.values():
            f.close()

        # Afficher un résumé
        print(f"\nFichier GFF split avec succès!")
        print(f"Nombre de familles trouvées: {len(family_files)}")
        print("\nRésumé par famille:")
        for family in sorted(family_counts.keys()):
            output_file = os.path.join(output_dir, f"{input_basename}_{family}.gff")
            print(f"  {family}: {family_counts[family]} lignes -> {output_file}")

        # Retourner un dictionnaire avec les chemins
        result = {family: os.path.join(output_dir, f"{input_basename}_{family}.gff")
                  for family in family_files.keys()}
        return result

    except FileNotFoundError:
        print(f"Erreur: Le fichier {gff_file} n'existe pas")
        return {}
    except Exception as e:
        print(f"Erreur lors du traitement: {e}")
        return {}


# Exemple d'utilisation
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python split_gff_by_family.py <fichier_gff> [répertoire_sortie]")
        print("\nExemple:")
        print("  python split_gff_by_family.py mon_fichier.gff")
        print("  python split_gff_by_family.py mon_fichier.gff ./resultats/")
        sys.exit(1)

    gff_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None

    # Filtrer d'abord les gènes
    filtered_file = filter_gff_genes(gff_file)

    # Puis split par famille uniquement sur les gènes filtrés
    if filtered_file:
        split_gff_by_family(filtered_file, output_dir)
