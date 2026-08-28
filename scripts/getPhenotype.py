'''python getPhenotype.py \
    --input_h5 /nas/users/duyen/pseudomonas/SNP_Clstr/test_SNP_pipeline/test_AB/snp_output/all.95.8.snp.h5 \
    --metadata ~/pseudomonas/PATRIC_genomes_AMR.txt \
    --antibiotic meropenem \
    --species "Acinetobacter baumannii" \
    --output_file /nas/users/duyen/pseudomonas/SNP_Clstr/test_SNP_pipeline/test_AB/snp_output/snp_{antibiotic}.h5'''
import os
import argparse
import tables
import pandas as pd


# ============================================================
# Default settings
# ============================================================

DEFAULT_H5 = "/nas/users/duyen/pseudomonas/SNP_Clstr/test_SNP_pipeline/snp_output/all_snp.h5"
DEFAULT_METADATA = "/nas/users/duyen/pseudomonas/PATRIC_genomes_AMR.txt"
DEFAULT_GENLIST = None
DEFAULT_SPECIES = "Pseudomonas aeruginosa"


# ============================================================
# 1. Read SNP HDF5
# ============================================================

def load_snp_h5(h5_path):
    print("\n==============================================")
    print("Reading SNP HDF5")
    print("==============================================")
    print(f"[INFO] HDF5: {h5_path}")

    if not os.path.exists(h5_path):
        raise FileNotFoundError(f"HDF5 file not found: {h5_path}")

    with tables.open_file(h5_path, mode="r") as h5:
        genomes = [x.decode() for x in h5.root.index[:]]
        print(f"[INFO] Genomes in HDF5: {len(genomes)}")

        cluster_dfs = []
        cluster_count = 0
        total_features = 0

        for node in h5.root.snp_matrix:
            cluster_name = node._v_name
            matrix = node[:]

            feature_node = h5.root.feature_names._f_get_child(cluster_name)
            feature_names = [x.decode() for x in feature_node[:]]

            if matrix.shape[0] != len(genomes):
                raise ValueError(f"Row mismatch in {cluster_name}: matrix has {matrix.shape[0]} rows, but HDF5 index has {len(genomes)} genomes.")

            if matrix.shape[1] != len(feature_names):
                raise ValueError(f"Feature mismatch in {cluster_name}: matrix has {matrix.shape[1]} columns, but there are {len(feature_names)} feature names.")

            cluster_df = pd.DataFrame(matrix, index=genomes, columns=feature_names)
            cluster_dfs.append(cluster_df)

            cluster_count += 1
            total_features += len(feature_names)

            print(f"[INFO] {cluster_name}: {matrix.shape[0]} × {matrix.shape[1]}")

    if not cluster_dfs:
        raise ValueError("No SNP clusters found in HDF5.")

    print("\n[INFO] Combining clusters...")

    df = pd.concat(cluster_dfs, axis=1)
    del cluster_dfs

    print(f"[INFO] Number of clusters: {cluster_count}")
    print(f"[INFO] Total SNP features: {total_features}")
    print(f"[INFO] Final SNP matrix: {df.shape[0]} genomes × {df.shape[1]} features")

    return df


# ============================================================
# 2. Read PATRIC phenotype information
# ============================================================

def load_antibiotic_phenotypes(metadata_path, antibiotic, genlist_path, species):
    print("\n==============================================")
    print("Reading phenotype information")
    print("==============================================")
    print(f"[INFO] Metadata file: {metadata_path}")
    print(f"[INFO] Species: {species}")
    print(f"[INFO] Antibiotic: {antibiotic}")
    print(f"[INFO] Genlist: {genlist_path}")

    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    
    # --------------------------------------------------------
    # Read phenotype table
    # --------------------------------------------------------

    cols = ["genome_id", "genome_name", "antibiotic", "resistant_phenotype", "laboratory_typing_method"]

    df = pd.read_csv(metadata_path, sep="\t", usecols=cols, dtype=str)

    # --------------------------------------------------------
    # Clean columns
    # --------------------------------------------------------

    for col in cols:
        df[col] = df[col].fillna("").str.strip()

    # --------------------------------------------------------
    # Keep requested species
    # --------------------------------------------------------

    df = df[df["genome_name"].str.contains(species, case=False, na=False)].copy()

    print(f"[INFO] Records for species: {len(df)}")

    # --------------------------------------------------------
    # Keep requested antibiotic
    # --------------------------------------------------------

    df = df[df["antibiotic"].str.lower() == antibiotic.lower()].copy()

    print(f"[INFO] Records for {antibiotic}: {len(df)}")

    # --------------------------------------------------------
    # Keep only Resistant / Susceptible
    # --------------------------------------------------------

    df["phenotype_lower"] = df["resistant_phenotype"].str.lower()

    df = df[df["phenotype_lower"].isin(["resistant", "susceptible"])].copy()

    print(f"[INFO] Valid R/S records: {len(df)}")

    # --------------------------------------------------------
    # Convert phenotype
    # Resistant = 1
    # Susceptible = 0
    # --------------------------------------------------------

    df["resistant_phenotype"] = df["phenotype_lower"].map({"resistant": 1, "susceptible": 0})

    # --------------------------------------------------------
    # Normalize genome IDs
    # --------------------------------------------------------

    df["genome_id"] = df["genome_id"].astype(str).str.strip()
    df["genome_id"] = df["genome_id"].str.replace(".fna", "", regex=False)

    # --------------------------------------------------------
    # Optional genlist filtering
    # --------------------------------------------------------

    if genlist_path is not None:

        if not os.path.exists(genlist_path):
            raise FileNotFoundError(
                f"Genlist file not found: {genlist_path}"
            )

        with open(genlist_path, "r", encoding="utf-8") as f:
            gen_ids = [
                line.strip().replace(".fna", "")
                for line in f
                if line.strip()
            ]

        gen_set = set(gen_ids)

        print(f"[INFO] Genomes in genlist: {len(gen_set)}")

        # Keep only genomes present in genlist
        df = df[df["genome_id"].isin(gen_set)].copy()

        print(
            f"[INFO] Phenotype records after genlist: {len(df)}"
        )

    else:

        print("[INFO] No genlist provided.")
        print("[INFO] Using all phenotype records.")

    # --------------------------------------------------------
    # Select preferred phenotype
    #
    # MIC > disk_diffusion > first available
    # --------------------------------------------------------

    def pick_preferred(row_group):
        methods = row_group["laboratory_typing_method"].str.lower().str.strip()

        mic_row = row_group[methods == "mic"]

        if not mic_row.empty:
            return mic_row.iloc[0]

        disk_row = row_group[methods == "disk_diffusion"]

        if not disk_row.empty:
            return disk_row.iloc[0]

        return row_group.iloc[0]

    # --------------------------------------------------------
    # One phenotype per genome
    # --------------------------------------------------------

    df_unique = df.groupby("genome_id", group_keys=False).apply(pick_preferred)

    # --------------------------------------------------------
    # Create phenotype dictionary
    # --------------------------------------------------------

    phenotype = df_unique["resistant_phenotype"].astype(int).to_dict()

    resistant = sum(value == 1 for value in phenotype.values())
    susceptible = sum(value == 0 for value in phenotype.values())

    # --------------------------------------------------------
    # Count selected laboratory methods
    # --------------------------------------------------------

    selected_methods = df_unique["laboratory_typing_method"].str.lower().str.strip()

    mic_count = (selected_methods == "mic").sum()
    disk_count = (selected_methods == "disk_diffusion").sum()
    other_count = len(df_unique) - mic_count - disk_count

    print(f"[INFO] Phenotyped genomes: {len(phenotype)}")
    print(f"[INFO] Resistant: {resistant}")
    print(f"[INFO] Susceptible: {susceptible}")
    print(f"[INFO] Selected MIC records: {mic_count}")
    print(f"[INFO] Selected disk_diffusion records: {disk_count}")
    print(f"[INFO] Selected other/fallback records: {other_count}")

    return phenotype


# ============================================================
# 3. Create antibiotic-specific SNP DataFrame
# ============================================================

def create_antibiotic_dataframe(snp_df, phenotype):
    print("\n==============================================")
    print("Selecting genomes")
    print("==============================================")

    selected_genomes = [genome for genome in snp_df.index if genome in phenotype]

    print(f"[INFO] SNP genomes: {len(snp_df)}")
    print(f"[INFO] Phenotyped genomes: {len(phenotype)}")
    print(f"[INFO] Genomes in both: {len(selected_genomes)}")

    if not selected_genomes:
        raise ValueError("No genomes are shared between SNP HDF5 and phenotype table.")

    # --------------------------------------------------------
    # Select SNP rows
    # --------------------------------------------------------

    drug_df = snp_df.loc[selected_genomes].copy()

    # --------------------------------------------------------
    # Add phenotype
    # --------------------------------------------------------

    drug_df.insert(0, "resistant_phenotype", [phenotype[genome] for genome in drug_df.index])

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    resistant = (drug_df["resistant_phenotype"] == 1).sum()
    susceptible = (drug_df["resistant_phenotype"] == 0).sum()

    print(f"[INFO] Final matrix: {drug_df.shape[0]} genomes × {drug_df.shape[1]} columns")
    print(f"[INFO] Resistant: {resistant}")
    print(f"[INFO] Susceptible: {susceptible}")

    return drug_df


# ============================================================
# 4. Save antibiotic DataFrame
# ============================================================

def save_antibiotic_dataframe(df, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    print("\n==============================================")
    print("Saving HDF5")
    print("==============================================")
    print(f"[INFO] Output: {output_path}")

    df.to_hdf(output_path, key="snp_matrix", mode="w", format="fixed")

    print("[DONE] HDF5 saved.")


# ============================================================
# 5. Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Extract antibiotic-specific SNP matrices from an SNP HDF5 file using PATRIC antibiotic phenotype data.")

    parser.add_argument("--input_h5", default=DEFAULT_H5, help="Input SNP HDF5 file created by the SNP extraction pipeline.")
    parser.add_argument("--metadata", default=DEFAULT_METADATA, help="PATRIC_genomes_AMR.txt file.")
    parser.add_argument("--genlist", default=None, help="Optional genome list containing genomes to keep.")
    parser.add_argument("--species", default=DEFAULT_SPECIES, help="Species name to keep, e.g. 'Pseudomonas aeruginosa'.")
    parser.add_argument("--antibiotic", required=True, help="Antibiotic to extract, e.g. meropenem.")
    parser.add_argument("--output_file", required=True, help="Full path of the output HDF5 files.")

    args = parser.parse_args()

    # ========================================================
    # 1. Load complete SNP matrix
    # ========================================================

    snp_df = load_snp_h5(args.input_h5)

    # ========================================================
    # 2. Load phenotype
    # ========================================================

    phenotype = load_antibiotic_phenotypes(args.metadata, args.antibiotic, args.genlist, args.species)

    # ========================================================
    # 3. Select antibiotic-specific genomes
    # ========================================================

    drug_df = create_antibiotic_dataframe(snp_df, phenotype)

    # ========================================================
    # 4. Output file
    # ========================================================

    output_path = args.output_file.replace(
        "{antibiotic}",
        args.antibiotic.replace("/", "_").replace(" ", "_")
    )

    # ========================================================
    # 5. Save
    # ========================================================

    save_antibiotic_dataframe(drug_df, output_path)

    # ========================================================
    # 6. Final summary
    # ========================================================

    print("\n==============================================")
    print("[DONE] Antibiotic SNP extraction complete.")
    print(f"[INFO] Species: {args.species}")
    print(f"[INFO] Antibiotic: {args.antibiotic}")
    print(f"[INFO] Final shape: {drug_df.shape}")
    print(f"[INFO] Output: {output_path}")
    print("==============================================")


# ============================================================
# Run
# ============================================================

if __name__ == "__main__":
    main()
