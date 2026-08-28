###python
###This is the 7 version of extract_snp.py
###Run with only one input: .cdhit.clstr, .aligned.faa (Mafft aligned file), and all.ffn (all nucleotide seqs
###python
###This is the 7 version of extract_snp.py
###Run with only one input: .cdhit.clstr, .aligned.faa (Mafft aligned file), and all.ffn (all nucleotide seqs
###python
###This is the 7 version of extract_snp.py
###Run with only one input: .cdhit.clstr, .aligned.faa (Mafft aligned file), and all.ffn (all nucleotide seqs
import os, tables
import re
import glob
import argparse
from pathlib import Path
import pandas as pd
from tqdm import tqdm
import numpy as np

# ============================================================
# Default settings
# ============================================================

DEFAULT_INPUT_DIR = "/nas/users/duyen/pseudomonas/SNP_Clstr/test_SNP_pipeline/CDHIT"
DEFAULT_FFN_FILE = "/nas/users/duyen/pseudomonas/ffn/all.ffn"
DEFAULT_GENLIST = None

# ============================================================
# 1. Parse FASTA
# ============================================================

def parse_aligned_fasta(filepath):
    """
    Parse the combined MAFFT-aligned FASTA file.

    Example header:

        >287.999|fig|287.999.peg.1|Q003_00001
        MKT---A...

    Returns:

        {
            "287.999|fig|287.999.peg.1|Q003_00001":
                "MKT---A...",
            ...
        }
    """

    sequences = {}
    current_id = None
    seq_parts = []

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:

        for raw_line in f:

            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):

                # Save previous sequence
                if current_id is not None:
                    sequences[current_id] = "".join(seq_parts)

                header = line[1:].strip()

                current_id = header.split()[0]

                seq_parts = []

            else:
                seq_parts.append(line)

        if current_id is not None:
            sequences[current_id] = "".join(seq_parts)

    return sequences


# ============================================================
# 2. Parse protein header
# ============================================================

def parse_protein_header(header):

    """
    Handles both structural types:
    Type A (Old): >470.7359|fig|470.7359.peg.3568|RQ11_12830
    Type B (New): >470.7359__fig|470.7359.peg.3568|RQ11_12830

    Example protein header:

    >287.999|fig|287.999.peg.1|Q003_00001 Methyl-accepting chemotaxis sensor ...

    Returns:

        genome_id = 287.999
        full_id    = fig|287.999.peg.1|Q003_00001
    """
    header = header[1:].strip()
    first_field = header.split()[0]

    # Rule 1: Check for the double underscore format first
    if "__" in first_field:
        parts = first_field.split("__", 1)
        genome_id = parts[0]
        full_id = parts[1]
        return genome_id, full_id

    # Rule 2: Fallback to standard pipe delimiter mapping
    elif "|" in first_field:
        parts = first_field.split("|")
        if len(parts) >= 3:
            genome_id = parts[0]
            full_id = "|".join(parts[1:])
            return genome_id, full_id

    # Fallback default if completely unparseable
    print(f"[WARNING] Cannot parse header structure: {first_field}")
    return None, None

# ============================================================
# 3. Load combined FFN
# ============================================================

def load_ffn(ffn_path):

    sequences = {}
    current_id = None
    seq_parts = []

    with open(ffn_path, "r", encoding="utf-8", errors="ignore") as f:

        for raw_line in f:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(seq_parts)

                _, full_id = parse_protein_header(line)
                current_id = full_id
                seq_parts = []

            else:
                seq_parts.append(line)

        if current_id is not None:
            sequences[current_id] = "".join(seq_parts)

    return sequences


# ============================================================
# 5. Convert protein alignment to nucleotide alignment
# ============================================================

def protein_alignment_to_nt(protein_alignment, nucleotide_sequence):

    nucleotide_sequence = nucleotide_sequence.upper()

    protein_ungapped = protein_alignment.replace("-", "")

    expected_nt_length = len(protein_ungapped) * 3

    if (
        len(nucleotide_sequence) == expected_nt_length + 3
        and nucleotide_sequence[-3:] in {"TAA", "TAG", "TGA"}
    ):

        nucleotide_sequence = nucleotide_sequence[:-3]

    if len(nucleotide_sequence) % 3 != 0:

        print(
            f"        [MAP ERROR] CDS length is not divisible by 3: "
            f"{len(nucleotide_sequence)}"
        )

        return None

    actual_nt_length = len(nucleotide_sequence)

    if expected_nt_length != actual_nt_length:

        print(
            f"        [MAP ERROR] Protein/CDS length mismatch:"
        )

        print(
            f"            Protein aligned length: "
            f"{len(protein_alignment)}"
        )

        print(
            f"            Protein ungapped length: "
            f"{len(protein_ungapped)}"
        )

        print(
            f"            Expected CDS length: "
            f"{expected_nt_length}"
        )

        print(
            f"            Actual CDS length: "
            f"{actual_nt_length}"
        )

        print(
            f"            Terminal codon: "
            f"{nucleotide_sequence[-3:]}"
        )

        return None

    result = []

    nt_index = 0

    for amino_acid in protein_alignment:

        if amino_acid == "-":

            result.extend(["-", "-", "-"])

        else:

            codon = nucleotide_sequence[nt_index:nt_index + 3]
            result.extend(list(codon))
            nt_index += 3

    if nt_index != len(nucleotide_sequence):

        print("[MAP ERROR] Not all nucleotide sequence was consumed.")

        print(f"Consumed: {nt_index}")

        print(f"Available: {len(nucleotide_sequence)}")

        return None

    return result


# ============================================================
# 6. Encode nucleotide
# ============================================================

def encode_base(base):

    base = base.upper()

    if base == "A":

        return 2

    elif base == "T":

        return 3

    elif base == "C":

        return 5

    elif base == "G":

        return 4

    else:

        return 0


# ============================================================
# 7. Parse CD-HIT clusters
# ============================================================
''' Nested list, each list is a cluster
clusters = [
    [
        "287.999|fig|287.999.peg.1|Q003_00001",
        "288.100|fig|288.100.peg.5|Q003_00002",
        "300.200|fig|300.200.peg.10|Q003_00003"
    ],
    [
        "287.999|fig|287.999.peg.50|Q003_00050",
        "288.100|fig|288.100.peg.60|Q003_00060"
    ]
]
]'''
def parse_cdhit_clusters(clstr_file):

    clusters = []
    current_cluster = []

    with open(clstr_file, "r", encoding="utf-8", errors="ignore") as f:

        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">Cluster"):

                if current_cluster:
                    clusters.append(current_cluster)

                current_cluster = []

            else:

                if ">" not in line:
                    continue

                seq_id = line.split(">", 1)[1]

                if "..." in seq_id:
                    seq_id = seq_id.split("...", 1)[0]

                seq_id = seq_id.split()[0]
                current_cluster.append(seq_id)

        if current_cluster:
            clusters.append(current_cluster)

    return clusters


# ============================================================
# 8. Process one aligned CD-HIT cluster
# ============================================================

def process_one_cluster(
    cluster_number,
    cluster_ids, #cluster ids in each cluster in clusters list (.cdhit.clstr)
    protein_dict, # {"protein_id": "Protein sequence", ...} from Mafft aligned.faa
    ffn_dict, #{"fig|287.999.peg.1|Q003_00001": "ATGAAGACT...",...}
    genome_ids
):

    nucleotide_alignments = {}

    for sequence_id in cluster_ids:

        if sequence_id not in protein_dict:
            tqdm.write(
                f"[WARNING] Aligned sequence not found: {sequence_id}"
            )
            continue

        protein_sequence = protein_dict[sequence_id]
        genome_id, full_id = parse_protein_header(">" + sequence_id)

        if genome_id is None or full_id is None:
            print(f"[WARNING] Cannot parse header: {sequence_id}")
            continue

        if genome_ids is not None:
            if genome_id not in genome_ids:
                continue

        if full_id not in ffn_dict:
            print(f"[WARNING] {full_id} not found in all.ffn")
            continue

        nucleotide_sequence = ffn_dict[full_id]

        nt_alignment = protein_alignment_to_nt(
            protein_sequence,
            nucleotide_sequence
            )

        if nt_alignment is None:
            print(f"[WARNING] Cannot map nucleotide sequence for {sequence_id}")

            print(f"    Protein alignment length: {len(protein_sequence)}")

            print(
                f"    Protein ungapped length: "
                f"{len(protein_sequence.replace('-', ''))}")

            print(
                f"    Expected CDS length: "
                f"{len(protein_sequence.replace('-', '')) * 3}")

            print(f"    Actual CDS length: {len(nucleotide_sequence)}")

            print(f"    CDS length % 3: {len(nucleotide_sequence) % 3}")

            print(f"    Protein ID: {full_id}")

            continue

        nucleotide_alignments[sequence_id] = nt_alignment

    if len(nucleotide_alignments) < 2:
        return None

    alignment_lengths = {
        len(seq)
        for seq in nucleotide_alignments.values()
    }

    if len(alignment_lengths) != 1:
        print(
            f"[WARNING] Different alignment lengths "
            f"in Cluster {cluster_number}"
        )
        return None

    alignment_length = next(iter(alignment_lengths))
    feature_data = {}
    cluster_name = f"cluster{cluster_number}"

    for position in range(alignment_length):

        column = []

        for sequence in nucleotide_alignments.values():
            base = sequence[position]
            column.append(base.upper())

        # ----------------------------------------------------
        # Same SNP definition as your current script
        #
        # Ignore gaps when determining variability.
        # Keep bi-/tri-/multi-allelic sites.
        # ----------------------------------------------------

        valid_bases = {"A", "T", "C", "G"}

        valid_column = [
            base
            for base in column
            if base in valid_bases
        ]

        if len(set(valid_column)) < 2:
            continue

        feature_name = (f"{cluster_name}_{position}")

        if feature_name not in feature_data:
            feature_data[feature_name] = {}

        for protein_id, sequence in nucleotide_alignments.items():
            #genome_id = protein_id.split("|")[0]
            #feature_data[feature_name][genome_id] = encode_base(sequence[position])

        
            genome_id, full_id = parse_protein_header(
                ">" + protein_id
            )

            if genome_id is None:
                continue

            feature_data[feature_name][genome_id] = encode_base(
                sequence[position]
            )


    if not feature_data:
        return None

    return feature_data


# ============================================================
# 9. Load genome list
# ============================================================

def load_genome_list(genlist_path):

    genome_ids = set()

    with open(
        genlist_path, "r", encoding="utf-8") as f:

        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.lower().endswith(".fna"):
                genome_id = line[:-4]

            elif line.lower().endswith(".faa"):
                genome_id = line[:-4]

            else:
                genome_id = line

            genome_ids.add(genome_id)

    return genome_ids


# ============================================================
# 10. Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Extract nucleotide SNP features from one "
            "combined MAFFT alignment using CD-HIT clusters."
        )
    )

    parser.add_argument(
        "--aligned_fasta",
        default="all.95.0.aligned.faa",
        help="Combined MAFFT alignment"
    )

    parser.add_argument(
        "--clstr",
        default="all.95.0.cdhit.clstr",
        help="CD-HIT cluster file"
    )

    parser.add_argument(
        "--ffn",
        default=DEFAULT_FFN_FILE,
        help="Combined all.ffn file"
    )

    parser.add_argument(
        "--genlist",
        default=None,
        help="Optional genome list"
    )

    parser.add_argument(
        "--output_file",
        default="all_snp.h5",
        help="Full path and filename for the output HDF5 file"
    )

    args = parser.parse_args()

    # ========================================================
    # 1. Genome list
    # ========================================================

    genome_ids = None

    if args.genlist is not None:

        if not os.path.exists(args.genlist):

            print(
                f"[ERROR] Genome list not found: "
                f"{args.genlist}"
            )

            return

        genome_ids = load_genome_list(args.genlist)

        print("[INFO] Genome filtering enabled.")

        print(f"[INFO] Genomes in genlist: {len(genome_ids)}")

    else:

        print("[INFO] Genome filtering disabled.")

    # ========================================================
    # 2. Create output directory
    # ========================================================

    output_dir = os.path.dirname(args.output_file)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # ========================================================
    # 3. Input files
    # ========================================================

    aligned_fasta = args.aligned_fasta

    clstr_file = args.clstr

    ffn_file = args.ffn

    if not os.path.exists(aligned_fasta):
        print(
            f"[ERROR] Alignment file not found: "
            f"{aligned_fasta}"
        )

        return

    if not os.path.exists(clstr_file):
        print(
            f"[ERROR] CD-HIT cluster file not found: "
            f"{clstr_file}"
        )

        return

    if not os.path.exists(ffn_file):
        print(f"[ERROR] Combined FFN not found: {ffn_file}")

        return

    # ========================================================
    # 4. Read combined aligned FASTA
    # ========================================================

    print("\n========================================")

    print("Reading combined MAFFT alignment")

    print("========================================")

    protein_dict = parse_aligned_fasta(aligned_fasta)

    print(f"[INFO] Aligned proteins: {len(protein_dict)}")

    # ========================================================
    # 5. Read CD-HIT clusters
    # ========================================================

    print("\n========================================")

    print("Reading CD-HIT clusters")

    print("========================================")

    clusters = parse_cdhit_clusters(clstr_file)

    print(
        f"[INFO] Total clusters: "
        f"{len(clusters)}"
    )

    # ========================================================
    # 6. Load combined FFN
    # ========================================================

    print("\n========================================")

    print("Loading combined FFN")

    print("========================================")

    ffn_cache = {}

    if "all" not in ffn_cache:
        ffn_cache["all"] = load_ffn(ffn_file)

    ffn_dict = ffn_cache["all"]

    print(f"[INFO] FFN sequences: {len(ffn_dict)}")

    # ========================================================
    # 7. Collect genome IDs
    # ========================================================

    all_sequence_ids = set()

    print("\n[INFO] Collecting genome IDs...")

    for sequence_id in protein_dict:

        genome_id, full_id = parse_protein_header(">" + sequence_id)

        if genome_id is None:
            continue

        if genome_ids is not None:
            if genome_id not in genome_ids:
                continue

        all_sequence_ids.add(genome_id)

    all_sequence_ids = sorted(all_sequence_ids)

    if not all_sequence_ids:
        print("[ERROR] No genome IDs found.")

        return

    print(f"[INFO] Total genomes: {len(all_sequence_ids)}")

    # ========================================================
    # 8. Create HDF5
    # ========================================================

    output_path = args.output_file

    if os.path.exists(output_path):
        print(f"[INFO] Removing existing HDF5: {output_path}")
        os.remove(output_path)

    with tables.open_file(output_path, mode="w") as h5:

        h5.create_array(
            "/",
            "index",
            obj=np.array(
                all_sequence_ids,
                dtype="S"
            )
        )

        matrix_group = h5.create_group(
            "/",
            "snp_matrix"
        )

        feature_group = h5.create_group(
            "/",
            "feature_names"
        )

        total_features = 0
        processed_clusters = 0
        failed_clusters = 0

        # ====================================================
        # 9. Process each cluster
        # ====================================================

        for cluster_number, cluster_ids in enumerate(
            tqdm(
                clusters,
                desc="Processing clusters"
            )
        ):

            try:

                feature_data = process_one_cluster(
                    cluster_number,
                    cluster_ids,
                    protein_dict,
                    ffn_dict,
                    genome_ids
                )

                if not feature_data:
                    continue

                cluster_df = pd.DataFrame.from_dict(
                    feature_data,
                    orient="columns"
                )

                if cluster_df.empty:

                    del feature_data
                    del cluster_df

                    continue

                cluster_df = cluster_df.fillna(0)

                cluster_df = cluster_df.astype("uint8")

                cluster_df = cluster_df.reindex(
                    all_sequence_ids,
                    fill_value=0
                )

                cluster_array = cluster_df.to_numpy(dtype=np.uint8)

                n_features = (cluster_array.shape[1])

                cluster_name = (f"cluster{cluster_number}")

                h5.create_array(
                    matrix_group,
                    cluster_name,
                    obj=cluster_array
                )

                h5.create_array(
                    feature_group,
                    cluster_name,
                    obj=np.array(
                        cluster_df.columns.astype(str),
                        dtype="S"
                    )
                )

                processed_clusters += 1
                total_features += n_features

                print(f"\n[SAVED] {cluster_name}")

                print(f"        Rows: {cluster_array.shape[0]}")

                print(f"        Features: {cluster_array.shape[1]}")

                del cluster_array
                del cluster_df
                del feature_data

            except Exception as e:
                failed_clusters += 1

                print(f"\n[ERROR] Cluster {cluster_number}: {e}")

    # ========================================================
    # 10. Summary
    # ========================================================

    print(
        "\n=============================================="
    )

    print(
        "[DONE] Final HDF5 created:"
    )

    print(
        f"       {output_path}"
    )

    print(
        f"[INFO] Number of genomes: "
        f"{len(all_sequence_ids)}"
    )

    print(
        f"[INFO] Number of SNP features: "
        f"{total_features}"
    )

    print(
        f"[INFO] Number of processed clusters: "
        f"{processed_clusters}"
    )

    print(
        f"[INFO] Number of failed clusters: "
        f"{failed_clusters}"
    )

    print(
        f"[INFO] FFN sequences loaded: "
        f"{len(ffn_dict)}"
    )

    print(
        "=============================================="
    )


if __name__ == "__main__":
    main()
