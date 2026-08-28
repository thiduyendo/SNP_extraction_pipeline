import argparse
import os
import subprocess


def read_fasta(fasta_file):
    """Read FASTA sequences into a dictionary."""

    sequences = {}
    current_id = None  #1402491.3|fig|1402491.3.peg.1991|Q004_01960
    current_header = None #>1402491.3|fig|1402491.3.peg.1991|Q004_01960 Protein Annotation
    sequence = []

    with open(fasta_file, "r", encoding="utf-8", errors="ignore") as infile:

        for line in infile:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if current_id is not None:
                    sequences[current_id] = (
                        current_header,
                        "".join(sequence)
                    )

                current_header = line[1:]
                current_id = current_header.split()[0]
                sequence = []

            else:

                sequence.append(line)

        if current_id is not None:
            sequences[current_id] = (
                current_header,
                "".join(sequence)
            )
    return sequences


def read_clusters(clstr_file):
    """Read CD-HIT clusters."""

    clusters = []
    current_cluster = []

    with open(clstr_file, "r", encoding="utf-8", errors="ignore") as infile:

        for line in infile:
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

                seq_id = line.split(">", 1)[1] #287.8491|fig|287.8491.peg.2059...*

                if "..." in seq_id:
                    seq_id = seq_id.split("...", 1)[0] #287.8491|fig|287.8491.peg.2059"

                seq_id = seq_id.split()[0]

                current_cluster.append(seq_id)

        if current_cluster:
            clusters.append(current_cluster)

    return clusters


def run_mafft(cluster_ids, sequences, mafft):

    fasta_lines = []
    valid_sequences = 0

    for seq_id in cluster_ids:

        if seq_id not in sequences:
            print(f"[WARNING] Sequence not found: {seq_id}")
            continue

        header, sequence = sequences[seq_id]
        fasta_lines.append(f">{header}")
        fasta_lines.append(sequence)
        valid_sequences += 1

    # MAFFT requires at least 2 valid sequences
    if valid_sequences < 2:
        return None

    fasta_text = "\n".join(fasta_lines) + "\n"

    cmd = [mafft, "--auto", "--anysymbol", "-"]

    result = subprocess.run(
        cmd,
        input=fasta_text,
        text=True,
        capture_output=True,
        check=True
    )

    return result.stdout

def main():

    parser = argparse.ArgumentParser(
        description="Run MAFFT on CD-HIT clusters and create one aligned FASTA file."
    )

    parser.add_argument(
        "--fasta",
        default="all.95.0.faa",
        help="Protein faa filename"
    )

    parser.add_argument(
        "--clstr",
        default="all.95.0.cdhit.clstr",
        help="CD-HIT cluster file"
    )

    parser.add_argument(
        "--output",
        default="all.95.0.aligned.faa",
        help="Output aligned FASTA filename"
    )

    parser.add_argument(
        "--mafft",
        default="mafft",
        help="MAFFT executable"
    )

    args = parser.parse_args()

    # ========================================================
    # Input files
    # ========================================================

    fasta = args.fasta
    clstr_file = args.clstr
    output_file = args.output

    # ========================================================
    # Check files
    # ========================================================

    if not os.path.exists(fasta):

        raise FileNotFoundError(
            f"Protein faa FASTA not found: {fasta}"
        )

    if not os.path.exists(clstr_file):

        raise FileNotFoundError(
            f"CD-HIT cluster file not found: {clstr_file}"
        )

    # ========================================================
    # Read Protein faa FASTA
    # ========================================================

    print("\n========================================")
    print("Reading Protein faa FASTA")
    print("========================================")

    sequences = read_fasta(fasta)

    print(
        f"Protein sequences: {len(sequences)}"
    )

    # ========================================================
    # Read clusters
    # ========================================================

    print("\n========================================")
    print("Reading CD-HIT clusters")
    print("========================================")

    clusters = read_clusters(clstr_file)

    print(f"Clusters: {len(clusters)}")

    # ========================================================
    # Remove old output
    # ========================================================

    if os.path.exists(output_file):

        print(
            f"\nRemoving existing output: {output_file}")

        os.remove(output_file)

    # ========================================================
    # Process clusters
    # ========================================================

    print("\n========================================")
    print("Running MAFFT")
    print("========================================")

    processed = 0
    skipped = 0

    with open(output_file, "w", encoding="utf-8") as outfile:

        for cluster_number, cluster_ids in enumerate(clusters):

            print(
                f"\nCluster {cluster_number}: "
                f"{len(cluster_ids)} sequences"
            )

            # ------------------------------------------------
            # MAFFT needs at least two sequences
            # ------------------------------------------------

            if len(cluster_ids) < 2:

                skipped += 1

                print(
                    "  Skipped: only one sequence"
                )

                continue

            try:

                aligned = run_mafft(
                    cluster_ids,
                    sequences,
                    args.mafft
                )

                if aligned is None:

                    skipped += 1

                    print(
                        "  Skipped: MAFFT returned no alignment"
                    )

                    continue

                outfile.write(
                    aligned
                )

                outfile.flush()

                processed += 1

                print("  Alignment written")

            except subprocess.CalledProcessError as e:

                skipped += 1

                print(
                    f"  [MAFFT ERROR] Cluster {cluster_number}"
                )

                if e.stderr:
                    print(e.stderr)

    # ========================================================
    # Summary
    # ========================================================

    print("\n========================================")
    print("MAFFT completed!")
    print("========================================")

    print(
        f"Input FASTA      : {fasta}"
    )

    print(
        f"Cluster file     : {clstr_file}"
    )

    print(
        f"Output alignment : {output_file}"
    )

    print(
        f"Total clusters   : {len(clusters)}"
    )

    print(
        f"Processed        : {processed}"
    )

    print(
        f"Skipped          : {skipped}"
    )

    print("========================================")


if __name__ == "__main__":
    main()
