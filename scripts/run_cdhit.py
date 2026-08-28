import argparse
import glob
import os
import subprocess

def read_genlist(path):
    """Read genome IDs from genlist.txt."""
    genome_ids = set()

    with open(path) as f:
        for line in f:
            genome_id = line.strip()

            if not genome_id:
                continue

            # Remove .faa if present
            if genome_id.endswith(".fna"):
                genome_id = genome_id[:-4]
            genome_ids.add(genome_id)
    return genome_ids

def combine_fasta(faa_files, output_file):
    """Combine FASTA files and prepend genome ID to each header."""
    with open(output_file, "w") as out:
        for faa in faa_files:
            filename = os.path.basename(faa)
            if ".PATRIC.faa" in filename:
                genome_id = filename.replace(".PATRIC.faa", "")
            else:
                genome_id = filename.replace(".faa", "")

            print(f"[INFO] Processing: {filename}")
            print(f"       Genome ID: {genome_id}")
            
            with open(faa) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith(">"):
                        line = f">{genome_id}|{line[1:]}"

                    out.write(line + "\n")


def count_sequences(fasta_file):
    """Count FASTA sequences."""
    with open(fasta_file) as f:
        return sum(line.startswith(">") for line in f)


def main():
    parser = argparse.ArgumentParser(
        description="Combine PATRIC FASTA files and run CD-HIT."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Directory containing *.PATRIC.faa files"
    )

    parser.add_argument(
        "--output",
        default="./CDHIT",
        help="Output directory (default: ./CDHIT)"
    )

    parser.add_argument(
        "--genlist",
        default=None,
        help="Optional file containing genome IDs to include"
    )

    parser.add_argument(
        "--identity",
        type=float,
        default=0.95,
        help="CD-HIT identity threshold (default: 0.95)"
    )

    parser.add_argument(
        "--memory",
        type=int,
        default=10000,
        help="CD-HIT memory in MB (default: 10000)"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="CD-HIT threads (default: 8)"
    )
    
    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix for output files, e.g. all.95.8"
    )

    parser.add_argument(
        "--cd-hit",
        default="cd-hit",
        help="CD-HIT executable (default: cd-hit)"
    )

    args = parser.parse_args()

    # ========================================================
    # Find FASTA files
    # ========================================================

    faa_files = sorted(
        glob.glob(
            os.path.join(args.input, "*.faa")
        )
    )

    if not faa_files:
        raise FileNotFoundError(
            f"No *.faa files found in {args.input}"
        )

    # ========================================================
    # Optional genome filtering
    # ========================================================

    if args.genlist:

        genome_ids = read_genlist(args.genlist)

        faa_files = [
            faa for faa in faa_files
            if os.path.basename(faa).replace(
                ".PATRIC.faa", ""
            ) in genome_ids
        ]

        if not faa_files:
            raise ValueError(
                "No FASTA files matched the genome IDs in genlist."
            )

    print(f"FASTA files selected: {len(faa_files)}")

    # ========================================================
    # Output files
    # ========================================================

    os.makedirs(args.output, exist_ok=True)

    prefix = args.prefix

    combined_fasta = os.path.join(
        args.output,
        f"{prefix}.faa"
    )

    cdhit_prefix = os.path.join(
        args.output,
        f"{prefix}.cdhit"
    )

    cdhit_output = os.path.join(
        args.output,
        f"{prefix}.cdhit.faa"
    )

    clstr_file = os.path.join(
        args.output,
        f"{prefix}.cdhit.clstr"
    )

    # ========================================================
    # Combine FASTA
    # ========================================================

    print("\nCombining FASTA files...")

    combine_fasta(
        faa_files,
        combined_fasta
    )

    n_sequences = count_sequences(combined_fasta)

    print(f"Protein sequences: {n_sequences}")
    print(f"Combined FASTA: {combined_fasta}")

    # ========================================================
    # Run CD-HIT
    # ========================================================

    cmd = [
        args.cd_hit,
        "-i", combined_fasta,
        "-o", cdhit_prefix,
        "-c", str(args.identity),
        "-d", "0",
        "-M", str(args.memory),
        "-T", str(args.threads),
    ]

    print("\nRunning CD-HIT:")
    print(" ".join(cmd))
    print()

    subprocess.run(cmd, check=True)

    # ========================================================
    # Check output
    # ========================================================

    if not os.path.exists(cdhit_prefix):
        raise FileNotFoundError(
            f"CD-HIT output not found: {cdhit_prefix}"
        )

    if not os.path.exists(clstr_file):
        raise FileNotFoundError(
            f"CD-HIT cluster file not found: {clstr_file}"
        )

    # Rename CD-HIT output to the final .faa filename
    os.rename(cdhit_prefix, cdhit_output)

    print("\n========================================")
    print("Completed successfully!")
    print("========================================")
    print(f"Genomes          : {len(faa_files)}")
    print(f"Proteins         : {n_sequences}")
    print(f"Combined FASTA   : {combined_fasta}")
    print(f"CD-HIT FASTA     : {cdhit_output}")
    print(f"Cluster file     : {clstr_file}")
    print("Cluster FASTAs   : not created")


if __name__ == "__main__":
    main()
