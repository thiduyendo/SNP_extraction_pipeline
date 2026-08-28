import os
import glob
import argparse
import re


def find_ffn_files(input_dir):
    """
    Find FFN files in the input directory.

    Supports:
        *.PATRIC.ffn
        *.ffn
    """

    ffn_files = sorted(
        glob.glob(
            os.path.join(input_dir, "*.ffn")
        )
    )

    if not ffn_files:
        ffn_files = sorted(
            glob.glob(
                os.path.join(input_dir, "*.ffn")
            )
        )

    return ffn_files


def add_genome_id_to_header(header):
    """
    Convert:

        >fig|1402486.3.peg.1|P999_00001 description...

    to:

        >1402486.3|fig|1402486.3.peg.1|P999_00001 description...

    The genome ID is extracted from the FIG ID.
    """

    # Remove ">"
    header_without_gt = header[1:].strip()

    # Get first field before whitespace
    first_field = header_without_gt.split()[0]

    # Already has genome ID
    if first_field.count("|") >= 2:
        first_part = first_field.split("|", 1)[0]

        if re.match(r"^\d+\.\d+$", first_part):
            return header

    # Expected:
    # fig|1402486.3.peg.1|P999_00001
    match = re.match(
        r"^(fig\|(\d+\.\d+)\.peg\.[^|]+.*)$",
        first_field
    )

    if not match:
        print(
            f"[WARNING] Cannot extract genome ID from header: "
            f"{header}"
        )
        return header

    full_id = match.group(1)
    genome_id = match.group(2)

    # Replace the first field while keeping description
    remainder = header_without_gt[len(first_field):]

    return f">{genome_id}|{full_id}{remainder}"


def combine_ffn_files(ffn_files, output_file):

    output_dir = os.path.dirname(output_file)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    total_files = len(ffn_files)
    total_sequences = 0

    print("========================================")
    print("Combining FFN files")
    print("========================================")

    print(f"[INFO] Input files: {total_files}")

    with open(
        output_file,
        "w",
        encoding="utf-8"
    ) as outfile:

        for i, ffn_file in enumerate(ffn_files, start=1):

            print(
                f"[{i}/{total_files}] Processing: "
                f"{os.path.basename(ffn_file)}"
            )

            with open(
                ffn_file,
                "r",
                encoding="utf-8",
                errors="ignore"
            ) as infile:

                for line in infile:

                    if line.startswith(">"):

                        total_sequences += 1

                        line = add_genome_id_to_header(
                            line.rstrip("\n")
                        ) + "\n"

                    outfile.write(line)

    print()
    print("========================================")
    print("[DONE] Combined FFN created:")
    print(f"       {output_file}")
    print(f"[INFO] Input files: {total_files}")
    print(f"[INFO] Total sequences: {total_sequences}")
    print("========================================")


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Combine multiple FFN files into one all.ffn file "
            "and add genome IDs to FASTA headers."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Directory containing FFN files"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output combined FFN file"
    )

    args = parser.parse_args()

    if not os.path.isdir(args.input):
        raise FileNotFoundError(
            f"Input directory not found: {args.input}"
        )

    ffn_files = find_ffn_files(args.input)

    if not ffn_files:
        raise FileNotFoundError(
            f"No .ffn files found in: {args.input}"
        )

    combine_ffn_files(
        ffn_files,
        args.output
    )


if __name__ == "__main__":
    main()
