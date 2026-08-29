import glob

configfile: "config_PA.yaml"

PREFIX = config["prefix"]

RESULTS = "results"

CDHIT_DIR = f"{RESULTS}/cdhit"
MAFFT_DIR = f"{RESULTS}/mafft"
FFN_DIR = f"{RESULTS}/ffn"
SNP_DIR = f"{RESULTS}/snp_output"

INPUT_FAA = config["input_faa"]
FFN_DIR_INPUT = config["ffn_dir"]

METADATA = config["metadata"]
GENLIST = config.get("genlist", None)

SPECIES = config["species"]
ANTIBIOTICS = config["antibiotics"]


rule all:
    input:
        expand(
            f"{SNP_DIR}/snp_{{antibiotic}}.h5",
            antibiotic=ANTIBIOTICS
        )


rule run_cdhit:
    output:
        combined_faa = f"{CDHIT_DIR}/{PREFIX}.faa",
        cdhit_faa = f"{CDHIT_DIR}/{PREFIX}.cdhit.faa",
        clstr = f"{CDHIT_DIR}/{PREFIX}.cdhit.clstr"

    params:
        faa_dir = INPUT_FAA,
        output_dir = CDHIT_DIR,
        identity = config["cdhit"]["identity"],
        memory = config["cdhit"]["memory"],
        threads = config["cdhit"]["threads"],
        prefix = PREFIX,
        genlist = f"--genlist {GENLIST}" if GENLIST else ""

    threads:
        config["cdhit"]["threads"]

    shell:
        """
        python scripts/run_cdhit.py \
            --input {params.faa_dir} \
            --output {params.output_dir} \
            --prefix {params.prefix} \
            --identity {params.identity} \
            --memory {params.memory} \
            --threads {params.threads} \
            {params.genlist}
        """


rule run_mafft:
    input:
        fasta = f"{CDHIT_DIR}/{PREFIX}.faa",
        clstr = f"{CDHIT_DIR}/{PREFIX}.cdhit.clstr"

    output:
        aligned = f"{MAFFT_DIR}/{PREFIX}.aligned.faa"

    params:
        mafft = config["mafft"]["executable"]

    threads:
        config["cdhit"]["threads"]

    shell:
        """
        mkdir -p {MAFFT_DIR}

        python scripts/run_mafft.py \
            --fasta {input.fasta} \
            --clstr {input.clstr} \
            --output {output.aligned} \
            --mafft {params.mafft}
        """


rule combine_ffn:
    input:
        ffn_files = lambda wildcards: sorted(
            glob.glob(f"{FFN_DIR_INPUT}/*.ffn")
        )

    output:
        ffn = f"{FFN_DIR}/all.ffn"

    shell:
        """
        mkdir -p {FFN_DIR}

        python scripts/combine_ffn.py \
            --input {FFN_DIR_INPUT} \
            --output {output.ffn}
        """


rule extract_snp:
    input:
        aligned_fasta = f"{MAFFT_DIR}/{PREFIX}.aligned.faa",
        clstr = f"{CDHIT_DIR}/{PREFIX}.cdhit.clstr",
        ffn = f"{FFN_DIR}/all.ffn"

    output:
        snp_h5 = f"{SNP_DIR}/{PREFIX}.snp.h5"

    params:
        genlist = f"--genlist {GENLIST}" if GENLIST else ""

    threads:
        1

    shell:
        """
        mkdir -p {SNP_DIR}

        python scripts/extract_snp.py \
            --aligned_fasta {input.aligned_fasta} \
            --clstr {input.clstr} \
            --ffn {input.ffn} \
            --output_file {output.snp_h5} \
            {params.genlist}
        """


rule get_phenotype:
    input:
        snp_h5 = f"{SNP_DIR}/{PREFIX}.snp.h5"

    output:
        h5 = f"{SNP_DIR}/snp_{{antibiotic}}.h5"

    params:
        metadata = METADATA,
        species = SPECIES,
        genlist = f"--genlist {GENLIST}" if GENLIST else ""

    shell:
        """
        python scripts/getPhenotype.py \
            --input_h5 {input.snp_h5} \
            --metadata {params.metadata} \
            --antibiotic {wildcards.antibiotic} \
            --species "{params.species}" \
            --output_file {output.h5} \
            {params.genlist}
        """
