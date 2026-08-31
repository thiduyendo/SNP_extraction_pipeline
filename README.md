# SNP Extraction Pipeline for Bacteria

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![Snakemake](https://img.shields.io/badge/Snakemake-workflow-blue.svg)](https://snakemake.readthedocs.io/)
[![GitHub Downloads](https://img.shields.io/github/downloads/thiduyendo/SNP_extraction_pipeline/total.svg)](https://github.com/thiduyendo/SNP_extraction_pipeline/releases)

A reproducible **Snakemake workflow for reference-free SNP feature extraction from bacterial coding sequences**.

The pipeline identifies homologous protein sequences using CD-HIT, aligns them with MAFFT, reconstructs codon-preserving nucleotide alignments from corresponding FFN sequences, performs sequence-consistency checks, and extracts genome-level SNP features in HDF5 format.

The workflow was developed for comparative bacterial genomics and is particularly suitable for generating SNP-based feature matrices for **population-level genomic analysis and machine-learning applications, including antimicrobial-resistance (AMR) prediction**.

---

## 🧬 Why this pipeline?

Whole-genome bacterial datasets contain multiple complementary sources of genetic variation. Differences between genomes can arise from **gene-content variation**, such as genes that are present in one genome but absent in another, as well as from **nucleotide-level variation within homologous genes**.

This pipeline focuses primarily on the extraction of **nucleotide-level SNP features within homologous coding sequences**.

Protein sequences from all input genomes are first grouped into homologous clusters using **CD-HIT**. These clusters provide the homology framework for identifying corresponding coding sequences across genomes. Protein sequences within each cluster are then aligned using **MAFFT**, and the corresponding nucleotide coding sequences are mapped onto the protein alignments in a **codon-preserving, reference-free manner**.

Variable nucleotide positions are subsequently extracted as genome-level SNP features and stored as an HDF5 matrix.

### What type of genomic variation does it capture?

The resulting SNP representation can capture:

* single-nucleotide substitutions within homologous genes;
* synonymous nucleotide variation;
* nonsynonymous nucleotide variation;
* fine-scale sequence differences among bacterial genomes; and
* nucleotide variation within homologous protein clusters.

### Relationship to gene-content variation

Gene-content and SNP variation represent complementary aspects of bacterial genomic diversity:

| Representation | Captures                                                |
| -------------- | ------------------------------------------------------- |
| Gene content   | Which homologous genes are present or absent            |
| SNP variation  | How the nucleotide sequences of homologous genes differ |

A gene-content presence/absence matrix can therefore complement the SNP representation when both types of variation are required.

> **Note:** The current version of this pipeline primarily focuses on **SNP extraction and SNP-matrix generation**. It does not currently generate a separate genome × gene-cluster presence/absence matrix.

---

## 🧬 Reference-free SNP extraction

A key feature of this pipeline is its **reference-free, homologous-cluster-based strategy**.

Unlike conventional reference-based SNP calling, the workflow does not require one genome to be designated as the reference genome.

Instead, the workflow follows:

```text
Protein FASTA files
        │
        ▼
   CD-HIT clustering
        │
        ▼
 Homologous protein clusters
        │
        ▼
    MAFFT alignment
        │
        ▼
 Corresponding FFN sequences
        │
        ▼
Codon-preserving nucleotide mapping
        │
        ▼
 Sequence consistency checks
        │
        ▼
    SNP extraction
        │
        ▼
 Genome × SNP matrix
```

For each homologous cluster, nucleotide sequences are reconstructed according to the protein-level alignment. A nucleotide position is retained as a SNP feature when at least two valid nucleotide alleles (`A`, `T`, `C`, or `G`) are observed among the genomes represented in that alignment.

No single genome is required to provide the reference allele.

### Why is a reference-free strategy useful?

Reference-based approaches provide a convenient coordinate system, but using a single reference genome can introduce limitations when analyzing genetically diverse bacterial populations.

For example, a reference genome may contain genes or genomic regions that are absent from other genomes. Highly divergent sequences may also be poorly represented or difficult to map to the selected reference.

The reference-free strategy used here instead derives variation from **homologous sequences directly observed across the input genomes**.

This provides several advantages:

### 1. Reduced dependence on a single reference genome

No individual strain is assumed to represent the genomic structure of the entire population.

This is useful when analyzing collections containing substantial genetic diversity or multiple lineages.

### 2. Better compatibility with genomic diversity

Homologous clusters are identified from the input dataset itself. Nucleotide variation can therefore be analyzed within coding sequences represented among the analyzed genomes without requiring every genome to conform to the structure of one reference genome.

### 3. Reduced reference-allele dependence

SNP features are defined from nucleotide alleles observed within each homologous alignment rather than from differences relative to one predetermined reference allele.

The goal is therefore to describe **variation within the analyzed population**, rather than variation relative to an arbitrarily selected strain.

### 4. Natural integration with protein homology

Protein clustering establishes homologous relationships before nucleotide variation is extracted.

This provides a biologically informed framework in which:

```text
protein homology
       ↓
sequence alignment
       ↓
codon-preserving nucleotide alignment
       ↓
SNP variation
```

is used to define comparable nucleotide positions.

### 5. Useful for comparative and machine-learning analyses

The resulting SNP features represent nucleotide variation directly observed across the analyzed genomes.

They can therefore be used as genome-level predictors in downstream applications such as:

* comparative genomics;
* bacterial population analysis;
* genotype–phenotype association;
* antimicrobial-resistance prediction; and
* machine-learning models using genomic features.

> **Important:** This pipeline is **not a conventional reference-based SNP caller**. It is designed to extract comparative SNP features from homologous coding sequences across a collection of bacterial genomes.

---

## 🛡️ Sequence consistency and mapping quality control

Because nucleotide sequences are reconstructed from protein alignments, the correctness of the **protein–FFN correspondence** is critical.

An incorrect protein/FFN pairing, truncated CDS, or inconsistent annotation could otherwise cause nucleotide positions to be incorrectly assigned to amino-acid alignment positions.

The pipeline therefore performs explicit **sequence-consistency checks during protein-to-nucleotide mapping** before SNP features are extracted.

For each protein/FFN pair, the expected coding-sequence length is calculated from the ungapped protein sequence:

```text
Expected CDS length = ungapped protein length × 3
```

The mapping procedure checks that:

* the CDS length is compatible with codon-based mapping;
* the CDS length is divisible by three;
* an additional terminal stop codon is handled appropriately when present;
* the actual CDS length is consistent with the expected length derived from the protein sequence;
* the complete nucleotide sequence can be consumed during codon-by-codon mapping; and
* nucleotide alignments generated within a homologous cluster have consistent alignment lengths.

If these checks fail, the affected sequence or cluster is not used for SNP extraction and diagnostic information is reported.

### What problems can these checks detect?

These checks can help identify:

* inconsistent genome annotations;
* incorrect protein–FFN correspondence;
* truncated coding sequences;
* unexpected terminal stop-codon structure;
* inconsistent sequence identifiers;
* incomplete CDS records; and
* other input-data inconsistencies.

This provides an additional layer of protection against generating biologically incorrect SNP features from mismatched protein and nucleotide sequences.

### Codon-preserving nucleotide mapping

Each aligned amino acid corresponds to three nucleotides.

For example:

```text
Protein:       A    L    -    G
               ↓    ↓    ↓    ↓
Nucleotide:   ATG  CTG  ---  GGT
```

A protein alignment gap (`-`) is therefore represented by three nucleotide gaps (`---`), preserving the reading frame.

This allows the workflow to transfer the homologous relationship established at the protein level into nucleotide space while maintaining codon structure.

The resulting nucleotide alignment is therefore constrained by both:

```text
Protein-level homology
          +
Coding-sequence consistency
          ↓
Reliable nucleotide alignment
```

before SNP features are extracted.

---

## 🔬 How SNPs are detected

For each homologous protein cluster:

1. Protein sequences are aligned using MAFFT.
2. The corresponding FFN coding sequences are retrieved.
3. Protein alignments are converted into nucleotide alignments in a codon-preserving manner.
4. Protein/CDS consistency checks are performed during mapping.
5. Each nucleotide alignment column is examined across genomes.
6. Variable nucleotide positions are retained as SNP features.

Only standard nucleotide bases are considered valid:

```text
A
T
C
G
```

Gaps (`-`) and ambiguous/non-standard bases such as `N` are **not considered valid SNP alleles**.

### Example of a retained SNP

```text
Genome 1    A
Genome 2    A
Genome 3    G
Genome 4    A
```

This position is retained because two valid nucleotide alleles, `A` and `G`, are observed.

### Example of a non-SNP position

```text
Genome 1    A
Genome 2    A
Genome 3    -
Genome 4    N
```

This position is **not** considered a SNP because the gap (`-`) and ambiguous nucleotide (`N`) are excluded from allele determination. Among the valid nucleotides, only `A` is observed.

Both bi-allelic and multi-allelic variable positions can be retained.

For example:

```text
Genome 1    A
Genome 2    G
Genome 3    C
Genome 4    A
```

contains three valid alleles (`A`, `G`, and `C`) and is therefore considered a variable position.

### SNP feature naming

The resulting feature name has the form:

```text
cluster3_757
```

where:

* `cluster3` identifies the homologous protein cluster;
* `757` identifies the nucleotide position within the cluster alignment.

> **Note:** The position is an alignment position within the homologous cluster and is **not a genomic coordinate on a reference chromosome**.

---

## ⚙️ Pipeline at a glance

| Step                | Tool       | Main purpose                                 |
| ------------------- | ---------- | -------------------------------------------- |
| Protein clustering  | **CD-HIT** | Identify homologous protein groups           |
| Protein alignment   | **MAFFT**  | Establish homologous sequence positions      |
| Nucleotide mapping  | **Python** | Map coding sequences onto protein alignments |
| Sequence QC         | **Python** | Validate protein/CDS consistency             |
| SNP extraction      | **Python** | Identify variable nucleotide positions       |
| Phenotype filtering | **Python** | Generate antibiotic-specific SNP matrices    |

The primary output is a **genome × SNP feature matrix stored in HDF5 format**.

> 💡 **Tip:** The workflow is managed by Snakemake, so interrupted runs can be resumed and completed intermediate results can be reused when they are compatible with the current inputs and parameters.

---

## 🚀 Quick start

```bash
git clone https://github.com/thiduyendo/SNP_extraction_pipeline.git
cd SNP_extraction_pipeline

conda env create -f environment.yml
conda activate snp_pipeline

# Check the workflow without running it
snakemake -s Snakefile -n

# Run the example workflow
snakemake -s Snakefile --cores 4
```

After a successful run, the main SNP matrix is generated at:

```text
results/snp_output/all.95.0.snp.h5
```

See **[Troubleshooting](#7-troubleshooting-and-resuming-the-workflow)** for how to resume failed runs, reuse intermediate results, run specific targets, or force a rule to rerun.

---

# 1. Requirements

The pipeline requires:

* Python 3.9+
* Snakemake
* CD-HIT
* MAFFT
* Python packages listed in `environment.yml`

Recommended installation using Conda/Mamba:

```bash
conda env create -f environment.yml
conda activate snp_pipeline
```

Alternatively:

```bash
mamba env create -f environment.yml
mamba activate snp_pipeline
```

Check that the required external programs are available:

```bash
cd-hit -h
mafft --version
snakemake --version
```

---

# 2. Repository structure

After downloading the repository, the directory should look approximately like:

```text
SNP_extraction_pipeline/
├── config_PA.yaml
├── environment.yml
├── Snakefile
├── README.md
├── CITATION.cff
├── .gitignore
├── scripts/
│   ├── combine_ffn.py
│   ├── extract_snp.py
│   ├── getPhenotype.py
│   ├── run_cdhit.py
│   └── run_mafft.py
└── data/
    └── example/
        ├── faa/
        ├── ffn/
        ├── genlist.txt
        └── PATRIC_genomes_AMR.txt
```

The `data/example/` directory contains a small example dataset for testing the workflow.

Large experimental datasets and generated results should **not** be committed to GitHub.

---

# 3. Input data

The pipeline requires protein FASTA files and corresponding nucleotide coding-sequence (FFN) files.

Genome lists and phenotype metadata can optionally be provided for selecting genomes and generating antibiotic-specific SNP datasets.

## 3.1 Protein FASTA files

Each genome should have a protein FASTA file.

For example:

```text
data/example/faa/
├── 1163395.3.PATRIC.faa
├── 1402486.3.PATRIC.faa
├── 1402487.3.PATRIC.faa
└── ...
```

The filename is used to determine the genome ID.

For example:

```text
1402486.3.PATRIC.faa
```

corresponds to:

```text
1402486.3
```

---

## 3.2 Nucleotide FFN files

Each genome should have a corresponding nucleotide coding-sequence file:

```text
data/example/ffn/
├── 1163395.3.PATRIC.ffn
├── 1402486.3.PATRIC.ffn
├── 1402487.3.PATRIC.ffn
└── ...
```

The protein and nucleotide files must correspond to the same genome and compatible gene annotations.

The pipeline combines these files into:

```text
results/ffn/all.ffn
```

The combined FFN file is then used to map protein alignments back to nucleotide coding sequences.

---

## 3.3 Genome list

`genlist.txt` contains the genome IDs that should be included in the analysis.

Example:

```text
1402486.3.fna
1402487.3.fna
1402488.3.fna
1402489.3.fna
```

The genome list is optional.

If `genlist` is specified in the configuration file, only those genomes are processed.

If it is omitted, all available genomes are used.

---

## 3.4 Phenotype metadata

The phenotype metadata file is used by `getPhenotype.py` to generate antibiotic-specific SNP datasets.

The metadata should contain genome identifiers and antibiotic susceptibility information required by the pipeline.

The exact column structure depends on the metadata source and analysis.

An example metadata file is provided in:

```text
data/example/PATRIC_genomes_AMR.txt
```

---

# 4. Configure the pipeline

Edit:

```text
config_PA.yaml
```

A typical configuration is:

```yaml
species: "Pseudomonas aeruginosa"

prefix: "all.95.0"

input_faa: "data/example/faa"

ffn_dir: "data/example/ffn"

metadata: "data/example/PATRIC_genomes_AMR.txt"

genlist: "data/example/genlist.txt"

cdhit:
  identity: 0.95
  memory: 10000
  threads: 30

mafft:
  executable: "mafft"

antibiotics:
  - tobramycin
```

### Important

The paths above are **relative paths** and are intended to make the example workflow portable.

For your own dataset, replace the paths with the locations of your files.

---

# 5. CD-HIT parameters

The following parameters control protein clustering:

```yaml
cdhit:
  identity: 0.95
  memory: 10000
  threads: 30
```

For example:

```yaml
identity: 0.95
```

means that proteins are clustered using a 95% sequence identity threshold.

The output prefix is controlled by:

```yaml
prefix: "all.95.0"
```

This produces files such as:

```text
results/cdhit/all.95.0.faa
results/cdhit/all.95.0.cdhit.faa
results/cdhit/all.95.0.cdhit.clstr
```

The CD-HIT cluster definitions are subsequently used to identify homologous protein sequences for alignment and SNP extraction.

---

# 6. Run the pipeline

First activate the Conda environment:

```bash
conda activate snp_pipeline
```

From the repository root:

```bash
snakemake -s Snakefile --cores 30
```

For a dry run before executing:

```bash
snakemake -s Snakefile -n
```

For a more detailed dry run:

```bash
snakemake -s Snakefile -n -p
```

The `-p` option prints the commands that Snakemake would execute.

---

# 7. Troubleshooting and resuming the workflow

## 7.1 Resume an interrupted run

Snakemake can resume an interrupted or failed workflow without rerunning completed steps.

If the pipeline stops because of an error, fix the underlying problem and rerun:

```bash
snakemake -s Snakefile --cores 30
```

Snakemake checks the existing output files and executes only the jobs that are missing or need to be updated.

For example, if CD-HIT and FFN preparation completed successfully but MAFFT failed:

```text
Combine FASTA       ✓
CD-HIT               ✓
MAFFT                ✗
Combine FFN          ✓
SNP extraction       ✗
Phenotype filtering  ✗
```

After fixing the MAFFT problem, rerunning:

```bash
snakemake -s Snakefile --cores 30
```

will reuse completed outputs and continue from the missing steps.

> **Recommendation:** Do not delete completed intermediate files just to resume a failed run. Keep them in place and let Snakemake determine which jobs still need to run.

---

## 7.2 Run only a downstream step

A specific output can be requested directly. Snakemake will automatically determine which dependencies are required.

For example:

```bash
snakemake -s Snakefile \
    results/snp_output/all.95.0.snp.h5 \
    --cores 30
```

If upstream outputs are available and considered up-to-date, they will not be regenerated.

Similarly, an antibiotic-specific SNP matrix can be generated with:

```bash
snakemake -s Snakefile \
    results/snp_output/snp_tobramycin.h5 \
    --cores 30
```

---

## 7.3 Reuse existing intermediate files

Intermediate results can be reused to avoid repeating computationally expensive steps.

For example, if the CD-HIT and MAFFT results already exist:

```text
results/cdhit/all.95.0.cdhit.faa
results/cdhit/all.95.0.cdhit.clstr
results/mafft/all.95.0.aligned.faa
```

you can request the downstream SNP output:

```bash
snakemake -s Snakefile \
    results/snp_output/all.95.0.snp.h5 \
    --cores 30
```

Snakemake will use existing outputs when they are considered up-to-date.

**Important:** Existing intermediate files should only be reused if they were generated from compatible input data and pipeline parameters. For example, CD-HIT results generated with a different clustering threshold or a different genome dataset should not be reused.

---

## 7.4 Force a rule to run again

Use this only when you intentionally want to regenerate an existing output.

If the MAFFT rule is named `run_mafft`:

```bash
snakemake -s Snakefile \
    --forcerun run_mafft \
    --cores 30
```

A dry run is recommended before forcing a rule:

```bash
snakemake -s Snakefile \
    --forcerun run_mafft \
    -n -p
```

---

## 7.5 Check what Snakemake will run

Before executing the workflow:

```bash
snakemake -s Snakefile -n
```

For detailed commands:

```bash
snakemake -s Snakefile -n -p
```

This is particularly useful when resuming a failed workflow or reusing intermediate files.

---

## 7.6 Starting from an existing intermediate result

The pipeline can be continued from an existing intermediate result if that result matches the expected output of the corresponding workflow step.

For example, if a valid aligned protein FASTA already exists at:

```text
results/mafft/all.95.0.aligned.faa
```

the downstream SNP extraction can be generated without repeating CD-HIT or MAFFT:

```bash
snakemake -s Snakefile \
    results/snp_output/all.95.0.snp.h5 \
    --cores 30
```

The existing alignment must correspond to the same genome dataset, cluster definitions, and expected file format used by the downstream SNP extraction step.

---

# 8. Output files

The pipeline creates a `results/` directory automatically.

The main output structure is:

```text
results/
├── cdhit/
│   ├── all.95.0.faa
│   ├── all.95.0.cdhit.faa
│   └── all.95.0.cdhit.clstr
│
├── ffn/
│   └── all.ffn
│
├── mafft/
│   └── all.95.0.aligned.faa
│
└── snp_output/
    ├── all.95.0.snp.h5
    └── snp_tobramycin.h5
```

The primary output is:

```text
results/snp_output/all.95.0.snp.h5
```

Phenotype-specific SNP matrices are generated when antibiotic information is configured.

Generated results are excluded from Git using `.gitignore`.

---

# 9. HDF5 output structure

The main SNP file:

```text
results/snp_output/all.95.0.snp.h5
```

contains:

```text
/
├── index
├── snp_matrix/
│   ├── cluster0
│   ├── cluster1
│   ├── cluster2
│   └── ...
│
└── feature_names/
    ├── cluster0
    ├── cluster1
    ├── cluster2
    └── ...
```

### `index`

Contains the genome IDs.

For example:

```text
1402486.3
1402487.3
1402488.3
...
```

### `snp_matrix`

Contains the encoded SNP matrices for individual homologous clusters.

### `feature_names`

Contains the feature names corresponding to the columns of each cluster matrix.

---

# 10. SNP encoding

The nucleotide bases are encoded numerically:

```text
A → 2
T → 3
G → 4
C → 5
other/invalid → 0
```

For example:

```text
             cluster3_757
1402486.3          2
1402487.3          2
1402488.3          4
1402489.3          2
```

Here:

```text
2 = A
4 = G
```

so the position contains an A/G polymorphism.

---

# 11. Reading the HDF5 file with Python

The HDF5 file can be inspected using PyTables.

```python
import tables

h5 = tables.open_file(
    "results/snp_output/all.95.0.snp.h5",
    mode="r"
)

print(h5)
```

To inspect the genome IDs:

```python
genomes = h5.root.index[:]

print(genomes)
```

To list available clusters:

```python
print(h5.root.snp_matrix._v_children.keys())
```

To inspect feature names for a cluster:

```python
features = h5.root.feature_names.cluster3[:]

print(features)
```

To inspect the SNP matrix for the same cluster:

```python
matrix = h5.root.snp_matrix.cluster3[:]

print(matrix.shape)
print(matrix)
```

The number of rows corresponds to genomes, while the number of columns corresponds to SNP features in that cluster.

---

# 12. Reading an antibiotic-specific SNP matrix

The phenotype-specific HDF5 files generated by the pipeline can be read using pandas.

For example, for meropenem:

```python
import pandas as pd

h5_file = "results/snp_output/snp_meropenem.h5"

# Check available HDF5 keys
with pd.HDFStore(h5_file, mode="r") as store:
    print(store.keys())
```

The available key should be checked before loading the matrix.

For example, if the file contains the appropriate SNP matrix key:

```python
df = pd.read_hdf(h5_file, key="snp_matrix")

print("Shape:", df.shape)
print(df.head())
```

The resulting DataFrame contains the SNP features used for the antibiotic-specific analysis.

To inspect feature names:

```python
print(df.columns.tolist())
```

To inspect genome IDs:

```python
print(df.index.tolist())
```

For example:

```text
             cluster3_757  cluster3_758  cluster7_102
1402486.3               2              2              4
1402487.3               2              4              4
1402488.3               4              2              4
```

The values represent encoded nucleotide alleles:

```text
A → 2
T → 3
G → 4
C → 5
other/invalid → 0
```

The HDF5 key should always be checked first because the exact key depends on how the HDF5 file was written.

For example:

```python
with pd.HDFStore(h5_file, mode="r") as store:
    print(store.keys())

    for key in store.keys():
        print(key)
        print(store[key].shape)
```

This provides a simple way to inspect the contents of any `snp_<antibiotic>.h5` file.

---

# 13. Reproducibility

For reproducible analyses, we recommend recording:

* the Git commit or release used for the analysis;
* the configuration file;
* the input dataset;
* the Conda environment;
* the versions of CD-HIT and MAFFT;
* the Snakemake version; and
* the final HDF5 output.

The repository provides:

```text
environment.yml
CITATION.cff
```

to help document the computational environment and citation information.

For published analyses, we recommend citing the exact release or archived version used to generate the results.

---

# 14. Example workflow

A complete example using the included test data is:

```bash
git clone https://github.com/thiduyendo/SNP_extraction_pipeline.git

cd SNP_extraction_pipeline

conda env create -f environment.yml

conda activate snp_extraction_pipeline

# Check the workflow
snakemake -s Snakefile -n

# Run the example
snakemake -s Snakefile --cores 4
```

After successful completion:

```text
results/snp_output/all.95.0.snp.h5
```

and the configured antibiotic-specific files should be available.

---

# 15. Using your own dataset

To run the pipeline on your own bacterial dataset:

### Step 1 — Prepare protein FASTA files

```text
data/faa/
├── genome1.PATRIC.faa
├── genome2.PATRIC.faa
└── ...
```

### Step 2 — Prepare corresponding FFN files

```text
data/ffn/
├── genome1.PATRIC.ffn
├── genome2.PATRIC.ffn
└── ...
```

### Step 3 — Create a genome list if desired

```text
data/genlist.txt
```

### Step 4 — Prepare phenotype metadata

Only required if antibiotic-specific SNP datasets are needed:

```text
data/PATRIC_genomes_AMR.txt
```

### Step 5 — Update `config_PA.yaml`

```yaml
input_faa: "data/faa"
ffn_dir: "data/ffn"
metadata: "data/PATRIC_genomes_AMR.txt"
genlist: "data/genlist.txt"
```

### Step 6 — Select the antibiotic

```yaml
antibiotics:
  - tobramycin
```

### Step 7 — Run

```bash
snakemake -s Snakefile --cores 30
```

---

# 16. Important considerations

The nucleotide mapping assumes that the protein sequences and corresponding FFN coding sequences are consistent.

In particular:

* protein and FFN records should correspond to the same gene;
* protein and nucleotide annotations should originate from compatible genome annotations;
* coding sequences should normally have lengths divisible by three;
* coding sequences should represent the expected CDS;
* genome IDs must be consistent between protein FASTA, FFN, genome list, and phenotype metadata; and
* protein sequences within a homologous cluster should be sufficiently related for meaningful alignment.

The pipeline performs sequence-consistency and mapping-integrity checks during nucleotide reconstruction and reports warnings when sequences cannot be successfully mapped.

Because SNP feature positions are defined within homologous cluster alignments, they should **not** be interpreted as genomic coordinates on a particular reference genome.

---

# 17. Citation

If you use this pipeline in a publication, please cite the software repository and the specific release used for your analysis.

Citation information is provided in:

```text
CITATION.cff
```

For archived releases, please use the corresponding Zenodo DOI to cite the specific version used in your analysis.

For reproducibility, we recommend citing the exact release or version used in your analysis rather than the current state of the repository.

---

