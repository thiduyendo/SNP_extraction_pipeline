# SNP Extraction Pipeline for Bacteria

A reproducible Snakemake pipeline for extracting nucleotide-level SNP features from *Pseudomonas aeruginosa* protein sequences.

The pipeline combines protein sequences, clusters homologous proteins with CD-HIT, aligns clustered proteins with MAFFT, maps protein alignments back to nucleotide coding sequences (FFN), and extracts variable nucleotide positions as SNP features.

The final SNP features are stored in HDF5 format and can optionally be filtered by antibiotic phenotype.

---

## Overview

The pipeline performs the following steps:

```text
Protein FASTA files
        │
        ▼
   Combine FASTA
        │
        ▼
      CD-HIT
   protein clustering
        │
        ├── .cdhit.faa
        └── .cdhit.clstr
        │
        ▼
      MAFFT
 protein alignment
        │
        ▼
   Combine FFN files
        │
        ▼
 Map protein alignment
 to nucleotide sequences
        │
        ▼
    SNP extraction
        │
        ▼
 all.<prefix>.snp.h5
        │
        ▼
 Phenotype filtering
        │
        ▼
 snp_<antibiotic>.h5
```

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

The pipeline requires four types of input.

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

The protein and nucleotide files must correspond to the same genome IDs.

The pipeline combines these files into:

```text
results/ffn/all.ffn
```

The combined FFN file is then used to map the protein alignment back to nucleotide sequences.

---

## 3.3 Genome list

`genlist.txt` contains the genome IDs that should be included in the analysis.

Example:

```text
1402486.3
1402487.3
1402488.3
1402489.3
```

The genome list is optional.

If `genlist` is specified in the configuration file, only those genomes are processed.

If it is omitted, all available genomes are used.

---

## 3.4 Phenotype metadata

The phenotype metadata file should contain genome IDs and antibiotic susceptibility information required by `getPhenotype.py`.

The exact column structure depends on the metadata used for the analysis.

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
  - meropenem
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

# 7. Run only a specific output

For example, to generate the complete SNP HDF5 file:

```bash
snakemake -s Snakefile \
    results/snp_output/all.95.0.snp.h5 \
    --cores 1
```

To generate an antibiotic-specific SNP matrix:

```bash
snakemake -s Snakefile \
    results/snp_output/snp_meropenem.h5 \
    --cores 1
```

Snakemake automatically determines which previous steps are required.

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
    └── snp_meropenem.h5
```

Generated results are excluded from Git using `.gitignore`.

---

# 9. How SNPs are detected

For each CD-HIT protein cluster:

1. Protein sequences are aligned using MAFFT.
2. The protein alignment is converted to a nucleotide alignment using the corresponding FFN sequences.
3. Each nucleotide alignment column is examined.
4. Only standard nucleotide bases are considered:

```text
A
T
C
G
```

Gaps (`-`) and ambiguous/non-standard bases such as `N` are not considered valid SNP alleles.

A position is classified as a SNP when at least two different valid nucleotide bases are observed among the genomes.

For example:

```text
Genome 1    A
Genome 2    A
Genome 3    G
Genome 4    A
```

is a SNP because both `A` and `G` occur.

A column such as:

```text
Genome 1    A
Genome 2    A
Genome 3    -
Genome 4    A
```

is not considered a SNP because the only valid nucleotide allele is `A`.

The resulting feature name has the form:

```text
cluster3_757
```

where:

* `cluster3` = CD-HIT cluster
* `757` = nucleotide alignment position

---

# 10. HDF5 output structure

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

Contains the encoded SNP matrices for individual CD-HIT clusters.

### `feature_names`

Contains the feature names corresponding to the columns of each cluster matrix.

---

# 11. SNP encoding

The nucleotide bases are encoded numerically:

```text
A → 2
T → 3
G → 4
C → 5
other/invalid → 0
```

Therefore, an example SNP feature could look like:

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

# 12. Reading the HDF5 file with Python

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

# 13. Reading an antibiotic-specific SNP matrix

The phenotype-specific HDF5 files generated by the pipeline can be read using pandas.

For example, for meropenem:

```python
import pandas as pd

h5_file = "results/snp_output/snp_meropenem.h5"

# Check available HDF5 keys
with pd.HDFStore(h5_file, mode="r") as store:
    print(store.keys())
```

The available key should be used to load the SNP matrix. For example, if the file contains the key `data`:

```python
df = pd.read_hdf(h5_file, key="snp_matrix")

print("Shape:", df.shape)
print(df.head())
```

The resulting DataFrame contains the SNP features used for the antibiotic-specific analysis.

To inspect the feature names:

```python
print(df.columns.tolist())
```

To inspect the genome IDs:

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

# 14. Reproducibility

For reproducible analyses, we recommend recording:

* the Git commit or release used for the analysis
* the configuration file
* the input dataset
* the Conda environment
* the versions of CD-HIT and MAFFT
* the final HDF5 output

The repository provides:

```text
environment.yml
CITATION.cff
```

to help document the computational environment and citation information.

---

# 15. Example workflow

A complete example using the included test data is:

```bash
git clone <REPOSITORY-URL>

cd SNP_extraction_pipeline

conda env create -f environment.yml

conda activate snp_pipeline

snakemake -s Snakefile -n

snakemake -s Snakefile --cores 4
```

After successful completion:

```text
results/snp_output/all.95.0.snp.h5
```

and the configured antibiotic-specific files should be available.

---

# 16. Using your own dataset

To run the pipeline on your own *Pseudomonas aeruginosa* dataset:

### Step 1

Prepare protein FASTA files:

```text
my_data/faa/
├── genome1.PATRIC.faa
├── genome2.PATRIC.faa
└── ...
```

### Step 2

Prepare corresponding FFN files:

```text
my_data/ffn/
├── genome1.PATRIC.ffn
├── genome2.PATRIC.ffn
└── ...
```

### Step 3

Create a genome list if desired:

```text
my_data/genlist.txt
```

### Step 4

Prepare phenotype metadata:

```text
my_data/PATRIC_genomes_AMR.txt
```

### Step 5

Update `config_PA.yaml`:

```yaml
input_faa: "my_data/faa"
ffn_dir: "my_data/ffn"
metadata: "my_data/PATRIC_genomes_AMR.txt"
genlist: "my_data/genlist.txt"
```

### Step 6

Select the antibiotic:

```yaml
antibiotics:
  - meropenem
```

### Step 7

Run:

```bash
snakemake -s Snakefile --cores 30
```

---

# 17. Important considerations

The nucleotide mapping assumes that the protein sequences and corresponding FFN coding sequences are consistent.

In particular:

* protein and FFN records must correspond to the same gene
* coding sequences should represent the expected CDS
* CDS sequences should normally have lengths divisible by three
* protein and nucleotide records should originate from the same genome annotation
* genome IDs must be consistent between protein FASTA, FFN, genome list, and phenotype metadata

The pipeline reports mapping warnings when a protein sequence cannot be successfully mapped to its nucleotide sequence.

---

# 18. Citation

If you use this pipeline in a publication, please cite the repository and the associated archived release/DOI.

Citation information is provided in:

```text
CITATION.cff
```

Once the repository is archived through Zenodo, the DOI for the specific release should be used when citing the version of the pipeline used for an analysis.

---

# License

See the repository license for terms of use.

