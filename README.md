# IBDMap

IBDMap is a multithreaded scalable application for IBD (Identity By Descent) mapping of binary traits. IBDMap implements a test statistic developed by Browning & Thompson et al., 2012 [cite] that detects regions where IBD sharing is statistically enriched among case-case pairs than among case-control pairs.

IBDMap addresses key challenges in IBD mapping by using a map-reduce approach to analyze genome-wide relative IBD segment enrichment, making it tractable for biobank-scale analysis. IBDMap’s computational flexibility enables execution in a wide array of compute environments.

### Map-reduce approach

Biobank scale datasets can require substantial computation time to complete enough permutations to reach genome-wide significance. For researchers working in a high performance computing environment we recommend a map-reduce strategy to enable analysis in reasonable periods of time.

The genome is conveniently broken into 22 autosomal components. Further, each job can be submitted with a subset of the total desired permutations. For example, if the user wishes to conduct a genome-wide analysis with 100,000 permutations on each breakpoint, then they can submit ten, 10,000 permutation runs for each chromosome, for a total of 100,000 permutations across all chromosomes. Each job can specify the seed used to initialize the random number generator used for permutation. Users should ensure that the seed is specified, and different, for each job. Jobs submitted without specifying the seed will generate a random seed from std::random_device as necessary.

The seed should be specified and different for each job on a chromosome. However, the seeds should be identical between chromosomes. so if you have five jobs for chromosome one, and your seeds are {1, 2, 3, 4, 5} for each job, then they should also be the same for chromosomes 2-22 as well. Results can be combined using IBDReduce, a complementary Python script designed to use a small amount of memory to collapse all the results from across chromosomes, and across permutation sets within a chromosome.

## Key features of IBDMap

- Implements flexible multiple testing correction options (FWER and FDR approaches)
- Includes supplementary statistical assessments, such as:
    - Likelihood ratio test for signal localization
    - Covariate adjustment for controlling environmental factors
- Compatible with output from common IBD detection tools (GERMLINE, iLASH, hap-IBD)
- Memory and CPU-efficient multithreading framework
- Fully containerized with local (via Dockerfile) and web-based Dockerhub implementations

## Compilation and Installation

### Submodule initialization

---

IBDMap has submodules that need to be initialized after cloning the directory from GitHub.

```bash
git submodule init
git submodule update
```

### C++ dependency installation

---

Dependencies to build dynamically:

- C++ compiler supporting C++17
- Armadillo: version >= 8.6
- Boost C++ Library

Dependencies to do a mostly static build on Linux:

- OpenBLAS or Intel MKL
- gfortran if using OpenBLAS

Full build environment for Ubuntu 20.04: 

```
apt install cmake build-essential libarmadillo-dev libboost-iostreams-dev libopenblas-dev gfortran
```

### IBDMap compilation

---

Create a build directory and run cmake and make from within the build directory:

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

If you're working in a high performance computing environment you can tell cmake where to look for both Armadillo and Boost.

If cmake fails to detect Armadillo, but you're sure it is available, you may need to direct cmake to the library. For example, this might happen when compiling on a cluster with packages in non-standard locations. 

In that case the following should work:

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DARMADILLO_INCLUDE_DIR=<path_to_armadillo>/include/ -DARMADILLO_LIBRARY=<path_to_armadillo>/lib64/libarmadillo.so
make
```

If you're having trouble with cmake detecting the correct compiler, try the following:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=<path_to_compiler> ..
make
```

The Boost location may also need to be specified if it isn't installed in a typical location:

```bash
cmake -DBOOST_ROOT=<path_to_boost> ..
```

You can combine the above as necessary. Earlier versions of the Armadillo library may work, but haven't been tested. 

If you need to change the compiler used from the one automatically detected to another, try a newer compiler:

```bash
cmake -DCMAKE_CXX_COMPILER=<path_to_executable> ..
```

If using OpenBLAS with cmake older than 3.21, it will not be able to find `lapack.a` for static linking. You can tell cmake where to find it. On Ubuntu 20.04:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DLAPACK_LIBRARIES=/usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.a ..
```

### IBDreduce dependency installation

---

The "reduce" step of IBDMap is done by the `ibdreduce.py` script. This script is built for Python >=3.9. 

IBDReduce requires the following Python packages:

* `numpy` (>=1.26.4 and <2.0.0)
* `scipy` (>=1.13.0 and <2.0.0)
* `zstandard` (>=0.22.0 and <1.0.0)

These dependencies are listed in a requirements.txt file with the IBDReduce sub directory. It is recommend that you install the dependencies into a virtualenv such as conda or venv. The requirements.txt file was made from a python3.11 environment to be consistent with the docker environment so it is recommend that you use python3.11

```bash
# virtual environment creation (python3.11 must be install)
python3.11 -m venv venv
# activate environment
source venv/bin/activate
# use pip to install required packages
pip install --upgrade pip && pip install -r IBDReduce/requirements.txt
```
Before running IBDReduce, make sure that this environment is activate otherwise you will get an error saying the dependencies are not installed.

*note on virtual environments:*
We recommend installing dependencies and running IBDReduce inside a [venv](https://docs.python.org/3/library/venv.html) or [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html) environment. All Python dependencies are also installed automatically if using the Docker build steps.

## Execution

### IBDMap

---

IBDMap requires several input files:

1. **Phenotype information** for each sample (case/control/exclusion status). This statuses are represented as 1/0/NA, respectively. Samples lacking phenotype information are automatically dropped from analysis. This file is expected to be tab separated and the first line is expected to have a header. first column of this file is expected to be titled “GRID”. Flag: `--pheno <FILE>`
2. **IBD segments** identified using an outside program (e.g., GERMLINE, iLASH, or hap-IBD) processed into IBDMap's input format. The output from this segment detection programs must first be converted to the proper format using the “[convert_ilash.py](https://github.com/rjbohlender/ibdmap/blob/master/tools/convert_ilash.py)” python script. Flag: `--input <FILE>`
3. **Genetic map file** for positioning breakpoints. Ensure that this is in the same build as your IBD segments. Flag: `--gmap <FILE>`
4. **Optional:**
    1. **DASH-processed data**. Only needed if DASH was run on the pairwise IBD data.  Flag(s): `--info <FILE>` (and `--dash`)
    2. **Covariates file.** Flag: `--cov <FILE>`   

Please use `ibdmap --help` to view all runtime options. 

### IBDReduce

---

IBDReduce is a companion Python tool (`IBDReduce/ibdreduce.py` ) that collapses IBDMap’s permutation output across chromosomes to generate an extreme value distribution for deriving adjusted p-values. This script is run as the second step of the IBD Mapping pipeline.

IBDReduce requires several inputs:

1. The directory (`—-prefix <FILE>`) and filename suffix (`--suffix <FILE>`) of the corresponding IBDMap results
2. The number of permutations used per IBDMap command executed, represented by `--nperm <NUM>`  and `--nruns <NUM>` . For example, if you used 500,000 total permutations and split it across 5 commands executed for each chromosome (for parallelization or otherwise), you would specify `--nperm 100000 --nruns 5`. 
3. The same genetic map location as for IBDMap. Flag: `--gmap <FILE>`
4. Output directory for IBDReduce results. Flag: `--output <FILE>`
5. Optional:
    1. Likelihood ratio test: use the same phenotype file provided to IBDMap. Flag: `--lrt <FILE>` 

#### Selecting breakpoints in IBDReduce

Use `--breakpoints <FILE>` with `IBDReduce/ibdreduce_v3.py` to analyze a subset
of breakpoints from a full IBDMap run. The file contains whitespace-separated
chromosome and base-pair position columns, using the same genome build as the
input results:

```text
CHROM POS
chr1 1000000
1 2500000
2 1800000
```

The `CHROM POS` header is optional. Blank lines and lines beginning with `#`
are ignored, extra columns are ignored, and duplicate locations are used once.
Locations can be listed in any order. Chromosomes must be autosomes 1–22, with
an optional `chr` prefix. Plain text, gzip, and Zstandard location files are
supported. The first two columns of a previous IBDReduce output can also be used
directly, including its metadata and header.

```bash
python3 IBDReduce/ibdreduce_v3.py \
  --prefix results/ --suffix 'chr{i}.run{j}.txt.zst' \
  --nperm 10000 --nruns 5 --gmap maps/ \
  --breakpoints selected_breakpoints.txt --output selected_results.tsv
```

Only chromosomes represented in the list are read. If `--single` is supplied,
the list must stay within that chromosome scope (including `--null`, if supplied).
Every requested location must exist in the input, and the selected rows must
have matching locations and ordering across permutation runs; mismatches fail
with an error.

Selection applies while reading all supported result formats, before allocating
the statistic matrix. Genetic-distance weights, the genome-wide baseline, and
multiple-testing adjustments are calculated using the selected data. Weights
are recomputed between consecutive selected locations, so results match running
IBDReduce on input files containing only those rows. Selected all-zero rows
still undergo the usual inactive-breakpoint removal.

Reproducing another run requires matching observed statistics and permutations
at the selected locations, including consistent permutation identities across
chromosomes. Without `--breakpoints`, IBDReduce continues to use all input rows.

### Docker Implementation

---

We have provided a Dockerfile in the top level IBDMap folder for streamlined building and execution. The same Docker image is additionally available at Dockerhub here: https://hub.docker.com/r/jtb114/ibdmap. This image includes both the ibdmap executable and the ibdreduce script

#### Pulling from Dockerhub

```bash
docker pull jtb114/ibdmap:latest
```

If you are in a compute environment that prevents you from running containers as sudo then you should use singularity instead of Docker. Singularity can pull the image using the following command:

```bash
singularity pull ibdmap.sif docker:jtb114/ibdmap:latest
```

This command will create a new singularity image called `ibdmap.sif`. The tool can then be run using the following command to view the help menu:

```bash
# Docker image
docker run --rm {image ID/tag} ibdmap -h
docker run --rm {image ID/tag} python3 /app/ibdreduce_v3.py -h

# Singularity image
singularity exec ibdmap.sif ibdmap -h
singularity exec ibdmap.sif python3 /app/ibdreduce_v3.py -h
```

*note:* the ibdreduce script is kept in the /app/ folder of the image so it is required to specify the full path as shown above otherwise you will get an error

## Data Formatting

### Inputs

---

We support segment files from GERMLINE, hap-IBD, iLASH, RaPID, and custom column
layouts. Use `tools/ibd_to_binary.py` to generate IBDF v3 input for the
`binary_format` branch, or `tools/ibd_input.py` to generate legacy text input
with the following columns:

```
chr pos add del
```

The `add` column represents a list of sample pairs with segments starting at the breakpoint. The `del` column is a list of sample pairs with segments terminating at the breakpoint.

#### Converting raw segments into reusable inputs

Both `ibd_to_binary.py` and `ibd_input.py` discover sample IDs directly from
the raw segment file. They do not read phenotype files or exclude samples
based on phenotype status. IBDMap applies the phenotype exclusions separately
for each analysis. The default minimum segment length is zero; use `--min`
only when you explicitly want a length cutoff.

```bash
# Binary input for the binary_format branch
python3 tools/ibd_to_binary.py \
  --input chr1.ilash.gz --format ilash --output reusable --chrom 1

# Legacy text input
python3 tools/ibd_input.py \
  --input chr1.ilash.gz --format ilash --output reusable --chrom 1
```

These commands create `reusable.chr1.ibdf` with `reusable.chr1.samples`, and
`reusable.chr1.small.txt.gz`, respectively. Each chromosome's binary mapping
uses sample IDs in first-appearance order. Always supply that chromosome's
companion `.samples` file to IBDMap. The raw converters no longer write a shared
`<prefix>.samples` file.

For compatibility with older commands, `--pheno` / `-p` is still accepted by
these two raw-segment converters, but is explicitly ignored—even if the file
contains `NA` values or omits samples. No phenotype file is needed. Existing
binary or text files that already lost segments through phenotype filtering
must be regenerated from the original segment data to recover those segments.

#### Converting existing IBDMap text input to binary

The `binary_format` branch's IBDMap executable requires IBDF v3 binary input.
If you already have a `chr`, `pos`, `add`, `del` text file from the older
converters, use `tools/ibd_text_to_binary.py`:

```bash
python3 tools/ibd_text_to_binary.py \
  --input old_results.chr1.small.txt.gz \
  --output converted \
  --chrom 1
```

This creates `converted.chr1.ibdf` and `converted.chr1.samples`. Use them with IBDMap:

```bash
./ibdmap --input converted.chr1.ibdf --samples converted.chr1.samples \
  --chromosome 1 --pheno phenotypes.txt --gmap chr1.gmap.gz \
  --seed 123 --output results
```

The text input must have a `chr pos add del` header and four tab-separated
columns. Each event column contains space-separated `cM:sample1-sample2`
entries or `NA`. Plain, gzip, and Zstandard input are supported; Zstandard
requires the Python `zstandard` package. DASH cluster input and IBDMap
statistic/permutation output are not supported by this converter.

By default, sample names are inferred from the text and assigned numeric IDs
in first-appearance order. Each chromosome has its own `.chrN.samples` mapping
because sample presence and ordering can differ across chromosomes. Always
pass the companion mapping for the chromosome being analyzed to `--samples`;
converting another chromosome with the same output prefix preserves existing
chromosomes' mappings. IBDMap resolves these IDs to names before matching the
phenotype file, so numeric IDs do not need to agree across chromosomes.

A pair must have exactly one possible split into
two nonempty sample names. Ambiguous pairs such as `a-b-c` cause an error with
the input line number and a suggestion to provide `--pheno`.

Optionally supply `--pheno phenotypes.txt` to resolve pairs containing
hyphenated sample IDs. This file must have a header and sample IDs in its first
column. When supplied, every sample referenced in the text must appear in that
file, and each pair must resolve uniquely. The converter retains all phenotype
statuses and leaves phenotype exclusions to IBDMap.

Conversion preserves every input breakpoint, including empty breakpoints, and
segment multiplicities. Same-position additions and deletions cancel. Input
must contain one chromosome with strictly increasing, unique positions and
enough preceding additions for each deletion. Conversion reads the file in
order, keeping active segment counts and the output index in memory.

Optional settings:

- `--min 3.0`: retain only segments at least 3 cM long. The default is `0.0`,
  retaining all lengths. Breakpoint locations remain present after filtering.
- `--checkpoint-interval 1000` (or `-N 1000`): write a full active-set checkpoint
  every 1,000 breakpoints; this is the default.

IBDF v3 stores lengths as float32, so decimal lengths may round slightly and
can cross a later IBDMap length cutoff. Set `--min` to the intended cutoff
during conversion to filter using the original text lengths.

If the data have been processed with DASH then there are two additional expected files. 

The first is the input file (`--input <FILE>`), which should have the following columns:

```
 chr pos n.cluster n.cluster.haplotype n.cluster.pair n.singleton.pair cluster.add cluster.del singleton.add singleton.del
```

The second is an info file (`--info <FILE>`) that contains information about haplotype clusters found by DASH used by IBDMap to filter segments. We expect the following columns:

```
chr clusterID start end count freq cM
```

`freq` is the frequency of a cluster in the sample. `cM` is segment length.

Reasonable defaults are provided for the corresponding `--cm <NUM>` and `--af <NUM>` options. Keep in mind any similar filtering parameters that might have been applied prior to running DASH.

### Outputs

---

#### IBDMap

IBDMap returns zstandard-compressed files, one per chromosome, containing breakpoint-specific test statistic and permutation results. Each file contains the following tab-separated columns (no header):

```bash
chr pos cscs-rate cscn-rate cncn-rate test-statistic perm1 perm2 ... permN
```

The first six columns are present for all runs, and the remaining columns represent the results of N permutations. 

#### IBDReduce

IBDReduce returns one file per phenotype with adjusted p-value computations per breakpoint. Each file contains the following tab-separated columns (info lines and header shown):

```
# your IBDReduce command string that was executed goes here
# Genome-wide Average: XXX
# Total breakpoints: YYY
CHROM	POS	cM	PVal	PValCI	PAdj	PAdjCutoff	Success	Permutation	Delta	LLik
```

The IBDReduce


## Citation

If you use IBDMap in your research, please cite:

```
Bohlender RJ, Evans GF, Baker JT, Landman JM, Frankel EG, Morrow AR, Petty LE, Petty AS, Bastarache L,
Chen HH, Zawistowski M, Samuels DC, Below JE, Huff CD. IBDMap: biobank scale identity-by-descent
mapping software for binary traits. [Publication details forthcoming]
```
