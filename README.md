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

We support several different file formats for IBD data. IBDMap is compatible with output from GERMLINE, hap-IBD, and iLASH. After shared segments have been detected using one of these programs, the data must be reformatted using the provided Python script `tools/convert_ilash.py`. Once the data are reformatted it should have the following columns:

```
chr pos add del
```

The `add` column represents a list of sample pairs with segments starting at the breakpoint. The `del` column is a list of sample pairs with segments terminating at the breakpoint.

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
