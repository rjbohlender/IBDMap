# IBDMap

IBDMap is a multithreaded scalable application for IBD mapping in Biobank
scale datasets.

## Compiling

This project has submodules that need to be initialized after cloning the directory from GitHub. 

```bash
git submodule init
git submodule update
```

Dependencies to build dynamically:
    - C++ compiler supporting C++17
    - Armadillo: version >= 8.6
    - Boost C++ Library

Dependencies to do a mostly static build on Linux:
    - OpenBLAS or Intel MKL
    - gfortran if using OpenBLAS
    - Full build environment for Ubuntu 20.04: `apt install cmake build-essential libarmadillo-dev libboost-iostreams-dev libopenblas-dev gfortran`
    
Create a build directory and run CMAKE and make from within the build directory:

```bash
mkdir build && cd build

cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

If you're working in a high performance computing environment you can tell CMake where to look for both Armadillo and
Boost.

If cmake fails to detect Armadillo, but you're sure it is available, you may need to direct cmake to the library. For
example, when compiling on a cluster, with packages in non-standard locations. In that case the following should work:

```bash
mkdir build && cd build

cmake -DCMAKE_BUILD_TYPE=Release -DARMADILLO_INCLUDE_DIR=<path_to_armadillo>/include/ -DARMADILLO_LIBRARY=<path_to_armadillo>/lib64/libarmadillo.so
make
```

If you're having trouble with cmake detecting the correct compiler, try the following.

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=<path_to_compiler> ..
make

```

The location for Boost may need to be specified if it isn't installed in a typical location.

```bash
cmake -DBOOST_ROOT=<path_to_boost> ..
```

You can combine the above as necessary. Earlier versions of the Armadillo library may work, but haven't been tested. If
you need to change the compiler used from the one automatically detected to another, perhaps newer compiler:

```bash
cmake -DCMAKE_CXX_COMPILER=<path_to_executable> ..
```

If using OpenBLAS with CMake older than 3.21, it will not be able to
find lapack.a for static linking. You can tell cmake where to find it.
On Ubuntu 20.04:
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DLAPACK_LIBRARIES=/usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.a ..
```

## Program Execution

Our IBD Mapping process requires several input files. Phenotypes information should be provided for each sample. If any
samples lack phenotype information, they are automatically dropped from analysis. IBD segments identified using
an outside program, e.g., GERMLINE or iLASH, should be processed into our input format. Each breakpoint is a segment
transition point where we either add new segment carriers or end IBD segments and delete carriers of those segments.
Each breakpoint is the unit being tested in permutation. Also required are a genetic map file, and if the user is
providing DASH processed data, an additional info file containing information about the identified segments.

We read each breakpoint, update the data vector according to the additions and deletions for the breakpoint, and then
submit the breakpoint to a process queue to await a thread. We calculate the original value of the test statistic, and
store the permuted test statistics, eventually outputing for all permutations.

## Map-Reduce

Biobank scale datasets can require substantial computation time to complete enough permutations to reach genome-wide
significance. For researchers working in a high performance computing environment we recommend a map-reduce strategy to
enable analysis in reasonable periods of time.

The genome is conveniently broken into 22 autosomal components. Further, each job can be submitted with a subset of the
total desired permutations. For example, if the user wishes to conduct a genome-wide analysis with 100,000 permutations
on each breakpoint, then they can submit ten, 10,000 permutation runs for each chromosome, for a total of 100,000
permutations across all chromosomes. Each job can specify the seed used to initialize the random number generator used
for permutation. Users should ensure that the seed is specified, and different, for each job. Jobs submitted without
specifying the seed will generate a random seed from std::random_device as necessary.

To be clear regarding the seeds used in jobs, the seed should be specified and different for each job on a chromosome.
However, the seeds should be identical between chromosomes. so if you have five jobs for chromosome one, and your seeds
are {1, 2, 3, 4, 5} for each job, then they should also be the same for chromosomes 2-22 as well.

Results can be combined using a python package provided along with carvaIBD. The python package IBDreduce is designed to
use a small amount of memory to collapse all the results from across chromosomes, and across permutation sets within a
chromosome.

## Data Formatting

We support a couple different file formats for IBD data. If the data has been formatted from GERMLINE, hapibd, or iLASH,
then we require the following columns:

```tsv
chr pos add del
```

The add column represents a list of sample pairs with segments starting at the breakpoint. The del column is a list of
sample pairs with segments terminating at the breakpoint.

If the data has been passed through DASH then there are two expected files. The data file should have the following
columns:

```tsv
 chr pos n.cluster n.cluster.haplotype n.cluster.pair n.singleton.pair cluster.add cluster.del singleton.add singleton.del
```

The second file is an info file that contains information about the clusters found by DASH. We expect the following
columns:

```tsv
chr clusterID start end count freq cM
```

freq is the frequency of the cluster in the sample, and cM is the length of the segment. Both are used to filter out
segments. Reasonable defaults are provided for the corresponding --cm and --af options.

Both file formats can be created with scripts provided in the tools directory. They take data from GERMLINE, hap-ibd,
iLASH, or any of the three that have been passed through DASH, and generate our input format.

## Running a fully static build
A script exists at `build-scripts/static-build-on-focal.sh` that should be able to prepare
the environment, and do a statically linked build on Ubuntu 20.04.
