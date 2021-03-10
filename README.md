# carvaIBD

carvaIBD is a multithreaded scalable application for IBD mapping in Biobank
scale datasets.

## Compiling

This project has submodules that need to be initialized after cloning the directory from github. 

```bash
git submodule init
git submodule update
```

Dependencies:
    - C++ compiler supporting C++17
    - Armadillo: version >= 8.6
    - Boost C++ Library
    
Create a build directory and run CMAKE and make from within the build directory:

```bash
mkdir build && cd build

cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

If you're working in a high performance computing environment or another environment
where packages may not be installed in standard places, you can tell CMake where to
look for both Armadillo and Boost.

If cmake fails to detect Armadillo, but you're sure it is available,
you may need to direct cmake to the library, e.g., when compiling on a
cluster, with packages in non-standard locations. In that case the
following should work:

```bash
mkdir build && cd build

cmake -DCMAKE_BUILD_TYPE=Release -DARMADILLO_INCLUDE_DIR=<path_to_armadillo>/include/ -DARMADILLO_LIBRARY=<path_to_armadillo>/lib64/libarmadillo.so
make
```

If you're having trouble with cmake detecting the correct compiler.

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=<path_to_compiler> ..
make

```

The location for Boost may need to be specified if it isn't installed in a
typical location.

```bash
cmake -DBOOST_ROOT=<path_to_boost> ..
```

You can combine the above as necessary. Earlier versions of the
Armadillo library may work, but haven't been tested. If you need to
change the compiler used from the one automatically detected to
another, perhaps newer compiler:

```bash
cmake -DCMAKE_CXX_COMPILER=<path_to_executable> ..
```

## Map-Reduce

Biobank scale datasets can require substantial computation time to complete
enough permutations to reach genome-wide significance. For researchers working
in a high performance computing environment we recommend a map-reduce strategy
to enable analysis in reasonable periods of time. 

The genome is conveniently broken into 22 autosomal components. Further, each
job can be submitted with a subset of the total desired permutations. For
example, if the user wishes to conduct a genome-wide analysis with 100,000
permutations on each breakpoint, then they can submit ten, 10,000 permutation
runs for each chromosome, for a total of 100,000 permutations across all
chromosomes. Each job can specify the seed used to initialize the random number
generator used for permutation. Users should ensure that the seed is specified,
and different, for each job. Jobs submitted without specifying the seed will
generate a random seed from std::random_device as necessary.

Results can be combined using a python package provided along with carvaIBD. The
python package IBDreduce is designed to use a small amount of memory to collapse
all the results from across chromosomes, and across permutation sets within a
chromosome.

## Data Formatting

We support a couple different file formats for IBD data. If the data has been formatted from GERMLINE, hapibd, or iLASH, then the following columns are required:

```tsv
chr pos segments pairs add del
```

If the data has been passed through DASH then there are two expected files. The data file should have the following columns:

```tsv
 chr pos n.cluster n.cluster.haplotype n.cluster.pair n.singleton.pair cluster.add cluster.del singleton.add singleton.del
```

The second file is an info file that contains information about the clusters found by DASH. The following columns are expected:

```tsv
chr clusterID start end count
```

Additionally there are two optional columns that are required if the corresponding filtering arguments are used. The
filtering arguments are -af which filters on cluster frequency, and -cM which filters on cluster length. The columns
are:

```tsv
freq cM
```

