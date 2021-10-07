# carvaIBD

carvaIBD is a multithreaded scalable application for IBD mapping in Biobank
scale datasets.

## Compiling

This project has submodules that need to be initialized after cloning the directory from github. 

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

If using OpenBLAS with CMake older than 3.21, it will not be able to 
find lapack.a for static linking. You can tell cmake where to find it. 
On Ubuntu 20.04:
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DLAPACK_LIBRARIES=/usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.a ..
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

To be clear regarding the seeds used in jobs, the seed should be specified and different
for each job on a chromosome. However, the seeds should be identical between chromosomes.
so if you have five jobs for chromosome one, and your seeds are {1, 2, 3, 4, 5} for
each job, then they should also be the same for chromosomes 2-22 as well.

Results can be combined using a python package provided along with carvaIBD. The
python package IBDreduce is designed to use a small amount of memory to collapse
all the results from across chromosomes, and across permutation sets within a
chromosome.

## Running a fully static build
A script exists at `build-scripts/static-build-on-focal.sh` that should be able to prepare
the environment, and do a statically linked build on Ubuntu 20.04.
