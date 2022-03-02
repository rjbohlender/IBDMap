# This build container image ends up about 3GB
FROM ubuntu:20.04 AS build-container

# libboost-python on Ubuntu 20.04 is python 3.8

RUN apt-get update &&  \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    cmake \
    git \
    intel-mkl \
    libarmadillo-dev \
    libboost-iostreams-dev \
    libboost-numpy-dev \
    libboost-python-dev \
    python3.8-dev


WORKDIR /app
COPY . .
RUN git clean -xfd
RUN git submodule update --init

RUN mkdir build
WORKDIR /app/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
# It'd be smart to detect system threads, but alas
RUN make -j8

# Make a leaner run container without build dependencies
# How much leaner is it? It's about 2.3GB. Most of this size is Intel MKL
FROM ubuntu:20.04 as run-container

# We seem to not need libboost-numpy at runtime?!
RUN apt-get update &&  \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    intel-mkl \
    libboost-iostreams-dev \
    libboost-python-dev \
    python3.8-dev

WORKDIR /app
COPY --from=build-container /app/build/carvaIBD carvaIBD

ENTRYPOINT ["./carvaIBD"]
