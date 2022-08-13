# This build container image ends up about 3GB
FROM ubuntu:22.04 AS build-container

# libboost-python on Ubuntu 22.04 is python 3.10
# liboost-python-dev includes python3.10-dev with it
# We don't need intel-mkl or any other specific BLAS lib because carvaIBD only uses Armadillo headers
# No actual BLAS functions are called
RUN apt-get update &&  \
    apt-get dist-upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    cmake \
    libarmadillo-dev \
    libboost-iostreams-dev \
    libboost-numpy-dev \
    libboost-python-dev

WORKDIR /app
COPY . .
RUN rm -rf build

RUN mkdir build
WORKDIR /app/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
# It'd be smart to detect system threads, but alas
RUN make -j4

# Make a leaner run container without build dependencies
# How much leaner is it? It's about 170MB
FROM ubuntu:22.04 as run-container

# We seem to not need libboost-numpy at runtime?!
RUN apt-get update &&  \
    apt-get -y dist-upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    libboost-iostreams1.74.0 \
    libboost-numpy1.74.0 \
    libboost-python1.74.0 \
    libpython3.10

WORKDIR /app
COPY --from=build-container /app/build/carvaIBD carvaIBD
COPY --from=build-container /app/build/*.so .

ENTRYPOINT ["./carvaIBD"]
