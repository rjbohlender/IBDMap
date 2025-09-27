# This build container image ends up about 3GB
FROM ubuntu:22.04 AS ibdmap-build-container

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# libboost-python on Ubuntu 22.04 is python 3.10
# liboost-python-dev includes python3.10-dev with it
# We don't need intel-mkl or any other specific BLAS lib because carvaIBD only uses Armadillo headers
# No actual BLAS functions are called
RUN apt-get update &&  \
    apt-get dist-upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    cmake \
    git \
    libarmadillo-dev \
    libboost-iostreams-dev \
    libboost-numpy-dev \
    libboost-python-dev

WORKDIR /app
COPY . .
# Add this step to fetch the submodule content

RUN git submodule update --init --recursive
# Make sure their is no previous build directory
RUN rm -rf build

RUN mkdir build
WORKDIR /app/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
# It'd be smart to detect system threads, but alas
RUN make -j4

# We need to install the necessary python dependencies
FROM ubuntu:22.04 AS ibdreduce-build-container

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# We need to install curl
RUN apt-get update && \
    apt-get install -y curl python3.11 python3.11-venv git && \
    rm -rf /var/lib/apt/lists/* 

WORKDIR /app

COPY IBDReduce/requirements.txt .

RUN python3.11 -m venv ibdreduce_venv
RUN . ibdreduce_venv/bin/activate

ENV PATH="/app/ibdreduce_venv/bin:$PATH"
RUN pip3 install --upgrade pip && pip3 install -r requirements.txt

# Make a leaner run container without build dependencies
# How much leaner is it? It's about 170MB
FROM ubuntu:22.04 AS run-container

RUN rm /bin/sh && ln -s /bin/bash /bin/sh
# We seem to not need libboost-numpy at runtime?!
RUN apt-get update &&  \
    apt-get -y dist-upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    libboost-iostreams1.74.0 \
    libboost-numpy1.74.0 \
    libboost-python1.74.0 \
    liblapack3 \
    libpython3.10 \
    python3.11 \
    python3.11-venv

WORKDIR /app
COPY --from=ibdmap-build-container /app/build/ibdmap ibdmap
COPY --from=ibdmap-build-container /app/build/*.so .
COPY --from=ibdreduce-build-container /app/ibdreduce_venv ./ibdreduce_venv
# need to copy the scripts required for ibdreduce
COPY IBDReduce/geneticmap.py .
COPY IBDReduce/ibdreduce_v3.py .

ENV PATH="/app/ibdreduce_venv/bin:$PATH"

ENV PATH="/app:$PATH"
