# This build container image ends up about 3GB
FROM ubuntu:20.04 AS build-container

# libboost-python on Ubuntu 20.04 is python 3.8
# liboost-python-dev includes python3.8-dev with it
RUN apt-get update &&  \
    apt-get dist-upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    git \
    intel-mkl \
    libarmadillo-dev \
    libboost-iostreams-dev \
    libboost-numpy-dev \
    libboost-python-dev \
    software-properties-common \
    wget

# Add Kitware Cmake repository to get a newer cmake version
RUN wget -v -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null
RUN apt-add-repository "deb https://apt.kitware.com/ubuntu/ focal main"
RUN apt-get update && \
    apt-get install -y kitware-archive-keyring && \
    rm /etc/apt/trusted.gpg.d/kitware.gpg && \
    apt-get install -y cmake


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
# Maybe it's opened with dlopen and this won't work
# Even weirder, we seem not to need intel mkl. wtf is going on.
RUN apt-get update &&  \
    apt-get -y dist-upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    libboost-iostreams-dev \
    libboost-python-dev

WORKDIR /app
COPY --from=build-container /app/build/carvaIBD carvaIBD
COPY --from=build-container /app/build/ibdlib.so ibdlib.so

ENTRYPOINT ["./carvaIBD"]
