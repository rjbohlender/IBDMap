# This build container image ends up about 3GB
FROM ubuntu:20.04 AS build-container

# libboost-python on Ubuntu 20.04 is python 3.8
# liboost-python-dev includes python3.8-dev with it
# We don't need intel-mkl because carvaIBD only uses Armadillo headers
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

## Add Kitware Cmake repository to get a newer cmake version
#RUN apt-get install -y \
#    software-properties-common \
#    wget
#RUN wget -v -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null
#RUN apt-add-repository "deb https://apt.kitware.com/ubuntu/ focal main"
#RUN apt-get install -y kitware-archive-keyring && \
#    rm /etc/apt/trusted.gpg.d/kitware.gpg && \
#    apt-get install -y cmake


WORKDIR /app
COPY . .
RUN rm -rf build
# RUN git submodule update --init

RUN mkdir build
WORKDIR /app/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
# It'd be smart to detect system threads, but alas
RUN make -j4

# Make a leaner run container without build dependencies
# How much leaner is it? It's about 170MB
FROM ubuntu:20.04 as run-container

# We seem to not need libboost-numpy at runtime?!
# Maybe it's opened with dlopen and this won't work
# Even weirder, we seem not to need intel mkl. wtf is going on.
RUN apt-get update &&  \
    apt-get -y dist-upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    libboost-iostreams1.71.0 \
    libboost-numpy1.71.0 \
    libboost-python1.71.0 \
    libpython3.8

WORKDIR /app
COPY --from=build-container /app/build/carvaIBD carvaIBD
COPY --from=build-container /app/build/*.so .

## Section to build a tar of executable, library, and all libs they need
#RUN mkdir libs
#RUN ldd carvaIBD | grep "=> /" | awk '{print $3}' | xargs -I '{}' cp -nv '{}' libs
#RUN ldd ibdlib.so | grep "=> /" | awk '{print $3}' | xargs -I '{}' cp -nv '{}' libs
#
#WORKDIR /
#RUN tar -czf ibdmap.tar.gz app/

ENTRYPOINT ["./carvaIBD"]
