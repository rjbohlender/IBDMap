#!/bin/sh

apt-get update
apt-get -y upgrade
apt-get -y cmake build-essential libarmadillo-dev libboost-iostreams-dev intel-mkl
git submodule update --init

cd static-build
mkdir "build"
cd build

cmake -DCMAKE_BUILD_TYPE=Release ..
make -j`nproc`
