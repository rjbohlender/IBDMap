# This build container image ends up about 3GB
FROM ubuntu:22.04 AS carva-build-container

# libboost-python on Ubuntu 22.04 is python 3.10
# liboost-python-dev includes python3.10-dev with it
RUN apt-get update &&  \
    apt-get dist-upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    cmake \
    libarmadillo-dev \
    libboost-iostreams-dev \
    libboost-numpy-dev \
    libboost-python-dev \
    liblapack-dev

WORKDIR /app
COPY . .
RUN rm -rf build

RUN mkdir build
WORKDIR /app/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
# It'd be smart to detect system threads, but alas
RUN make -j4

# Now we can create the python venv with ibdreduce
FROM ubuntu:22.04 as python-build-container
# changing the working directory to be app
WORKDIR /app/

# We need to install curl to download poetry
RUN apt-get update && \ 
    apt-get -y dist-upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    curl \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/* 

# Now we can configure our environment to install things with poetry
ARG BUILD_ENV

ENV BUILD_ENV=${BUILD_ENV} \
PIP_NO_CACHE_DIR=1 \
POETRY_VERSION=1.8.3 \
POETRY_NO_INTERACTION=1 \
POETRY_VIRTUALENVS_CREATE=false \
POETRY_CACHE_DIR='/var/cache/pypoetry' \
POETRY_HOME='/usr/local'

# Install poetry
RUN curl -sSL https://install.python-poetry.org | python3 - \
    && poetry self add poetry-plugin-bundle 



# Copy the requirements file into the container
COPY ./IBDReduce/ibdreduce ./ibdreduce
COPY ./IBDReduce/pyproject.toml /app/pyproject.toml
COPY ./IBDReduce/poetry.lock /app/poetry.lock
COPY ./IBDReduce/README.md /app/README.md

# Use poetry to install dependencies into the virtualenv
RUN poetry bundle venv /opt/venv/ 

# Make a leaner run container without build dependencies
# How much leaner is it? It's about 170MB
FROM ubuntu:22.04 as run-container

LABEL maintainer="belowlab"

# We seem to not need libboost-numpy at runtime?!
RUN apt-get update &&  \
    apt-get -y dist-upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    libboost-iostreams1.74.0 \
    libboost-numpy1.74.0 \
    libboost-python1.74.0 \
    libpython3.10 \
    liblapack3 \
    python3 \
    && rm -rf /var/lib/apt/lists/* 

WORKDIR /app

COPY --from=carva-build-container /app/build/ibdmap ibdmap
# COPY --from=carva-build-container /app/build/*.so .
COPY --from=python-build-container /opt/venv /opt/venv

# lets add IBDMap and /opt/venv/ to the path
ENV PATH="/app:$PATH"
ENV PATH="/opt/venv/bin:$PATH"


