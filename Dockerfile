FROM ubuntu:20.04

LABEL org.opencontainers.image.source=https://github.com/precimed/simu
LABEL org.opencontainers.image.description="Simu-Linux"
LABEL org.opencontainers.image.licenses=GPL-3.0

ENV TZ=Europe
ENV DEBIAN_FRONTEND noninteractive

WORKDIR /tmp

RUN apt-get update && apt-get install -y ca-certificates && \
    update-ca-certificates && \
    apt-get install -y git build-essential libboost-all-dev cmake && \
    rm -rf /var/lib/apt/lists/*

RUN git clone --recurse-submodules https://github.com/precimed/simu && \
    cd simu && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install && \
    cd /tmp && \
    rm -rf simu

WORKDIR /
ENTRYPOINT ["/usr/local/bin/simu"]