FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
ARG EMSDK_VERSION=3.1.71
ARG VCPKG_COMMIT=84bab45d415d22042bd0b9081aea57f362da3f35
ARG DUCKDB_WASM_NPM_VERSION=1.29.0

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    ccache \
    cmake \
    curl \
    git \
    make \
    ninja-build \
    pkg-config \
    python3 \
    python3-venv \
    nodejs \
    npm \
    unzip \
    zip \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/emsdk \
    && cd /opt/emsdk \
    && curl -sSL "https://github.com/emscripten-core/emsdk/archive/${EMSDK_VERSION}.tar.gz" | tar -xz --strip-components=1 \
    && ./emsdk install ${EMSDK_VERSION} \
    && ./emsdk activate ${EMSDK_VERSION}

RUN mkdir -p /opt/vcpkg \
    && cd /opt/vcpkg \
    && git init \
    && git remote add origin https://github.com/microsoft/vcpkg.git \
    && git fetch --depth=1 origin ${VCPKG_COMMIT} \
    && git checkout FETCH_HEAD \
    && ./bootstrap-vcpkg.sh -disableMetrics

RUN mkdir -p /opt/duckdb-wasm-runtime \
    && cd /opt/duckdb-wasm-runtime \
    && curl -fsSL "https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@${DUCKDB_WASM_NPM_VERSION}/dist/duckdb-browser.mjs" -o duckdb-browser.mjs \
    && curl -fsSL "https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@${DUCKDB_WASM_NPM_VERSION}/dist/duckdb-browser-eh.worker.js" -o duckdb-browser-eh.worker.js \
    && curl -fsSL "https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@${DUCKDB_WASM_NPM_VERSION}/dist/duckdb-eh.wasm" -o duckdb-eh.wasm

ENV EMSDK=/opt/emsdk
ENV VCPKG_ROOT=/opt/vcpkg
ENV VCPKG_TOOLCHAIN_PATH=/opt/vcpkg/scripts/buildsystems/vcpkg.cmake
ENV PATH="/opt/vcpkg:/opt/emsdk:/opt/emsdk/upstream/emscripten:${PATH}"

WORKDIR /work
