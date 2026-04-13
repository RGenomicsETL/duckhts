FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
ARG EMSDK_VERSION=3.1.68
ARG VCPKG_COMMIT=84bab45d415d22042bd0b9081aea57f362da3f35

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

# Ubuntu 24.04 ships ninja 1.11.1; newer vcpkg requires >= 1.13.1 and tries to
# download it at cmake-configure time (inside the container, where network may
# be restricted).  Pre-install the binary so vcpkg finds it without downloading.
RUN curl -fsSL https://github.com/ninja-build/ninja/releases/download/v1.13.1/ninja-linux.zip \
      -o /tmp/ninja.zip \
    && unzip /tmp/ninja.zip ninja -d /usr/local/bin \
    && rm /tmp/ninja.zip \
    && chmod +x /usr/local/bin/ninja \
    && ninja --version

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

ENV EMSDK=/opt/emsdk
ENV VCPKG_ROOT=/opt/vcpkg
ENV VCPKG_TOOLCHAIN_PATH=/opt/vcpkg/scripts/buildsystems/vcpkg.cmake
ENV PATH="/opt/vcpkg:/opt/emsdk:/opt/emsdk/upstream/emscripten:${PATH}"

WORKDIR /work
