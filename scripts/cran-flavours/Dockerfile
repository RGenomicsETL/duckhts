# syntax=docker/dockerfile:1
# Local CRAN-flavour images: Fedora 44 with R-devel built from a pinned revision.
#
#   TOOLCHAIN=clang23  LLVM 23.1.4 (release binary) with libc++, as CRAN's
#                      r-devel-linux-x86_64-fedora-clang flavour
#   TOOLCHAIN=gcc      Fedora's GCC 16
#   NOLD=1             R configured with --disable-long-double (CRAN's noLD)
#
# Configuration source: https://www.stats.ox.ac.uk/pub/bdr/Rconfig/
# Deviations from CRAN, recorded in the receipt of every run: gfortran supplies
# Fortran (CRAN uses flang 23 for the clang flavour), no recommended packages,
# no -mtune=native, no LTO, and a pinned Fedora image, not a rolling one.
ARG FEDORA=fedora@sha256:6c75d5bf57cb0fa5aa4b92c6a83c86c791644496d9ac230de7711f5b8ec3b898
FROM ${FEDORA}

ARG TOOLCHAIN=gcc
ARG NOLD=0
ARG R_REVISION=90760
ARG LLVM_URL=https://github.com/llvm/llvm-project/releases/download/llvmorg-23.1.4/LLVM-23.1.4-Linux-X64.tar.zst
ARG LLVM_SHA256=efdd1dbd089b99bfcfcde7b8eb1236163f2e88df5a981cca4f783ec9cb5e0378

RUN dnf install -y --setopt=install_weak_deps=False \
        bzip2-devel cmake curl diffutils file findutils gcc gcc-c++ gcc-gfortran \
        glibc-langpack-en libcurl-devel libdeflate-devel libicu-devel make \
        openssl-devel pandoc patch pcre2-devel perl pkgconf-pkg-config qpdf \
        readline-devel subversion tar texinfo unzip which xz xz-devel zlib-devel zstd \
    && dnf clean all

RUN if [ "$TOOLCHAIN" = "clang23" ]; then \
        curl --fail --location --retry 3 --silent --show-error -o /tmp/llvm.tar.zst "$LLVM_URL" \
        && echo "${LLVM_SHA256}  /tmp/llvm.tar.zst" | sha256sum -c - \
        && mkdir -p /opt/llvm && tar --use-compress-program="zstd -d --long=31" -xf /tmp/llvm.tar.zst -C /opt/llvm --strip-components=1 \
        && rm /tmp/llvm.tar.zst; \
    fi

ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

# R-devel at a pinned Subversion revision, configured like CRAN's flavour.
RUN set -eu; \
    svn checkout -q -r "$R_REVISION" https://svn.r-project.org/R/trunk /src/R; \
    cd /src/R; \
    flags=""; [ "$NOLD" = "1" ] && flags="--disable-long-double"; \
    if [ "$TOOLCHAIN" = "clang23" ]; then \
        llvm=/opt/llvm; \
        export CC="$llvm/bin/clang" CXX="$llvm/bin/clang++ -stdlib=libc++" OBJC="$llvm/bin/clang" FC=gfortran \
            CFLAGS="-O3 -Wall -pedantic -Wp,-D_FORTIFY_SOURCE=3" \
            C17FLAGS="-O3 -Wall -pedantic -Wno-strict-prototypes -Wp,-D_FORTIFY_SOURCE=3" \
            FFLAGS="-O2 -pedantic" \
            CXXFLAGS="-O3 -Wall -pedantic -frtti -Wp,-D_FORTIFY_SOURCE=3 -Wno-missing-template-arg-list-after-template-kw" \
            LDFLAGS="-L$llvm/lib/x86_64-unknown-linux-gnu -Wl,-rpath,$llvm/lib/x86_64-unknown-linux-gnu"; \
        lapack="--without-lapack"; \
    else \
        export CFLAGS="-g -O2 -Wall -pedantic -Werror=format-security -Wp,-D_FORTIFY_SOURCE=3 -fexceptions -fstack-protector-strong -fstack-clash-protection -fcf-protection -Werror=implicit-function-declaration -Wstrict-prototypes" \
            CXXFLAGS="-g -O2 -Wall -pedantic -Wno-ignored-attributes -Wno-parentheses -Werror=format-security -Wp,-D_FORTIFY_SOURCE=3 -fexceptions -fstack-protector-strong -fstack-clash-protection -fcf-protection" \
            FFLAGS="-g -O2 -Wall -pedantic"; \
        lapack=""; \
    fi; \
    ./configure --prefix=/opt/R --without-x --without-tcltk --without-recommended-packages \
        --without-cairo --disable-java $lapack $flags; \
    make -j"$(nproc)"; \
    make install; \
    rm -rf /src/R

# Late layer, so adding a tool does not rebuild R: R CMD check --as-cran warns
# "a complete check needs the 'checkbashisms' script" without it.
RUN dnf install -y --setopt=install_weak_deps=False devscripts-checkbashisms && dnf clean all

ENV PATH=/opt/R/bin:/opt/llvm/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/llvm/lib/x86_64-unknown-linux-gnu
RUN Rscript --vanilla -e 'cat(R.version.string, "\n"); print(capabilities()["long.double"])'
