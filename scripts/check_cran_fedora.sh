#!/usr/bin/env bash
# Reproduce CRAN-like R package checks on Fedora.
#
# Local use launches the same container/profile that GitHub Actions uses:
#   scripts/check_cran_fedora.sh          # R-hub Fedora/GCC 16
#   scripts/check_cran_fedora_clang.sh    # Fedora 44 / Clang 22 CRAN-like
#
# The Clang profile tracks the relevant public CRAN compiler configuration at
# https://www.stats.ox.ac.uk/pub/bdr/Rconfig/r-devel-linux-x86_64-fedora-clang.
# It deliberately uses Fedora's packaged R, not CRAN's upstream R-devel
# runtime. The job is therefore CRAN-like rather than an exact CRAN clone.
set -euo pipefail

PROFILE="${DUCKHTS_FEDORA_PROFILE:-gcc}"
case "$PROFILE" in
    gcc)
        DEFAULT_IMAGE="ghcr.io/r-hub/containers/gcc16:latest"
        ;;
    clang)
        # Fedora 44 image pinned on 2026-07-25. The exact Fedora R and Clang
        # packages used below are pinned and verified after installation.
        DEFAULT_IMAGE="fedora@sha256:6c75d5bf57cb0fa5aa4b92c6a83c86c791644496d9ac230de7711f5b8ec3b898"
        ;;
    *)
        echo "ERROR: DUCKHTS_FEDORA_PROFILE must be gcc or clang (got '$PROFILE')" >&2
        exit 1
        ;;
esac
IMAGE="${DUCKHTS_FEDORA_IMAGE:-$DEFAULT_IMAGE}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ "${DUCKHTS_FEDORA_CONTAINER:-0}" != "1" ]]; then
    RLIBS_CACHE="${DUCKHTS_FEDORA_RLIBS:-${HOME}/.cache/duckhts-fedora-${PROFILE}-r-libs}"
    mkdir -p "$RLIBS_CACHE"
    exec docker run --rm \
        -v "$REPO_ROOT:/repo" \
        -v "$RLIBS_CACHE:/root/R-libs" \
        -w /repo \
        -e DUCKHTS_FEDORA_CONTAINER=1 \
        -e DUCKHTS_FEDORA_PROFILE="$PROFILE" \
        "$IMAGE" \
        bash scripts/check_cran_fedora.sh
fi

DNF_PACKAGES=(
    bzip2-devel
    cmake
    curl
    libcurl-devel
    libdeflate-devel
    openssl-devel
    pandoc
    patch
    pkgconf-pkg-config
    xz-devel
    zlib-devel
)
dnf install -y --setopt=install_weak_deps=False "${DNF_PACKAGES[@]}"

stage_pinned_rpm() {
    local url="$1"
    local sha256="$2"
    local destination="$3/$(basename "$url")"

    curl --fail --location --retry 3 --connect-timeout 15 --silent --show-error \
        --output "$destination" "$url"
    echo "${sha256}  ${destination}" | sha256sum -c -
}

if [[ "$PROFILE" = "clang" ]]; then
    # Pin the complete R/Clang package closure to immutable Koji artifacts.
    # The Fedora image fixes the starting filesystem; the SHA-256 checks below
    # prevent live Fedora repository metadata from changing the R runtime or
    # compiler selected by this CRAN-like profile.
    FEDORA_R_DEVEL_NEVRA="R-devel-4.6.1-1.fc44"
    FEDORA_CLANG_NEVRA="clang-22.1.8-4.fc44"
    FEDORA_R_KOJI_BASE="https://kojipkgs.fedoraproject.org/packages/R/4.6.1/1.fc44/x86_64"
    # Fedora builds Clang 22.1.8 as part of the LLVM source package; its
    # immutable RPMs are therefore under /packages/llvm/, not /packages/clang/.
    FEDORA_LLVM_KOJI_BASE="https://kojipkgs.fedoraproject.org/packages/llvm/22.1.8/4.fc44/x86_64"
    PINNED_RPM_DIR="$(mktemp -d /tmp/duckhts-fedora-rpms.XXXXXX)"

    stage_pinned_rpm "$FEDORA_R_KOJI_BASE/R-4.6.1-1.fc44.x86_64.rpm" \
        "51dd70927b3dc40f71e02da165ff38e65ab7e465c80e5b753de624c14466dc32" "$PINNED_RPM_DIR"
    stage_pinned_rpm "$FEDORA_R_KOJI_BASE/R-core-4.6.1-1.fc44.x86_64.rpm" \
        "ef5274cb6dbcc0cb9e3f8c275243403363c93e2ff1efeaa9bbb0c5887a208908" "$PINNED_RPM_DIR"
    stage_pinned_rpm "$FEDORA_R_KOJI_BASE/R-core-devel-4.6.1-1.fc44.x86_64.rpm" \
        "6f27f80d97a3f3fb5ee018dc98cbfdec7cdbf5e67a17d30f7df1ce3d8c1eeacd" "$PINNED_RPM_DIR"
    stage_pinned_rpm "$FEDORA_R_KOJI_BASE/R-devel-4.6.1-1.fc44.x86_64.rpm" \
        "b635fd75d7f051a702bda2d96e89d9a5e1767ef5eb5eaa891f5dccbf7196fe2f" "$PINNED_RPM_DIR"
    stage_pinned_rpm "$FEDORA_R_KOJI_BASE/R-java-4.6.1-1.fc44.x86_64.rpm" \
        "0a4ba1745ce6a1f6ac663ad1ba776e07e690b17fe96ba5fc7171abce279e3fd7" "$PINNED_RPM_DIR"
    stage_pinned_rpm "$FEDORA_R_KOJI_BASE/R-java-devel-4.6.1-1.fc44.x86_64.rpm" \
        "ad5a20e060a80f6aeff9934f815c02ed3886b28c8ffb366df1cd9faf213ef9d0" "$PINNED_RPM_DIR"
    stage_pinned_rpm "$FEDORA_R_KOJI_BASE/libRmath-4.6.1-1.fc44.x86_64.rpm" \
        "07e7e105afe258065f7b3ee9b502d1f71bc92c381856b9f69a96e2db0c40c2f9" "$PINNED_RPM_DIR"
    stage_pinned_rpm "$FEDORA_R_KOJI_BASE/libRmath-devel-4.6.1-1.fc44.x86_64.rpm" \
        "7c8d08b156604f2e0bf1317988846b4f44663aa5a724e7ca78806018c643f227" "$PINNED_RPM_DIR"
    stage_pinned_rpm "$FEDORA_LLVM_KOJI_BASE/clang-22.1.8-4.fc44.x86_64.rpm" \
        "02e9a1f927a3b3e363035ae91978deda302fb81bc81da4e5113eca4ce60b07fa" "$PINNED_RPM_DIR"
    stage_pinned_rpm "$FEDORA_LLVM_KOJI_BASE/clang-libs-22.1.8-4.fc44.x86_64.rpm" \
        "cb61748b1a32d23b7f770c36b1301d89da775105712b9f05d12b16fd526e9da4" "$PINNED_RPM_DIR"
    stage_pinned_rpm "$FEDORA_LLVM_KOJI_BASE/clang-resource-filesystem-22.1.8-4.fc44.x86_64.rpm" \
        "22edf155d2ff9c4d24177a8c99e7066d12da6368e793d0f8bc075c6b10bf2b72" "$PINNED_RPM_DIR"
    stage_pinned_rpm "$FEDORA_LLVM_KOJI_BASE/llvm-filesystem-22.1.8-4.fc44.x86_64.rpm" \
        "a427bb7d241f20dbd1cc759c63786107f7ba308425d91971811a4d2474453f6d" "$PINNED_RPM_DIR"
    stage_pinned_rpm "$FEDORA_LLVM_KOJI_BASE/llvm-libs-22.1.8-4.fc44.x86_64.rpm" \
        "094f63a7991a92d9ccaeb3984e8609986277dc43e12a613c01f8c69dbb10ece5" "$PINNED_RPM_DIR"

    dnf install -y --setopt=install_weak_deps=False \
        "$PINNED_RPM_DIR"/*.rpm \
        devscripts-checkbashisms \
        glibc-langpack-en \
        make
    rm -rf "$PINNED_RPM_DIR"
    test "$(rpm -q --qf '%{NAME}-%{VERSION}-%{RELEASE}' R-devel)" = "$FEDORA_R_DEVEL_NEVRA"
    test "$(rpm -q --qf '%{NAME}-%{VERSION}-%{RELEASE}' clang)" = "$FEDORA_CLANG_NEVRA"
    export DUCKHTS_CC="${DUCKHTS_CC:-clang}"
    # Fedora's packaged R was built with GCC and supplies GCC-only -specs
    # options in R CMD config CFLAGS. Use the public CRAN Clang C flags when
    # compiling DuckHTS itself rather than treating those foreign options as a
    # project warning.
    export DUCKHTS_CFLAGS="${DUCKHTS_CFLAGS:--O3 -Wall -pedantic -Wp,-D_FORTIFY_SOURCE=3}"
    export DUCKHTS_LDFLAGS="${DUCKHTS_LDFLAGS-}"
    export DUCKHTS_STRICT_WARNINGS=1
    export LANG=en_US.UTF-8
    export LC_ALL=en_US.UTF-8
    # These are the CRAN Fedora-Clang compilation/check controls relevant to
    # an extension package.  --as-cran supplies the standard CRAN checks.
    export _R_CHECK_COMPILATION_FLAGS_=true
    export _R_CHECK_INSTALL_DEPENDS_=true
    export _R_CHECK_SUGGESTS_ONLY_=true
    export _R_CHECK_NO_RECOMMENDED_=true
    export _R_CHECK_DEPRECATED_DEFUNCT_=true
    export _R_CHECK_FF_CALLS_=registration
    export _R_CHECK_PRAGMAS_=true
    "${DUCKHTS_CC}" --version
fi
R --version | head -n 1
R CMD config CC

mkdir -p /root/R-libs
export R_LIBS_USER=/root/R-libs

DEPS_LOG=/tmp/rduckhts-fedora-dependencies.log
if ! R -q -e '
    packages <- c("DBI", "duckdb", "Rtinycc", "tinytest")
    installed <- rownames(installed.packages(lib.loc = Sys.getenv("R_LIBS_USER")))
    needed <- setdiff(packages, installed)
    if (length(needed)) {
        install.packages(
            needed,
            lib = Sys.getenv("R_LIBS_USER"),
            repos = "https://cloud.r-project.org"
        )
    }
' >"$DEPS_LOG" 2>&1; then
    tail -n 200 "$DEPS_LOG" >&2
    exit 1
fi
tail -n 20 "$DEPS_LOG"

CHECK_ROOT="$(mktemp -d /tmp/rduckhts-fedora-check.XXXXXX)"
cd "$CHECK_ROOT"
R CMD build --no-manual "$REPO_ROOT/r/Rduckhts"
PACKAGE_TARBALL="$(find "$CHECK_ROOT" -maxdepth 1 -name 'Rduckhts_*.tar.gz' -print -quit)"
test -n "$PACKAGE_TARBALL"
set +e
R CMD check --as-cran --no-manual "$PACKAGE_TARBALL"
check_status=$?
set -e

CHECK_LOG="$CHECK_ROOT/Rduckhts.Rcheck/00check.log"
test -f "$CHECK_LOG"
if [ "$check_status" -ne 0 ]; then
    INSTALL_LOG="$CHECK_ROOT/Rduckhts.Rcheck/00install.out"
    if [ -f "$INSTALL_LOG" ]; then
        echo "--- $INSTALL_LOG (last 240 lines) ---" >&2
        tail -n 240 "$INSTALL_LOG" >&2
    fi
fi
if grep -Eq '^Status:.*(ERROR|WARNING)' "$CHECK_LOG"; then
    echo "Fedora ${PROFILE} check produced an ERROR or WARNING" >&2
    exit 1
fi
exit "$check_status"
