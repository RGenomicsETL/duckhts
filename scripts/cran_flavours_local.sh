#!/usr/bin/env bash
# Run CRAN-flavour checks of one source tarball on this host, in cached images.
#
#   scripts/cran_flavours_local.sh [--git-rev REV | --tarball FILE] [--out DIR]
#                                  [--jobs N] LANE...
#
# Lanes (see scripts/cran-flavours/Dockerfile for what each image is):
#   fedora-gcc16    R-hub Fedora/GCC 16 with R-devel
#   fedora-clang23  Fedora 44, LLVM 23.1.4 with libc++, R-devel r90760
#   nold-fedora     Fedora 44, GCC 16, R-devel r90760 --disable-long-double
#   nold-ubuntu     R-hub noLD image (R-devel r90728)
#   all             every lane above
#
# --git-rev builds the tarball from that commit (for example the revision that
# was published to CRAN) so a known failure can be reproduced before a fix is
# trusted. Without either option the tarball is built from the working tree.
# Every lane writes receipt.txt, dependencies.log, check.out, status.txt and the
# .Rcheck directory under OUT/<lane>/ and the script prints a summary table.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CACHE="${DUCKHTS_FLAVOUR_CACHE:-$HOME/.cache/duckhts-cran-flavours}"
JOBS="${DUCKHTS_JOBS:-5}"
GIT_REV=""
TARBALL=""
OUT=""
LANES=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --git-rev) GIT_REV="$2"; shift 2 ;;
        --tarball) TARBALL="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        --jobs) JOBS="$2"; shift 2 ;;
        -h|--help) sed -n '2,21p' "$0"; exit 0 ;;
        *) LANES+=("$1"); shift ;;
    esac
done
[[ ${#LANES[@]} -gt 0 ]] || { echo "name at least one lane (or all)" >&2; exit 2; }
if [[ " ${LANES[*]} " == *" all "* ]]; then
    LANES=(fedora-gcc16 fedora-clang23 nold-fedora nold-ubuntu)
fi

R_REVISION=90760
declare -A IMAGE=(
    [fedora-clang23]="duckhts-cran:clang23-r${R_REVISION}"
    [nold-fedora]="duckhts-cran:nold-gcc-r${R_REVISION}"
)
declare -A BUILD_ARGS=(
    [fedora-clang23]="--build-arg TOOLCHAIN=clang23 --build-arg R_REVISION=${R_REVISION}"
    [nold-fedora]="--build-arg TOOLCHAIN=gcc --build-arg NOLD=1 --build-arg R_REVISION=${R_REVISION}"
)

ensure_image() {
    local lane="$1"
    if ! docker image inspect "${IMAGE[$lane]}" >/dev/null 2>&1; then
        echo "building ${IMAGE[$lane]} (cached afterwards)"
        # shellcheck disable=SC2086
        docker build -t "${IMAGE[$lane]}" ${BUILD_ARGS[$lane]} "$REPO_ROOT/scripts/cran-flavours"
    fi
}

# One tarball for every lane.
WORK="$(mktemp -d "${TMPDIR:-/tmp}/duckhts-flavours.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
if [[ -z "$TARBALL" ]]; then
    if [[ -n "$GIT_REV" ]]; then
        mkdir -p "$WORK/src"
        git -C "$REPO_ROOT" archive "$GIT_REV" r/Rduckhts | tar -x -C "$WORK/src"
        source_dir="$WORK/src/r/Rduckhts"
    else
        source_dir="$REPO_ROOT/r/Rduckhts"
    fi
    mkdir -p "$WORK/tarball"
    (cd "$WORK/tarball" && R CMD build --no-manual "$source_dir" >"$WORK/build.log" 2>&1) \
        || { cat "$WORK/build.log" >&2; exit 1; }
    TARBALL="$(find "$WORK/tarball" -name 'Rduckhts_*.tar.gz' -print -quit)"
fi
SHA256="$(sha256sum "$TARBALL" | cut -d' ' -f1)"
OUT="${OUT:-$CACHE/results/${SHA256:0:12}}"
mkdir -p "$OUT"
cp -n "$TARBALL" "$OUT/" 2>/dev/null || true
echo "tarball $(basename "$TARBALL") sha256 $SHA256"
echo "results $OUT"

run_lane() {
    local lane="$1"
    local lane_out="$OUT/$lane"
    rm -rf "$lane_out"
    mkdir -p "$lane_out" "$CACHE/rlibs/$lane"
    case "$lane" in
        fedora-clang23|nold-fedora)
            ensure_image "$lane"
            docker run --rm \
                -v "$REPO_ROOT:/repo:ro" -v "$(dirname "$TARBALL"):/in:ro" \
                -v "$lane_out:/out" -v "$CACHE/rlibs/$lane:/root/R-libs" \
                -e DUCKHTS_LANE="$lane" -e DUCKHTS_TARBALL="/in/$(basename "$TARBALL")" \
                -e DUCKHTS_TARBALL_SHA256="$SHA256" -e DUCKHTS_CHECK_OUT=/out \
                -e DUCKHTS_JOBS="$JOBS" -e DUCKHTS_IMAGE_REF="${IMAGE[$lane]}" \
                "${IMAGE[$lane]}" bash /repo/scripts/cran-flavours/lane.sh
            ;;
        fedora-gcc16)
            DUCKHTS_TARBALL="$TARBALL" DUCKHTS_TARBALL_SHA256="$SHA256" \
                DUCKHTS_CHECK_OUT="$lane_out" DUCKHTS_FEDORA_RLIBS="$CACHE/rlibs/$lane" \
                bash "$REPO_ROOT/scripts/check_cran_fedora.sh"
            ;;
        nold-ubuntu)
            DUCKHTS_TARBALL="$TARBALL" DUCKHTS_TARBALL_SHA256="$SHA256" \
                DUCKHTS_CHECK_OUT="$lane_out" DUCKHTS_NOLD_RLIBS="$CACHE/rlibs/$lane" \
                bash "$REPO_ROOT/scripts/check_cran_nold.sh"
            ;;
        *) echo "unknown lane $lane" >&2; return 2 ;;
    esac
}

# Images are built one at a time before the lanes start together.
for lane in "${LANES[@]}"; do
    [[ -n "${IMAGE[$lane]:-}" ]] && ensure_image "$lane"
done
pids=()
for lane in "${LANES[@]}"; do
    # set +e: a failing lane must still record its exit status.
    (set +e; run_lane "$lane" >"$OUT/$lane.log" 2>&1; echo "$?" >"$OUT/$lane.exit") &
    pids+=($!)
done
wait "${pids[@]}" || true

printf '\n%-16s %-6s %s\n' lane exit status
overall=0
for lane in "${LANES[@]}"; do
    code="$(cat "$OUT/$lane.exit" 2>/dev/null || echo '?')"
    status="$(cat "$OUT/$lane/status.txt" 2>/dev/null | tr '\n' ' ')"
    printf '%-16s %-6s %s\n' "$lane" "$code" "${status:-no status (see $OUT/$lane.log)}"
    [[ "$code" == 0 ]] || overall=1
done
exit "$overall"
