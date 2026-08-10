#!/usr/bin/env bash
# Shared cache and provenance helpers for explicit DuckHTS staging scripts.
# Source this file; it deliberately does not perform network I/O.

_duckhts_cache_default() {
    if [[ -n "${DUCKHTS_CACHE_DIR:-}" ]]; then
        printf '%s\n' "$DUCKHTS_CACHE_DIR"
    elif [[ -n "${XDG_CACHE_HOME:-}" ]]; then
        printf '%s/duckhts\n' "$XDG_CACHE_HOME"
    else
        printf '%s/.cache/duckhts\n' "${HOME:?HOME is required when DUCKHTS_CACHE_DIR is unset}"
    fi
}

DUCKHTS_CACHE_DIR="$(_duckhts_cache_default)"
export DUCKHTS_CACHE_DIR

# Return a directory below the common cache root, rejecting an absolute path or
# a parent traversal so staging callers cannot escape the user-selected cache.
duckhts_cache_subdir() {
    local relative="$1"
    if [[ -z "$relative" || "$relative" = /* || "$relative" = *".."* ]]; then
        echo "invalid DuckHTS cache-relative path: $relative" >&2
        return 2
    fi
    printf '%s/%s\n' "${DUCKHTS_CACHE_DIR%/}" "$relative"
}

# Write a small, human-readable provenance record next to a staged artifact.
# Arguments after PATH are FIELD=VALUE pairs. This records acquisition and
# transformation context; it is not a content-digest mechanism.
duckhts_write_provenance() {
    local path="$1"
    shift
    local parent partial pair field value
    parent="$(dirname "$path")"
    mkdir -p "$parent"
    partial="${path}.partial.$$"
    {
        printf 'field\tvalue\n'
        printf 'recorded_at_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        for pair in "$@"; do
            if [[ "$pair" != *=* ]]; then
                echo "provenance field must be FIELD=VALUE: $pair" >&2
                rm -f "$partial"
                return 2
            fi
            field="${pair%%=*}"
            value="${pair#*=}"
            if [[ -z "$field" || "$field" = *$'\t'* || "$field" = *$'\n'* ||
                "$value" = *$'\t'* || "$value" = *$'\n'* ]]; then
                echo "invalid tabular provenance field" >&2
                rm -f "$partial"
                return 2
            fi
            printf '%s\t%s\n' "$field" "$value"
        done
    } >"$partial"
    mv -f "$partial" "$path"
}
