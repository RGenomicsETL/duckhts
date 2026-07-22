cmake_minimum_required(VERSION 3.5...3.29)
option(DUCKDB_WASM_EXTENSION "Whether compiling for Wasm target" OFF)

###
# Configuration
###
if(NOT DEFINED EXTENSION_NAME)
    message(FATAL_ERROR "DuckDB extension name is required")
endif()

add_definitions(-DDUCKDB_EXTENSION_NAME=${EXTENSION_NAME})

# ------------------------------------------------------------------
# macOS target architecture (must be set before project()).
#
# Prefer the vcpkg target triplet when provided to avoid host-arch leakage
# during cross-compilation on Apple Silicon runners.
# ------------------------------------------------------------------
if(APPLE)
    if(NOT DEFINED OSX_BUILD_ARCH OR "${OSX_BUILD_ARCH}" STREQUAL "")
        if(DEFINED VCPKG_TARGET_TRIPLET)
            if(VCPKG_TARGET_TRIPLET MATCHES "^x64-osx")
                set(OSX_BUILD_ARCH "x86_64")
            elseif(VCPKG_TARGET_TRIPLET MATCHES "^arm64-osx")
                set(OSX_BUILD_ARCH "arm64")
            endif()
        endif()
        if((NOT DEFINED OSX_BUILD_ARCH OR "${OSX_BUILD_ARCH}" STREQUAL "") AND
           DEFINED CMAKE_OSX_ARCHITECTURES AND NOT "${CMAKE_OSX_ARCHITECTURES}" STREQUAL "")
            list(GET CMAKE_OSX_ARCHITECTURES 0 OSX_BUILD_ARCH)
        endif()
    endif()
    if(DEFINED OSX_BUILD_ARCH AND NOT "${OSX_BUILD_ARCH}" STREQUAL "")
        set(CMAKE_OSX_ARCHITECTURES "${OSX_BUILD_ARCH}" CACHE STRING "" FORCE)
        message(STATUS "macOS target arch: ${OSX_BUILD_ARCH}")
    endif()
endif()

# Avoid "manually-specified variables were not used" warnings when the
# CI passes TARGET_DUCKDB_VERSION_* for header fetching but the stable API
# path doesn't consume them.
if(DEFINED TARGET_DUCKDB_VERSION_MAJOR OR DEFINED TARGET_DUCKDB_VERSION_MINOR OR DEFINED TARGET_DUCKDB_VERSION_PATCH)
    set(_duckdb_target_release_version
        "${TARGET_DUCKDB_VERSION_MAJOR}.${TARGET_DUCKDB_VERSION_MINOR}.${TARGET_DUCKDB_VERSION_PATCH}"
        CACHE INTERNAL "DuckDB release version passed by CI"
    )
endif()

# Only override the C API version when building with the unstable API.
# For stable API builds, the header's own defaults (currently v1.2.0) are correct.
# TARGET_DUCKDB_VERSION_* is the DuckDB *release* version (for header fetching),
# which is different from the C extension API version.
if (DEFINED DUCKDB_EXTENSION_API_VERSION_UNSTABLE)
    add_definitions(-DDUCKDB_EXTENSION_API_VERSION_UNSTABLE=${DUCKDB_EXTENSION_API_VERSION_UNSTABLE})
    if (DEFINED TARGET_DUCKDB_VERSION_MAJOR)
        add_definitions(-DDUCKDB_EXTENSION_API_VERSION_MAJOR=${TARGET_DUCKDB_VERSION_MAJOR})
    endif()
    if (DEFINED TARGET_DUCKDB_VERSION_MINOR)
        add_definitions(-DDUCKDB_EXTENSION_API_VERSION_MINOR=${TARGET_DUCKDB_VERSION_MINOR})
    endif()
    if (DEFINED TARGET_DUCKDB_VERSION_PATCH)
        add_definitions(-DDUCKDB_EXTENSION_API_VERSION_PATCH=${TARGET_DUCKDB_VERSION_PATCH})
    endif()
endif()

###
# Build
###
project(${EXTENSION_NAME} LANGUAGES C)

# Create Extension library
set(EXTENSION_SOURCES
        src/duckhts.c
        src/simd/duckhts_simd_dispatch.c
        src/simd/duckhts_simd_scalar.c
        src/simd/duckhts_simd_avx2.c
        src/simd/duckhts_simd_avx512.c
        src/simd/duckhts_simd_neon.c
        src/simd/duckhts_simd_wasm_simd128.c
        src/bcf_reader.c
        src/bam_reader.c
        src/bam_pileup.c
        src/bgzip.c
        src/hts_index_builder.c
        src/liftover_udf.c
        src/bcftools_norm_udf.c
        src/munge_udf.c
        src/seq_reader.c
        src/fastq_qc.c
        src/interval_udf.c
        src/bigwig_reader.c
        src/libbigwig_hfile_io.c
        src/kmer_udf.c
        src/tabix_reader.c
        src/vep_parser.c
        src/duckvep/kernel/src/duckvep_kernel.c
        src/duckvep/kernel/src/duckvep_so.c
        src/duckvep/kernel/src/duckvep_sweep.c
        src/duckvep/kernel/src/duckvep_classify.c
        src/duckvep/kernel/src/duckvep_effect.c
        src/duckvep/kernel/src/duckvep_sv.c
        src/duckvep/kernel/src/duckvep_delta.c
        src/duckvep/kernel/src/duckvep_projection.c
        src/duckvep/kernel/src/duckvep_transcript_edit.c
        src/duckvep/kernel/src/duckvep_hgvs.c
        src/duckvep/kernel/src/duckvep_codon.c
        src/duckvep/kernel/src/duckvep_coding.c
        src/duckvep/kernel/src/duckvep_haplotype.c
        src/duckvep/duckvep_variant_tile.c
        src/duckvep/duckvep_model.c
        src/duckvep/duckvep_annotate.c
        src/duckvep/duckvep_ensembl.c
        src/duckvep/duckvep_sql.c
        src/hts_meta_reader.c
        src/quality_encoding.c
        src/quality_encoding_reader.c
        src/score_udf.c
        src/mosdepth_table.c
        src/bam_bin_counts.c
        src/samtools_idxstats_table.c
        src/bam_bed_coverage.c
        src/cgranges_api.c
        src/variantkey_udf.c
        third_party/cgranges/cgranges.c
        third_party/libBigWig/bwRead.c
        third_party/libBigWig/bwValues.c
        src/bcftools_filter.c
        src/bcftools_shim.c
        src/wasm_http_hfile.c
)

if (DUCKDB_WASM_EXTENSION)
	# Keep SIMD128 specific to the DuckDB-Wasm runtime.  webR and native builds
	# compile this translation unit as a stub unless their own build contract
	# explicitly opts in; a global -msimd128 would silently change those ABIs.
	set_property(SOURCE src/simd/duckhts_simd_wasm_simd128.c APPEND PROPERTY
		COMPILE_OPTIONS -msimd128)
	add_library(${EXTENSION_NAME} STATIC ${EXTENSION_SOURCES})
	# C files must use the same setjmp/longjmp mechanism as the DuckDB WASM
	# runtime (__wasm_longjmp / __wasm_setjmp).  Without -fwasm-exceptions the
	# C compiler falls back to the old emscripten_longjmp import, which DuckDB
	# WASM does not export, causing an "indirect call signature mismatch" trap
	# at extension load time.
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fwasm-exceptions")
	target_compile_definitions(${EXTENSION_NAME} PRIVATE DUCKHTS_WASM_DUCKDB_RUNTIME=1)
	target_compile_options(${EXTENSION_NAME} PRIVATE -include ${CMAKE_SOURCE_DIR}/src/include/wasm_socket_compat.h)
else()
	add_library(${EXTENSION_NAME} SHARED ${EXTENSION_SOURCES})
endif()

# Hide symbols
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(VISIBILITY_INLINES_HIDDEN ON)
# Strip symbols on Linux (not supported by AppleClang)
if(NOT APPLE)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -s")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")
endif()

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(${EXTENSION_NAME} PRIVATE -Wpedantic)
endif()

# Include own headers
target_include_directories(${EXTENSION_NAME} PRIVATE src/include)
target_include_directories(${EXTENSION_NAME} PRIVATE src/duckvep/kernel/include)
target_include_directories(${EXTENSION_NAME} PRIVATE src/duckvep/kernel/src)
target_include_directories(${EXTENSION_NAME} PRIVATE third_party/cgranges)
target_include_directories(${EXTENSION_NAME} PRIVATE third_party/libBigWig)
target_include_directories(${EXTENSION_NAME} PRIVATE third_party/variantkey/include)
target_compile_definitions(${EXTENSION_NAME} PRIVATE NOCURL=1)

# Include DuckDB C API headers
target_include_directories(${EXTENSION_NAME} PRIVATE duckdb_capi)

# ---- htslib ----
# Build vendored htslib using its own configure + make.  CMake drives the
# process via ExternalProject_Add so it runs once, in-tree, producing a
# static libhts.a with -fPIC.  htslib's configure auto-detects every
# available library (zlib, bz2, lzma, libdeflate, libcurl, openssl …).
include(ExternalProject)

set(HTSLIB_SRC_DIR "${CMAKE_SOURCE_DIR}/third_party/htslib")
set(HTSLIB_BUILD_DIR "${HTSLIB_SRC_DIR}")  # build in-tree like upstream expects

# ------------------------------------------------------------------
# Dependency detection.
#
# We use find_package() which works everywhere:
#   • In CI the vcpkg toolchain intercepts these and finds vcpkg ports
#     (declared in vcpkg.json).
#   • On Linux with system packages, CMake's standard modules find them.
#   • On macOS, Homebrew installs into /opt/homebrew (ARM) or /usr/local
#     (Intel) — both are in CMake's default search paths.
# ------------------------------------------------------------------

# zlib — always required
find_package(ZLIB REQUIRED)

# Optional compression libraries
find_package(BZip2 QUIET)
find_package(LibLZMA QUIET)

# libdeflate has no standard FindModule; use find_library + find_path
find_library(LIBDEFLATE_LIBRARY NAMES deflate)
find_path(LIBDEFLATE_INCLUDE_DIR NAMES libdeflate.h)
if(LIBDEFLATE_LIBRARY AND LIBDEFLATE_INCLUDE_DIR)
    set(LIBDEFLATE_FOUND TRUE)
else()
    set(LIBDEFLATE_FOUND FALSE)
endif()

# Optional network / crypto
find_package(CURL QUIET)
find_package(OpenSSL QUIET)
find_package(PkgConfig QUIET)

# Threads (portable: resolves to -lpthread on Linux, nothing on macOS)
find_package(Threads REQUIRED)

# MinGW regex (regcomp/regexec are in libgnurx)
if(WIN32 AND MINGW)
    find_library(GNUREGEX_LIBRARY NAMES gnurx regex)
    find_library(MINGW_LIBPSL_LIBRARY NAMES psl)
    find_library(MINGW_LIBIDN2_LIBRARY NAMES idn2)
    find_library(MINGW_LIBUNISTRING_LIBRARY NAMES unistring)
    find_library(MINGW_LIBICONV_LIBRARY NAMES iconv)
    find_library(MINGW_LIBSSH2_LIBRARY NAMES ssh2)
    find_library(MINGW_LIBGCRYPT_LIBRARY NAMES gcrypt)
    find_library(MINGW_LIBGPG_ERROR_LIBRARY NAMES gpg-error)
    find_library(MINGW_LIBZSTD_LIBRARY NAMES zstd)
endif()
if(WIN32 AND MINGW AND PKG_CONFIG_FOUND)
    pkg_check_modules(PC_LIBCURL QUIET libcurl)
endif()

# ------------------------------------------------------------------
# Build paths for htslib's Makefile.
#
# When vcpkg is active, libraries are installed under
#   ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/
# and we pass CPPFLAGS/LDFLAGS/LIBS directly to make. On plain system builds
# these are empty and harmless.
# ------------------------------------------------------------------
set(HTSLIB_EXTRA_CFLAGS "")
set(HTSLIB_EXTRA_CPPFLAGS "")
set(HTSLIB_BUILD_CPPFLAGS "")
set(HTSLIB_EXTRA_LDFLAGS "")

# duckdb-wasm uses -fwasm-exceptions throughout; htslib must match so its
# setjmp/longjmp uses __wasm_longjmp rather than emscripten_longjmp.
if(DUCKDB_WASM_EXTENSION)
    set(HTSLIB_EXTRA_CFLAGS "-fwasm-exceptions")
    # Keep the shared compatibility shim wasm-only.
    # The duckdb-wasm host exports native i64 file APIs as orig$lseek /
    # orig$ftruncate, so mark the top-level DuckDB wasm build explicitly while
    # leaving the webR/package wasm path on the default symbols.
    set(HTSLIB_BUILD_CPPFLAGS "-DDUCKHTS_WASM_DUCKDB_RUNTIME=1 -include ${CMAKE_SOURCE_DIR}/src/include/wasm_socket_compat.h")
endif()

# ------------------------------------------------------------------
# macOS cross-compilation support (htslib build flags).
# ------------------------------------------------------------------
if(APPLE AND DEFINED OSX_BUILD_ARCH AND NOT "${OSX_BUILD_ARCH}" STREQUAL "")
    set(HTSLIB_EXTRA_CFLAGS "-arch ${OSX_BUILD_ARCH}")
    set(HTSLIB_EXTRA_LDFLAGS "-arch ${OSX_BUILD_ARCH}")
    message(STATUS "macOS htslib arch: ${OSX_BUILD_ARCH}")
endif()

if(DEFINED VCPKG_INSTALLED_DIR AND DEFINED VCPKG_TARGET_TRIPLET)
    set(_vcpkg_prefix "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")
    if(EXISTS "${_vcpkg_prefix}")
        set(HTSLIB_EXTRA_CFLAGS "${HTSLIB_EXTRA_CFLAGS} -I${_vcpkg_prefix}/include")
        set(HTSLIB_EXTRA_LDFLAGS "${HTSLIB_EXTRA_LDFLAGS} -L${_vcpkg_prefix}/lib")
        set(HTSLIB_EXTRA_CPPFLAGS "${HTSLIB_EXTRA_CPPFLAGS} -I${_vcpkg_prefix}/include")
        set(HTSLIB_BUILD_CPPFLAGS "${HTSLIB_BUILD_CPPFLAGS} -I${_vcpkg_prefix}/include")
        message(STATUS "vcpkg prefix for htslib: ${_vcpkg_prefix}")
    endif()
endif()

# Strip leading whitespace that may appear when no arch flags were prepended
string(STRIP "${HTSLIB_EXTRA_CFLAGS}" HTSLIB_EXTRA_CFLAGS)
string(STRIP "${HTSLIB_EXTRA_CPPFLAGS}" HTSLIB_EXTRA_CPPFLAGS)
string(STRIP "${HTSLIB_BUILD_CPPFLAGS}" HTSLIB_BUILD_CPPFLAGS)
string(STRIP "${HTSLIB_EXTRA_LDFLAGS}" HTSLIB_EXTRA_LDFLAGS)

set(HTSLIB_MAKE_ARGS
    "CC=${CMAKE_C_COMPILER}"
    "CFLAGS=-O2 -fPIC ${HTSLIB_EXTRA_CFLAGS}"
    "LDFLAGS=${HTSLIB_EXTRA_LDFLAGS}"
)

if(NOT "${HTSLIB_BUILD_CPPFLAGS}" STREQUAL "")
    list(APPEND HTSLIB_MAKE_ARGS "CPPFLAGS=${HTSLIB_BUILD_CPPFLAGS}")
endif()

# Static-library transitive dependencies.
# When vcpkg (or similar) provides static .a archives, configure's link
# checks like `gcc test.c -lcurl` fail because libcurl.a itself needs
# -lssl -lcrypto -lz -lpthread etc.  Autotools' LIBS variable is
# appended to every link check, so we collect all transitive deps here.
set(HTSLIB_EXTRA_LIBS "")
if(OPENSSL_FOUND AND NOT DUCKDB_WASM_EXTENSION)
    list(APPEND HTSLIB_EXTRA_LIBS "-lssl" "-lcrypto")
endif()
list(APPEND HTSLIB_EXTRA_LIBS "-lz")
if(LIBLZMA_FOUND AND NOT DUCKDB_WASM_EXTENSION)
    list(APPEND HTSLIB_EXTRA_LIBS "-llzma")
endif()
if(BZIP2_FOUND)
    list(APPEND HTSLIB_EXTRA_LIBS "-lbz2")
endif()
if(APPLE)
    # macOS: pthreads are in libSystem (no -lpthread), and curl with
    # SecureTransport needs macOS frameworks for TLS.
    list(APPEND HTSLIB_EXTRA_LIBS "-lpthread")
    if(CURL_FOUND)
        list(APPEND HTSLIB_EXTRA_LIBS
            "-framework CoreFoundation"
            "-framework SystemConfiguration"
            "-framework Security"
        )
    endif()
else()
    list(APPEND HTSLIB_EXTRA_LIBS "-lpthread")
    if(NOT DUCKDB_WASM_EXTENSION)
        list(APPEND HTSLIB_EXTRA_LIBS "-ldl")
    endif()
endif()
if(WIN32 AND MINGW AND GNUREGEX_LIBRARY)
    list(APPEND HTSLIB_EXTRA_LIBS "-lgnurx")
endif()
list(JOIN HTSLIB_EXTRA_LIBS " " HTSLIB_EXTRA_LIBS_STR)
list(APPEND HTSLIB_MAKE_ARGS "LIBS=${HTSLIB_EXTRA_LIBS_STR}")

message(STATUS "htslib make args: ${HTSLIB_MAKE_ARGS}")

# htslib uses GNU Makefiles (not Ninja), so we must always use `make`
# even when the outer CMake build uses Ninja.
find_program(MAKE_COMMAND NAMES gmake make REQUIRED)
find_program(SH_COMMAND NAMES sh bash REQUIRED)

set(HTSLIB_CONFIGURE_FLAGS "")
if(NOT BZIP2_FOUND)
    list(APPEND HTSLIB_CONFIGURE_FLAGS --disable-bz2)
endif()
if(NOT LIBLZMA_FOUND)
    list(APPEND HTSLIB_CONFIGURE_FLAGS --disable-lzma)
endif()
if(NOT LIBDEFLATE_FOUND)
    list(APPEND HTSLIB_CONFIGURE_FLAGS --disable-libdeflate)
endif()
if(DUCKDB_WASM_EXTENSION)
    # Curl-based HTTP and plugin loading both fail in a WASM SIDE_MODULE context
    # (no socket access in Web Workers, no dlopen support).  HTTP/HTTPS access
    # is provided by wasm_http_hfile.c (synchronous XHR via EM_ASM) instead.
    # Curl's internal function pointer callbacks also produce i64 call_indirect
    # signatures that mismatch duckdb-wasm's function table ABI.
    # Match the stricter webR/browser path: no curl, no remote plugin loaders,
    # and no wasm-side lzma dependency.
    list(APPEND HTSLIB_CONFIGURE_FLAGS --disable-libcurl --disable-plugins --disable-s3 --disable-gcs --disable-lzma)
elseif(WIN32 AND MINGW)
    # MinGW/Rtools builds should not rely on htslib's checked-in config.mk.
    # Keep plugins off for now, but allow libcurl when the outer link can
    # provide the full dependency chain.
    list(APPEND HTSLIB_CONFIGURE_FLAGS --disable-plugins)
elseif(NOT CURL_FOUND)
    list(APPEND HTSLIB_CONFIGURE_FLAGS --disable-libcurl)
endif()

set(HTSLIB_CONFIGURE_ENV_VARS "")
if(DUCKDB_WASM_EXTENSION)
    # Match the package configure path: skip htslib's recv() library probe in
    # browser wasm builds, where socket libraries do not exist.
    list(APPEND HTSLIB_CONFIGURE_ENV_VARS "ac_cv_search_recv=none required")
endif()

set(HTSLIB_CONFIGURE_COMMAND
    ${CMAKE_COMMAND} -E env
    "CC=${CMAKE_C_COMPILER}"
    "CPPFLAGS=${HTSLIB_EXTRA_CPPFLAGS}"
    "CFLAGS=-O2 -fPIC ${HTSLIB_EXTRA_CFLAGS}"
    "LDFLAGS=${HTSLIB_EXTRA_LDFLAGS}"
    "LIBS=${HTSLIB_EXTRA_LIBS_STR}"
    ${HTSLIB_CONFIGURE_ENV_VARS}
    ${SH_COMMAND} ./configure ${HTSLIB_CONFIGURE_FLAGS}
)

ExternalProject_Add(htslib_build
    SOURCE_DIR        "${HTSLIB_SRC_DIR}"
    BUILD_IN_SOURCE   TRUE
    CONFIGURE_COMMAND ${HTSLIB_CONFIGURE_COMMAND}
    BUILD_COMMAND     ${MAKE_COMMAND} -j lib-static ${HTSLIB_MAKE_ARGS}
    BUILD_ALWAYS      TRUE
    INSTALL_COMMAND   ""   # no install step; we link libhts.a directly
    BUILD_BYPRODUCTS  "${HTSLIB_BUILD_DIR}/libhts.a"
    LOG_CONFIGURE     TRUE
    LOG_BUILD         TRUE
)

# Import the static library so CMake knows about it
add_library(hts STATIC IMPORTED GLOBAL)
set_target_properties(hts PROPERTIES
    IMPORTED_LOCATION "${HTSLIB_BUILD_DIR}/libhts.a"
    INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_SRC_DIR}"
)
add_dependencies(hts htslib_build)

# ------------------------------------------------------------------
# Link hts + transitive dependencies (only what was found).
# Uses CMake imported targets for correctness & portability.
# ------------------------------------------------------------------
set(HTSLIB_LINK_LIBS hts ZLIB::ZLIB)

if(BZIP2_FOUND)
    list(APPEND HTSLIB_LINK_LIBS BZip2::BZip2)
endif()
if(LIBLZMA_FOUND)
    list(APPEND HTSLIB_LINK_LIBS LibLZMA::LibLZMA)
endif()
if(LIBDEFLATE_FOUND)
    list(APPEND HTSLIB_LINK_LIBS ${LIBDEFLATE_LIBRARY})
endif()
if(CURL_FOUND AND NOT DUCKDB_WASM_EXTENSION)
    if(WIN32 AND MINGW)
        if(PC_LIBCURL_FOUND)
            string(STRIP "${PC_LIBCURL_LDFLAGS}" PC_LIBCURL_LDFLAGS_STRIPPED)
            separate_arguments(PC_LIBCURL_LDFLAGS_LIST NATIVE_COMMAND "${PC_LIBCURL_LDFLAGS_STRIPPED}")
            list(APPEND HTSLIB_LINK_LIBS ${PC_LIBCURL_LDFLAGS_LIST})
        else()
            list(APPEND HTSLIB_LINK_LIBS
                -lcurl
                -lbcrypt
                -lssl
                -lcrypto
                -lcrypt32
                -lwsock32
                -lwldap32
            )
            if(MINGW_LIBPSL_LIBRARY)
                list(APPEND HTSLIB_LINK_LIBS ${MINGW_LIBPSL_LIBRARY})
            endif()
            if(MINGW_LIBIDN2_LIBRARY)
                list(APPEND HTSLIB_LINK_LIBS ${MINGW_LIBIDN2_LIBRARY})
            endif()
            if(MINGW_LIBUNISTRING_LIBRARY)
                list(APPEND HTSLIB_LINK_LIBS ${MINGW_LIBUNISTRING_LIBRARY})
            endif()
            if(MINGW_LIBICONV_LIBRARY)
                list(APPEND HTSLIB_LINK_LIBS ${MINGW_LIBICONV_LIBRARY})
            endif()
            if(MINGW_LIBSSH2_LIBRARY)
                list(APPEND HTSLIB_LINK_LIBS ${MINGW_LIBSSH2_LIBRARY})
            endif()
            if(MINGW_LIBGCRYPT_LIBRARY)
                list(APPEND HTSLIB_LINK_LIBS ${MINGW_LIBGCRYPT_LIBRARY})
            endif()
            if(MINGW_LIBGPG_ERROR_LIBRARY)
                list(APPEND HTSLIB_LINK_LIBS ${MINGW_LIBGPG_ERROR_LIBRARY})
            endif()
            if(MINGW_LIBZSTD_LIBRARY)
                list(APPEND HTSLIB_LINK_LIBS ${MINGW_LIBZSTD_LIBRARY})
            endif()
        endif()
    else()
        list(APPEND HTSLIB_LINK_LIBS CURL::libcurl)
    endif()
endif()

if(OPENSSL_FOUND AND NOT DUCKDB_WASM_EXTENSION)
    list(APPEND HTSLIB_LINK_LIBS OpenSSL::Crypto)
endif()
if(WIN32 AND MINGW AND GNUREGEX_LIBRARY)
    list(APPEND HTSLIB_LINK_LIBS ${GNUREGEX_LIBRARY})
endif()
if(WIN32 AND MINGW)
    list(APPEND HTSLIB_LINK_LIBS ws2_32)
endif()

# Threads::Threads is portable (pthread on Linux, nothing on macOS)
list(APPEND HTSLIB_LINK_LIBS Threads::Threads)

# libm: needed on Linux/musl, not present on Windows
if(NOT APPLE AND NOT WIN32)
    list(APPEND HTSLIB_LINK_LIBS m)
endif()

message(STATUS "htslib link libs: ${HTSLIB_LINK_LIBS}")

target_link_libraries(${EXTENSION_NAME} PRIVATE
    ${HTSLIB_LINK_LIBS}
)

if(DUCKDB_WASM_EXTENSION)
    # The DuckDB wasm packaging step only sees libduckhts.a, so make that
    # archive self-contained by merging htslib and any static archive
    # dependencies directly into it after the normal static library build.
    set(DUCKHTS_WASM_FAT_ARCHIVE_INPUTS
        "${HTSLIB_BUILD_DIR}/libhts.a"
    )

    if(TARGET ZLIB::ZLIB)
        list(APPEND DUCKHTS_WASM_FAT_ARCHIVE_INPUTS "$<TARGET_FILE:ZLIB::ZLIB>")
    elseif(ZLIB_LIBRARIES)
        list(APPEND DUCKHTS_WASM_FAT_ARCHIVE_INPUTS ${ZLIB_LIBRARIES})
    endif()

    if(BZIP2_FOUND)
        if(TARGET BZip2::BZip2)
            list(APPEND DUCKHTS_WASM_FAT_ARCHIVE_INPUTS "$<TARGET_FILE:BZip2::BZip2>")
        elseif(BZIP2_LIBRARIES)
            list(APPEND DUCKHTS_WASM_FAT_ARCHIVE_INPUTS ${BZIP2_LIBRARIES})
        endif()
    endif()

    if(LIBLZMA_FOUND)
        if(TARGET LibLZMA::LibLZMA)
            list(APPEND DUCKHTS_WASM_FAT_ARCHIVE_INPUTS "$<TARGET_FILE:LibLZMA::LibLZMA>")
        elseif(LIBLZMA_LIBRARIES)
            list(APPEND DUCKHTS_WASM_FAT_ARCHIVE_INPUTS ${LIBLZMA_LIBRARIES})
        endif()
    endif()

    if(LIBDEFLATE_FOUND)
        list(APPEND DUCKHTS_WASM_FAT_ARCHIVE_INPUTS "${LIBDEFLATE_LIBRARY}")
    endif()

    if(CURL_FOUND AND NOT DUCKDB_WASM_EXTENSION)
        if(WIN32 AND MINGW)
            list(APPEND DUCKHTS_WASM_FAT_ARCHIVE_INPUTS ${PC_LIBCURL_LDFLAGS_LIST})
        elseif(TARGET CURL::libcurl)
            list(APPEND DUCKHTS_WASM_FAT_ARCHIVE_INPUTS "$<TARGET_FILE:CURL::libcurl>")
        endif()
    endif()

    if(OPENSSL_FOUND AND NOT DUCKDB_WASM_EXTENSION AND TARGET OpenSSL::Crypto)
        list(APPEND DUCKHTS_WASM_FAT_ARCHIVE_INPUTS "$<TARGET_FILE:OpenSSL::Crypto>")
    endif()

    if(WIN32 AND MINGW AND GNUREGEX_LIBRARY)
        list(APPEND DUCKHTS_WASM_FAT_ARCHIVE_INPUTS "${GNUREGEX_LIBRARY}")
    endif()

    list(REMOVE_DUPLICATES DUCKHTS_WASM_FAT_ARCHIVE_INPUTS)
    string(JOIN "\n" DUCKHTS_WASM_FAT_ARCHIVE_INPUTS_CONTENT ${DUCKHTS_WASM_FAT_ARCHIVE_INPUTS})
    file(GENERATE
        OUTPUT "${CMAKE_BINARY_DIR}/wasm_fat_archive_inputs.txt"
        CONTENT "${DUCKHTS_WASM_FAT_ARCHIVE_INPUTS_CONTENT}\n"
    )

    add_dependencies(${EXTENSION_NAME} htslib_build)
    add_custom_command(TARGET ${EXTENSION_NAME} POST_BUILD
        COMMAND
            ${CMAKE_COMMAND}
            -DARCHIVE=$<TARGET_FILE:${EXTENSION_NAME}>
            -DINPUT_LIST=${CMAKE_BINARY_DIR}/wasm_fat_archive_inputs.txt
            -DARCHIVER=${CMAKE_AR}
            -DRANLIB=${CMAKE_RANLIB}
            -P ${CMAKE_SOURCE_DIR}/cmake/merge_static_archives.cmake
        VERBATIM
    )
endif()
