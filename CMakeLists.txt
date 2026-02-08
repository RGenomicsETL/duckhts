cmake_minimum_required(VERSION 3.5...3.29)
option(DUCKDB_WASM_EXTENSION "Whether compiling for Wasm target" OFF)

###
# Configuration
###
if(NOT DEFINED EXTENSION_NAME)
    message(FATAL_ERROR "DuckDB extension name is required")
endif()

add_definitions(-DDUCKDB_EXTENSION_NAME=${EXTENSION_NAME})

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
        src/bcf_reader.c
        src/bam_reader.c
        src/seq_reader.c
        src/tabix_reader.c
        src/vep_parser.c
)

if (DUCKDB_WASM_EXTENSION)
	add_library(${EXTENSION_NAME} STATIC ${EXTENSION_SOURCES})
else()
	add_library(${EXTENSION_NAME} SHARED ${EXTENSION_SOURCES})
endif()

# Hide symbols
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(VISIBILITY_INLINES_HIDDEN ON)
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -s")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")

# Include own headers
target_include_directories(${EXTENSION_NAME} PRIVATE src/include)

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
# Probe optional dependencies and build configure flags dynamically.
# htslib's ./configure will fail when --enable-X is given but the lib
# is missing, so we only pass flags for libraries we actually find.
# ------------------------------------------------------------------

# Probe all dependencies with find_library (not find_package) so that
# vcpkg toolchain files don't intercept the search and fail when the
# libs are installed as system packages but not via vcpkg.

# zlib — required by htslib; we always link it.
find_library(LIBZ NAMES z)

# Optional compression libraries
find_library(LIBBZ2   NAMES bz2)
find_library(LIBLZMA  NAMES lzma)
find_library(LIBDEFLATE NAMES deflate)

# Optional network / crypto libraries
find_library(LIBCURL  NAMES curl)
find_library(LIBCRYPTO NAMES crypto)

# Build the htslib configure argument list
set(HTSLIB_CONFIGURE_ARGS
    "CFLAGS=-fPIC -O2"
    "--disable-plugins"
)

if(LIBCURL)
    list(APPEND HTSLIB_CONFIGURE_ARGS "--enable-libcurl")
else()
    list(APPEND HTSLIB_CONFIGURE_ARGS "--disable-libcurl")
endif()

if(LIBCURL AND LIBCRYPTO)
    list(APPEND HTSLIB_CONFIGURE_ARGS "--enable-gcs" "--enable-s3")
else()
    list(APPEND HTSLIB_CONFIGURE_ARGS "--disable-gcs" "--disable-s3")
endif()

if(LIBDEFLATE)
    list(APPEND HTSLIB_CONFIGURE_ARGS "--with-libdeflate")
else()
    list(APPEND HTSLIB_CONFIGURE_ARGS "--without-libdeflate")
endif()

message(STATUS "htslib configure args: ${HTSLIB_CONFIGURE_ARGS}")

ExternalProject_Add(htslib_build
    SOURCE_DIR        "${HTSLIB_SRC_DIR}"
    BUILD_IN_SOURCE   TRUE
    CONFIGURE_COMMAND "${HTSLIB_SRC_DIR}/configure"
                      ${HTSLIB_CONFIGURE_ARGS}
    BUILD_COMMAND     ${CMAKE_MAKE_PROGRAM} -j lib-static
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

# Link hts and only the transitive dependencies that are actually present.
# zlib is always needed; use the found path or fall back to -lz.
if(LIBZ)
    set(HTSLIB_LINK_LIBS hts ${LIBZ})
else()
    set(HTSLIB_LINK_LIBS hts z)
endif()

if(LIBBZ2)
    list(APPEND HTSLIB_LINK_LIBS ${LIBBZ2})
endif()
if(LIBLZMA)
    list(APPEND HTSLIB_LINK_LIBS ${LIBLZMA})
endif()
if(LIBDEFLATE)
    list(APPEND HTSLIB_LINK_LIBS ${LIBDEFLATE})
endif()
if(LIBCURL)
    list(APPEND HTSLIB_LINK_LIBS ${LIBCURL})
endif()
if(LIBCRYPTO)
    list(APPEND HTSLIB_LINK_LIBS ${LIBCRYPTO})
endif()

list(APPEND HTSLIB_LINK_LIBS pthread m)

message(STATUS "htslib link libs: ${HTSLIB_LINK_LIBS}")

target_link_libraries(${EXTENSION_NAME} PRIVATE
    ${HTSLIB_LINK_LIBS}
)
