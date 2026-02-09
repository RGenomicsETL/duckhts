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
# Strip symbols on Linux (not supported by AppleClang)
if(NOT APPLE)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -s")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")
endif()

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

# Threads (portable: resolves to -lpthread on Linux, nothing on macOS)
find_package(Threads REQUIRED)

# MinGW regex (regcomp/regexec are in libgnurx)
if(WIN32 AND MINGW)
    find_library(GNUREGEX_LIBRARY NAMES gnurx regex)
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
set(HTSLIB_EXTRA_LDFLAGS "")

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
        message(STATUS "vcpkg prefix for htslib: ${_vcpkg_prefix}")
    endif()
endif()

# Strip leading whitespace that may appear when no arch flags were prepended
string(STRIP "${HTSLIB_EXTRA_CFLAGS}" HTSLIB_EXTRA_CFLAGS)
string(STRIP "${HTSLIB_EXTRA_LDFLAGS}" HTSLIB_EXTRA_LDFLAGS)

set(HTSLIB_MAKE_ARGS
    "CC=${CMAKE_C_COMPILER}"
    "CFLAGS=-O2 -fPIC ${HTSLIB_EXTRA_CFLAGS}"
    "LDFLAGS=${HTSLIB_EXTRA_LDFLAGS}"
)

# Only pass -I flags in CPPFLAGS (not -arch or other compiler flags).
if(DEFINED _vcpkg_prefix AND EXISTS "${_vcpkg_prefix}")
    list(APPEND HTSLIB_MAKE_ARGS "CPPFLAGS=-I${_vcpkg_prefix}/include")
endif()

# Static-library transitive dependencies.
# When vcpkg (or similar) provides static .a archives, configure's link
# checks like `gcc test.c -lcurl` fail because libcurl.a itself needs
# -lssl -lcrypto -lz -lpthread etc.  Autotools' LIBS variable is
# appended to every link check, so we collect all transitive deps here.
set(HTSLIB_EXTRA_LIBS "")
if(OPENSSL_FOUND)
    list(APPEND HTSLIB_EXTRA_LIBS "-lssl" "-lcrypto")
endif()
list(APPEND HTSLIB_EXTRA_LIBS "-lz")
if(LIBLZMA_FOUND)
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
    list(APPEND HTSLIB_EXTRA_LIBS "-ldl")
endif()
if(WIN32 AND MINGW AND GNUREGEX_LIBRARY)
    list(APPEND HTSLIB_EXTRA_LIBS "-lgnurx")
endif()
list(JOIN HTSLIB_EXTRA_LIBS " " HTSLIB_EXTRA_LIBS_STR)
list(APPEND HTSLIB_MAKE_ARGS "LIBS=${HTSLIB_EXTRA_LIBS_STR}")

message(STATUS "htslib make args: ${HTSLIB_MAKE_ARGS}")

if(WIN32 AND NOT MINGW)
    # MSVC builds can't run autotools configure. Use vcpkg/system htslib.
    find_path(HTSLIB_INCLUDE_DIR NAMES htslib/hts.h)
    find_library(HTSLIB_LIBRARY NAMES htslib hts)
    if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
        message(FATAL_ERROR "htslib not found. On Windows MSVC, install htslib via vcpkg or provide HTSLIB_INCLUDE_DIR/HTSLIB_LIBRARY.")
    endif()
    add_library(hts UNKNOWN IMPORTED GLOBAL)
    set_target_properties(hts PROPERTIES
        IMPORTED_LOCATION "${HTSLIB_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIR}"
    )
else()
    # htslib uses GNU Makefiles (not Ninja), so we must always use `make`
    # even when the outer CMake build uses Ninja.
    find_program(MAKE_COMMAND NAMES gmake make REQUIRED)

    ExternalProject_Add(htslib_build
        SOURCE_DIR        "${HTSLIB_SRC_DIR}"
        BUILD_IN_SOURCE   TRUE
        CONFIGURE_COMMAND ""
        BUILD_COMMAND     ${MAKE_COMMAND} -j lib-static ${HTSLIB_MAKE_ARGS}
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
endif()

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
if(CURL_FOUND)
    list(APPEND HTSLIB_LINK_LIBS CURL::libcurl)
endif()
if(OPENSSL_FOUND)
    list(APPEND HTSLIB_LINK_LIBS OpenSSL::Crypto)
endif()
if(WIN32 AND MINGW AND GNUREGEX_LIBRARY)
    list(APPEND HTSLIB_LINK_LIBS ${GNUREGEX_LIBRARY})
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
