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

# ------------------------------------------------------------------
# Build paths for htslib's configure.
#
# When vcpkg is active, libraries are installed under
#   ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/
# and htslib's autotools configure needs CPPFLAGS/LDFLAGS/PKG_CONFIG_PATH
# pointing there.  On plain system builds these are empty and harmless.
# ------------------------------------------------------------------
set(HTSLIB_EXTRA_CFLAGS "")
set(HTSLIB_EXTRA_LDFLAGS "")
set(HTSLIB_PKG_CONFIG_PATH "")

# ------------------------------------------------------------------
# macOS cross-compilation support.
#
# The CI passes -DOSX_BUILD_ARCH=x86_64 (or arm64) for universal builds.
# When building for x86_64 on an arm64 Mac (or vice-versa) we need to
# tell both the C compiler and htslib's autotools configure about the
# target architecture.
# ------------------------------------------------------------------
set(HTSLIB_HOST_FLAG "")
if(APPLE)
    if(NOT DEFINED OSX_BUILD_ARCH OR "${OSX_BUILD_ARCH}" STREQUAL "")
        if(DEFINED CMAKE_OSX_ARCHITECTURES AND NOT "${CMAKE_OSX_ARCHITECTURES}" STREQUAL "")
            list(GET CMAKE_OSX_ARCHITECTURES 0 OSX_BUILD_ARCH)
        elseif(DEFINED VCPKG_TARGET_TRIPLET)
            if(VCPKG_TARGET_TRIPLET MATCHES "^x64-osx")
                set(OSX_BUILD_ARCH "x86_64")
            elseif(VCPKG_TARGET_TRIPLET MATCHES "^arm64-osx")
                set(OSX_BUILD_ARCH "arm64")
            endif()
        endif()
    endif()
    if(DEFINED OSX_BUILD_ARCH AND NOT "${OSX_BUILD_ARCH}" STREQUAL "")
        set(CMAKE_OSX_ARCHITECTURES "${OSX_BUILD_ARCH}" CACHE STRING "" FORCE)
        set(HTSLIB_EXTRA_CFLAGS "-arch ${OSX_BUILD_ARCH}")
        set(HTSLIB_EXTRA_LDFLAGS "-arch ${OSX_BUILD_ARCH}")
        if("${OSX_BUILD_ARCH}" STREQUAL "x86_64")
            set(HTSLIB_HOST_FLAG "--host=x86_64-apple-darwin")
        elseif("${OSX_BUILD_ARCH}" STREQUAL "arm64")
            set(HTSLIB_HOST_FLAG "--host=aarch64-apple-darwin")
        endif()
        message(STATUS "macOS cross-compile: arch=${OSX_BUILD_ARCH} host=${HTSLIB_HOST_FLAG}")
    endif()
endif()

if(DEFINED VCPKG_INSTALLED_DIR AND DEFINED VCPKG_TARGET_TRIPLET)
    set(_vcpkg_prefix "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")
    if(EXISTS "${_vcpkg_prefix}")
        set(HTSLIB_EXTRA_CFLAGS "${HTSLIB_EXTRA_CFLAGS} -I${_vcpkg_prefix}/include")
        set(HTSLIB_EXTRA_LDFLAGS "${HTSLIB_EXTRA_LDFLAGS} -L${_vcpkg_prefix}/lib")
        set(HTSLIB_PKG_CONFIG_PATH "${_vcpkg_prefix}/lib/pkgconfig")
        message(STATUS "vcpkg prefix for htslib: ${_vcpkg_prefix}")
    endif()
endif()

# Strip leading whitespace that may appear when no arch flags were prepended
string(STRIP "${HTSLIB_EXTRA_CFLAGS}" HTSLIB_EXTRA_CFLAGS)
string(STRIP "${HTSLIB_EXTRA_LDFLAGS}" HTSLIB_EXTRA_LDFLAGS)

set(HTSLIB_CONFIGURE_ARGS
    "CFLAGS=-fPIC -O2 ${HTSLIB_EXTRA_CFLAGS}"
    "LDFLAGS=${HTSLIB_EXTRA_LDFLAGS}"
    "--disable-plugins"
)

# --host for cross-compilation (must come before variable assignments)
if(NOT "${HTSLIB_HOST_FLAG}" STREQUAL "")
    list(INSERT HTSLIB_CONFIGURE_ARGS 0 "${HTSLIB_HOST_FLAG}")
endif()

# Autotools convention: CPPFLAGS for -I, CFLAGS for -O/-f etc.
# Some configure scripts only honour CPPFLAGS for include-path probes.
# Only pass -I flags in CPPFLAGS (not -arch or other compiler flags).
if(DEFINED _vcpkg_prefix AND EXISTS "${_vcpkg_prefix}")
    list(APPEND HTSLIB_CONFIGURE_ARGS "CPPFLAGS=-I${_vcpkg_prefix}/include")
endif()

# Pass PKG_CONFIG_PATH so configure's pkg-config checks find vcpkg libs
if(NOT "${HTSLIB_PKG_CONFIG_PATH}" STREQUAL "")
    list(APPEND HTSLIB_CONFIGURE_ARGS "PKG_CONFIG_PATH=${HTSLIB_PKG_CONFIG_PATH}")
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
list(JOIN HTSLIB_EXTRA_LIBS " " HTSLIB_EXTRA_LIBS_STR)
list(APPEND HTSLIB_CONFIGURE_ARGS "LIBS=${HTSLIB_EXTRA_LIBS_STR}")

if(BZIP2_FOUND)
    list(APPEND HTSLIB_CONFIGURE_ARGS "--enable-bz2")
else()
    list(APPEND HTSLIB_CONFIGURE_ARGS "--disable-bz2")
endif()

if(LIBLZMA_FOUND)
    list(APPEND HTSLIB_CONFIGURE_ARGS "--enable-lzma")
else()
    list(APPEND HTSLIB_CONFIGURE_ARGS "--disable-lzma")
endif()

if(CURL_FOUND)
    list(APPEND HTSLIB_CONFIGURE_ARGS "--enable-libcurl")
else()
    list(APPEND HTSLIB_CONFIGURE_ARGS "--disable-libcurl")
endif()

if(CURL_FOUND AND OPENSSL_FOUND)
    list(APPEND HTSLIB_CONFIGURE_ARGS "--enable-gcs" "--enable-s3")
else()
    list(APPEND HTSLIB_CONFIGURE_ARGS "--disable-gcs" "--disable-s3")
endif()

if(LIBDEFLATE_FOUND)
    list(APPEND HTSLIB_CONFIGURE_ARGS "--with-libdeflate")
else()
    list(APPEND HTSLIB_CONFIGURE_ARGS "--without-libdeflate")
endif()

message(STATUS "htslib configure args: ${HTSLIB_CONFIGURE_ARGS}")

# htslib uses GNU Makefiles (not Ninja), so we must always use `make`
# even when the outer CMake build uses Ninja.
find_program(MAKE_COMMAND NAMES gmake make REQUIRED)

ExternalProject_Add(htslib_build
    SOURCE_DIR        "${HTSLIB_SRC_DIR}"
    BUILD_IN_SOURCE   TRUE
    CONFIGURE_COMMAND "${HTSLIB_SRC_DIR}/configure"
                      ${HTSLIB_CONFIGURE_ARGS}
    BUILD_COMMAND     ${MAKE_COMMAND} -j lib-static
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
if(CURL_FOUND)
    list(APPEND HTSLIB_LINK_LIBS CURL::libcurl)
endif()
if(OPENSSL_FOUND)
    list(APPEND HTSLIB_LINK_LIBS OpenSSL::Crypto)
endif()

# Threads::Threads is portable (pthread on Linux, nothing on macOS)
list(APPEND HTSLIB_LINK_LIBS Threads::Threads)

# libm: needed on Linux/musl, harmless no-op on macOS
if(NOT APPLE)
    list(APPEND HTSLIB_LINK_LIBS m)
endif()

message(STATUS "htslib link libs: ${HTSLIB_LINK_LIBS}")

target_link_libraries(${EXTENSION_NAME} PRIVATE
    ${HTSLIB_LINK_LIBS}
)
