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
# Build htslib from vendored sources via its own CMakeLists.txt.
# This detects available system libraries (zlib, bz2, lzma, libdeflate,
# curl, openssl) and compiles htslib as a static library with -fPIC.
add_subdirectory(third_party/htslib)
target_link_libraries(${EXTENSION_NAME} PRIVATE hts)
