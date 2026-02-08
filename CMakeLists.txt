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

# ---- htslib integration ----
# Prefer vendored third_party/htslib if present, else fall back to system
set(HTSLIB_DIR "${CMAKE_SOURCE_DIR}/third_party/htslib")
if(EXISTS "${HTSLIB_DIR}/htslib/hts.h")
    message(STATUS "Using vendored htslib from ${HTSLIB_DIR}")
    target_include_directories(${EXTENSION_NAME} PRIVATE "${HTSLIB_DIR}")

    if(EXISTS "${HTSLIB_DIR}/libhts.a")
        target_link_libraries(${EXTENSION_NAME} PRIVATE "${HTSLIB_DIR}/libhts.a")

        # Parse htslib's config.mk to discover which libraries it was built with.
        # config.mk has a main LIBS line plus noplugin_LIBS (for statically-linked
        # plugin deps like curl, crypto).  We collect -l flags from both.
        if(EXISTS "${HTSLIB_DIR}/config.mk")
            file(READ "${HTSLIB_DIR}/config.mk" _htslib_config)

            # Collect all -l flags from LIBS, noplugin_LIBS, LIBCURL_LIBS, CRYPTO_LIBS lines
            string(REGEX MATCHALL "-l([^ \t\n]+)" _all_lib_flags "${_htslib_config}")
            set(_seen_libs "")
            foreach(_flag ${_all_lib_flags})
                string(REGEX REPLACE "^-l" "" _lib "${_flag}")
                if(NOT _lib IN_LIST _seen_libs)
                    list(APPEND _seen_libs "${_lib}")
                    target_link_libraries(${EXTENSION_NAME} PRIVATE ${_lib})
                endif()
            endforeach()
            message(STATUS "htslib link libs (from config.mk): ${_seen_libs}")

            # Always need pthread and m
            target_link_libraries(${EXTENSION_NAME} PRIVATE pthread m)
        else()
            # Fallback: common htslib dependencies
            target_link_libraries(${EXTENSION_NAME} PRIVATE z bz2 lzma deflate curl crypto pthread m)
        endif()
    elseif(EXISTS "${HTSLIB_DIR}/libhts.so")
        target_link_libraries(${EXTENSION_NAME} PRIVATE "${HTSLIB_DIR}/libhts.so")
    else()
        message(WARNING "Vendored htslib headers found but no library. Run 'make' in third_party/htslib/ first.")
    endif()
else()
    # Try pkg-config for system htslib
    find_package(PkgConfig)
    if(PkgConfig_FOUND)
        pkg_check_modules(HTSLIB htslib)
        if(HTSLIB_FOUND)
            message(STATUS "Using system htslib via pkg-config")
            target_include_directories(${EXTENSION_NAME} PRIVATE ${HTSLIB_INCLUDE_DIRS})
            target_link_libraries(${EXTENSION_NAME} PRIVATE ${HTSLIB_LINK_LIBRARIES})
        endif()
    endif()
    if(NOT HTSLIB_FOUND)
        message(STATUS "No vendored or system htslib found; trying -lhts")
        target_link_libraries(${EXTENSION_NAME} PRIVATE hts z)
    endif()
endif()
