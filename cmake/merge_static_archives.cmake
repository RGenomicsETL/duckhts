if(NOT DEFINED ARCHIVE OR NOT DEFINED INPUT_LIST OR NOT DEFINED ARCHIVER)
    message(FATAL_ERROR "ARCHIVE, INPUT_LIST, and ARCHIVER are required")
endif()

if(NOT EXISTS "${ARCHIVE}")
    message(FATAL_ERROR "Archive to merge does not exist: ${ARCHIVE}")
endif()

if(NOT EXISTS "${INPUT_LIST}")
    message(FATAL_ERROR "Archive input list does not exist: ${INPUT_LIST}")
endif()

if(NOT DEFINED RANLIB OR "${RANLIB}" STREQUAL "")
    set(RANLIB "")
endif()

set(MERGE_ROOT "${ARCHIVE}.merge")
file(REMOVE_RECURSE "${MERGE_ROOT}")
file(MAKE_DIRECTORY "${MERGE_ROOT}")

set(ARCHIVES_TO_PACK "${ARCHIVE}")
file(STRINGS "${INPUT_LIST}" EXTRA_ARCHIVES)
foreach(extra_archive IN LISTS EXTRA_ARCHIVES)
    if(NOT extra_archive STREQUAL "" AND EXISTS "${extra_archive}" AND extra_archive MATCHES [[\.a$]])
        list(APPEND ARCHIVES_TO_PACK "${extra_archive}")
    endif()
endforeach()

set(ALL_OBJECTS "")
set(ARCHIVE_INDEX 0)
foreach(current_archive IN LISTS ARCHIVES_TO_PACK)
    math(EXPR ARCHIVE_INDEX "${ARCHIVE_INDEX} + 1")
    set(EXTRACT_DIR "${MERGE_ROOT}/${ARCHIVE_INDEX}")
    file(MAKE_DIRECTORY "${EXTRACT_DIR}")

    execute_process(
        COMMAND "${ARCHIVER}" x "${current_archive}"
        WORKING_DIRECTORY "${EXTRACT_DIR}"
        RESULT_VARIABLE EXTRACT_RESULT
    )
    if(NOT EXTRACT_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to extract archive: ${current_archive}")
    endif()

    file(GLOB OBJECT_FILES "${EXTRACT_DIR}/*.o")
    list(APPEND ALL_OBJECTS ${OBJECT_FILES})
endforeach()

if(ALL_OBJECTS STREQUAL "")
    message(FATAL_ERROR "No object files were extracted while merging archives")
endif()

file(REMOVE "${ARCHIVE}")
execute_process(
    COMMAND "${ARCHIVER}" qc "${ARCHIVE}" ${ALL_OBJECTS}
    RESULT_VARIABLE ARCHIVE_RESULT
)
if(NOT ARCHIVE_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to create merged archive: ${ARCHIVE}")
endif()

if(NOT "${RANLIB}" STREQUAL "")
    execute_process(
        COMMAND "${RANLIB}" "${ARCHIVE}"
        RESULT_VARIABLE RANLIB_RESULT
    )
    if(NOT RANLIB_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to run ranlib on merged archive: ${ARCHIVE}")
    endif()
endif()
