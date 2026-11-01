# Check only the library objects, not transitive platform math/runtime code.
if(NOT DEFINED NM OR NOT DEFINED ARCHIVE)
    message(FATAL_ERROR "NM and ARCHIVE are required")
endif()
execute_process(COMMAND "${NM}" -u "${ARCHIVE}"
    RESULT_VARIABLE status OUTPUT_VARIABLE symbols ERROR_VARIABLE errors)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "Cannot inspect library symbols: ${errors}")
endif()
if(symbols MATCHES "(^|[ \t\r\n])_?(malloc|calloc|realloc|free|aligned_alloc|strdup)([ \t\r\n]|$)")
    message(FATAL_ERROR "Unexpected direct heap allocator dependency:\n${symbols}")
endif()
message(STATUS "No direct heap allocator dependencies")
