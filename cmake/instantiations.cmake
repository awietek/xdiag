
# Detects sources containing BEGIN_INSTANTIATION_GROUP markers
function(detect_grouped_sources INPUT_SOURCES OUT_GROUPED OUT_NORMAL)

  set(GROUPED)
  set(NORMAL)

  foreach(SRC ${INPUT_SOURCES})

    # Read file content
    file(READ ${CMAKE_SOURCE_DIR}/xdiag/${SRC} FILE_CONTENTS)

    # Check for marker
    string(REGEX MATCH
      "BEGIN_INSTANTIATION_GROUP\\([A-Za-z0-9_]+\\)"
      HAS_GROUP
      "${FILE_CONTENTS}")

    if(HAS_GROUP)
      list(APPEND GROUPED ${SRC})
    else()
      list(APPEND NORMAL ${SRC})
    endif()

  endforeach()

  # Return values to caller
  set(${OUT_GROUPED} ${GROUPED} PARENT_SCOPE)
  set(${OUT_NORMAL}  ${NORMAL}  PARENT_SCOPE)

endfunction()


# Configure-time: compute the output file paths that would be generated from a
# grouped source, without writing anything.
#
# ABSOLUTE_SOURCE_FILE : absolute path to the .cpp source file
# SOURCE_ROOT          : root used to compute the relative subdirectory
#                        (i.e. ${CMAKE_SOURCE_DIR}/xdiag)
# OUTPUT_DIR           : directory under which generated files are placed,
#                        preserving the relative subdirectory structure
# OUT_FILES            : output variable receiving the list of absolute paths
function(get_instantiation_outputs ABSOLUTE_SOURCE_FILE SOURCE_ROOT OUTPUT_DIR OUT_FILES)
  file(READ "${ABSOLUTE_SOURCE_FILE}" FILE_CONTENTS)
  file(RELATIVE_PATH REL_PATH "${SOURCE_ROOT}" "${ABSOLUTE_SOURCE_FILE}")
  get_filename_component(REL_DIR   "${REL_PATH}" DIRECTORY)
  get_filename_component(BASE_NAME "${REL_PATH}" NAME_WE)

  string(REGEX MATCHALL
    "BEGIN_INSTANTIATION_GROUP\\([A-Za-z0-9_]+\\)"
    MATCHES "${FILE_CONTENTS}")

  set(OUT_LIST)
  foreach(MATCH ${MATCHES})
    string(REGEX REPLACE
      "BEGIN_INSTANTIATION_GROUP\\(([A-Za-z0-9_]+)\\)" "\\1"
      GROUP_NAME "${MATCH}")
    list(APPEND OUT_LIST "${OUTPUT_DIR}/${REL_DIR}/${BASE_NAME}_${GROUP_NAME}.cpp")
  endforeach()

  set(${OUT_FILES} "${OUT_LIST}" PARENT_SCOPE)
endfunction()


# Build-time: parse a grouped source file and write one .cpp per instantiation
# group.  Intended to be called from script mode (cmake -P).
#
# ABSOLUTE_SOURCE_FILE : absolute path to the .cpp source file
# SOURCE_ROOT          : root for computing the relative subdirectory
#                        (i.e. ${CMAKE_SOURCE_DIR}/xdiag)
# INCLUDE_ROOT         : root for computing #include paths
#                        (i.e. ${CMAKE_SOURCE_DIR})
# OUTPUT_DIR           : directory under which generated files are placed
function(generate_instantiation_files ABSOLUTE_SOURCE_FILE SOURCE_ROOT INCLUDE_ROOT OUTPUT_DIR)

  file(READ "${ABSOLUTE_SOURCE_FILE}" FILE_CONTENTS)
  file(RELATIVE_PATH REL_PATH   "${SOURCE_ROOT}"  "${ABSOLUTE_SOURCE_FILE}")
  file(RELATIVE_PATH INCL_PATH  "${INCLUDE_ROOT}" "${ABSOLUTE_SOURCE_FILE}")

  get_filename_component(REL_DIR   "${REL_PATH}"  DIRECTORY)
  get_filename_component(BASE_NAME "${REL_PATH}"  NAME_WE)

  # Header info for the local #include -> angle-bracket substitution
  string(REPLACE ".cpp" ".hpp" HEADER_INCL_PATH "${INCL_PATH}")
  get_filename_component(HEADER_FILENAME "${HEADER_INCL_PATH}" NAME)

  ######################
  # Split into lines, working around CMake's semicolon-as-list-separator.
  # We also encode '[' and ']' because CMake's list parser treats an unmatched
  # '[' as a bracket that swallows subsequent ';' separators (e.g. the comment
  # "// Thread t writes to the contiguous range [rep_offset[t], rep_offset[t+1])."
  # contains an outer '[' closed by ')' not ']', which would merge all later
  # list elements into one giant element).
  set(_SEMICOLON     "__CMAKE_SEMICOLON__")
  set(_BACKSLASH     "__CMAKE_BACKSLASH__")
  set(_LINESEP       "__CMAKE_LINESEP__")
  set(_OPENBRACKET   "__CMAKE_OPENBRACKET__")
  set(_CLOSEBRACKET  "__CMAKE_CLOSEBRACKET__")

  string(REPLACE ";"    "${_SEMICOLON}"    FILE_CONTENTS "${FILE_CONTENTS}")
  string(REPLACE "\\"   "${_BACKSLASH}"    FILE_CONTENTS "${FILE_CONTENTS}")
  string(REPLACE "["    "${_OPENBRACKET}"  FILE_CONTENTS "${FILE_CONTENTS}")
  string(REPLACE "]"    "${_CLOSEBRACKET}" FILE_CONTENTS "${FILE_CONTENTS}")
  string(REPLACE "\r\n" "\n"               FILE_CONTENTS "${FILE_CONTENTS}")
  string(REPLACE "\n"   "${_LINESEP}"      FILE_CONTENTS "${FILE_CONTENTS}")
  string(REPLACE "${_LINESEP}" ";" LINES "${FILE_CONTENTS}")

  set(GROUP_NAME "")
  set(BLOCK_LINES "")
  set(PREAMBLE "")

  foreach(LINE IN LISTS LINES)
    string(REPLACE "${_BACKSLASH}"    "\\" LINE "${LINE}")
    string(REPLACE "${_SEMICOLON}"    ";"  LINE "${LINE}")
    string(REPLACE "${_OPENBRACKET}"  "["  LINE "${LINE}")
    string(REPLACE "${_CLOSEBRACKET}" "]"  LINE "${LINE}")

    # Detect BEGIN_INSTANTIATION_GROUP
    if(LINE MATCHES "//[ ]*BEGIN_INSTANTIATION_GROUP\\(([A-Za-z0-9_]+)\\)")
      string(REGEX REPLACE
        ".*BEGIN_INSTANTIATION_GROUP\\(([A-Za-z0-9_]+)\\).*" "\\1"
        GROUP_NAME "${LINE}")
      set(BLOCK_LINES "")
      continue()
    endif()

    # Detect END_INSTANTIATION_GROUP -> write the file
    if(LINE MATCHES "//[ ]*END_INSTANTIATION_GROUP")
      set(GENERATED_FILE "${OUTPUT_DIR}/${REL_DIR}/${BASE_NAME}_${GROUP_NAME}.cpp")
      file(MAKE_DIRECTORY "${OUTPUT_DIR}/${REL_DIR}")
      file(WRITE "${GENERATED_FILE}" "${PREAMBLE}\n${BLOCK_LINES}\n")
      message(STATUS "${REL_DIR}/${BASE_NAME}_${GROUP_NAME}.cpp")
      set(GROUP_NAME "")
      continue()
    endif()

    # Accumulate block lines while inside a group
    if(NOT GROUP_NAME STREQUAL "")
      string(APPEND BLOCK_LINES "${LINE}\n")
      continue()
    endif()

    # Everything outside a group goes into the shared preamble
    if(LINE MATCHES "^[ \t]*#include[ ]+\"${HEADER_FILENAME}\"")
      string(APPEND PREAMBLE "#include <${HEADER_INCL_PATH}>\n")
    else()
      string(APPEND PREAMBLE "${LINE}\n")
    endif()

  endforeach()

endfunction()
