# Script-mode entry point for build-time instantiation file generation.
#
# Invoked by add_custom_command as:
#   cmake -DSOURCE_FILE=<abs>
#         -DSOURCE_ROOT=<abs>
#         -DINCLUDE_ROOT=<abs>
#         -DOUTPUT_DIR=<abs>
#         -P generate_instantiation.cmake

include("${CMAKE_CURRENT_LIST_DIR}/instantiations.cmake")

generate_instantiation_files(
  "${SOURCE_FILE}"
  "${SOURCE_ROOT}"
  "${INCLUDE_ROOT}"
  "${OUTPUT_DIR}"
)
