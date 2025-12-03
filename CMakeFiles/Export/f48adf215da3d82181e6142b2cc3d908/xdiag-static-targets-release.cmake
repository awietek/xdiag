#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "xdiag::xdiag" for configuration "Release"
set_property(TARGET xdiag::xdiag APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(xdiag::xdiag PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libxdiag.a"
  )

list(APPEND _cmake_import_check_targets xdiag::xdiag )
list(APPEND _cmake_import_check_files_for_xdiag::xdiag "${_IMPORT_PREFIX}/lib64/libxdiag.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
