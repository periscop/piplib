# Try to find the piplib library

# PIPLIB_OSL_FOUND       - System has piplib lib
# PIPLIB_OSL_INCLUDE_DIR - The piplib include directory
# PIPLIB_OSL_LIBRARY     - Library needed to use piplib


if (PIPLIB_OSL_INCLUDE_DIR AND PIPLIB_OSL_LIBRARY)
	# Already in cache, be silent
	set(PIPLIB_OSL_FIND_QUIETLY TRUE)
endif()

find_path(PIPLIB_OSL_INCLUDE_DIR NAMES piplib/piplib_OSL.h)
find_library(PIPLIB_OSL_LIBRARY NAMES piplib_OSL)

if (PIPLIB_OSL_LIBRARY AND PIPLIB_OSL_INCLUDE_DIR)
	message(STATUS "Library piplib_osl found =) ${PIPLIB_OSL_LIBRARY}")
else()
	message(STATUS "Library piplib_osl not found =(")
endif()


include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PIPLIB_OSL DEFAULT_MSG PIPLIB_OSL_INCLUDE_DIR PIPLIB_OSL_LIBRARY)

mark_as_advanced(PIPLIB_OSL_INCLUDE_DIR PIPLIB_OSL_LIBRARY)
