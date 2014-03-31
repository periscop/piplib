# Try to find the piplib library

# PIPLIB_OSL_INT_FOUND       - System has piplib lib
# PIPLIB_OSL_INT_INCLUDE_DIR - The piplib include directory
# PIPLIB_OSL_INT_LIBRARY     - Library needed to use piplib


if (PIPLIB_OSL_INT_INCLUDE_DIR AND PIPLIB_OSL_INT_LIBRARY)
	# Already in cache, be silent
	set(PIPLIB_OSL_INT_FIND_QUIETLY TRUE)
endif()

find_path(PIPLIB_OSL_INT_INCLUDE_DIR NAMES piplib/piplib_osl_int.h)
find_library(PIPLIB_OSL_INT_LIBRARY NAMES piplib_osl_int)

if (PIPLIB_OSL_INT_LIBRARY AND PIPLIB_OSL_INT_INCLUDE_DIR)
	message(STATUS "Library piplib_osl_int found =) ${PIPLIB_OSL_INT_LIBRARY}")
else()
	message(STATUS "Library piplib_osl_int not found =(")
endif()


include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PIPLIB_OSL_INT DEFAULT_MSG PIPLIB_OSL_INT_INCLUDE_DIR PIPLIB_OSL_INT_LIBRARY)

mark_as_advanced(PIPLIB_OSL_INT_INCLUDE_DIR PIPLIB_OSL_INT_LIBRARY)
