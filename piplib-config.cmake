# Try to find the piplib library

# PIPLIB_FOUND       - System has piplib lib
# PIPLIB_INCLUDE_DIR - The piplib include directory
# PIPLIB_LIBRARY     - Library needed to use piplib


if (PIPLIB_INCLUDE_DIR AND PIPLIB_LIBRARY)
	# Already in cache, be silent
	set(PIPLIB_FIND_QUIETLY TRUE)
endif()

find_path(PIPLIB_INCLUDE_DIR NAMES piplib/PIPLIB.h)
find_library(PIPLIB_LIBRARY NAMES PIPLIB)

if (PIPLIB_LIBRARY AND PIPLIB_INCLUDE_DIR)
	message(STATUS "Library PIPLIB found =) ${PIPLIB_LIBRARY}")
else()
	message(STATUS "Library PIPLIB not found =(")
endif()


include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PIPLIB DEFAULT_MSG PIPLIB_INCLUDE_DIR PIPLIB_LIBRARY)

mark_as_advanced(PIPLIB_INCLUDE_DIR PIPLIB_LIBRARY)
