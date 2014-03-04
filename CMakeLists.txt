cmake_minimum_required(VERSION 2.8)


set(VERSION "1.4.0")


# User's settings - C Flags

# 	set(release "TRUE")
	set(release "FALSE")

	# Release
	if (release)
		set(CMAKE_C_FLAGS "-O3")
	# Debug # valgrind --show-reachable=yes --leak-check=full -v exe
	else()
		set(CMAKE_C_FLAGS "-O0 -g3")
	endif()

# User's settings - General C Flags
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -std=c99 -pedantic")


# Build documentation
	find_program(texi2pdf_exe texi2pdf)
	if(texi2pdf_exe)
		file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/doc DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
		configure_file("${CMAKE_CURRENT_BINARY_DIR}/doc/gitversion.texi.in" "${CMAKE_CURRENT_BINARY_DIR}/doc/gitversion.texi")
		add_custom_target(
			doc
			${texi2pdf_exe} piplib.texi --output=${CMAKE_CURRENT_BINARY_DIR}/piplib.pdf
			WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doc"
			COMMENT "Generating documentation (pdf) (with texi2pdf)" VERBATIM
		)
	else()
		message (STATUS "texi2pdf not found :( Documentation can not be built")
	endif()


# osl

	find_package(osl REQUIRED)

# GMP
	message(STATUS "---")
	find_library(gmp_LIB gmp)
	if (gmp_LIB)
		message (STATUS "Library gmp found =) ${gmp_LIB}")
		add_definitions(-DOSL_GMP_IS_HERE)
	else()
		message (STATUS "Library gmp not found :(")
	endif()

# Include directories (to use #include <> instead of #include "")

	# include/osl/scop.h
	configure_file("./include/piplib/version.h.in" "./include/piplib/version.h")
	include_directories("${CMAKE_CURRENT_BINARY_DIR}/include/piplib")
	# piplib
	set(pip_INCLUDE "./include")
	message(STATUS "Include pip = ${pip_INCLUDE}")
	include_directories("${pip_INCLUDE}")


# Compiler log
	message(STATUS "---")
	message(STATUS "C compiler = ${CMAKE_C_COMPILER}")
	if (release)
		message(STATUS "Mode Release")
	else()
		message(STATUS "Mode Debug")
	endif()
	message(STATUS "C flags    = ${CMAKE_C_FLAGS}")


# Library

	# files .c
	file(
		GLOB_RECURSE
		sources
		source/*
	)

	message(STATUS "---")

	# Dynamic
	message(STATUS "Add library piplib (shared)")
	add_library(
		piplib
		SHARED
		${sources}
	)

	# Static
	message(STATUS "Add library piplib (static)")
	add_library(
		piplib_static
		STATIC
		${sources}
	)
	set_target_properties(piplib_static PROPERTIES OUTPUT_NAME piplib)
	target_link_libraries(piplib_static ${OSL_LIBRARY})


# Executables & tests

	message(STATUS "---")

	set(exe "pip")
	message(STATUS "Add executable ${exe}")
	add_executable(${exe} "source/maind.c")
	target_link_libraries(${exe} piplib_static)

	set(exe "exemple")
	message(STATUS "Add executable ${exe}")
	add_executable(${exe} "example/example.c")
	target_link_libraries(${exe} piplib_static)

	# Test

	message(STATUS "---")
	
	find_package(PythonInterp REQUIRED)
	
		if (PYTHONINTERP_FOUND)

		enable_testing()

		file(
			GLOB_RECURSE
			examples
			example/*.pip
		)

		foreach(test ${examples}) # Tests in test directory are ignored...
		
			message(STATUS "Add Test ${test}")
			STRING(REPLACE ".pip" ".ll" test_output "${test}")
			
			add_test(
				"${test}_pip"
				"${CMAKE_CURRENT_SOURCE_DIR}/example/check_source_result.py"
				"${test}"
				"${test_output}"
				"${CMAKE_CURRENT_BINARY_DIR}/exemple"
			)
			
		endforeach()

	endif()


# Install

	install(TARGETS piplib LIBRARY DESTINATION lib)
	install(TARGETS piplib_static ARCHIVE DESTINATION lib)
	install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
	install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include/" DESTINATION include FILES_MATCHING PATTERN "*.h")
	install(FILES piplib-config.cmake DESTINATION lib/piplib)


# Little help

	message(STATUS "---")
	message(STATUS "You can execute:")
	message(STATUS "    make         # To compile osl library")
	if (PYTHONINTERP_FOUND)
		message(STATUS "    make test    # To execute tests")
	endif()
	message(STATUS "    make install # To install library, include and CMake module")
	message(STATUS "                 # If you need root access:")
	message(STATUS "                 #     sudo make install")
	message(STATUS "                 #     su -c \"make install\"")
	if(texi2pdf_exe)
		message(STATUS "    make doc     # To generate the documentation")
	endif()

	message(STATUS "---")
