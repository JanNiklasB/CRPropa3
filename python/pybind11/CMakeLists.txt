cmake_minimum_required(VERSION 3.14)

project(crpropa-pybind11-bindings)

add_library(crpropa SHARED IMPORTED)
set_target_properties(
	crpropa 
	PROPERTIES 
	IMPORTED_LOCATION "${CRPropa_BINARY_DIR}/libcrpropa.so"
)

set(Python_FIND_VIRTUALENV "FIRST")
set(Python_FIND_STRATEGY "LOCATION")
find_package(Python 3.0 COMPONENTS Interpreter Development NumPy)

# build it on pybind11:
find_package(pybind11)
if(NOT pybind11_FOUND)
	include(FetchContent)
	if(NOT pybind11_POPULATED)
		FetchContent_Populate(
			pybind11
			GIT_REPOSITORY https://github.com/RosettaCommons/pybind11.git
		)
	endif()
		
	# set CMAKE_ARGS before add_subdirectory per set:
	set(PYBIND11_TEST OFF CACHE BOOL "")
	
	add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/createBindings.sh ${CMAKE_CURRENT_BINARY_DIR}/createBindings.sh @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/binder.cfg ${CMAKE_CURRENT_BINARY_DIR}/binder.cfg @ONLY)

#generate bindings:
execute_process(
	COMMAND sh ${CMAKE_CURRENT_BINARY_DIR}/createBindings.sh
)

file(READ ${CMAKE_CURRENT_BINARY_DIR}/bindings/crpropa.sources SOURCES)

#found this technique to build the cmake list here:
#http://public.kitware.com/pipermail/cmake/2007-May/014236.html
#maybe there is a better way...
STRING(REGEX REPLACE ";" "\\\\;" SOURCES "${SOURCES}")
STRING(REGEX REPLACE "\n" ";" SOURCES "${SOURCES}")
list(TRANSFORM SOURCES PREPEND ${CMAKE_CURRENT_BINARY_DIR}/bindings/)

# add pybind modules:
pybind11_add_module(
	crpropa_pybind SHARED 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Vector3.cpp
	
	${SOURCES}
)

target_include_directories(
	crpropa_pybind PUBLIC
	pybind11::headers
	${CRPropa_SOURCE_DIR}/include
	${CRPROPA_EXTRA_INCLUDES}
)
target_link_libraries(
	crpropa_pybind PUBLIC
	${Python_LIBRARIES}
	pybind11::module
	crpropa
)
set_target_properties(
	crpropa_pybind 
	PROPERTIES 
	LIBRARY_OUTPUT_DIRECTORY "${CRPropa_BINARY_DIR}/crpropa-pybind"
	OUTPUT_NAME "crpropa"
	PREFIX "" SUFFIX ".so"
)

# move __init__.py and other needed files to build folder for tests and usage without installation
file(
	COPY
	"${CRPropa_SOURCE_DIR}/python/crpropa/"
	DESTINATION
	"${CRPropa_BINARY_DIR}/crpropa-pybind"
)