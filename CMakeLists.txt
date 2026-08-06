cmake_minimum_required(VERSION 3.20)

################################################################################
# Global configuration types
################################################################################
set(CMAKE_CONFIGURATION_TYPES
    "Debug"
    "Release"
    CACHE STRING "" FORCE
)

################################################################################
# Global compiler options
################################################################################
if(MSVC)
    # remove default flags provided with CMake for MSVC
    set(CMAKE_CXX_FLAGS "")
    set(CMAKE_CXX_FLAGS_DEBUG "")
    set(CMAKE_CXX_FLAGS_RELEASE "")
endif()

################################################################################
# Global linker options
################################################################################
if(MSVC)
    # remove default flags provided with CMake for MSVC
    set(CMAKE_EXE_LINKER_FLAGS "")
    set(CMAKE_MODULE_LINKER_FLAGS "")
    set(CMAKE_SHARED_LINKER_FLAGS "")
    set(CMAKE_STATIC_LINKER_FLAGS "")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS_DEBUG "${CMAKE_MODULE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_STATIC_LINKER_FLAGS_DEBUG "${CMAKE_STATIC_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS_RELEASE "${CMAKE_MODULE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_STATIC_LINKER_FLAGS_RELEASE "${CMAKE_STATIC_LINKER_FLAGS}")
endif()

## Change default oneMKL interface
set(MKL_INTERFACE lp64)

## specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(PROJECT_NAME Giraffe)
## Project name and version
project(Giraffe VERSION 2.0.137 LANGUAGES CXX)

## Change startup project to Giraffe
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Giraffe)


## If true, the computer must have oneMKL installed, 
# Otherwise, it is assumed that MKL DLLs will be used in .exe 
# directory or its folder is already added to system path 
option(USE_oneMKL "MKL -> oneMKL from Intel oneAPI" ON)


#set(EIGEN_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/eigen/)
#include_directories(${EIGEN_ROOT_DIR})



#-- Search for math/eigensolver dependencies (platform specific)
if(WIN32)
	#-- Search for oneMKL
	find_package(MKL CONFIG REQUIRED)
else()
	#-- macOS/Linux: use OpenBLAS (BLAS+LAPACK) and arpack-ng from Homebrew instead of oneMKL/vcpkg
	if(APPLE)
		# openblas/arpack/libomp are keg-only on macOS (Accelerate.framework owns the
		# default BLAS/LAPACK slot), so their .pc files aren't on the default
		# pkg-config search path. Add the Homebrew keg paths for both Apple Silicon
		# and Intel prefixes before looking them up.
		set(_giraffe_brew_pkgconfig_dirs "")
		foreach(_prefix /opt/homebrew /usr/local)
			foreach(_pkg openblas arpack)
				if(EXISTS "${_prefix}/opt/${_pkg}/lib/pkgconfig")
					list(APPEND _giraffe_brew_pkgconfig_dirs "${_prefix}/opt/${_pkg}/lib/pkgconfig")
				endif()
			endforeach()
		endforeach()
		if(_giraffe_brew_pkgconfig_dirs)
			string(REPLACE ";" ":" _giraffe_brew_pkgconfig_path "${_giraffe_brew_pkgconfig_dirs}")
			set(ENV{PKG_CONFIG_PATH} "${_giraffe_brew_pkgconfig_path}:$ENV{PKG_CONFIG_PATH}")
		endif()
	endif()

	find_package(PkgConfig REQUIRED)
	pkg_check_modules(OPENBLAS REQUIRED IMPORTED_TARGET openblas)
	pkg_check_modules(ARPACK REQUIRED IMPORTED_TARGET arpack)

	# Homebrew's libomp is keg-only (no pkg-config/cmake package), so locate it by hand.
	find_path(LIBOMP_INCLUDE_DIR omp.h
		HINTS /opt/homebrew/opt/libomp/include /usr/local/opt/libomp/include)
	find_library(LIBOMP_LIBRARY omp
		HINTS /opt/homebrew/opt/libomp/lib /usr/local/opt/libomp/lib)
	if(NOT LIBOMP_INCLUDE_DIR OR NOT LIBOMP_LIBRARY)
		message(FATAL_ERROR "libomp not found. Install it with 'brew install libomp'.")
	endif()
endif()

## Organize files in solution explorer with a "filter tree"
FILE(GLOB_RECURSE files
	${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
	)
#-- Create the tree in the VS solution explorer
set(TREE_FILES ${files})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${TREE_FILES})

## Configure a header file to pass some of the CMake settings to the source code
#configure_file(${PROJECT_NAME}Config.h.in ${PROJECT_NAME}Config.h)


## Other files to include in Visual Studio solution explorer
set(FILES_TO_ADD README.md install_dependencies.bat build.bat)

## Create the executable


################################################################################
# Target
################################################################################
add_executable(${PROJECT_NAME} ${TREE_FILES} ${FILES_TO_ADD})
#use_props(${PROJECT_NAME} "${CMAKE_CONFIGURATION_TYPES}" "${DEFAULT_CXX_PROPS}")
set(ROOT_NAMESPACE Giraffe)

################################################################################
# Output directory
################################################################################
set_target_properties(${PROJECT_NAME} PROPERTIES
    OUTPUT_DIRECTORY_DEBUG   "${CMAKE_SOURCE_DIR}/bin/${CMAKE_VS_PLATFORM_NAME}/$<CONFIG>/"
    OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/bin/${CMAKE_VS_PLATFORM_NAME}/$<CONFIG>/"
)
set_target_properties(${PROJECT_NAME} PROPERTIES
    INTERPROCEDURAL_OPTIMIZATION_RELEASE "TRUE"
)
################################################################################
# Include directories
################################################################################
target_include_directories(${PROJECT_NAME} PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/dependencies/exprtk
	${CMAKE_CURRENT_SOURCE_DIR}/dependencies/eigen
)

if(WIN32)
	target_include_directories(${PROJECT_NAME} PUBLIC
		$<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES> # oneMKL
		${CMAKE_CURRENT_SOURCE_DIR}/dependencies/vcpkg/packages/arpack-ng_x64-windows/include/arpack
	)
	target_link_libraries(${PROJECT_NAME}
	PUBLIC
		$<LINK_ONLY:MKL::MKL>
		${CMAKE_CURRENT_SOURCE_DIR}/dependencies/vcpkg/packages/arpack-ng_x64-windows/lib/libarpack.lib
	)
	target_compile_definitions(${PROJECT_NAME} PRIVATE "GIRAFFE_USE_MKL")
else()
	target_include_directories(${PROJECT_NAME} PUBLIC ${LIBOMP_INCLUDE_DIR})
	target_link_libraries(${PROJECT_NAME}
	PUBLIC
		PkgConfig::ARPACK
		PkgConfig::OPENBLAS
		${LIBOMP_LIBRARY}
	)
endif()

################################################################################
# Compile definitions
################################################################################
if(MSVC)
target_compile_definitions(${PROJECT_NAME} PRIVATE
    "_CRT_SECURE_NO_WARNINGS;"
    "_MBCS"
)
endif()

##
## Define different options for oneMKL or DLLs
##
#set(ADDITIONAL_LIBRARY_DEPENDENCIES "arpack_win64")
	
if(MSVC)
	target_compile_options(${PROJECT_NAME}
	  PUBLIC
		$<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>
	  PRIVATE
		/W3 # '/WX' to treat warnings as erros in VS
		$<$<CONFIG:Debug>:
			/Od;
			/sdl;
			/bigobj
		>
		$<$<CONFIG:Release>:
			/O2;
			/Oi;
			/sdl-;
			/Gy
		>
		/MP;
		${DEFAULT_CXX_DEBUG_INFORMATION_FORMAT};
		${DEFAULT_CXX_EXCEPTION_HANDLING};
		/openmp
		)
else()
	target_compile_options(${PROJECT_NAME}
	  PRIVATE
		-Wall -Wextra -pedantic # '-Werror' to treat warnings as errors
		# Legacy codebase was only ever built with MSVC, which is far more lax
		# about these; silence the categories that just add noise without
		# flagging real bugs, keep the rest.
		-Wno-unused-variable
		-Wno-unused-but-set-variable
		-Wno-unused-parameter
		-Wno-sign-compare
		-Wno-format-security
		-Wno-deprecated-declarations
		-Wno-invalid-utf8
		-Xpreprocessor -fopenmp
		# The codebase was only ever built with MSVC, whose standard headers
		# transitively pull in <cstdlib>/<cfloat>/etc. libc++/libstdc++ don't,
		# so force-include the handful of C headers many files rely on without
		# including directly (atoi/atof, FLT_MAX/DBL_EPSILON, ...).
		"SHELL:-include cstdlib"
		"SHELL:-include cfloat"
		"SHELL:-include cstring"
		"SHELL:-include cstdio"
		$<$<CONFIG:Debug>:-Og;-g>
		$<$<CONFIG:Release>:-O2>
		)
endif()


## Compiler flags
if (MSVC)
	# Precompiled header 
    #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /YuPCH.h /FIPCH.h")
	#set_source_files_properties(src/PCH.cpp	PROPERTIES COMPILE_FLAGS "/YcPCH.h")
	# Change the Visual Studio working directory
	set_target_properties(${PROJECT_NAME} PROPERTIES
		#VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/inputs")
		VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()

