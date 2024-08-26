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



#-- Search for oneMKL
find_package(MKL CONFIG REQUIRED)
#package arpack-ng
#find_package(arpackng CONFIG REQUIRED)

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
	$<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES> # oneMKL
	${CMAKE_CURRENT_SOURCE_DIR}/dependencies/exprtk
	${CMAKE_CURRENT_SOURCE_DIR}/dependencies/eigen
	${CMAKE_CURRENT_SOURCE_DIR}/dependencies/vcpkg/packages/arpack-ng_x64-windows/include/arpack
)

target_link_libraries(${PROJECT_NAME} 
PUBLIC 
	$<LINK_ONLY:MKL::MKL> 
	${CMAKE_CURRENT_SOURCE_DIR}/dependencies/vcpkg/packages/arpack-ng_x64-windows/lib/libarpack.lib
	)

################################################################################
# Compile definitions
################################################################################
target_compile_definitions(${PROJECT_NAME} PRIVATE
    "_CRT_SECURE_NO_WARNINGS;"
    "_MBCS"
)

##
## Define different options for oneMKL or DLLs
##
#set(ADDITIONAL_LIBRARY_DEPENDENCIES "arpack_win64")
	
target_compile_options(${PROJECT_NAME} 
  PUBLIC 
	$<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>
  PRIVATE
	$<$<CXX_COMPILER_ID:MSVC>:/W3> # '/WX' to treat warnings as erros in VS
	$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -pedantic> # '-Werror' to treat warnings as errors in other compilers
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
	/W3;
	${DEFAULT_CXX_DEBUG_INFORMATION_FORMAT};
	${DEFAULT_CXX_EXCEPTION_HANDLING};
	/openmp
	)


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

