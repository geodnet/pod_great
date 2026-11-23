#Minimum requirement of CMake version : 3.5.0
cmake_minimum_required(VERSION 3.5)

#Project name and version number
project(great VERSION 0.1.0)

set(CMAKE_EXPORT_COMPILE_COMMANDS True)

# C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/std:c++14)
endif()

option(USE_OPENMP  "using openmp for parallle" OFF)

if (USE_OPENMP)
  add_definitions(-DUSE_OPENMP)
endif()

# important setting path
if(CMAKE_SYSTEM_NAME MATCHES "Windows")
	message(STATUS "==================    Windows  ヾ(ｏ･ω･)ﾉ  nice ======================")
	message(STATUS "                                                                     ")
	find_path(zlib_ROOT        HINTS "${zlib_ROOT}"         "$ENV{zlib_ROOT}             ")
	find_path(Eigen_ROOT       HINTS "${Eigen_ROOT}"        "$ENV{Eigen_ROOT}            ")
else()
	message(STATUS "==================    Not Windows o(▼皿▼メ;)o   ======================")
	message(STATUS "                                                                      ")
	set(zlib_ROOT        ../third-party/zlib/zlib-1.2.11/zlib)
	set(Eigen_ROOT       ../third-party/Eigen)
endif()

include_directories(${zlib_ROOT})
include_directories(${Eigen_ROOT})

# global value
# set lib name
set( libGNut    libGNut       )
set( libGREAT   libGREAT      )

# set app name
set( podflt     great_podflt              )

# set lib path
set( libGNutSrc  ${PROJECT_SOURCE_DIR}/libGNut)
set( libGREATSrc ${PROJECT_SOURCE_DIR}/libGREAT)

message(STATUS "====================================================")
message(STATUS "   (<ゝω・)☆ check the important path here          ")
message(STATUS "====================================================")
message(STATUS "The path of  zlib        is : ${zlib_ROOT}          ")
message(STATUS "The path of  Eigen       is : ${Eigen_ROOT}         ")
message(STATUS "                                                    ")
message(STATUS "The ROOT path is  : ${PROJECT_SOURCE_DIR} ")
message(STATUS "The sorce file for lib GNut  is : ${libGNutSrc}     ")
message(STATUS "The sorce file for lib GREAT is : ${libGREATSrc}    ")
message(STATUS "                                                    ")
message(STATUS "====================================================")

# Compile default settings
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
	add_compile_options(-w)
	if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64")
		add_compile_options(-m64 -mavx -mavx2 -mfma)
	endif()
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
	add_compile_options(-pthread)
elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
	add_compile_options(/bigobj)
	add_compile_options(/MP)
	add_compile_options(/w)
endif()

#Set the ROOT and subdirectory
set(ROOT ${PROJECT_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/Bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Lib)

set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_RELEASE_POSTFIX "")
set(CMAKE_RELWITHDEBINFO_POSTFIX "rd")
set(CMAKE_MINSIZEREL_POSTFIX "s")

if(CMAKE_BUILD_TYPE MATCHES "Release")
	set(CMAKE_BUILD_POSTFIX "${CMAKE_RELEASE_POSTFIX}")
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
	set(CMAKE_BUILD_POSTFIX "${CMAKE_DEBUG_POSTFIX}")
elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
	set(CMAKE_BUILD_POSTFIX "${CMAKE_RELWITHDEBINFO_POSTFIX}")
elseif(CMAKE_BUILD_TYPE MATCHES "MinSizeRel")
	set(CMAKE_BUILD_POSTFIX "${CMAKE_MINSIZEREL_POSTFIX}")
else()
	set(CMAKE_BUILD_POSTFIX "")
endif()

# For lib source files (libGREAT first for circular dependency resolution)
add_subdirectory(${libGREATSrc}     ${CMAKE_BINARY_DIR}/${libGREAT})
add_subdirectory(${libGNutSrc}      ${CMAKE_BINARY_DIR}/${libGNut})

# podflt application
if(EXISTS "${ROOT}/app/great_podflt")
    message(STATUS "Found podflt application directory")
    add_subdirectory(${ROOT}/app/great_podflt ${CMAKE_BINARY_DIR}/${podflt})
else()
    message(WARNING "podflt application directory not found at: ${ROOT}/app/great_podflt")
endif()

# group
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
SET_PROPERTY(TARGET ${libGNut}      PROPERTY FOLDER "LIB")
SET_PROPERTY(TARGET ${libGREAT}     PROPERTY FOLDER "LIB")
if(TARGET ${podflt})
    SET_PROPERTY(TARGET ${podflt}        PROPERTY FOLDER "POD_APP")
endif()

#Output Messages for debug the Cmake
message(STATUS "operation system is : ${CMAKE_SYSTEM}")
message(STATUS "current platform is : ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMake version is    : ${CMAKE_SYSTEM_VERSION}")
message(STATUS "C compiler is : ${CMAKE_C_COMPILER}")
message(STATUS "C++ compiler is : ${CMAKE_CXX_COMPILER}")
message(STATUS "The program main directory is : ${ROOT}")