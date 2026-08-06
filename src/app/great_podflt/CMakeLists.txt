#Cmake 最低版本号要求
cmake_minimum_required(VERSION 3.5)

#项目信息
project(${podflt})

file(GLOB header_files 	*.h)
file(GLOB source_files	*.cpp)

source_group("CMake Files" FILES CMakeLists.txt)
source_group("Header Files" FILES header_files)
source_group("Soruce Files" FILES source_files)


set(include_path
	${libGNutSrc}
    ${libGREATSrc}
	${zlib_ROOT})
include_directories(${include_path})

add_executable(${PROJECT_NAME} ${header_files} ${source_files})

if(CMAKE_SYSTEM_NAME MATCHES "Windows")
	set(link_path
		${CMAKE_BINARY_DIR}/Lib/Debug
		${CMAKE_BINARY_DIR}/Lib/Release
		${CMAKE_BINARY_DIR}/Lib/RelWithDebInfo
		${CMAKE_BINARY_DIR}/Lib/MinSizeRel)
	link_directories(${link_path} )
else()
	set(link_path
		${CMAKE_BINARY_DIR}/Lib)
	link_directories(${link_path} )
endif()


set(lib_list
	${libGNut}
	${libGREAT})
target_link_libraries(${PROJECT_NAME} ${lib_list})

add_dependencies(${PROJECT_NAME} ${lib_list})
