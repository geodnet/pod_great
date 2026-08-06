cmake_minimum_required(VERSION 3.5)
project(libGNut VERSION 0.1.0)

# Gather source files by subdirectory
file(GLOB gall_h       gall/*.h)
file(GLOB gall_cpp     gall/*.cpp)
file(GLOB gcoders_h    gcoders/*.h)
file(GLOB gcoders_cpp  gcoders/*.cpp)
file(GLOB gdata_h      gdata/*.h)
file(GLOB gdata_cpp    gdata/*.cpp)
file(GLOB gio_h        gio/*.h)
file(GLOB gio_cpp      gio/*.cpp)
file(GLOB gmodels_h    gmodels/*.h)
file(GLOB gmodels_cpp  gmodels/*.cpp)
file(GLOB gproc_h      gproc/*.h)
file(GLOB gproc_cpp    gproc/*.cpp)
file(GLOB gset_h       gset/*.h)
file(GLOB gset_cpp     gset/*.cpp)
file(GLOB gutils_h     gutils/*.h)
file(GLOB gutils_cpp   gutils/*.cpp)
file(GLOB newmat_h     newmat/*.h)
file(GLOB newmat_cpp   newmat/*.cpp)
file(GLOB pugixml_h    pugixml/src/*.hpp)
file(GLOB pugixml_cpp  pugixml/src/*.cpp)

# Collect all
set(header_files ${gall_h} ${gcoders_h} ${gdata_h} ${gio_h} ${gmodels_h} ${gproc_h} ${gset_h} ${gutils_h} ${newmat_h} ${pugixml_h})
set(source_files ${gall_cpp} ${gcoders_cpp} ${gdata_cpp} ${gio_cpp} ${gmodels_cpp} ${gproc_cpp} ${gset_cpp} ${gutils_cpp} ${newmat_cpp} ${pugixml_cpp})

# Visual Studio source groups
source_group("Header Files\\gall"      FILES ${gall_h})
source_group("Header Files\\gcoders"   FILES ${gcoders_h})
source_group("Header Files\\gdata"     FILES ${gdata_h})
source_group("Header Files\\gio"       FILES ${gio_h})
source_group("Header Files\\gmodels"   FILES ${gmodels_h})
source_group("Header Files\\gproc"     FILES ${gproc_h})
source_group("Header Files\\gset"      FILES ${gset_h})
source_group("Header Files\\gutils"    FILES ${gutils_h})
source_group("Header Files\\newmat"    FILES ${newmat_h})
source_group("Header Files\\pugixml"   FILES ${pugixml_h})
source_group("Source Files\\gall"      FILES ${gall_cpp})
source_group("Source Files\\gcoders"   FILES ${gcoders_cpp})
source_group("Source Files\\gdata"     FILES ${gdata_cpp})
source_group("Source Files\\gio"       FILES ${gio_cpp})
source_group("Source Files\\gmodels"   FILES ${gmodels_cpp})
source_group("Source Files\\gproc"     FILES ${gproc_cpp})
source_group("Source Files\\gset"      FILES ${gset_cpp})
source_group("Source Files\\gutils"    FILES ${gutils_cpp})
source_group("Source Files\\newmat"    FILES ${newmat_cpp})
source_group("Source Files\\pugixml"   FILES ${pugixml_cpp})

# Static library
add_library(${PROJECT_NAME} STATIC ${header_files} ${source_files})

# Define DLL export macros (needed by pugixml, newmat, and G-Nut headers with BASE_LIBRARY_EXPORT)
target_compile_definitions(${PROJECT_NAME} PRIVATE
    GREAT_BASE_LIBRARY
    GREAT_MAT_LIBRARY)

# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC ${libGNutSrc})

# External dependencies
include_directories(${Eigen_ROOT})
include_directories(${zlib_ROOT})

# Depend on libGREAT for WHU extension headers (circular dep, valid for static libs)
target_link_libraries(${PROJECT_NAME} PUBLIC libGREAT)

# Platform-specific zlib linking
if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    link_directories(${CMAKE_BINARY_DIR}/Lib/Debug)
    link_directories(${CMAKE_BINARY_DIR}/Lib/Release)
    link_directories(${CMAKE_BINARY_DIR}/Lib/RelWithDebInfo)
    link_directories(${CMAKE_BINARY_DIR}/Lib/MinSizeRel)
    target_link_libraries(${PROJECT_NAME} PUBLIC
        $<$<CONFIG:Debug>:${zlib_ROOT}/zlibstaticd.lib>
        $<$<NOT:$<CONFIG:Debug>>:${zlib_ROOT}/zlibstatic.lib>)
else()
    link_directories(${CMAKE_BINARY_DIR}/Lib)
    link_directories(${zlib_ROOT})
    target_link_libraries(${PROJECT_NAME} PUBLIC z)
endif()

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
    endif()
endif()
