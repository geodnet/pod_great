cmake_minimum_required(VERSION 3.5)
project(libGREAT VERSION 0.1.0)

# Gather source files by subdirectory
file(GLOB gall_h       gall/*.h)
file(GLOB gall_cpp     gall/*.cpp)
file(GLOB gambfix_h    gambfix/*.h)
file(GLOB gambfix_cpp  gambfix/*.cpp)
file(GLOB gcoders_h    gcoders/*.h)
file(GLOB gcoders_cpp  gcoders/*.cpp)
file(GLOB gdata_h      gdata/*.h)
file(GLOB gdata_cpp    gdata/*.cpp)
file(GLOB gforce_h     gforce/*.h)
file(GLOB gforce_cpp   gforce/*.cpp)
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
file(GLOB exp_h        ExportHeaders/*.h)

# Collect all
set(header_files ${exp_h} ${gall_h} ${gambfix_h} ${gcoders_h} ${gdata_h} ${gforce_h} ${gio_h} ${gmodels_h} ${gproc_h} ${gset_h} ${gutils_h})
set(source_files ${gall_cpp} ${gambfix_cpp} ${gcoders_cpp} ${gdata_cpp} ${gforce_cpp} ${gio_cpp} ${gmodels_cpp} ${gproc_cpp} ${gset_cpp} ${gutils_cpp})

# Visual Studio source groups
source_group("ExportHeaders"           FILES ${exp_h})
source_group("Header Files\\gall"      FILES ${gall_h})
source_group("Header Files\\gambfix"   FILES ${gambfix_h})
source_group("Header Files\\gcoders"   FILES ${gcoders_h})
source_group("Header Files\\gdata"     FILES ${gdata_h})
source_group("Header Files\\gforce"    FILES ${gforce_h})
source_group("Header Files\\gio"       FILES ${gio_h})
source_group("Header Files\\gmodels"   FILES ${gmodels_h})
source_group("Header Files\\gproc"     FILES ${gproc_h})
source_group("Header Files\\gset"      FILES ${gset_h})
source_group("Header Files\\gutils"    FILES ${gutils_h})
source_group("Source Files\\gall"      FILES ${gall_cpp})
source_group("Source Files\\gambfix"   FILES ${gambfix_cpp})
source_group("Source Files\\gcoders"   FILES ${gcoders_cpp})
source_group("Source Files\\gdata"     FILES ${gdata_cpp})
source_group("Source Files\\gforce"    FILES ${gforce_cpp})
source_group("Source Files\\gio"       FILES ${gio_cpp})
source_group("Source Files\\gmodels"   FILES ${gmodels_cpp})
source_group("Source Files\\gproc"     FILES ${gproc_cpp})
source_group("Source Files\\gset"      FILES ${gset_cpp})
source_group("Source Files\\gutils"    FILES ${gutils_cpp})

# Static library
add_library(${PROJECT_NAME} STATIC ${header_files} ${source_files})

# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC ${libGREATSrc})

# Depend on libGNut for G-Nut base types (circular dep, valid for static libs)
target_link_libraries(${PROJECT_NAME} PUBLIC libGNut)

# Define all DLL export macros (must match old per-library definitions for MSVC template instantiation)
target_compile_definitions(${PROJECT_NAME} PRIVATE
    GREAT_BASE_LIBRARY
    GREAT_OI_LIBRARY
    LSQ_BASE_LIBRARY
    GREAT_MAT_LIBRARY
    GREAT_AMB_LIBRARY)

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
    endif()
endif()
