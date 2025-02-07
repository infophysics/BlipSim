cmake_minimum_required(VERSION 3.19 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 20)
add_definitions(-w)

project(Blip)

set(CMAKE_BUILD_TYPE Debug)

##############################################################
# Debug mode
option(DEBUG_MODE "Compile with debug options" OFF)
if(DEBUG_MODE)
    message(STATUS "Compiling in Debug mode.")
    add_compile_options(-g)
endif()
##############################################################

##############################################################
# profiling mode
option(PROFILING_MODE "Compile with profiling turned on" ON)
if(PROFILING_MODE)
    add_definitions(-DBLIP_PROFILING)
    message(STATUS "Compiling with profiling mode turned on")
endif()
##############################################################

##############################################################
# Multithreading mode
option(MULTITHREADING "Compile with multithreading mode" OFF)
if(MULTITHREADING)
    if(DEFINED ENV{G4MULTITHREADED})
        add_definitions(-DBLIP_MULTITHREADED)
        message(STATUS "Compiling with multithreading enabled.")
    else()
        message(STATUS 
            "Attempt to use multithreading failed - 
            G4MULTITHREADED is not defined! 
            Maybe you didn't compile G4 with multithreading."
        )
    endif()
endif()
##############################################################

##############################################################
# Geant version
option(GEANT4_10 "Compile with Geant4 version 10" OFF)
if(GEANT4_10)
    add_definitions(-DBLIP_GEANT_10)
    message(STATUS "Compiling against Geant4 version 10.")
    find_package(Geant4 10.07.02 EXACT REQUIRED ui_all vis_all)
    include(${Geant4_USE_FILE})
else()
    message(STATUS "Compiling against Geant4 version 11.")
    find_package(Geant4 11 REQUIRED ui_all vis_all)
    include(${Geant4_USE_FILE})
endif()
message(STATUS "Geant4 installation directory -- ${Geant4_DIR}")
##############################################################

##############################################################
# Whether to use ROOT for old BLIP-I distribution
option(USE_ROOT "Use ROOT for Blip-I analysis" OFF)
if(USE_ROOT)
    message(STATUS "Compiling with ROOT.")
    find_package(ROOT 6.26 CONFIG REQUIRED COMPONENTS Core)
    message(STATUS "ROOT installation directory -- ${ROOT_DIR}")
    message(STATUS "ROOT include directories -- ${ROOT_INCLUDE_DIRS}")
    #message("ROOT Libraries: ${ROOT_LIBRARIES}")
    add_definitions(-DBLIP_ROOT)
endif()
##############################################################

##############################################################
# Whether to use yaml-cpp
option(USE_YAML "Use YAML for configuration files" ON)
if(USE_YAML)
    message(STATUS "Compiling with yaml-cpp.")
    add_definitions(-DBLIP_YAML)
    set(YAML_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/vendor/yaml-cpp/include)
    add_subdirectory(${PROJECT_SOURCE_DIR}/vendor/yaml-cpp)
    set(YAML_LIBS yaml-cpp)
else()
    set(YAML_INCLUDE_DIR "")
    set(YAML_LIBS "")
endif()
##############################################################

##############################################################
# use the visualization manager
option(USE_VIS "Whether to use the visualization." OFF)
if(USE_VIS)
    add_definitions(-DBLIP_USE_VIS)
    message(STATUS "Compiling with the visualization enabled.")
endif()
##############################################################

##############################################################
# use the ui executive/manager
option(USE_UI "Whether to use the ui interface." OFF)
if(USE_UI)
    add_definitions(-DBLIP_USE_UI)
    message(STATUS "Compiling with the ui enabled.")
endif()
##############################################################

##############################################################
# Get all header directories
MACRO(HEADER_DIRECTORIES return_list)
    FILE(GLOB_RECURSE new_list *.hh)
    SET(dir_list "")
    FOREACH(file_path ${new_list})
        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
        SET(dir_list ${dir_list} ${dir_path})
    ENDFOREACH()
    LIST(REMOVE_DUPLICATES dir_list)
    SET(${return_list} ${dir_list})
ENDMACRO()
HEADER_DIRECTORIES(header_list)
##############################################################

##############################################################
# generate sources and headers
include_directories(
    ${header_list}
    ${ROOT_INCLUDE_DIRS}
    ${YAML_INCLUDE_DIR}
)

file(GLOB_RECURSE sources 
    ${PROJECT_SOURCE_DIR}/src/*.cxx
)
file(GLOB_RECURSE headers 
    ${PROJECT_SOURCE_DIR}/include/*.hh
)
##############################################################

##############################################################
# Create library
add_library(blip SHARED ${sources} ${headers})
target_link_libraries(blip PUBLIC ${Geant4_LIBRARIES})

if(USE_ROOT)
    target_link_libraries(
        blip PUBLIC ROOT::Core ${ROOT_LIBRARIES}
    )
endif()
if(USE_YAML)
    target_link_libraries(
        blip PUBLIC ${YAML_LIBS}
    )
endif()
##############################################################

##############################################################
# create programs
add_executable(blip_run 
    src/programs/Blip.cc 
)
target_link_libraries(blip_run blip)
##############################################################

##############################################################
# create directories
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/macros/)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data/)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/analysis/)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/configs/)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/plots/)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/gdml/)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/)

##############################################################
# copy macro and data files
file(GLOB MACRO_FILES
    "macros/*.mac"
)
file(GLOB DATA_FILES
    "data/*.*"
)
file(GLOB ANALYSIS_FILES
    "analysis/*.py"
)
file(GLOB_RECURSE CONFIG_FILES
    "configs/*.yaml"
)
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR}/macros/)
file(COPY ${DATA_FILES} DESTINATION ${PROJECT_BINARY_DIR}/data/)
file(COPY ${ANALYSIS_FILES} DESTINATION ${PROJECT_BINARY_DIR}/analysis/)
file(COPY ${CONFIG_FILES} DESTINATION ${PROJECT_BINARY_DIR}/configs/)
##############################################################

##############################################################
add_custom_target(Blip DEPENDS blip)
##############################################################