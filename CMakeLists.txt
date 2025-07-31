cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(BUILD_PYBINDINGS True CACHE BOOL "Build pyMCVSPEC")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(mcvspec VERSION 1.0.0)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Generate WD mass radius relationship header
set(MR_TXT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/data/mass_radius.txt)
set(MR_HEADER_PATH ${CMAKE_CURRENT_BINARY_DIR}/mass_radius.hh)

add_custom_command(
    OUTPUT ${MR_HEADER_PATH}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/data/txt_to_header.py -i ${MR_TXT_PATH} -o ${MR_HEADER_PATH}
    DEPENDS ${MR_TXT_PATH}
    COMMENT "Generating header from ${MR_TXT_PATH}"
)

add_custom_target(embed_wd_mass_radius DEPENDS ${MR_HEADER_PATH})

if(NOT DEFINED ENV{HEADAS})
    message(FATAL_ERROR "HEADAS must be defined and point to your installation of HEASoft")
endif()

# gather info about HEASOFT env, reproduce logic of build naming
set(HEADAS "$ENV{HEADAS}/")
get_filename_component(HEASOFT_DIR "${HEADAS}" DIRECTORY)
execute_process (COMMAND bash -c
    ". ${HEASOFT_DIR}/BUILD_DIR/config.guess"
    OUTPUT_VARIABLE HEA_SYS_GUESS)
execute_process (COMMAND bash -c
    ". ${HEASOFT_DIR}/BUILD_DIR/config.sub ${HEA_SYS_GUESS}"
    OUTPUT_VARIABLE HEA_SYSTEM_NAME)
string(STRIP "${HEA_SYSTEM_NAME}" HEA_SYSTEM_NAME)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX
    "${CMAKE_BINARY_DIR}" CACHE PATH "MCVSPEC install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

if(BUILD_PYBINDINGS)
    find_package(pybind11 CONFIG REQUIRED)

    pybind11_add_module(_pymcvspec pymcvspec/bindings.cpp
                            source/Cataclysmic_Variable.cpp
                            source/integration.cpp)

    add_dependencies(_pymcvspec embed_wd_mass_radius)

    # MCVSPEC Headers
    target_include_directories(_pymcvspec PUBLIC "${PROJECT_SOURCE_DIR}/include" ${CMAKE_CURRENT_BINARY_DIR})

    # HEASOFT Headers
    target_include_directories(_pymcvspec PUBLIC
        "${HEASOFT_DIR}/Xspec/${HEA_SYSTEM_NAME}/include"
        "${HEASOFT_DIR}/Xspec/${HEA_SYSTEM_NAME}/include/XSFunctions"
        "${HEASOFT_DIR}/Xspec/${HEA_SYSTEM_NAME}/include/XSFunctions/Utilities"
        "${HEASOFT_DIR}/heacore/${HEA_SYSTEM_NAME}/include"
        "${HEASOFT_DIR}/tcltk/${HEA_SYSTEM_NAME}/include"
    )
    target_link_libraries(_pymcvspec PRIVATE "${HEASOFT_DIR}/Xspec/${HEA_SYSTEM_NAME}/lib/libXSFunctions.dylib")
    set(PYMCVSPEC_DIR "${CMAKE_INSTALL_PREFIX}/pymcvspec")
    install(TARGETS _pymcvspec DESTINATION ${PYMCVSPEC_DIR})
    target_compile_definitions(_pymcvspec PRIVATE)
endif()

# compile xspec model
set(XSPEC_MODEL_DIR "${CMAKE_INSTALL_PREFIX}/xspec")
file(MAKE_DIRECTORY "${XSPEC_MODEL_DIR}")

set(SOURCE_FILES ${CMAKE_SOURCE_DIR}/source/Cataclysmic_Variable.cpp;
                 ${CMAKE_SOURCE_DIR}/source/integration.cpp)
set(MODEL_FILES ${CMAKE_SOURCE_DIR}/model_funcs/polarspec.cpp;
                ${CMAKE_SOURCE_DIR}/model_funcs/ipsepc.cpp;
                ${CMAKE_SOURCE_DIR}/model_funcs/deqspec.cpp;
                ${CMAKE_SOURCE_DIR}/model_funcs/mcvspec.dat)
set(HEADER_FILES ${CMAKE_SOURCE_DIR}/include/Cataclysmic_Variable.hh;
                 ${CMAKE_SOURCE_DIR}/include/integration.hh;
                 ${CMAKE_SOURCE_DIR}/include/constants.hh;
                 ${CMAKE_SOURCE_DIR}/include/tableau.hh;
                 ${MR_HEADER_PATH})

foreach(FILE IN LISTS SOURCE_FILES MODEL_FILES HEADER_FILES)
    get_filename_component(FNAME "${FILE}" NAME)
    file(CREATE_LINK "${FILE}" "${XSPEC_MODEL_DIR}/${FNAME}" SYMBOLIC)
endforeach()

# clean up xspec dir if still has old files
if(EXISTS "${XSPEC_MODEL_DIR}/lpack_mcvspec.cxx")
    file(REMOVE "${XSPEC_MODEL_DIR}/lpack_mcvspec.cxx")
endif()

execute_process (COMMAND bash -c "initpackage mcvspec ${XSPEC_MODEL_DIR}/mcvspec.dat ${XSPEC_MODEL_DIR}")

add_custom_command(
    OUTPUT "${XSPEC_MODEL_DIR}/hmake.log"
    COMMAND bash -c "hmake > hmake.log"
    WORKING_DIRECTORY ${XSPEC_MODEL_DIR}
    COMMENT "Building XSPEC Model"
)
add_custom_target(mcvspec_xspec ALL DEPENDS "${XSPEC_MODEL_DIR}/hmake.log" VERBATIM)
