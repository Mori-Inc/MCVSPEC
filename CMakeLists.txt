cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(PYBIND11_FINDPYTHON ON)

find_package(Python 3.13 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

project(mcvspec VERSION 1.0.0)
#add_library(mcvspec STATIC source/Cataclysmic_Variable.cpp
#                           include/Cataclysmic_Variable.hh)
#
#pybind11_add_module(mcvspec source/Cataclysmic_Variable.cpp source/integration.cpp)
pybind11_add_module(mcvspec pymcvspec/bindings.cpp)

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

# MCVSPEC Headers
target_include_directories(mcvspec PUBLIC "${PROJECT_SOURCE_DIR}/include")

# HEASOFT Headers
target_include_directories(mcvspec PUBLIC
    "${HEASOFT_DIR}/Xspec/${HEA_SYSTEM_NAME}/include"
    "${HEASOFT_DIR}/Xspec/${HEA_SYSTEM_NAME}/include/XSFunctions"
    "${HEASOFT_DIR}/Xspec/${HEA_SYSTEM_NAME}/include/XSFunctions/Utilities"
    "${HEASOFT_DIR}/heacore/${HEA_SYSTEM_NAME}/include"
    "${HEASOFT_DIR}/tcltk/${HEA_SYSTEM_NAME}/include"
)
target_link_libraries(mcvspec PRIVATE "${HEASOFT_DIR}/Xspec/${HEA_SYSTEM_NAME}/lib/libXSFunctions.dylib")

target_compile_definitions(mcvspec PRIVATE)
