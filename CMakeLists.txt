cmake_minimum_required(VERSION 3.18)

project(xyz2smiles VERSION 0.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Eigen3 REQUIRED NO_MODULE)
# find_package(pybind11 REQUIRED)
find_package(RDKit REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)

set(SOURCE_DIR "src/${PROJECT_NAME}")
set(TARGET_H 
    "${SOURCE_DIR}/AdjMatrix.h"
    "${SOURCE_DIR}/Fileio.h"
    "${SOURCE_DIR}/Mol.h"
    "${SOURCE_DIR}/VecUtils.h"
    "${SOURCE_DIR}/xyz2smiles.h"
)

set(TARGET_SRC 
    "${SOURCE_DIR}/AdjMatrix.cpp"
    "${SOURCE_DIR}/Fileio.cpp"
    "${SOURCE_DIR}/Mol.cpp"
    "${SOURCE_DIR}/VecUtils.cpp"
    "${SOURCE_DIR}/xyz2smiles.cpp"
)

install(DIRECTORY DESTINATION "${CMAKE_SOURCE_DIR}/bin")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin")

add_executable("${PROJECT_NAME}" "${TARGET_SRC}"  "${TARGET_H}" "${SOURCE_DIR}/main.cpp")
target_link_libraries(
    "${PROJECT_NAME}"
    Eigen3::Eigen 
    RDKit::SmilesParse RDKit::GraphMol
    Boost::program_options
)

#pybind11_add_module(xyz2smiles ${SOURCES})
#target_link_libraries (xyz2smiles PRIVATE Eigen3::Eigen RDKit::SmilesParse RDKit::GraphMol Boost::program_options)

