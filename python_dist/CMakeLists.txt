cmake_minimum_required(VERSION 3.15)
project(deContamination LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(OpenMP REQUIRED)

pybind11_add_module(deContamination
    src/python_bindings.cpp
    src/Contamination.cpp
)

target_link_libraries(deContamination PRIVATE OpenMP::OpenMP_CXX)
target_include_directories(deContamination PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)


install(TARGETS deContamination LIBRARY DESTINATION .)