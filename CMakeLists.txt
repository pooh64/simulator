cmake_minimum_required(VERSION 3.10)

project(sim C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

find_package(Qt5 COMPONENTS Widgets Core Gui REQUIRED)

add_executable(sim
    cpu.cpp
    display.cpp
    instr.cpp
    life.cpp
    locus.cpp
    main.cpp
    mainwindow.cpp
    mapgenerator.cpp
    relaxation.cpp
    relaxationtor.cpp
)

set_target_properties(sim PROPERTIES AUTOMOC TRUE)
set_target_properties(sim PROPERTIES AUTOUIC TRUE)
target_link_libraries(sim Qt5::Widgets Qt5::Core Qt5::Gui)
target_include_directories(sim PUBLIC "${PROJECT_SOURCE_DIR}")
