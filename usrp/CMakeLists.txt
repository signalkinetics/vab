cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
set(TARGET_NAME sk-uhd)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(${TARGET_NAME})

find_package(UHD)
find_package(Boost COMPONENTS program_options filesystem REQUIRED)

include_directories(${UHD_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

# file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
foreach(srcfile ${SOURCES} )
    get_filename_component(srcname ${srcfile} NAME_WE)
    message("Adding ${srcname} executable from ${srcfile} file...")
    add_executable(${srcname} ${srcfile} )

    target_link_libraries(${srcname} ${Boost_LIBRARIES})
    target_link_libraries(${srcname} ${UHD_LIBRARIES})
endforeach(srcfile ${SOURCES})
