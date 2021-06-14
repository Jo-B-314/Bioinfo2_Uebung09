cmake_minimum_required(VERSION 2.8.12)

project(Uebung01)

# Add C++ compile flags
set(CXX_FLAGS
	-std=c++14               # Use C++14
    -pedantic                
    -Wall                    
)

# Set CIP Pool specific include paths
set(BALL_DIR "/opt/BALL/install/lib/cmake/BALL/")
set(INSTALLER_PATH /installer/import/debian/c++)
set(CIP_STD_LIBS c++ c++abi m c gcc_s gcc)
set(GTEST_SRC_DIR ${INSTALLER_PATH}/googletest-1.7.0)

link_directories(${INSTALLER_PATH}/libcxx/lib)

# Find packages
find_package(BALL REQUIRED)
find_package(Qt5Core REQUIRED)
find_package(Qt5Xml REQUIRED)
find_package(Qt5Network REQUIRED)

include_directories(
	${BALL_INCLUDE_DIRS}
	"${CMAKE_SOURCE_DIR}"
	"${CMAKE_BINARY_DIR}"
)
add_definitions(${CXX_FLAGS} ${BALL_ADDCXX_FLAGS})

# Enable testing
enable_testing()

# Find Boost
find_package(Boost 1.55.0 COMPONENTS iostreams filesystem system program_options REQUIRED)

if(NOT Boost_FOUND)
	message(FATAL_ERROR "Boost version 1.55.0 or higher is required")
endif()

SET(Boost_LIBRARIES
	${Boost_IOSTREAMS_LIBRARY}
	${Boost_FILESYSTEM_LIBRARY}
	${Boost_SYSTEM_LIBRARY}
	${Boost_PROGRAM_OPTIONS_LIBRARY}
)

include_directories(${Boost_INCLUDE_DIRS})

# Add a new library, which is used for the executable below and
# for the tests in test/CMakeLists.txt
# Here, you can also add new files to your project
add_library(Uebung01 SHARED
        Test.hpp
        Test.cpp
		getSecStruct.h
		Matrix.h
)

target_link_libraries(Uebung01 BALL ${Boost_LIBRARIES})


# Add new program example1
# Here, you can also add new programs
add_executable(example1 example.C)
target_link_libraries(example1 Uebung01 BALL)
add_executable(getSecStruct getSecStruct.C)
target_link_libraries(getSecStruct Uebung01 BALL)

add_subdirectory(tests)
