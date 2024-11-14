cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(sarcomere)


add_library(my_library INTERFACE)
target_include_directories(my_library INTERFACE include)


find_package(GSL REQUIRED)
# Find the HDF5 package
find_package(HDF5 REQUIRED COMPONENTS C CXX)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})



include_directories("${PROJECT_SOURCE_DIR}/autodiff")

# Set compiler flags for debugging symbols
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "-g")  # Add debug symbols




# add_executable(${PROJECT_NAME} main.cpp)
# target_link_libraries(${PROJECT_NAME} PRIVATE my_library)
# add_executable(minimal minimal.cpp)
# target_link_libraries(minimal PRIVATE my_library)
add_executable(${PROJECT_NAME} langevin.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE my_library)
add_executable(tests tests.cpp)
target_link_libraries(tests PRIVATE my_library gtest_main)

add_executable(autodiff autodiff.cpp)
target_link_libraries(autodiff PRIVATE my_library)



# Include GSL headers
target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS})
target_include_directories(my_library INTERFACE ${GSL_INCLUDE_DIRS})

# Link GSL libraries
target_link_libraries(${PROJECT_NAME} PRIVATE ${GSL_LIBRARIES})
target_link_libraries(my_library INTERFACE ${GSL_LIBRARIES})

# Include HDF5 headers
target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIRS})
target_include_directories(my_library INTERFACE ${HDF5_INCLUDE_DIRS})

# Link HDF5 libraries
target_link_libraries(${PROJECT_NAME} PRIVATE ${HDF5_LIBRARIES})
target_link_libraries(my_library INTERFACE ${HDF5_LIBRARIES})



# Add Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)
gtest_discover_tests(tests)

# Define custom target to run Valgrind
add_custom_target(valgrind
    COMMAND valgrind --leak-check=full ./${PROJECT_NAME}
    DEPENDS ${PROJECT_NAME}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running Valgrind on my_executable..."
)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)


