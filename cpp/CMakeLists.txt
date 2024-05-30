cmake_minimum_required(VERSION 3.0.0)
project(sarcomere)


add_library(my_library INTERFACE)
target_include_directories(my_library INTERFACE include)


find_package(GSL REQUIRED)
# Find the HDF5 package
find_package(HDF5 REQUIRED COMPONENTS C CXX)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})


# Set compiler flags for debugging symbols
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")




add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE my_library)
add_executable(tests tests.cpp)
target_link_libraries(tests PRIVATE my_library gtest_main)


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
