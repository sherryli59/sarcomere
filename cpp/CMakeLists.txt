cmake_minimum_required(VERSION 3.5)

# Force CMake to use the compiler from the Conda environment
if(DEFINED ENV{CONDA_PREFIX})
    set(CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX} CACHE PATH "Path to Conda environment")
    set(CMAKE_C_COMPILER $ENV{CONDA_PREFIX}/bin/gcc CACHE FILEPATH "C compiler")
    set(CMAKE_CXX_COMPILER $ENV{CONDA_PREFIX}/bin/g++ CACHE FILEPATH "C++ compiler")
    set(CMAKE_LIBRARY_PATH $ENV{CONDA_PREFIX}/lib CACHE PATH "Library path")
    set(CMAKE_INCLUDE_PATH $ENV{CONDA_PREFIX}/include CACHE PATH "Include path")
endif()

project(sarcomere)

# Set source and include directories
set(SRC_DIR "${PROJECT_SOURCE_DIR}/src")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

# Include the header files
include_directories(${INCLUDE_DIR})

# Find required libraries
find_package(GSL REQUIRED)
#find_package(GTest REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(OpenMP REQUIRED)
include_directories("${PROJECT_SOURCE_DIR}/autodiff")

# Default production builds to optimized code while retaining symbols for profiling.
# Callers can still explicitly select Debug with -DCMAKE_BUILD_TYPE=Debug.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Source files
set(SOURCES
    ${SRC_DIR}/sarcomere.cpp
    ${SRC_DIR}/geometry.cpp
    ${SRC_DIR}/interaction.cpp
    ${SRC_DIR}/utils.cpp
    ${SRC_DIR}/neighborlist.cpp
    ${SRC_DIR}/h5_utils.cpp
    ${SRC_DIR}/components.cpp
    ${SRC_DIR}/langevin.cpp
    ${SRC_DIR}/initialization.cpp


)



# FetchContent is built into CMake from version 3.11 onward.
include(FetchContent)

# Option to build Google Benchmark; disabled by default to avoid extra
# external dependencies (curl/OpenSSL) on clusters with mixed libs.
option(BUILD_BENCHMARK "Build Google Benchmark targets" OFF)
if(BUILD_BENCHMARK)
    # Optionally disable Google Benchmark's internal tests.
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable benchmark tests" FORCE)

    # Declare and download the Google Benchmark repository.
    FetchContent_Declare(
        benchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG v1.7.0  # Specify the desired tag/commit.
    )

    # Make the benchmark content available (it will be built as part of your project).
    FetchContent_MakeAvailable(benchmark)
endif()


# Create the library (STATIC or SHARED)
add_library(my_library STATIC ${SOURCES})  # Use SHARED for a dynamic library
target_include_directories(my_library PUBLIC ${INCLUDE_DIR})

# Link required libraries
target_link_libraries(my_library PUBLIC ${GSL_LIBRARIES} ${HDF5_LIBRARIES} Eigen3::Eigen OpenMP::OpenMP_CXX)

# Executables
add_executable(${PROJECT_NAME} ${SRC_DIR}/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE my_library)

add_executable(rod_mobility_calculator ${SRC_DIR}/rod_mobility_calculator.cpp)
target_include_directories(rod_mobility_calculator PRIVATE ${INCLUDE_DIR})

# Include directories for required libraries
target_include_directories(my_library PUBLIC ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

#check if the benchmark.cpp file exists
if(BUILD_BENCHMARK AND EXISTS ${SRC_DIR}/benchmark.cpp)
    # Now create your target and link it against Google Benchmark.
    add_executable(my_benchmark ${SRC_DIR}/benchmark.cpp)
    target_link_libraries(my_benchmark PRIVATE my_library benchmark::benchmark)
endif()

# Link GSL, HDF5, and OpenMP
target_link_libraries(${PROJECT_NAME} PRIVATE ${GSL_LIBRARIES} ${HDF5_LIBRARIES} OpenMP::OpenMP_CXX)

include(CTest)
if(BUILD_TESTING)
    add_executable(catch_bond_kinetics_test tests/catch_bond_kinetics_test.cpp)
    target_include_directories(catch_bond_kinetics_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME catch_bond_kinetics COMMAND catch_bond_kinetics_test)

    add_executable(am_motor_kinetics_test tests/am_motor_kinetics_test.cpp)
    target_include_directories(am_motor_kinetics_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME am_motor_kinetics COMMAND am_motor_kinetics_test)

    add_executable(actn2_load_kinetics_test tests/actn2_load_kinetics_test.cpp)
    target_include_directories(actn2_load_kinetics_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME actn2_load_kinetics COMMAND actn2_load_kinetics_test)

    add_executable(actn2_arm_kinetics_test tests/actn2_arm_kinetics_test.cpp)
    target_include_directories(actn2_arm_kinetics_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME actn2_arm_kinetics COMMAND actn2_arm_kinetics_test)

    add_executable(zdisc_tip_gap_test tests/zdisc_tip_gap_test.cpp)
    target_link_libraries(zdisc_tip_gap_test PRIVATE my_library)
    target_include_directories(zdisc_tip_gap_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME zdisc_tip_gap COMMAND zdisc_tip_gap_test)

    add_executable(pbc_round_test tests/pbc_round_test.cpp)
    target_link_libraries(pbc_round_test PRIVATE my_library)
    target_include_directories(pbc_round_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME pbc_round COMMAND pbc_round_test)

    add_executable(am_prefilter_test tests/am_prefilter_test.cpp)
    target_link_libraries(am_prefilter_test PRIVATE my_library)
    target_include_directories(am_prefilter_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME am_prefilter COMMAND am_prefilter_test)

    add_executable(barbed_lateral_test tests/barbed_lateral_test.cpp)
    target_link_libraries(barbed_lateral_test PRIVATE my_library)
    target_include_directories(barbed_lateral_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME barbed_lateral COMMAND barbed_lateral_test)

    add_executable(actn2_myosin_contact_test
                   tests/actn2_myosin_contact_test.cpp)
    target_link_libraries(actn2_myosin_contact_test PRIVATE my_library)
    add_test(NAME actn2_myosin_contact COMMAND actn2_myosin_contact_test)

    add_executable(zdisc_myosin_tip_contact_test
                   tests/zdisc_myosin_tip_contact_test.cpp)
    target_link_libraries(zdisc_myosin_tip_contact_test PRIVATE my_library)
    add_test(NAME zdisc_myosin_tip_contact COMMAND zdisc_myosin_tip_contact_test)

    add_executable(actn2_arm_state_integration_test
                   tests/actn2_arm_state_integration_test.cpp)
    target_link_libraries(actn2_arm_state_integration_test PRIVATE my_library)
    add_test(NAME actn2_arm_state_integration
             COMMAND actn2_arm_state_integration_test)

    # Manual single-link timing harness; intentionally not registered in CTest.
    add_executable(actn2_arm_benchmark tests/actn2_arm_benchmark.cpp)
    target_link_libraries(actn2_arm_benchmark PRIVATE my_library)

    # Dense ACTN2-tip solver A/B parity test and scalable timing harness.
    add_executable(actn2_tip_projection_benchmark
                   tests/actn2_tip_projection_benchmark.cpp)
    target_link_libraries(actn2_tip_projection_benchmark PRIVATE my_library)
    target_compile_options(actn2_tip_projection_benchmark PRIVATE -UNDEBUG)
    add_test(NAME actn2_tip_projection_parity
             COMMAND actn2_tip_projection_benchmark 200 2)
    set_tests_properties(actn2_tip_projection_parity PROPERTIES
                         ENVIRONMENT "OMP_NUM_THREADS=1")

    add_executable(explicit_crossbridge_mechanics_test
                   tests/explicit_crossbridge_mechanics_test.cpp)
    target_link_libraries(explicit_crossbridge_mechanics_test PRIVATE my_library)
    add_test(NAME explicit_crossbridge_mechanics
             COMMAND explicit_crossbridge_mechanics_test)

    add_executable(am_head_ownership_test tests/am_head_ownership_test.cpp)
    target_link_libraries(am_head_ownership_test PRIVATE my_library)
    add_test(NAME am_head_ownership COMMAND am_head_ownership_test)

    add_executable(lattice_spacing_test tests/lattice_spacing_test.cpp)
    target_include_directories(lattice_spacing_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME lattice_spacing COMMAND lattice_spacing_test)

    add_executable(rod_mobility_test tests/rod_mobility_test.cpp)
    target_include_directories(rod_mobility_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME rod_mobility COMMAND rod_mobility_test)

    add_executable(rod_mobility_statistics_test tests/rod_mobility_statistics_test.cpp)
    target_include_directories(rod_mobility_statistics_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME rod_mobility_statistics COMMAND rod_mobility_statistics_test)

    add_executable(aa_refractory_test tests/aa_refractory_test.cpp)
    target_include_directories(aa_refractory_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME aa_refractory COMMAND aa_refractory_test)

    add_executable(aa_association_kinetics_test tests/aa_association_kinetics_test.cpp)
    target_include_directories(aa_association_kinetics_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME aa_association_kinetics COMMAND aa_association_kinetics_test)

    add_executable(aa_neighbor_identity_test tests/aa_neighbor_identity_test.cpp)
    target_include_directories(aa_neighbor_identity_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME aa_neighbor_identity COMMAND aa_neighbor_identity_test)

    add_executable(counter_rng_test tests/counter_rng_test.cpp)
    target_include_directories(counter_rng_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME counter_rng COMMAND counter_rng_test)

    add_executable(titin_kinetics_test tests/titin_kinetics_test.cpp)
    target_include_directories(titin_kinetics_test PRIVATE ${INCLUDE_DIR})
    add_test(NAME titin_kinetics COMMAND titin_kinetics_test)

    add_executable(explicit_titin_test tests/explicit_titin_test.cpp)
    target_link_libraries(explicit_titin_test PRIVATE my_library)
    add_test(NAME explicit_titin COMMAND explicit_titin_test)

    add_executable(explicit_titin_balanced_restoring_test
                   tests/explicit_titin_balanced_restoring_test.cpp)
    target_link_libraries(explicit_titin_balanced_restoring_test
                          PRIVATE my_library)
    add_test(NAME explicit_titin_balanced_restoring
             COMMAND explicit_titin_balanced_restoring_test)

    add_executable(titin_actn2_capture_gate_test
                   tests/titin_actn2_capture_gate_test.cpp)
    target_link_libraries(titin_actn2_capture_gate_test PRIVATE my_library)
    add_test(NAME titin_actn2_capture_gate
             COMMAND titin_actn2_capture_gate_test)

    add_executable(motor_validation src/motor_validation.cpp)
    target_link_libraries(motor_validation PRIVATE my_library)
    target_compile_options(motor_validation PRIVATE -UNDEBUG)

    add_executable(myomesin_mband_test tests/myomesin_mband_test.cpp)
    target_link_libraries(myomesin_mband_test PRIVATE my_library)
    add_test(NAME myomesin_mband COMMAND myomesin_mband_test)

    add_executable(aa_pair_symmetry_test tests/aa_pair_symmetry_test.cpp)
    target_link_libraries(aa_pair_symmetry_test PRIVATE my_library)
    add_test(NAME aa_pair_symmetry COMMAND aa_pair_symmetry_test)

    add_executable(lattice_mechanics_test tests/lattice_mechanics_test.cpp)
    target_link_libraries(lattice_mechanics_test PRIVATE my_library)
    add_test(NAME lattice_mechanics COMMAND lattice_mechanics_test)

    add_executable(local_aa_link_test tests/local_aa_link_test.cpp)
    target_link_libraries(local_aa_link_test PRIVATE my_library)
    add_test(NAME local_aa_link COMMAND local_aa_link_test)

    add_executable(actin_sterics_test tests/actin_sterics_test.cpp)
    target_link_libraries(actin_sterics_test PRIVATE my_library)
    add_test(NAME actin_sterics COMMAND actin_sterics_test)

    add_executable(neighborlist_two_cell_test tests/neighborlist_two_cell_test.cpp)
    target_link_libraries(neighborlist_two_cell_test PRIVATE my_library)
    add_test(NAME neighborlist_two_cell COMMAND neighborlist_two_cell_test)

    add_executable(am_geometry_regression_test
                   tests/am_geometry_regression_test.cpp)
    target_link_libraries(am_geometry_regression_test PRIVATE my_library)
    add_test(NAME am_geometry_regression COMMAND am_geometry_regression_test)

    add_executable(am_mobility_direction_test
                   tests/am_mobility_direction_test.cpp)
    target_link_libraries(am_mobility_direction_test PRIVATE my_library)
    add_test(NAME am_mobility_direction COMMAND am_mobility_direction_test)

    add_executable(actn2_motor_removal_test
                   tests/actn2_motor_removal_test.cpp)
    target_link_libraries(actn2_motor_removal_test PRIVATE my_library)
    add_test(NAME actn2_motor_removal COMMAND actn2_motor_removal_test)

    add_executable(actn2_am_occlusion_test
                   tests/actn2_am_occlusion_test.cpp)
    target_link_libraries(actn2_am_occlusion_test PRIVATE my_library)
    add_test(NAME actn2_am_occlusion COMMAND actn2_am_occlusion_test)

    add_executable(minimal_balanced_initialization_test
                   tests/minimal_balanced_initialization_test.cpp)
    target_link_libraries(minimal_balanced_initialization_test PRIVATE my_library)
    add_test(NAME minimal_balanced_initialization
             COMMAND minimal_balanced_initialization_test)

    add_executable(minimal_double_myosin_balanced_initialization_test
                   tests/minimal_double_myosin_balanced_initialization_test.cpp)
    target_link_libraries(minimal_double_myosin_balanced_initialization_test
                          PRIVATE my_library)
    add_test(NAME minimal_double_myosin_balanced_initialization
             COMMAND minimal_double_myosin_balanced_initialization_test)

    add_executable(shared_lattice_initialization_test
                   tests/shared_lattice_initialization_test.cpp)
    target_link_libraries(shared_lattice_initialization_test PRIVATE my_library)
    add_test(NAME shared_lattice_initialization
             COMMAND shared_lattice_initialization_test)

    set_tests_properties(
        am_geometry_regression am_mobility_direction actn2_motor_removal
        actn2_am_occlusion
        minimal_balanced_initialization
        minimal_double_myosin_balanced_initialization
        shared_lattice_initialization
        titin_actn2_capture_gate explicit_titin_balanced_restoring
        myomesin_mband
        PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;HDF5_USE_FILE_LOCKING=FALSE")

    # The regression executables use C/C++ assert(). Production builds define
    # NDEBUG, so explicitly retain assertions for every test target even when
    # the library and simulator are compiled in Release/RelWithDebInfo mode.
    set(SARCOMERE_ASSERT_TEST_TARGETS
        catch_bond_kinetics_test
        am_motor_kinetics_test
        actn2_load_kinetics_test
        actn2_arm_kinetics_test
        actn2_myosin_contact_test
        am_prefilter_test
        pbc_round_test
        zdisc_tip_gap_test
        barbed_lateral_test
        zdisc_myosin_tip_contact_test
        actn2_arm_state_integration_test
        explicit_crossbridge_mechanics_test
        am_head_ownership_test
        lattice_spacing_test
        rod_mobility_test
        rod_mobility_statistics_test
        aa_refractory_test
        aa_association_kinetics_test
        aa_neighbor_identity_test
        counter_rng_test
        titin_kinetics_test
        explicit_titin_test
        explicit_titin_balanced_restoring_test
        titin_actn2_capture_gate_test
        myomesin_mband_test
        aa_pair_symmetry_test
        lattice_mechanics_test
        local_aa_link_test
        actin_sterics_test
        neighborlist_two_cell_test
        am_geometry_regression_test
        am_mobility_direction_test
        actn2_motor_removal_test
        actn2_am_occlusion_test
        minimal_balanced_initialization_test
        minimal_double_myosin_balanced_initialization_test
        shared_lattice_initialization_test)
    foreach(test_target IN LISTS SARCOMERE_ASSERT_TEST_TARGETS)
        target_compile_options(${test_target} PRIVATE -UNDEBUG)
    endforeach()
endif()

# Google Test (Uncomment if needed)
# include_directories(${GTEST_INCLUDE_DIRS})
# add_executable(tests ${SRC_DIR}/tests.cpp)
# target_link_libraries(tests PRIVATE my_library gtest_main)
# include(GoogleTest)
# gtest_discover_tests(tests)

# Define custom target to run Valgrind
add_custom_target(valgrind
    COMMAND valgrind --leak-check=full ./${PROJECT_NAME}
    DEPENDS ${PROJECT_NAME}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running Valgrind on ${PROJECT_NAME}..."
)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
