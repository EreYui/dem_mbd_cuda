cmake_minimum_required(VERSION 3.24)

# This must be set before enabling CUDA; otherwise CMake seeds its compiler
# default architecture first and the intended native default is never used.
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native CACHE STRING
        "CUDA architectures (defaults to the locally installed GPU)")
endif()

project(dem_mbd LANGUAGES CXX CUDA)

set(PROJECT_TARGET dem_mbd)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_EXTENSIONS OFF)

set(DEM_MBD_SOURCES
    src/async_output.cpp
    src/contact_geometry.cpp
    src/cuda_backend_cuda.cu
    src/cuda_simulation_integrator.cpp
    src/dyneq.cpp
    src/force.cpp
    src/grid.cpp
    src/main.cpp
    src/matrixtrans.cpp
    src/multibody.cpp
    src/particle.cpp
    src/simulation.cpp
    src/wall.cpp
    src/cuda/cuda_body_backend.cu
    src/cuda/cuda_contact_backend.cu
    src/cuda/cuda_grid_backend.cu
    src/cuda/cuda_particle_backend.cu
)

# Host lifecycle glue uses legacy host-only types and calls the CUDA C API.
set_source_files_properties(src/cuda_backend_cuda.cu PROPERTIES LANGUAGE CXX)

add_executable(${PROJECT_TARGET} ${DEM_MBD_SOURCES})
target_include_directories(${PROJECT_TARGET} PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/src"
)

find_package(CUDAToolkit REQUIRED)
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_TARGET} PRIVATE CUDA::cudart Threads::Threads)
target_compile_definitions(${PROJECT_TARGET} PRIVATE DEM_MBD_USE_CUDA=1)
set_target_properties(${PROJECT_TARGET} PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
    CUDA_RUNTIME_LIBRARY None
)

if(MSVC)
    target_compile_options(${PROJECT_TARGET} PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:/O2>
        $<$<COMPILE_LANGUAGE:CXX>:/EHsc>
        $<$<COMPILE_LANGUAGE:CXX>:/utf-8>
        $<$<COMPILE_LANGUAGE:CXX>:/Zc:preprocessor>
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/O2>
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/utf-8>
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor>
    )
else()
    target_compile_options(${PROJECT_TARGET} PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3>
        $<$<COMPILE_LANGUAGE:CUDA>:-O3>
    )
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(contact_geometry_tests
        tests/contact_geometry_test.cpp
        src/contact_geometry.cpp
    )
    target_include_directories(contact_geometry_tests PRIVATE include)
    if(MSVC)
        target_compile_options(contact_geometry_tests PRIVATE /utf-8)
    endif()
    add_test(NAME contact_geometry_tests COMMAND contact_geometry_tests)

    add_executable(file_utils_tests tests/file_utils_test.cpp)
    target_include_directories(file_utils_tests PRIVATE include)
    if(MSVC)
        target_compile_options(file_utils_tests PRIVATE /utf-8)
    endif()
    add_test(NAME file_utils_tests COMMAND file_utils_tests)
    set_tests_properties(file_utils_tests PROPERTIES
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    )

    add_executable(particle_schema_tests
        tests/particle_schema_test.cpp
        src/particle.cpp
    )
    target_include_directories(particle_schema_tests PRIVATE include)
    if(MSVC)
        target_compile_options(particle_schema_tests PRIVATE /utf-8)
    endif()
    add_test(NAME particle_schema_tests COMMAND particle_schema_tests)

    add_executable(output_compare tests/output_compare.cpp)
    target_compile_features(output_compare PRIVATE cxx_std_17)
    if(MSVC)
        target_compile_options(output_compare PRIVATE /utf-8)
    endif()

    add_executable(state_diagnostics tests/state_diagnostics.cpp)
    target_compile_features(state_diagnostics PRIVATE cxx_std_17)
    if(MSVC)
        target_compile_options(state_diagnostics PRIVATE /utf-8)
    endif()

    add_executable(scenario_validator tests/scenario_validator.cpp)
    target_compile_features(scenario_validator PRIVATE cxx_std_17)
    if(MSVC)
        target_compile_options(scenario_validator PRIVATE /utf-8)
    endif()

    file(COPY "${CMAKE_SOURCE_DIR}/validation/cases/"
         DESTINATION "${CMAKE_BINARY_DIR}/validation_cases")
    foreach(scenario free_fall wall_contact two_particle triangle_contact body_states)
        add_test(NAME cuda_${scenario}
            COMMAND ${PROJECT_TARGET}
                "${CMAKE_BINARY_DIR}/validation_cases/${scenario}")
        set_tests_properties(cuda_${scenario} PROPERTIES RESOURCE_LOCK cuda_gpu)
        add_test(NAME validate_${scenario}
            COMMAND scenario_validator ${scenario}
                "${CMAKE_BINARY_DIR}/validation_cases/${scenario}/OutputFile")
        set_tests_properties(validate_${scenario} PROPERTIES
            DEPENDS cuda_${scenario})
    endforeach()

    add_test(NAME cuda_device_test COMMAND ${PROJECT_TARGET} --cuda-info)
    add_test(NAME cuda_pure_dem_short
        COMMAND ${PROJECT_TARGET} "${CMAKE_SOURCE_DIR}/validation/ctest_data1_pure_dem")
    set_tests_properties(cuda_pure_dem_short PROPERTIES
        PASS_REGULAR_EXPRESSION "CUDA DEM solver"
        RESOURCE_LOCK cuda_gpu)
    add_test(NAME cuda_dem_mbd_short
        COMMAND ${PROJECT_TARGET} "${CMAKE_SOURCE_DIR}/validation/ctest_data1_mbd")
    set_tests_properties(cuda_dem_mbd_short PROPERTIES
        PASS_REGULAR_EXPRESSION "CUDA DEM-MBD solver"
        RESOURCE_LOCK cuda_gpu)
endif()
