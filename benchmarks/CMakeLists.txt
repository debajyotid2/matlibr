# Configure benchmarking.
set(BENCHMARK_BINS matlibr.benchmark)
file(GLOB_RECURSE BENCHMARK_SOURCE_FILES ${PROJECT_ROOT_DIR}/benchmarks/*.c)

include(AddCatch2)

message(STATUS "Building benchmarking binaries ...")
add_executable(${BENCHMARK_BINS} ${BENCHMARK_SOURCE_FILES})

target_include_directories(${BENCHMARK_BINS} PUBLIC ${PROJECT_ROOT_DIR}/src)
target_include_directories(${BENCHMARK_BINS} PUBLIC ${PROJECT_ROOT_DIR}/third_party/Catch2/include)
target_link_libraries(${BENCHMARK_BINS} PRIVATE 
    ${PROJECT_LIBRARY_NAME}
    Catch2::Catch2WithMain
)
