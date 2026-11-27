file(GLOB_RECURSE SOURCE_FILES ${PROJECT_ROOT_DIR}/src/*.c ${PROJECT_ROOT_DIR}/src/*.h)
add_library(${PROJECT_LIBRARY_NAME} STATIC ${SOURCE_FILES})

find_package(openblas REQUIRED)

target_include_directories(${PROJECT_LIBRARY_NAME} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_ROOT_DIR}/src>
    $<BUILD_INTERFACE:${PROJECT_ROOT_DIR}/third_party/OpenBLAS-0.3.28>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(${PROJECT_LIBRARY_NAME} PUBLIC m dl ${OPENBLAS_LIBRARIES})

target_compile_options(${PROJECT_LIBRARY_NAME} PRIVATE
    $<$<C_COMPILER_ID:GNU,CLang>:-O3 -march=native -Wall 
                -Wextra -pedantic -Wno-missing-field-initializers -Wno-unused-parameter>
    $<$<C_COMPILER_ID:MSVC>:/W4 /permissive->
)
