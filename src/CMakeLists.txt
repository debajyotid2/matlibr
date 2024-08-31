
file(GLOB_RECURSE SOURCE_FILES ${PROJECT_ROOT_DIR}/src/*.c)
add_library(${PROJECT_LIBRARY_NAME} STATIC ${SOURCE_FILES})

find_package(openblas REQUIRED)

# target_include_directories(${PROJECT_LIBRARY_NAME} PUBLIC ${PROJECT_ROOT_DIR}/src)
target_include_directories(${PROJECT_LIBRARY_NAME} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_ROOT_DIR}/src>
    $<INSTALL_INTERFACE:include>
)
target_include_directories(${PROJECT_LIBRARY_NAME} PRIVATE ${PROJECT_ROOT_DIR}/third_party/OpenBLAS-0.3.28)
target_link_libraries(${PROJECT_LIBRARY_NAME} PRIVATE m dl ${OPENBLAS_LIBRARIES})
