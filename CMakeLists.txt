cmake_minimum_required(VERSION 3.27.2)

project(matlibr)

set(CMAKE_C_COMPILER g++)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_PREFIX_PATH ${CMAKE_SOURCE_DIR}/third_party)
set(PROJECT_ROOT_DIR ${CMAKE_SOURCE_DIR})
set(PROJECT_LIBRARY_NAME matlibr)
set(PROJECT_VERSION 0.0.1)

add_subdirectory(src)

enable_testing()

add_subdirectory(tests)

# Set the version of the library
set_target_properties(${PROJECT_LIBRARY_NAME} PROPERTIES VERSION ${PROJECT_VERSION})

# Install the library and headers
install(TARGETS ${PROJECT_LIBRARY_NAME}
    EXPORT ${PROJECT_LIBRARY_NAME}Targets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)
install(FILES ${PROJECT_ROOT_DIR}/src/matrix.h DESTINATION include)

# Create a package configuration file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_LIBRARY_NAME}ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(EXPORT ${PROJECT_LIBRARY_NAME}Targets
    FILE ${PROJECT_LIBRARY_NAME}Targets.cmake
    NAMESPACE ${PROJECT_LIBRARY_NAME}::
    DESTINATION lib/cmake/${PROJECT_LIBRARY_NAME}
)

install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_LIBRARY_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_LIBRARY_NAME}ConfigVersion.cmake"
    DESTINATION lib/cmake/${PROJECT_LIBRARY_NAME}
)
