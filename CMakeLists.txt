# Minimum CMake version.
cmake_minimum_required (VERSION 2.8.5)

# Project and version numbers.
project (polyamri)
set (POLYAMRI_MAJOR_VERSION 0)
set (POLYAMRI_MINOR_VERSION 7)
set (POLYAMRI_PATCH_VERSION 0)
set (POLYAMRI_VERSION "${POLYAMRI_MAJOR_VERSION}.${POLYAMRI_MINOR_VERSION}.${POLYAMRI_PATCH_VERSION}")
message("-- Building polyamri (v${POLYAMRI_VERSION})")

# Use our own CMake stuff.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Look for headers in the top directory first.
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

# Look for Polymec and use its CMake settings.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${POLYMEC_PREFIX}/share/polymec/")
include(polymec)
if (NOT DEFINED POLYMEC_VERSION)
  message(FATAL_ERROR "Could not find polymec library. Please use the polymec=/path/to/polymec option.")
endif()
message("-- Found polymec library in ${POLYMEC_PREFIX} (v${POLYMEC_VERSION})")

# Build everything the way Polymec does.
set(BUILD_SHARED_LIBS ${POLYMEC_BUILD_SHARED_LIBS})
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${POLYMEC_PREFIX}/lib")

# Generate polyamri_version.h with the right version numbers.
add_custom_target(update_version_h ALL
                  python ${POLYMEC_PREFIX}/share/polymec/update_version_h.py polyamri ${POLYAMRI_VERSION} ${PROJECT_BINARY_DIR}/polyamri/polyamri_version.h)

# Generate polyamri.h at the top level.
add_custom_target(generate_top_header ALL
                  python ${PROJECT_SOURCE_DIR}/tools/generate_top_header.py ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)
                  
enable_testing()

add_subdirectory(polyamri)

# Generate a polyamri.cmake file that contains installation information.
set(CMAKE_INSTALL_PREFIX ${POLYMEC_PREFIX})
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Templates/polyamri.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/polyamri.cmake"
)
install(FILES ${PROJECT_BINARY_DIR}/include/polyamri.h DESTINATION include/polymec)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/polyamri.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polyamri_executable.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polyamri_test.cmake DESTINATION share/polymec)

