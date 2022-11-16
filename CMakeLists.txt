cmake_minimum_required (VERSION 3.1)

if(DEFINED ENV{VCPKG_DEFAULT_TRIPLET})
  message(STATUS "Setting default vcpkg target triplet to $ENV{VCPKG_DEFAULT_TRIPLET}")
  set(VCPKG_TARGET_TRIPLET $ENV{VCPKG_DEFAULT_TRIPLET})
endif()

if(ENABLE_TEST)
  list(APPEND VCPKG_MANIFEST_FEATURES "test")
endif()

if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
  message(STATUS "VCPKG found: $ENV{VCPKG_ROOT}")
  message(STATUS "Using VCPKG integration")
  message(STATUS "VCPKG_MANIFEST_FEATURES: ${VCPKG_MANIFEST_FEATURES}")
elseif(DEFINED CMAKE_TOOLCHAIN_FILE)
  message(STATUS "Using toolchain: ${CMAKE_TOOLCHAIN_FILE}")
  if(CMAKE_TOOLCHAIN_FILE MATCHES "vcpkg.cmake")
    message(STATUS "Toolchain uses VCPKG integration")
    message(STATUS "VCPKG_MANIFEST_FEATURES: ${VCPKG_MANIFEST_FEATURES}")
  endif()
else()
  message(WARNING "vcpkg not found, toolchain not defined, build might fail")
endif()

project (bbmutils)

include(CheckLibraryExists)

CHECK_LIBRARY_EXISTS(m sin "" HAVE_LIB_M)

if (HAVE_LIB_M)
  list(APPEND EXTRA_LIBS m)
endif (HAVE_LIB_M)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}" CACHE PATH "Install prefix" FORCE)
endif()

set(INSTALL_BIN_DIR      "bin"                     CACHE PATH "Path where exe and dll will be installed")
set(INSTALL_LIB_DIR      "lib"                     CACHE PATH "Path where lib will be installed")
set(INSTALL_INCLUDE_DIR  "include"                 CACHE PATH "Path where headers will be installed")
set(INSTALL_CMAKE_DIR    "share/${PROJECT_NAME}"   CACHE PATH "Path where cmake configs will be installed")

add_library(bbmutils ${CMAKE_CURRENT_LIST_DIR}/src/bbmutils.h ${CMAKE_CURRENT_LIST_DIR}/src/bbmutils.c)
target_include_directories(bbmutils PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src> $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>)
set_target_properties(bbmutils PROPERTIES PUBLIC_HEADER "${CMAKE_CURRENT_LIST_DIR}/src/bbmutils.h")

install(TARGETS bbmutils EXPORT BBMUtilsTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}"
  COMPONENT dev
)

if(ENABLE_TEST)
  add_executable (test_2d "${CMAKE_CURRENT_LIST_DIR}/test/test_2d.c")
  add_executable (test_3d "${CMAKE_CURRENT_LIST_DIR}/test/test_3d.c")
  add_executable (test_6d "${CMAKE_CURRENT_LIST_DIR}/test/test_6d.c")
  add_executable (test_es "${CMAKE_CURRENT_LIST_DIR}/test/test_es.c")

  target_link_libraries(test_2d PRIVATE bbmutils ${EXTRA_LIBS})
  target_link_libraries(test_3d PRIVATE bbmutils ${EXTRA_LIBS})
  target_link_libraries(test_6d PRIVATE bbmutils ${EXTRA_LIBS})
  target_link_libraries(test_es PRIVATE bbmutils ${EXTRA_LIBS})

  if(MSVC)
    target_compile_definitions(test_2d PRIVATE _USE_MATH_DEFINES)
    target_compile_definitions(test_2d PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_definitions(test_3d PRIVATE _USE_MATH_DEFINES)
    target_compile_definitions(test_3d PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_definitions(test_6d PRIVATE _USE_MATH_DEFINES)
    target_compile_definitions(test_6d PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_definitions(test_es PRIVATE _USE_MATH_DEFINES)
    target_compile_definitions(test_es PRIVATE _CRT_SECURE_NO_WARNINGS)
  endif()
  install (TARGETS test_2d test_3d test_6d test_es DESTINATION "${INSTALL_BIN_DIR}")
endif()

install(EXPORT BBMUtilsTargets
  FILE BBMUtilsConfig.cmake
  NAMESPACE BBMUtils::
  DESTINATION "${INSTALL_CMAKE_DIR}"
)
