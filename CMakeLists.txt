cmake_minimum_required (VERSION 3.1)
project (bbmutils)

add_definitions(-D_SCL_SECURE_NO_WARNINGS)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)
add_definitions(-D_USE_MATH_DEFINES)

include_directories(${CMAKE_SOURCE_DIR}/..)

add_executable (test_2d "${CMAKE_SOURCE_DIR}/src/test_2d")
add_executable (test_3d "${CMAKE_SOURCE_DIR}/src/test_3d")
add_executable (test_6d "${CMAKE_SOURCE_DIR}/src/test_6d")
add_executable (test_es "${CMAKE_SOURCE_DIR}/src/test_es")

add_library(bbmutils ${CMAKE_SOURCE_DIR}/src/bbmutils.c)

target_link_libraries(test_2d 
	bbmutils
	m
) 
target_link_libraries(test_3d 
	bbmutils
	m
) 
target_link_libraries(test_6d 
	bbmutils
	m
) 
target_link_libraries(test_es 
	bbmutils
	m
) 

install (TARGETS test_2d test_3d test_6d test_es
	DESTINATION ${CMAKE_SOURCE_DIR}/bin/)


