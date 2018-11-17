cmake_minimum_required(VERSION 3.0.0)

project(SCIPOptSuite VERSION 6.0.0)

set(ZIMPL_DIR ${CMAKE_BINARY_DIR})
set(SOPLEX_DIR ${CMAKE_BINARY_DIR})
set(SCIP_DIR ${CMAKE_BINARY_DIR})

# path to e.g. findGMP module
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scip/cmake/Modules)

option(GMP "should GMP be linked" ON)
if(GMP)
  find_package(BISON)
  find_package(FLEX)
  find_package(GMP)
  # only add ZIMPL subdirectory if required packages could be found
  if(${BISON_FOUND} AND ${FLEX_FOUND} AND ${GMP_FOUND})
    add_subdirectory(zimpl)
  endif()
endif()

add_subdirectory(soplex)
add_subdirectory(scip)

if(EXISTS ${PROJECT_SOURCE_DIR}/gcg/CMakeLists.txt)
  add_subdirectory(gcg)
endif()
if(EXISTS ${PROJECT_SOURCE_DIR}/ug/CMakeLists.txt)
  add_subdirectory(ug)
endif()

set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${SCIPOptSuite_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${SCIPOptSuite_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${SCIPOptSuite_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "Zuse Institute Berlin")
set(CPACK_NSIS_MODIFY_PATH on)
set(CPACK_PACKAGE_EXECUTABLES scip;SCIP soplex;SoPlex)
set(CPACK_PACKAGE_CONTACT "SCIP <scip@zib.de>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The SCIP Optimization Suite is a toolbox for generating and solving mixed integer nonlinear programs, in particular mixed integer linear programs, and constraint integer programs.")
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "http://scip.zib.de")
# autogenerate dependency information
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS on)
include(CPack)

enable_testing()

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  include(FeatureSummary)
  feature_summary(WHAT ALL)
endif()

