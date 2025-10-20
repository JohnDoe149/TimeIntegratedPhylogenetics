# Install script for directory: /Users/work/Desktop/JohnPhylo/lib/ncl/ncl

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/ncl" TYPE SHARED_LIBRARY FILES "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/libncl.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ncl/libncl.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ncl/libncl.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ncl/libncl.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/CMakeFiles/ncl_shared.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/ncl" TYPE STATIC_LIBRARY FILES "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/libncl.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ncl/libncl.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ncl/libncl.a")
    execute_process(COMMAND "/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ncl/libncl.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/CMakeFiles/ncl_static.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ncl" TYPE FILE FILES
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/ncl.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsallocatematrix.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsassumptionsblock.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsblock.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxscdiscretematrix.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxscharactersblock.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxscxxdiscretematrix.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsdatablock.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsdefs.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsdiscretedatum.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsdistancedatum.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsdistancesblock.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsexception.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsmultiformat.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxspublicblocks.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsreader.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxssetreader.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsstring.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxstaxaassociationblock.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxstaxablock.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxstoken.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxstreesblock.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsunalignedblock.h"
    "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/nxsutilcopy.h"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/work/Desktop/JohnPhylo/lib/ncl/ncl/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
