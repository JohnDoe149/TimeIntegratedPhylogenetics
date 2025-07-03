# Install script for directory: /workspaces/JohnPhylo/lib/ncl/ncl

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

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ncl/libncl.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ncl/libncl.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ncl/libncl.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/ncl" TYPE SHARED_LIBRARY FILES "/workspaces/JohnPhylo/build/lib/ncl/ncl/libncl.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ncl/libncl.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ncl/libncl.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ncl/libncl.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/ncl" TYPE STATIC_LIBRARY FILES "/workspaces/JohnPhylo/build/lib/ncl/ncl/libncl.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ncl" TYPE FILE FILES
    "/workspaces/JohnPhylo/lib/ncl/ncl/ncl.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsallocatematrix.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsassumptionsblock.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsblock.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxscdiscretematrix.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxscharactersblock.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxscxxdiscretematrix.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsdatablock.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsdefs.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsdiscretedatum.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsdistancedatum.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsdistancesblock.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsexception.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsmultiformat.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxspublicblocks.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsreader.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxssetreader.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsstring.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxstaxaassociationblock.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxstaxablock.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxstoken.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxstreesblock.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsunalignedblock.h"
    "/workspaces/JohnPhylo/lib/ncl/ncl/nxsutilcopy.h"
    )
endif()

