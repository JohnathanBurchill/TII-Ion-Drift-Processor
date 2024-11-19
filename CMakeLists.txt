# TII Cross-Track Ion Drift Processor: CMakeLists.txt

# Copyright (C) 2024  Johnathan K Burchill

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


CMAKE_MINIMUM_REQUIRED(VERSION 3.20)

PROJECT(TIICT)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

ADD_DEFINITIONS(-DHAVE_GETOPT_H)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    message( "-- Debugging enabled")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg -std=gnu99 ")
else (CMAKE_BUILD_TYPE STREQUAL Debug)
    message( "-- Optimization enabled")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -std=gnu99 ")
endif(CMAKE_BUILD_TYPE STREQUAL Debug)

INCLUDE(GNUInstallDirs)

# GSL
FIND_PACKAGE(GSL REQUIRED)
# SDL3
find_library(SDL3 SDL3 HINTS /usr/local/lib REQUIRED)

# From ChatGPT 18 Nov 2024
# Set RPATH for build and install
set(CMAKE_BUILD_RPATH "/usr/local/lib")
set(CMAKE_INSTALL_RPATH "/usr/local/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Prevent clearing RPATH during installation
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_SKIP_INSTALL_RPATH FALSE)

# Graphics stuff
FIND_PACKAGE(PkgConfig REQUIRED)
pkg_check_modules(LIBAV REQUIRED IMPORTED_TARGET
    libavformat
    libavcodec
    libswresample
    libswscale
    libavutil
)

# requires -lm on linux
FIND_LIBRARY(MATH m REQUIRED)
FIND_LIBRARY(TII tii REQUIRED)
FIND_LIBRARY(TIIGRAPHICS tiigraphics REQUIRED)
FIND_LIBRARY(CDF cdf HINTS /Applications/cdf/cdf/lib REQUIRED)
SET(LIBS ${LIBS} ${MATH} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${CDF} ${TIIGRAPHICS} ${TII})

INCLUDE_DIRECTORIES(include ${CMAKE_INSTALL_INCUDEDIR} ${INCLUDE_DIRS}
${GSL_INCLUDE_DIRS} /Applications/cdf/cdf/include /usr/local/include)

ADD_LIBRARY(Utilities utilities.c)
ADD_LIBRARY(Processor processing.c export.c loadData.c utilities.c visualize.c)
TARGET_INCLUDE_DIRECTORIES(Processor PRIVATE)
TARGET_LINK_LIBRARIES(Processor ${LIBS} PkgConfig::LIBAV)

ADD_EXECUTABLE(tiict0401 tiict.c)
TARGET_INCLUDE_DIRECTORIES(tiict0401 PRIVATE)
TARGET_LINK_LIBRARIES(tiict0401 PRIVATE Processor ${LIBS} ${CDF})

ADD_EXECUTABLE(tiictinteractive0401 tiictinteractive.c)
TARGET_INCLUDE_DIRECTORIES(tiictinteractive0401 PRIVATE)
TARGET_LINK_LIBRARIES(tiictinteractive0401 PRIVATE ${SDL3} Processor ${LIBS})

ADD_EXECUTABLE(tiictinspect tiictinspect.c)
TARGET_INCLUDE_DIRECTORIES(tiictinspect PRIVATE)
TARGET_LINK_LIBRARIES(tiictinspect Utilities ${LIBS})

ADD_EXECUTABLE(tiictqualityflags tiictqualityflags.c)
TARGET_INCLUDE_DIRECTORIES(tiictqualityflags PRIVATE )
TARGET_LINK_LIBRARIES(tiictqualityflags Utilities ${LIBS})

ADD_EXECUTABLE(tiictvariability tiictvariability.c)
TARGET_INCLUDE_DIRECTORIES(tiictvariability PRIVATE )
TARGET_LINK_LIBRARIES(tiictvariability Utilities ${LIBS})

ADD_SUBDIRECTORY(util)

INSTALL(TARGETS tiict0401 DESTINATION $ENV{HOME}/bin)
INSTALL(TARGETS tiictinteractive0401 DESTINATION $ENV{HOME}/bin)
INSTALL(TARGETS tiictinspect DESTINATION $ENV{HOME}/bin)
INSTALL(TARGETS tiictqualityflags DESTINATION $ENV{HOME}/bin)


