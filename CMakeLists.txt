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


CMAKE_MINIMUM_REQUIRED(VERSION 3.30)

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

# Graphics stuff
FIND_PACKAGE(PkgConfig REQUIRED)
pkg_check_modules(LIBAV REQUIRED IMPORTED_TARGET
    libavformat
    libavcodec
    libswresample
    libswscale
    libavutil
)

SET(CMAKE_SKIP_RPATH FALSE)

# requires -lm on linux
FIND_LIBRARY(MATH m REQUIRED)
FIND_LIBRARY(TII tii REQUIRED)
FIND_LIBRARY(TIIGRAPHICS tiigraphics REQUIRED)
FIND_LIBRARY(CDF cdf HINTS /Applications/cdf/cdf/lib REQUIRED)
SET(LIBS ${LIBS} ${MATH} ${GSL_LIBRARY} ${TII})

INCLUDE_DIRECTORIES(include ${CMAKE_INSTALL_INCUDEDIR} ${INCLUDE_DIRS}
${GSL_INCLUDE_DIRS} /Applications/cdf/cdf/include /usr/local/include)

ADD_LIBRARY(Utilities utilities.c)
ADD_LIBRARY(Processor processing.c export.c loadData.c utilities.c visualize.c)
TARGET_INCLUDE_DIRECTORIES(Processor PRIVATE)
TARGET_LINK_LIBRARIES(Processor ${LIBS} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${TIIGRAPHICS} PkgConfig::LIBAV)

ADD_EXECUTABLE(tiict0401 tiict.c)
TARGET_INCLUDE_DIRECTORIES(tiict0401 PRIVATE)
TARGET_LINK_LIBRARIES(tiict0401 PRIVATE Processor ${LIBS} ${CDF})

#ADD_EXECUTABLE(tiictinteractive0401 tiictinteractive.c)
#TARGET_INCLUDE_DIRECTORIES(tiictinteractive0401 PRIVATE)
#TARGET_LINK_LIBRARIES(tiictinteractive0401 PRIVATE Processor ${LIBS} ${CDF})

ADD_EXECUTABLE(tiictinspect tiictinspect.c)
TARGET_INCLUDE_DIRECTORIES(tiictinspect PRIVATE)
TARGET_LINK_LIBRARIES(tiictinspect Utilities ${LIBS} ${CDF})

ADD_EXECUTABLE(tiictqualityflags tiictqualityflags.c)
TARGET_INCLUDE_DIRECTORIES(tiictqualityflags PRIVATE )
TARGET_LINK_LIBRARIES(tiictqualityflags Utilities ${LIBS} ${CDF})

ADD_EXECUTABLE(tiictvariability tiictvariability.c)
TARGET_INCLUDE_DIRECTORIES(tiictvariability PRIVATE )
TARGET_LINK_LIBRARIES(tiictvariability Utilities ${LIBS} ${CDF})

ADD_SUBDIRECTORY(util)

INSTALL(TARGETS tiict0401 DESTINATION $ENV{HOME}/bin)
INSTALL(TARGETS tiictinspect DESTINATION $ENV{HOME}/bin)
INSTALL(TARGETS tiictqualityflags DESTINATION $ENV{HOME}/bin)


