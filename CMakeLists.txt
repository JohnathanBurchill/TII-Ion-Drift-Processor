# TII Cross-Track Ion Drift Processor: CMakeLists.txt

# Copyright (C) 2022  Johnathan K Burchill

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


PROJECT(TIICT)

CMAKE_MINIMUM_REQUIRED(VERSION 3.3)

ADD_DEFINITIONS(-DHAVE_GETOPT_H)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    message( "-- Debugging enabled")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -std=gnu99 ")
else (CMAKE_BUILD_TYPE STREQUAL Debug)
    message( "-- Optimization enabled")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -std=gnu99 ")
endif(CMAKE_BUILD_TYPE STREQUAL Debug)

# GSL
FIND_PACKAGE(GSL REQUIRED)

# requires -lm on linux
FIND_LIBRARY(MATH m)
SET(LIBS ${LIBS} ${MATH} ${GSL_LIBRARY} )


INCLUDE_DIRECTORIES(${INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})

ADD_EXECUTABLE(tiict tiict.c utilities.c processing.c)
TARGET_INCLUDE_DIRECTORIES(tiict PRIVATE)
TARGET_LINK_LIBRARIES(tiict ${LIBS} -lgslcblas -lcdf)

ADD_EXECUTABLE(tiictinspect tiictinspect.c utilities.c)
TARGET_INCLUDE_DIRECTORIES(tiictinspect PRIVATE)
TARGET_LINK_LIBRARIES(tiictinspect ${LIBS} -lgslcblas -lcdf)

ADD_EXECUTABLE(tiictqualityflags tiictqualityflags.c utilities.c)
TARGET_INCLUDE_DIRECTORIES(tiictqualityflags PRIVATE )
TARGET_LINK_LIBRARIES(tiictqualityflags ${MATH} -lcdf)

ADD_EXECUTABLE(tiictvariability tiictvariability.c utilities.c)
TARGET_INCLUDE_DIRECTORIES(tiictvariability PRIVATE )
TARGET_LINK_LIBRARIES(tiictvariability ${LIBS} -lgslcblas -lcdf)

ADD_EXECUTABLE(tiictbin tiictbin.c utilities.c)
TARGET_INCLUDE_DIRECTORIES(tiictbin PRIVATE )
TARGET_LINK_LIBRARIES(tiictbin ${LIBS} -lgslcblas -lcdf)

# ADD_EXECUTABLE(tiictfitstatistics tiictfitstatistics.c)
# TARGET_INCLUDE_DIRECTORIES(tiictfitstatistics PRIVATE )
# TARGET_LINK_LIBRARIES(tiictfitstatistics ${MATH})

install(TARGETS tiict DESTINATION $ENV{HOME}/bin)
install(TARGETS tiictinspect DESTINATION $ENV{HOME}/bin)
install(TARGETS tiictqualityflags DESTINATION $ENV{HOME}/bin)
install(TARGETS tiictbin DESTINATION $ENV{HOME}/bin)
# install(TARGETS tiictfitstatistics DESTINATION $ENV{HOME}/bin)


