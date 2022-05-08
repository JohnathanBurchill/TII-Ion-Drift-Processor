# CMakeLists.txt file for processor

PROJECT(CrossTrackCalibration)

CMAKE_MINIMUM_REQUIRED(VERSION 3.3)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/CMake")

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


INCLUDE_DIRECTORIES(${INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} /sw/include/gsl)

ADD_EXECUTABLE(CrossTrackCalibration CrossTrackCalibration.c utilities.c)
TARGET_INCLUDE_DIRECTORIES(CrossTrackCalibration PRIVATE /Applications/cdf/cdf/include)
TARGET_LINK_LIBRARIES(CrossTrackCalibration ${LIBS})

ADD_EXECUTABLE(CrossTrackValidation CrossTrackValidation.c utilities.c)
TARGET_INCLUDE_DIRECTORIES(CrossTrackValidation PRIVATE /Applications/cdf/cdf/include)
TARGET_LINK_LIBRARIES(CrossTrackValidation ${LIBS})

ADD_EXECUTABLE(CrossTrackQualityFlagDetermination CrossTrackQualityFlagDetermination.c utilities.c)
TARGET_INCLUDE_DIRECTORIES(CrossTrackQualityFlagDetermination PRIVATE /Applications/cdf/cdf/include)
TARGET_LINK_LIBRARIES(CrossTrackQualityFlagDetermination ${LIBS})

if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    message("-- Darwin linking")
    TARGET_LINK_LIBRARIES(CrossTrackCalibration -lgslcblas /Applications/cdf/cdf/lib/libcdf.a)
    TARGET_LINK_LIBRARIES(CrossTrackValidation -lgslcblas /Applications/cdf/cdf/lib/libcdf.a)
    TARGET_LINK_LIBRARIES(CrossTrackQualityFlagDetermination -lgslcblas /Applications/cdf/cdf/lib/libcdf.a)
endif()
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    message("-- Linux linking")
    TARGET_LINK_LIBRARIES(CrossTrackCalibration -lgslcblas -lcdf)
    TARGET_LINK_LIBRARIES(CrossTrackValidation -lgslcblas -lcdf)
    TARGET_LINK_LIBRARIES(CrossTrackQualityFlagDetermination -lgslcblas -lcdf)
endif ()

install(TARGETS CrossTrackCalibration DESTINATION $ENV{HOME}/bin)
install(TARGETS CrossTrackValidation DESTINATION $ENV{HOME}/bin)
install(TARGETS CrossTrackQualityFlagDetermination DESTINATION $ENV{HOME}/bin)
install(PROGRAMS process_TCT.sh parallel_process_TCT.sh compress_TCT.sh parallel_compress_TCT.sh cdf_files_to_compress.sh DESTINATION $ENV{HOME}/bin)

