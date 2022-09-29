/*

    TII Cross-Track Ion Drift Processor: tiiinspect.h

    Copyright (C) 2022  Johnathan K Burchill

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef TIICTINSPECT_H
#define TIICTINSPECT_H

#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include <cdf.h>

#define SOFTWARE_VERSION_STRING "tiictinspect 2022-05-19"

#define NUM_DATA_VARIABLES 6

void loadCrossTrackData(const char * filename, uint8_t **dataBuffers, long *numberOfRecords, bool *fourByteCalFlag);


#define ADDR(n, m) (((float*)dataBuffers[(n)]+(timeIndex + m)))
#define MEAS(n, m) ((float)(*(ADDR(n, m))))
#define TIME() ((double)*((double*)dataBuffers[0]+(timeIndex)))
#define VIY() (MEAS(1, 0))
#define VIYERROR() (MEAS(2, 0))
#define QDLAT() (MEAS(3, 0))
#define FLAG() ((uint16_t)*((uint16_t*)dataBuffers[4]+(timeIndex)))
#define CALFLAG() ((uint32_t)*((uint32_t*)dataBuffers[5]+(timeIndex)))

void analyzeMAD(const char satellite, uint8_t majorVersion, uint8_t minorVersion, bool fourByteCalFlag, uint8_t * dataBuffers[], long nRecs);

void analyzeTimeGaps(const char *filename, uint8_t * dataBuffers[], long nRecs, long *totalCount, long *gapCount);
void summarizeTimeGaps(const char *filename, uint8_t * dataBuffers[], long nRecs, long *totalCount, long *gapCount, long *largeGapCount);

void analyzeCalibrationFlag(const char satellite, uint8_t majorVersion, uint8_t minorVersion, bool fourByteCalFlag, uint8_t * dataBuffers[], long nRecs, bool numericOutput);

uint8_t getMajorVersion(const char *filename);
uint8_t getMinorVersion(const char *filename);

void getFileDate(const char *filename, char *startDate);

typedef enum AnalysisType {
    NO_ANALYSIS = 0,
    MAD_ANALYSIS = 1,
    CALIBRATION_ANALYSIS_NUMERIC = 2,
    CALIBRATION_ANALYSIS_GRAPHIC = 3,
    CALIBRATION_ANALYSIS_TIMEGAP = 4,
    CALIBRATION_ANALYSIS_TIMEGAPSUMMARY = 5
} AnalysisType;

#endif // TIICTINSPECT_H
