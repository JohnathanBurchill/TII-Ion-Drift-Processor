/*

    TII Cross-Track Ion Drift Processor: tiict.h

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

#ifndef TIICT_H
#define TIICT_H

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>

#include <cdf.h>

void closeCdf(CDFid id);

typedef enum DayType {
    PREVIOUS_DAY = -1,
    REQUESTED_DAY = 0,
    NEXT_DAY = 1
} DayType;

void loadCalData(const char *calDir, const char *calVersion, const char *satellite, const int year, const int month, const int day, uint8_t **dataBuffers, long *nRecs);

void setCalibrationFileName(const char *satellite, const int year, const int month, const int day, const char *calDir, const char *calVersion, char *calibrationFileName);

void loadCalDataFromDate(const DayType dayType, const char *calDir, const char *calVersion, const char *satellite, struct tm timestructure, uint8_t **dataBuffers, long *nRecs, long *calibrationMemorySize);

CDFstatus createVarFrom1DVar(CDFid id, char *name, long dataType, long startIndex, long stopIndex, void *buffer);
CDFstatus createVarFrom2DVar(CDFid id, char *name, long dataType, long startIndex, long stopIndex, void *buffer, uint8_t index, uint8_t dimSize);


CDFstatus addgEntry(CDFid id, long attrNum, long entryNum, const char *entry);


typedef struct varAttr {
    char * name;
    char * type;
    char * units;
    char * desc;
    double validMin;
    double validMax;
} varAttr;

CDFstatus addVariableAttributes(CDFid id, varAttr attr);

void addAttributes(CDFid id, const char *calVersion, const char *satellite, const char *version, double minTime, double maxTime);

long numberOfAvailableRecordsForDate(const char *satellite, const int year, const int month, const int day, const char *calDir, const char *calVersion);

void exportTCT16Cdfs(double startTime, double stopTime, const char *exportDir, const char *exportVersion, const char *calVersion, const char *satellite, long startIndex, long stopIndex, uint8_t **dataBuffers, float *ectFieldH, float *ectFieldV, float* bctField, float *viErrors, uint16_t *flags, uint32_t *fitInfo);

void exportTCT02Cdfs(double startTime, double stopTime, const char *exportDir, const char *exportVersion, const char *calVersion, const char *satellite, long startIndex, long stopIndex, uint8_t **dataBuffers, float *ectFieldH, float *ectFieldV, float* bctField, float *viErrors, uint16_t *flags, uint32_t *fitInfo);

// Returns true if a full 8 samples were downsampled to 1 sample, false otherwise
bool downSampleHalfSecond(long *index, long storageIndex, double t0, long maxIndex, uint8_t **dataBuffers, float *ectFieldH, float *ectFieldV, float *bctField, float *viErrors, uint16_t *flags, uint32_t *fitInfo);


#endif // TIICT_H
