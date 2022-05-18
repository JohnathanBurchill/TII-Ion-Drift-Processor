/*

    TII Cross-Track Ion Drift Processor: CrossTrackCalibration.h

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

#ifndef CROSSTRACKCALIBRATION_H
#define CROSSTRACKCALIBRATION_H

#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include <cdf.h>

#define SOFTWARE_VERSION_STRING "CrossTrackCalibration 2020-07-31"
#define SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING 300 // 5 minutes
#define SECONDS_OF_DATA_REQUIRED_FOR_EXPORTING 60 // 1 minutes
#define SECONDS_OF_BOUNDARY_DATA_REQUIRED_FOR_PROCESSING 1500 // 25 minutes to ensure coverage of each fit region

#define FLAGS_MAXIMUM_DRIFT_VALUE 8000.0 // 8 km/s maximum drift for flagging
#define DEFAULT_VI_ERROR -42.0
#define MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT 8
enum FITINFO_BIT_MASKS {
    FITINFO_OFFSET_NOT_REMOVED = (1 << 0),
    FITINFO_INCOMPLETE_REGION = (1 << 1),
    FITINFO_GSL_FIT_ERROR = (1 <<2),
    FITINFO_MAD_EXCEEDED = (1 << 3),
    FITINFO_DRIFT_MAGNITUDE_EXCEEDED = (1 << 4)
};
#define MAX_ALLOWED_CDF_GAP_SECONDS 600.0 // CDF export split into separate files at gaps exceeding 10 minutes
#define NUM_CAL_VARIABLES 15
#define NUM_EXPORT_VARIABLES 31
#define NUM_BUFFER_VARIABLES 30

#define GSL_FIT_MAXIMUM_ITERATIONS 500

void closeCdf(CDFid id);

#define ADDR(n, m, d) (((float*)dataBuffers[(n)]+(d*timeIndex + m)))
#define MEAS(n, m, d) ((float)(*(ADDR(n, m, d))))
#define TIME() ((double)*((double*)dataBuffers[0]+(timeIndex)))
#define MXH() (MEAS(1, 0, 2))
#define MYH() (MEAS(1, 1, 2))
#define MXV() (MEAS(2, 0, 2))
#define MYV() (MEAS(2, 1, 2))
#define VSATX() (MEAS(3, 0, 3))
#define VSATY() (MEAS(3, 1, 3))
#define VSATZ() (MEAS(3, 2, 3))
#define MLT() (MEAS(4, 0, 1))
#define QDLAT() (MEAS(5, 0, 1))
#define QDLON() (MEAS(6, 0, 1))
#define LAT() (MEAS(7, 0, 1))
#define LON() (MEAS(8, 0, 1))
#define RADIUS() (MEAS(9, 0, 1))
#define VCORX() (MEAS(10, 0, 3))
#define VCORY() (MEAS(10, 1, 3))
#define VCORZ() (MEAS(10, 2, 3))
#define VSATN() (MEAS(11, 0, 3))
#define VSATE() (MEAS(11, 1, 3))
#define VSATC() (MEAS(11, 2, 3))
#define BN() (MEAS(12, 0, 3))
#define BE() (MEAS(12, 1, 3))
#define BC() (MEAS(12, 2, 3))
#define VBIAS() (MEAS(13, 0, 1))
#define VFP() (MEAS(14, 0, 1))

typedef struct offset_model_fit_arguments {
    const uint8_t regionNumber;
    const char *regionName;
    const float lat1;
    const float lat2;
    const float lat3;
    const float lat4;
} offset_model_fit_arguments;

typedef enum DayType {
    PREVIOUS_DAY = -1,
    REQUESTED_DAY = 0,
    NEXT_DAY = 1
} DayType;

void loadCalData(const char *calDir, const char *calVersion, const char *satellite, const int year, const int month, const int day, uint8_t **dataBuffers, long *nRecs);

void setCalibrationFileName(const char *satellite, const int year, const int month, const int day, const char *calDir, const char *calVersion, char *calibrationFileName);

void loadCalDataFromDate(const DayType dayType, const char *calDir, const char *calVersion, const char *satellite, struct tm timestructure, uint8_t **dataBuffers, long *nRecs, long *calibrationMemorySize);

void removeOffsetsAndSetFlags(const char* satellite, offset_model_fit_arguments fitargs, long nRecs, uint8_t **dataBuffers, float *viErrors, uint16_t *flags, uint32_t *fitInfo, FILE* fitFile);

void updateDataQualityFlags(const char *satellite, uint8_t sensorIndex, uint8_t regionNumber, float driftValue, float mad, long timeIndex, uint16_t *flags, uint32_t *fitInfo);

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

float madThreshold(char satellite, int sensorIndex);


#endif // CROSSTRACKCALIBRATION_H
