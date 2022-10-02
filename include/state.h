/*

    TII Cross-Track Ion Drift Processor: state.h

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

#ifndef _STATE_H
#define _STATE_H

#include "settings.h"
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <cdf.h>

typedef struct offset_model_fit_arguments {
    const uint8_t regionNumber;
    const char *regionName;
    const float lat1;
    const float lat2;
    const float lat3;
    const float lat4;
} offset_model_fit_arguments;

typedef struct Arguments {
    const char* satellite;
    int year;
    int month;
    int day;
    const char* calVersion;
    const char* exportVersion;
    const char* calDir;
    const char* lpDir;
    const char* exportDir;
} Arguments;

typedef struct ProcessorState {

    int returnStatus;

    Arguments args;
    offset_model_fit_arguments fitargs[4];
    char fitLogFilename[FILENAME_MAX];
    FILE *fitFile;
    char processingLogFilename[FILENAME_MAX];
    FILE *processingLogFile;

    // Calibration CDF data
    char calibrationFileName[CDF_PATHNAME_LEN];
    uint8_t * dataBuffers[NUM_CAL_VARIABLES];
    long nRecs;
    long memoryAllocated;

    // LP data for floating potential
    double *lpTimes;
    float *lpPhiScHighGain;
    float *lpPhiScLowGain;
    float *lpPhiSc;
    size_t nLpRecs;
    float *potentials;
    bool usePotentials;

    // Additional data
    float *viErrors;
    uint16_t *flags;
    uint32_t *fitInfo;

    float *xhat;
    float *yhat;
    float *zhat;
    float *ectFieldH;
    float *ectFieldV;
    float *bctField;

} ProcessorState;

#endif // _STATE_H
