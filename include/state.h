/*

    TII Cross-Track Ion Drift Processor: state.h

    Copyright (C) 2024  Johnathan K Burchill

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
#include <gsl/gsl_multifit.h>

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

typedef struct BackgroundRemovalWorkspace
{
    long long beginIndex0;
    long long beginIndex1;
    long long endIndex0;
    long long endIndex1;
    long long modelDataIndex;
    long long modelDataMidPoint;
    gsl_vector *modelValues;
    gsl_vector *model1Values;
    gsl_vector *model2Values;
    long long numModelPoints; // data for fit, and for median calculations at each end
    long long numModel1Points;
    long long numModel2Points;
    gsl_vector *fitCoefficients;
    gsl_vector *work1;
    gsl_vector *work2;
    gsl_matrix *modelTimesMatrix;
    gsl_matrix *modelTimes1Matrix;
    gsl_matrix *modelTimes2Matrix;
    double epoch0;
    double tregion11;
    double tregion12;
    double tregion21;
    double tregion22;
    char startString[EPOCH_STRING_LEN+1];
    char stopString[EPOCH_STRING_LEN+1];
    size_t fitDegree; // fit degree (e.g., linear is 2) 

} BackgroundRemovalWorkspace_t;

typedef struct ProcessorState {

    int returnStatus;

    Arguments args;
    int nOptions;
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
    uint8_t *region;

    float *xhat;
    float *yhat;
    float *zhat;
    float *ectFieldH;
    float *ectFieldV;
    float *bctField;

    float *geoPotential;
    float *maxAbsGeopotentialSlope;
    float *geoPotentialDifference;
    float *geoPotentialDetrended;
    float *maxAbsGeopotentialDetrendedSlope;
    float *exAdjusted;
    float *exAdjustmentParameter;
    
    // Offset removal options
    uint8_t interval;
    bool setFlags;
    BackgroundRemovalWorkspace_t bgws;

} ProcessorState;

#endif // _STATE_H
