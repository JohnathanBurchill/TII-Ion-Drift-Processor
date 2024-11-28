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

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <cdf.h>
#include <gsl/gsl_multifit.h>

#include <tiigraphics/draw.h>

typedef enum MeasurementType {
    VELOCITY_MEASUREMENT = 0,
    ENERGY_MEASUREMENT,
} MeasurementType_enum;


typedef struct offset_model_fit_arguments {
    const uint8_t regionNumber;
    const char *regionName;
    const float lat1;
    const float lat2;
    const float lat3;
    const float lat4;
} offset_model_fit_arguments;

typedef struct Arguments {
    int argc;
    char **argv;
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

typedef struct ProcessorVariables {

    size_t memoryAllocated;
    long nRecs;

    // Exported variables
    // If changing this list, check freeing of memory and downsampling
    // algorithm
    double *timestamp;
    float *latitude;
    float *longitude;
    float *radius;
    float *qdlat;
    float *mlt;
    float *mxh;
    float *myh;
    float *mxv;
    float *myv;
    float *dxh;
    float *dyh;
    float *dxv;
    float *dyv;
    float *rh;
    float *rv;
    float *vmcph;
    float *vmcpv;
    float *vbiash;
    float *vbiasv;
    float *vfp;
    float *vsatx;
    float *vsaty;
    float *vsatz;
    float *enhRaw;
    float *envRaw;
    float *enh;
    float *env;
    float *vixh;
    float *vixherror;
    float *vixv;
    float *vixverror;
    float *viy;
    float *viyerror;
    float *viz;
    float *vizerror;
    float *vsatn;
    float *vsate;
    float *vsatc;
    float *bn;
    float *be;
    float *bc;
    float *vicrx;
    float *vicry;
    float *vicrz;
    float *geoelectricPotential;
    float *geoelectricPotentialDifference;
    float *maxAbsGeoelectricPotentialBaselineSlope;
    float *ehxAdjusted;
    float *ehxAdjustmentParameter;
    float *geoelectricPotentialDetrended;
    float *maxAbsGeoelectricPotentialDetrendedBaselineSlope;
    uint8_t *orbitRegion;
    uint16_t *flags;
    uint32_t *fitInfo;

    // LP
    int lpPotentialSource;
    size_t nLpRecs;
    double *lpTimes;
    float *lpPhiScHighGain;
    float *lpPhiScLowGain;
    float *lpPhiSc;
    float *potentials;

    // Internal variables
    float *xhat;
    float *yhat;
    float *zhat;
    float *ectxh;
    float *ectyh;
    float *ectzh;
    float *ectxv;
    float *ectyv;
    float *ectzv;
    float *bctx;
    float *bcty;
    float *bctz;


} ProcessorVariables_t;

typedef struct ProcessorState {

    int returnStatus;

    Arguments args;
    int nOptions;
    offset_model_fit_arguments fitargs[4];
    char fitLogFilename[FILENAME_MAX];
    FILE *fitFile;
    char processingLogFilename[FILENAME_MAX];
    FILE *processingLogFile;
    bool writeLogFiles;
    char processingDateString[32];

    // Calibration CDF data
    char tiiCalibrationFileName[CDF_PATHNAME_LEN];
    char lpCalibrationFileName[CDF_PATHNAME_LEN];
    // 16 Hz
    ProcessorVariables_t vars16hz;
    // 2 Hz
    ProcessorVariables_t vars2hz;
    // Pointer to variables in use
    ProcessorVariables_t *vars;

    bool useEofR;
    bool usePotentials;

    // Offset removal options
    uint8_t interval;
    bool setFlags;
    BackgroundRemovalWorkspace_t bgws;
    MeasurementType_enum measurementType;

    bool export2Hz;
    bool export16Hz;
    bool exportZip;
    bool exportVideo;

    bool visualizeResults;
    char *plotCommand;
    int defaultPlotHeight;
    int maxPlotsPerScreen;

    char *videoOutputDir;
    char videoFilename[FILENAME_MAX];
    bool printVideoFilename;
    double plotT0;
    double plotT1;
    Image *frames;
    int frameWidth;
    int frameHeight;
    int framesPerSecond;
    int nVideoFrames;

    volatile bool processorRunning;
    bool keepFrames;

} ProcessorState;

#endif // _STATE_H
