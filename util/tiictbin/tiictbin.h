/*

    TIICT Processor: util/tiictbin/tiictbin.h

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

// Based on SLIDEM processor util slidembin

#ifndef _TIICTBIN_H
#define _TIICTBIN_H

#include "statistics.h"

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>

#include <fts.h>
#include <cdf.h>

#define SOFTWARE_VERSION "1.0"
#define SOFTWARE_VERSION_STRING "slidembin 2023-03-21"

#define NUM_DATA_VARIABLES 7

typedef struct processingParameters
{
    int nOptions;
    char satelliteLetter;

    bool verbose;
    bool showFileProgress;

    char *cdfDirectory;
    char *inputFile;

    char *parameter;
    char *statistic;

    long nRecords;
    double *time;
    float *qdlat;
    float *mlt;
    float *values;
    uint16_t *flags;

    BinningState binningState;

    char *firstTimeString;
    char *lastTimeString;
    double firstTime;
    double lastTime;
    bool processAllSpaceSeries;

    long nFiles;

    int64_t flagIgnoreMask;
    uint64_t positiveFlagMask;
    bool flagMaskIsAnd;
    bool flagRaisedIsGood;

} ProcessingParameters;

void usage(char *name);
void about(void);
void parseCommandLine(ProcessingParameters *params, int argc, char *argv[]);
bool fileMatch(FTSENT *e, ProcessingParameters *params);
int processFile(ProcessingParameters *params);
int loadTiictData(ProcessingParameters *params);
CDFstatus loadCdfVariable(CDFid cdfId, char *variable, void **mem, long *nRecords);

void printQualityFlagTable(void);

#endif // _TIICTBIN_H
