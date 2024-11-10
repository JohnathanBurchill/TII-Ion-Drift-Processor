/*

    Science: statistics.h

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

// Adapted from TII-Ion-Drift-Processor
// https://github.com/JohnathanBurchill/TII-Ion-Drift-Processor/tree/tracis_flagging

#ifndef _STATISTICS_H
#define _STATISTICS_H

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define BIN_STORAGE_BLOCK_SIZE 10240 // Number of elements to increment bin storage by at a time

#define NSTATISTICS 7

enum STATISTICS_STATUS
{
    STATISTICS_OK = 0,
    STATISTICS_MEM,
    STATISTICS_NO_DATA,
    STATISTICS_POINTER,
    STATISTICS_UNSUPPORTED_STATISTIC    
};

enum BIN_STATUS
{
    BIN_OK = 0,
    BIN_VALUE_NOT_FINITE,
    BIN_QDLAT_OUTOFRANGE,
    BIN_MLT_OUTOFRANGE,
    BIN_SPECIFICATION,
    BIN_MEMORY
};

typedef struct binningState
{
    bool equalArea;

    double qdlatmin;
    double qdlatmax;
    double deltaqdlat;

    double mltmin;
    double mltmax;
    double deltamlt;

    int nQDLats;
    int nMLTs;
    double *nMltsVsLatitude;
    double *cumulativeMltsVsLatitude;

    double **binStorage;
    size_t *binSizes;
    size_t *binValidSizes;
    size_t *binMaxSizes;

    size_t nBins;
    long nValsRead;
    long nValsWithinBinLimits;
    long nValsBinned;

    bool flipParamWhenDawnward;
    bool flipParamWhenSunward;
    bool flipParamWhenDescending;

    double solidAngleUnit;
    double ringSolidAngle;

} BinningState;

int initBinningState(BinningState *binningState);

int allocateBinStorage(double ***bins, size_t **binSizes, size_t **binValidSizes, size_t **binMaxSizes, size_t nBins, size_t sizePerBin, bool useEqualArea);

int adjustBinStorage(double **bins, size_t *binMaxSizes, int mltQdLatIndex, long numberOfElementsToAdd);

void freeBinStorage(BinningState *binningState);

int binData(BinningState *binningState, double qdlat, double mlt, double value, bool includeValue);
void printBinningResults(BinningState *binningState, char *parameter, char *statistic);

void printAvailableStatistics(FILE *dest);
bool validStatistic(const char *statistic);

int calculateStatistic(const char *statistic, double **bins, size_t *binSizes, size_t mltQdIndex, void *returnValue);


#endif // _STATISTICS_H
