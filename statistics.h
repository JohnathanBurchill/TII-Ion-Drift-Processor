/*

    TII Cross-Track Ion Drift Processor: statistics.h

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

#ifndef _TIICT_STATISTICS_H
#define _TIICT_STATISTICS_H

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define BIN_STORAGE_BLOCK_SIZE 1024 // Number of elements to increment bin storage by

#define NSTATISTICS 6

enum STATISTICS_STATUS
{
    STATISTICS_OK = 0,
    STATISTICS_MEM,
    STATISTICS_NO_DATA,
    STATISTICS_POINTER,
    STATISTICS_UNSUPPORTED_STATISTIC    
};

int allocateBinStorage(float ***bins, size_t **binSizes, size_t **binMaxSizes, size_t nMLTs, size_t nQDLats, size_t sizePerBin);

int adjustBinStorage(float **bins, size_t *binMaxSizes, int mltQdLatIndex, long numberOfElementsToAdd);

void freeBinStorage(float **bins, size_t *binSizes, size_t *binMaxSizes, int nMLTs, int nQDLats);

void printAvailableStatistics(FILE *dest);
bool validStatistic(const char *statistic);

int calculateStatistic(const char *statistic, float **bins, size_t *binSizes, size_t mltQdIndex, void *returnValue);


#endif // _TIICT_STATISTICS_H
