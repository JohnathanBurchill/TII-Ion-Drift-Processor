/*

    TII Cross-Track Ion Drift Processor: lpData.h

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

#ifndef _LPDATA_H
#define _LPDATA_H

#include <stdint.h>
#include <stddef.h>

int getLpData(const char *lpDir, const char *satellite, const int year, const int month, const int day, uint8_t **dataBuffers, long nRecs, float **lpPhiScHighGain, float **lpPhiScLowGain, float **lpPhiSc, size_t *nLpRecs);

int getLpInputFilename(const char satelliteLetter, long year, long month, long day, const char *path, char *filename);

void loadLpInputs(const char *cdfFile, double **lpTime, double **lpPhiScHighGain, double **lpPhiScLowGain, double **lpPhiSc, long *numberOfRecords);

void interpolate(double *times, double *values, size_t nVals, double *requestedTimes, long nRequestedValues, float *newValues);

#endif // _LPDATA_H
