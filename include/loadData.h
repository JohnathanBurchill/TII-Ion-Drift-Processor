/*

    TII Cross-Track Ion Drift Processor: loadData.h

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

#ifndef _LOADDATA_H
#define _LOADDATA_H

#include "state.h"

#include <stdint.h>
#include <stddef.h>

typedef enum DayType {
    PREVIOUS_DAY = -1,
    REQUESTED_DAY = 0,
    NEXT_DAY = 1
} DayType;

int getLpInputFilename(const char satelliteLetter, long year, long month, long day, const char *path, char *filename);
int loadLpCalData(ProcessorState *state);
int getLpData(ProcessorState *state);
int loadLpInputs(const char *cdfFile, double **lpTime, double **lpPhiScHighGain, double **lpPhiScLowGain, double **lpPhiSc, double **ionDensity, uint32_t **flagbits, long *numberOfRecords);

int loadTiiData(ProcessorState *state);
void loadTiiCalDataFromDate(const DayType dayType, ProcessorState *state);
void loadTracisDataFromDate(const DayType dayType, ProcessorState *state);

int loadData(ProcessorState *state);

void loadCdf(const DayType dayType, const char * filename, const char *variables[], int nVars, ProcessorState *state, float sampleRate, uint8_t **dataBuffers, size_t *numberOfRecords, size_t *memoryAllocated);

void setCalibrationFileName(ProcessorState *state, DayType dayType, int year, int month, int day);
void setTracisFileName(ProcessorState *state, DayType dayType, int year, int month, int day);

int checkCalDataAvailability(ProcessorState *state, DayType dayType);

#endif // _LOADDATA_H
