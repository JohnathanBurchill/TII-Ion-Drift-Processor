/*

    TII Cross-Track Ion Drift Processor: loadData.h

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

#ifndef _LOADDATA_H
#define _LOADDATA_H

#include "state.h"

typedef enum CalibrationFileType {
    CAL_FILE_TII = 0,
    CAL_FILE_LP = 1,
} CalibrationFileType_enum;

typedef enum DayType {
    PREVIOUS_DAY = -1,
    REQUESTED_DAY = 0,
    NEXT_DAY = 1
} DayType;

enum LpPotentialSource {
    LP_POTENTIAL_NONE = 0,
    LP_POTENTIAL_HIGHGAIN = 1,
    LP_POTENTIAL_LOWGAIN = 2,
    LP_POTENTIAL_U_SC = 3,
    LP_POTENTIAL_UNKNOWN = 4,
};

typedef struct CalibrationVariable {
    char *name;
    void **memoryPointer;
    long recordValueOffset;
} CalibrationVariable_t;


int getLpInputFilename(const char satelliteLetter, long year, long month, long day, const char *path, char *filename);
int loadLpCalData(ProcessorState *state);
int getLpData(ProcessorState *state);
int loadLpInputs(const char *cdfFile, double **lpTime, double **lpPhiScHighGain, double **lpPhiScLowGain, double **lpPhiSc, size_t *numberOfRecords);

int loadCalData(ProcessorState *state);
void loadCalDataFromDate(const DayType dayType, ProcessorState *state, CalibrationVariable_t *variables, int nVariables);
int loadCdfVariable(ProcessorState *state, CDFid calCdfId, CalibrationVariable_t *variable, long startRecord, long stopRecord, long *calibrationMemorySize);
void freeVariable(void *var);
void freeVariables(ProcessorVariables_t *vars);
int reallocVariable(void **var, size_t newSize);
int reallocVariables(ProcessorVariables_t *vars, size_t newSize);

void setCalibrationFileName(ProcessorState *state, int year, int month, int day, CalibrationFileType_enum fileType);

int checkCalDataAvailability(ProcessorState *state);

#endif // _LOADDATA_H
