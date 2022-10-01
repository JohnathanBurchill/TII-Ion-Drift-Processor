/*

    TII Cross-Track Ion Drift Processor: export.h

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

#ifndef _TIICT_EXPORT_H
#define _TIICT_EXPORT_H

#include "state.h"

typedef struct varAttr {
    char * name;
    char * type;
    char * units;
    char * desc;
    double validMin;
    double validMax;
} varAttr;

int exportCdfs(ProcessorState *state);

CDFstatus createVarFrom1DVar(CDFid id, char *name, long dataType, long startIndex, long stopIndex, void *buffer);
CDFstatus createVarFrom2DVar(CDFid id, char *name, long dataType, long startIndex, long stopIndex, void *buffer, uint8_t index, uint8_t dimSize);


CDFstatus addgEntry(CDFid id, long attrNum, long entryNum, const char *entry);

CDFstatus addVariableAttributes(CDFid id, varAttr attr);

void addAttributes(CDFid id, const char *calVersion, const char *satellite, const char *version, double minTime, double maxTime);

int exportTCT16Cdfs(ProcessorState *state, double startTime, double stopTime, long startIndex, long stopIndex);
int exportTCT02Cdfs(ProcessorState *state, double startTime, double stopTime, long startIndex, long stopIndex);

int zipCdfFile(char *cdfFileName);

int initDirectories(Arguments *args);

void initLogFile(char *fitLogFileName, Arguments *a);

int initFitFiles(ProcessorState *state);

#endif // _TIICT_EXPORT_H
