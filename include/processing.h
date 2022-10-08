/*

    TII Cross-Track Ion Drift Processor: processing.h

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

#ifndef PROCESSING_H
#define PROCESSING_H

#include "state.h"

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>

enum FITINFO_BIT_MASKS {
    FiTINFO_OK = 0,
    FITINFO_OFFSET_NOT_REMOVED = (1 << 0),
    FITINFO_INCOMPLETE_REGION = (1 << 1),
    FITINFO_GSL_FIT_ERROR = (1 <<2 ),
    FITINFO_MAD_EXCEEDED = (1 << 3),
    FITINFO_DRIFT_MAGNITUDE_EXCEEDED = (1 << 4),
    FITINFO_TII_IMAGE_QUALITY = (1 << 5)
};

int initQualityData(ProcessorState *state);
int calibrateFlows(ProcessorState *state);

int removeOffsetsAndSetFlags(ProcessorState *state, bool setFlags);
int removeOffsetsAndSetFlagsForInterval(ProcessorState *state, uint8_t interval, bool setFlags);

void updateDataQualityFlags(ProcessorState *state, uint8_t sensorIndex, uint8_t regionNumber, float driftValue, float mad, long timeIndex);

float madThreshold(char satellite, int sensorIndex);

int initFields(ProcessorState *state);
int calculateFields(ProcessorState *state);

int runProcessor(int argc, char *argv[]);
int initProcessor(int argc, char **argv, ProcessorState *state);
int parseArguments(int argc, char **argv, ProcessorState *state);
void initHeader(ProcessorState *state);

int checkResult(int status, ProcessorState *state);
int shutdown(int status, ProcessorState *state);

#endif // PROCESSING_H
