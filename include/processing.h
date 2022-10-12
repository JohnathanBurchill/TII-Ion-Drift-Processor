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

enum CALIBRATION_FLAGS_BIT_MASKS {
    CALIBRATION_FLAGS_OK = 0,
    CALIBRATION_FLAGS_OFFSET_NOT_REMOVED = (1 << 0),
    CALIBRATION_FLAGS_INCOMPLETE_REGION = (1 << 1),
    CALIBRATION_FLAGS_GSL_FIT_ERROR = (1 <<2 ),
    CALIBRATION_FLAGS_MAD_EXCEEDED = (1 << 3),
    CALIBRATION_FLAGS_DRIFT_MAGNITUDE_EXCEEDED = (1 << 4),
    CALIBRATION_FLAGS_TII_IMAGE_QUALITY = (1 << 5),
    CALIBRATION_FLAGS_ION_DENSITY = (1 << 6),
    CALIBRATION_FLAGS_U_SC_QUALITY = (1 << 7)
};

int initQualityData(ProcessorState *state);
int calibrateFlows(ProcessorState *state);

int removeOffsetsAndSetFlags(ProcessorState *state, bool setFlags);
int removeOffsetsAndSetFlagsForInterval(ProcessorState *state, uint8_t interval, bool setFlags);

void setCalFlags(ProcessorState *state, long t1, long t2, int flagIndex, int flagValue);
void unsetCalFlags(ProcessorState *state, long t1, long t2, int flagIndex, int flagValue);

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
