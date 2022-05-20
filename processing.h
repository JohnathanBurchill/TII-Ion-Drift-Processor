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

#include <stdint.h>
#include <stdio.h>

enum FITINFO_BIT_MASKS {
    FITINFO_OFFSET_NOT_REMOVED = (1 << 0),
    FITINFO_INCOMPLETE_REGION = (1 << 1),
    FITINFO_GSL_FIT_ERROR = (1 <<2),
    FITINFO_MAD_EXCEEDED = (1 << 3),
    FITINFO_DRIFT_MAGNITUDE_EXCEEDED = (1 << 4)
};

typedef struct offset_model_fit_arguments {
    const uint8_t regionNumber;
    const char *regionName;
    const float lat1;
    const float lat2;
    const float lat3;
    const float lat4;
} offset_model_fit_arguments;


void removeOffsetsAndSetFlags(const char* satellite, offset_model_fit_arguments fitargs, long nRecs, uint8_t **dataBuffers, float *viErrors, uint16_t *flags, uint32_t *fitInfo, FILE* fitFile);

void updateDataQualityFlags(const char *satellite, uint8_t sensorIndex, uint8_t regionNumber, float driftValue, float mad, long timeIndex, uint16_t *flags, uint32_t *fitInfo);

float madThreshold(char satellite, int sensorIndex);


#endif // PROCESSING_H
