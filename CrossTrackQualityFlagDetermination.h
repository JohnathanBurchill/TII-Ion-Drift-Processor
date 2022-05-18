/*

    TII Cross-Track Ion Drift Processor: CrossTrackQualityFlagDetermination.h

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

#ifndef CROSSTRACKVALIDATION
#define CROSSTRACKVALIDATION

#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include <cdf.h>

#define SOFTWARE_VERSION_STRING "CrossTrackValidation 2020-05-11"

#define NUM_DATA_VARIABLES 6

void closeCdf(CDFid id);

void loadCrossTrackData(const char * filename, uint8_t **dataBuffers, long *numberOfRecords, bool *fourByteCalFlag);


#define ADDR(n, m) (((float*)dataBuffers[(n)]+(timeIndex + m)))
#define MEAS(n, m) ((float)(*(ADDR(n, m))))
#define TIME() ((double)*((double*)dataBuffers[0]+(timeIndex)))
#define QDLAT() (MEAS(1, 0))
#define VIY() (MEAS(2, 0))
#define VIYERROR() (MEAS(3, 0))
#define FLAG() ((uint16_t)*((uint16_t*)dataBuffers[4]+(timeIndex)))
#define CALFLAG() ((uint32_t)*((uint32_t*)dataBuffers[5]+(timeIndex)))


uint8_t getMinorVersion(const char *filename);


#endif // CROSSTRACKVALIDATION
