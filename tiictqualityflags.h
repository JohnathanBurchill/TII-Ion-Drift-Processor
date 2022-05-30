/*

    TII Cross-Track Ion Drift Processor: tiictqualityflags.h

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

#ifndef TIICTQUALITYFLAGS_H
#define TIICTQUALITYFLAGS_H

#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include <cdf.h>

#define SOFTWARE_VERSION "1.0"
#define SOFTWARE_VERSION_STRING "tiictqualityflags 2022-05-30"

#define NUM_DATA_VARIABLES 3

#define STATUS_INTERVAL_FRACTION 0.05

void closeCdf(CDFid id);

void loadCrossTrackData(const char * filename, uint8_t **dataBuffers, long *numberOfRecords, bool *fourByteCalFlag);


#define ADDR(n, m) (((float*)dataBuffers[(n)]+(timeIndex + m)))
#define MEAS(n, m) ((float)(*(ADDR(n, m))))
#define TIME() ((double)*((double*)dataBuffers[0]+(timeIndex)))
#define QDLAT() (MEAS(1, 0))
#define FLAG() ((uint16_t)*((uint16_t*)dataBuffers[2]+(timeIndex)))


uint8_t getMinorVersion(const char *filename);


#endif // TIICTQUALITYFLAGS_H
