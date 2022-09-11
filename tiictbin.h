/*

    TII Cross-Track Ion Drift Processor: tiictbin.h

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

#ifndef _TIICTBIN_H
#define _TIICTBIN_H

#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include <cdf.h>

#define SOFTWARE_VERSION "1.0"
#define SOFTWARE_VERSION_STRING "tiictbin 2022-09-10"

#define NUM_DATA_VARIABLES 7

void closeCdf(CDFid id);

CDFstatus loadCrossTrackData(const char *filename, uint8_t **dataBuffers, long *numberOfRecords, bool *fourByteCalFlag, const char *parameterName);

#define ADDR(n, m, d) (((float*)dataBuffers[(n)]+(d*timeIndex + m)))
#define MEAS(n, m, d) ((float)(*(ADDR(n, m, d))))
#define TIME() ((double)*((double*)dataBuffers[0]+(timeIndex)))
#define MLT() (MEAS(1, 0, 1))
#define QDLAT() (MEAS(2, 0, 1))
#define FLAG() ((uint16_t)*((uint16_t*)dataBuffers[3]+(timeIndex)))
#define CALFLAG() ((uint32_t)*((uint32_t*)dataBuffers[4]+(timeIndex)))

#define VSATN() (MEAS(5, 0, 1))

#define PARAMETER() (MEAS(6, 0, 1))

#define PARAMETERVEC2X() (MEAS(6, 0, 2))
#define PARAMETERVEC2Y() (MEAS(6, 1, 2))

#define PARAMETERVEC3X() (MEAS(6, 0, 3))
#define PARAMETERVEC3Y() (MEAS(6, 1, 3))
#define PARAMETERVEC3Z() (MEAS(6, 2, 3))


uint8_t getMinorVersion(const char *filename);

void usage(char *name);

#endif // _TIICTBIN_H
