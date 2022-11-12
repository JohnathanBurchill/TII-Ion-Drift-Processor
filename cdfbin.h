/*

    TII Cross-Track Ion Drift Processor: cdfbin.h

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

#ifndef _CDFBIN_H
#define _CDFBIN_H

#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include <cdf.h>

#define SOFTWARE_VERSION "1.0"
#define SOFTWARE_VERSION_STRING "tiictbin 2022-09-10"

#define NUM_DATA_VARIABLES 5
#define TIME_INDEX 0
#define QDLAT_INDEX 1
#define MLT_INDEX 2
#define FLAG_INDEX 3
#define PARAMETER_INDEX 4

CDFstatus loadCdfData(const char *filename, uint8_t **dataBuffers, long *numberOfRecords, const char *parameterName, const char *flagVarName, long *variableTypes);

#define ADDR(n, m, d) (((float*)dataBuffers[(n)]+(d*timeIndex + m)))
#define MEAS(n, m, d) ((float)(*(ADDR(n, m, d))))
#define TIME() ((double)*((double*)dataBuffers[0]+(timeIndex)))
#define MLT() (MEAS(1, 0, 1))
#define QDLAT() (MEAS(2, 0, 1))

#define PARAMETER() (MEAS(4, 0, 1))

#define PARAMETERVEC2X() (MEAS(4, 0, 2))
#define PARAMETERVEC2Y() (MEAS(4, 1, 2))

#define PARAMETERVEC3X() (MEAS(4, 0, 3))
#define PARAMETERVEC3Y() (MEAS(4, 1, 3))
#define PARAMETERVEC3Z() (MEAS(4, 2, 3))

#define ADDR8(n, m, d) (((double*)dataBuffers[(n)]+(d*timeIndex + m)))
#define MEAS8(n, m, d) ((*(ADDR8(n, m, d))))
#define MLT8() (MEAS8(1, 0, 1))
#define QDLAT8() (MEAS8(2, 0, 1))

#define PARAMETER8() (MEAS8(4, 0, 1))

#define PARAMETERVEC2X8() (MEAS8(4, 0, 2))
#define PARAMETERVEC2Y8() (MEAS8(4, 1, 2))

#define PARAMETERVEC3X8() (MEAS8(4, 0, 3))
#define PARAMETERVEC3Y8() (MEAS8(4, 1, 3))
#define PARAMETERVEC3Z8() (MEAS8(4, 2, 3))



uint8_t getMinorVersion(const char *filename);

void usage(char *name);

#endif // _CDFBIN_H
