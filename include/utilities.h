/*

    TII Cross-Track Ion Drift Processor: include/utilities.h

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

#ifndef UTILITIES_H
#define UTILITIES_H

#include <cdf.h>

// Constructs a full path to the export CDF file in the argument constructExportFileName.
void constructExportFileName(const char *dataset, double startTime, double stopTime, const char *exportDir, const char *exportVersion, const char *satellite, char *cdfFileName);

// Prints an error message from the CDFstatus
void printErrorMessageToFile(FILE *file, CDFstatus status);
void printErrorMessage(CDFstatus status);

void closeCdf(CDFid id);

int makeSureDirExists(const char *exportDir, const char *exportVersion, const char *subdir);

#endif // UTILITIES_H
