/*

    TII Cross-Track Ion Drift Processor: utilities.c

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

#include <stdio.h>
#include <cdf.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/errno.h>
#include <sys/stat.h>

#include "utilities.h"

// Prefix for all fprintf messages
extern char infoHeader[50];

#define DIRECTORY_PERMISSIONS 0755

// Generates the filename for exported CDF file, with full path
void constructExportFileName(const char *dataset, double startTime, double stopTime, const char *exportDir, const char *exportVersion, const char *satellite, char *cdfFileName) 
{
 
    long startyear, startmonth, startday, starthour, startminute, startsecond, startmillisecond;
    long stopyear, stopmonth, stopday, stophour, stopminute, stopsecond, stopmillisecond;

    EPOCHbreakdown(startTime, &startyear, &startmonth, &startday, &starthour, &startminute, &startsecond, &startmillisecond);
    EPOCHbreakdown(stopTime, &stopyear, &stopmonth, &stopday, &stophour, &stopminute, &stopsecond, &stopmillisecond);

    sprintf(cdfFileName, "%s/%s/%s/SW_EXPT_EFI%s_%s_%04ld%02ld%02ldT%02ld%02ld%02ld_%04ld%02ld%02ldT%02ld%02ld%02ld_%s", exportDir, exportVersion, dataset, satellite, dataset, startyear, startmonth, startday, starthour, startminute, startsecond, stopyear, stopmonth, stopday, stophour, stopminute, stopsecond, exportVersion);

}

void printErrorMessage(CDFstatus status)
{
    char errorMessage[CDF_STATUSTEXT_LEN + 1];
    CDFgetStatusText(status, errorMessage);
    fprintf(stdout, "%s%s\n", infoHeader, errorMessage);
}

void closeCdf(CDFid id)
{
    CDFstatus status;
    status = CDFcloseCDF(id);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
    }

}

int makeSureDirExists(const char *exportDir, const char *exportVersion, const char *subdir)
{
    int status = 0;
    char outDir[FILENAME_MAX] = {0};
    sprintf(outDir, "%s/%s", exportDir, exportVersion);

    errno = 0;
    status = access(outDir, F_OK);
    if (status != 0 && errno == ENOENT)
    {
        // Try to create the directory
        errno = 0;
        fprintf(stderr, "%sCreating %s\n", infoHeader, outDir);
        status = mkdir(outDir, DIRECTORY_PERMISSIONS);
    }

    if (status != 0)
        return status;


    sprintf(outDir, "%s/%s/%s", exportDir, exportVersion, subdir);
    status = access(outDir, F_OK);
    if (status != 0 && errno == ENOENT)
    {
        // Try to create the directory
        errno = 0;
        fprintf(stderr, "%sCreating %s\n", infoHeader, outDir);
        status = mkdir(outDir, DIRECTORY_PERMISSIONS);
    }

    return status;


}