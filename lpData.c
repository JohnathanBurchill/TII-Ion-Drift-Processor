/*

    TII Cross-Track Ion Drift Processor: lpData.c

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

#include "lpData.h"

#include "indexing.h"
#include "utilities.h"

#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <cdf.h>
#include <fts.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

extern char infoHeader[50];

int getLpData(const char *lpDir, const char *satellite, const int year, const int month, const int day, uint8_t **dataBuffers, long nRecs, float **lpPhiScHighGain, float **lpPhiScLowGain, float **lpPhiSc, size_t *nLpRecs)
{

    int status = 0;

    // Get data from previous day, requested date, and next day    
    int y = year;
    int m = month;
    int d = day;
    
    struct tm date = {0};
    date.tm_year = y - 1900;
    date.tm_mon = m - 1;
    date.tm_mday = d - 1;
    char lpFile[FILENAME_MAX];
    int res = 0;

    double *lpTimes2Hz = NULL;
    double *lpVsHg = NULL;
    double *lpVsLg = NULL;
    double *lpVs = NULL;

    for (int i = 0; i < 3; i++)
    {
        timegm(&date);
        bzero(lpFile, FILENAME_MAX);
        res = getLpInputFilename(satellite[0], date.tm_year+1900, date.tm_mon+1, date.tm_mday, lpDir, lpFile);
        if (res != 0)
        {
            fprintf(stderr, "%sNo LP data for %4d%02d%02d\n", infoHeader, date.tm_year+1900, date.tm_mon+1, date.tm_mday);
            continue;
        }
        fprintf(stderr, "%sLoading LP data from %s\n", infoHeader, lpFile);
        
        loadLpInputs(lpFile, &lpTimes2Hz, &lpVsHg, &lpVsLg, &lpVs, nLpRecs);

        date.tm_mday = date.tm_mday + 1;
        
    }

    // interpolate LP data to TII times
    *lpPhiScHighGain = malloc(sizeof(float) * nRecs);
    *lpPhiScLowGain = malloc(sizeof(float) * nRecs);
    *lpPhiSc = malloc(sizeof(float) * nRecs);
    if (*lpPhiScHighGain == NULL || *lpPhiScLowGain == NULL || *lpPhiSc == NULL)
    {
        return -1;
    }
    bzero(*lpPhiScHighGain, sizeof(float) * nRecs);
    bzero(*lpPhiScLowGain, sizeof(float) * nRecs);
    bzero(*lpPhiSc, sizeof(float) * nRecs);

    double *tiiTime = (double*)dataBuffers[0];

    interpolate(lpTimes2Hz, lpVsHg, *nLpRecs, tiiTime, nRecs, *lpPhiScHighGain);
    interpolate(lpTimes2Hz, lpVsLg, *nLpRecs, tiiTime, nRecs, *lpPhiScLowGain);
    interpolate(lpTimes2Hz, lpVsLg, *nLpRecs, tiiTime, nRecs, *lpPhiSc);

    return status;

}

// Adapted from SLIDEM
int getLpInputFilename(const char satelliteLetter, long year, long month, long day, const char *path, char *filename)
{
	char *searchPath[2] = {NULL, NULL};
    searchPath[0] = (char *)path;

	FTS * fts = fts_open(searchPath, FTS_PHYSICAL | FTS_NOCHDIR, NULL);	
	if (fts == NULL)
	{
		printf("Could not open directory %s for reading.", path);
		return -1;
	}
	FTSENT * f = fts_read(fts);

    bool gotHmFile = false;
    long fileYear;
    long fileMonth;
    long fileDay;
    long lastVersion = -1;
    long fileVersion;
	while(f != NULL)
	{
		if ((strlen(f->fts_name) == 59 || strlen(f->fts_name) == 70) && *(f->fts_name+11) == satelliteLetter && strncmp(f->fts_name+13, "LP_HM", 5) == 0)
		{
            char fyear[5] = { 0 };
            char fmonth[3] = { 0 };
            char fday[3] = { 0 };
            char version[5] = { 0 };
            strncpy(fyear, f->fts_name + 19, 4);
            fileYear = atol(fyear);
            strncpy(fmonth, f->fts_name + 23, 2);
            fileMonth = atol(fmonth);
            strncpy(fday, f->fts_name + 25, 2);
            fileDay = atol(fday);
            strncpy(version, f->fts_name + 51, 4);
            fileVersion = atol(version);
            if (fileYear == year && fileMonth == month && fileDay == day && fileVersion > lastVersion)
            {
                lastVersion = fileVersion;
                sprintf(filename, "%s", f->fts_path);
                gotHmFile = true;
            }
		}
		f = fts_read(fts);
	}

	fts_close(fts);

    if (gotHmFile)
    {
        return 0;
    }
    else
    {
        return -1;
    }

}



// Adapted from SLIDEM
void loadLpInputs(const char *cdfFile, double **lpTime, double **lpPhiScHighGain, double **lpPhiScLowGain, double **lpPhiSc, long *numberOfRecords)
{
    // Open the CDF file with validation
    CDFsetValidate(VALIDATEFILEoff);
    CDFid cdfId;
    CDFstatus status;
    // Attributes
    long attrN;
    long entryN;
    char attrName[CDF_ATTR_NAME_LEN256+1];
    long attrScope, maxEntry;

    // Check CDF info
    long decoding, encoding, majority, maxrRec, numrVars, maxzRec, numzVars, numAttrs, format;

    long numBytesToAdd, numVarBytes, numValues, totalBytes;
    long varNum, dataType, numElems, numRecs, numDims, recVary;
    long dimSizes[CDF_MAX_DIMS], dimVarys[CDF_MAX_DIMS];
    CDFdata data;

    status = CDFopenCDF(cdfFile, &cdfId);
    if (status != CDF_OK) 
    {
        printErrorMessage(status);
        fprintf(stdout, "%s Could not open CDF file. Skipping this date.\n", infoHeader);
        return;
    }

    status = CDFgetFormat(cdfId, &format);
    status = CDFgetDecoding(cdfId, &decoding);
    status = CDFinquireCDF(cdfId, &numDims, dimSizes, &encoding, &majority, &maxrRec, &numrVars, &maxzRec, &numzVars, &numAttrs);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        fprintf(stdout, "\n%s Problem with CDF file. Skipping this date.\n", infoHeader);
        closeCdf(cdfId);
        return;
    }
    int nVariables = 4;
    char * variables[4] = {"Timestamp", "Vs_hgn", "Vs_lgn", "U_SC"};

    for (uint8_t i = 0; i<nVariables; i++)
    {
        status = CDFconfirmzVarExistence(cdfId, variables[i]);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            fprintf(stdout, "\n%s Error reading variable %s from CDF file. Skipping this date.\n", infoHeader, variables[i]);
            closeCdf(cdfId);
                return;
        }
        else
        {
            // fprintf(stdout, "%s OK\n", infoHeader);
        }
    }
    
    for (uint8_t i = 0; i < nVariables; i++)
    {
        varNum = CDFgetVarNum(cdfId, variables[i]);
        status = CDFreadzVarAllByVarID(cdfId, varNum, &numRecs, &dataType, &numElems, &numDims, dimSizes, &recVary, dimVarys, &data);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            fprintf(stdout, "%s Error loading data for %s. Skipping this date.\n", infoHeader, variables[i]);
            closeCdf(cdfId);
            CDFdataFree(data);
            return;
        }
        // Calculate new size of memory to allocate
        status = CDFgetDataTypeSize(dataType, &numVarBytes);
        numValues = 1;
        for (uint8_t j = 0; j < numDims; j++)
        {
            numValues *= dimSizes[j];
        }
        numBytesToAdd = numValues * numRecs * numVarBytes;
        totalBytes = sizeof(double)* (*numberOfRecords) + numBytesToAdd;
        switch(i)
        {
            case 0: // time
                *lpTime = (double*) realloc(*lpTime, totalBytes);
                memcpy(*lpTime + (*numberOfRecords), data, numBytesToAdd);
                break;
            case 1: // Vs_hgn
                *lpPhiScHighGain = (double*) realloc(*lpPhiScHighGain, totalBytes);
                memcpy(*lpPhiScHighGain + (*numberOfRecords), data, numBytesToAdd);
                break;
            case 2: // Vs_lgn
                *lpPhiScLowGain = (double*) realloc(*lpPhiScLowGain, totalBytes);
                memcpy(*lpPhiScLowGain + (*numberOfRecords), data, numBytesToAdd);
                break;
            case 3: // U_SC
                *lpPhiSc = (double*) realloc(*lpPhiSc, totalBytes);
                memcpy(*lpPhiSc + (*numberOfRecords), data, numBytesToAdd);
                break;
            default:
                return;
        }
        CDFdataFree(data);
    }
    // close CDF
    closeCdf(cdfId);

    // Update number of records found and memory allocated
    *numberOfRecords += numRecs;
    // *totalMemoryAllocated = fpMemorySize;

}


// Copied and modified from TRACIS interpolate.c
void interpolate(double *times, double *values, size_t nVals, double *requestedTimes, long nRequestedValues, float *newValues)
{
    
    size_t lastIndex = 0;
    double thisTime = 0;
    double t1 = 0, t2 = 0, dt = 0;
    double fraction = 0;
    double v1 = 0, v2 = 0, dv = 0;
    double x = 0, y = 0, z = 0;
    double rad = 0;
    double lat = 0;
    double lon = 0;
    for (size_t i = 0; i < nRequestedValues; i++)
    {
        thisTime = requestedTimes[i];
        // 
        while (times[lastIndex] <= thisTime && lastIndex < nVals) 
        {
            lastIndex++;
        }
        // Extrapolate earlier or later, or the times are the same
        if (times[lastIndex] >= thisTime || lastIndex == nVals - 1)
        {
            newValues[i] = (float)values[lastIndex];
        }
        // Interpolate
        else
        {
            t1 = times[lastIndex];
            t2 = times[lastIndex+1];
            // Assumes t2 > t1
            dt = thisTime - t1;
            fraction = dt / (t2 - t1);
            newValues[i] = (float)(values[lastIndex] + (values[lastIndex+1] - values[lastIndex]) * fraction);
        }
    }

    return;

}
