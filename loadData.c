/*

    TII Cross-Track Ion Drift Processor: loadData.c

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

#include "loadData.h"
#include "errors.h"
#include "indexing.h"
#include "utilities.h"
#include "processing.h"

#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <cdf.h>
#include <fts.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

extern char infoHeader[50];

// Adapted from SLIDEM
int getLpInputFilename(const char satelliteLetter, long year, long month, long day, const char *path, char *filename)
{
	char *searchPath[2] = {NULL, NULL};
    searchPath[0] = (char *)path;

	FTS * fts = fts_open(searchPath, FTS_PHYSICAL | FTS_NOCHDIR, NULL);	
	if (fts == NULL)
	{
		printf("Could not open directory %s for reading.", path);
		return TIICT_DIRECTORY_READ;
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
        return TIICT_OK;
    else
        return TIICT_NO_LP_HM_DATA;

}

int loadLpCalData(ProcessorState *state)
{
    // Get LP floating potentials
    state->lpTimes = NULL;
    state->lpPhiScHighGain = NULL;
    state->lpPhiScLowGain = NULL;
    state->lpPhiSc = NULL;
    state->nLpRecs = 0;

    int status = TIICT_OK;
    state->usePotentials = strcmp(state->args.exportVersion, "0401") >= 0; 
    if (state->usePotentials)
    {
        status = getLpData(state);
        if (status == TIICT_OK)
        {
            if (state->nLpRecs < LP_MIN_NUMBER_OF_POTENTIALS)
            {
                fprintf(state->processingLogFile, "%sNot enough (%lu) LP potentials imported.\n", infoHeader, state->nLpRecs);
                status =  TIICT_NO_LP_HM_DATA;
            }
        }
    }

    return status;

}

int getLpData(ProcessorState *state)
{

    // Get data from previous day, requested date, and next day    
    int y = state->args.year;
    int m = state->args.month;
    int d = state->args.day;
    
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

    int status = TIICT_OK;

    for (int i = 0; i < 3; i++)
    {
        timegm(&date);
        bzero(lpFile, FILENAME_MAX);
        res = getLpInputFilename(state->args.satellite[0], date.tm_year+1900, date.tm_mon+1, date.tm_mday, state->args.lpDir, lpFile);
        if (res != 0)
        {
            date.tm_mday = date.tm_mday + 1;
            continue;
        }
        
        status = loadLpInputs(lpFile, &lpTimes2Hz, &lpVsHg, &lpVsLg, &lpVs, &state->nLpRecs);
        if (status == TIICT_MEMORY)
        {
            fprintf(state->processingLogFile, "%sUnable to allocate memory for LP data.", infoHeader);
            return status;
        }

        date.tm_mday = date.tm_mday + 1;
        
    }
    if (state->nLpRecs == 0)
    {
        fprintf(state->processingLogFile, "%sNo LP records found.", infoHeader);
        return TIICT_NO_LP_HM_DATA;
    }

    // interpolate LP data to TII times
    state->lpPhiScHighGain = malloc(sizeof(float) * state->nRecs);
    state->lpPhiScLowGain = malloc(sizeof(float) * state->nRecs);
    state->lpPhiSc = malloc(sizeof(float) * state->nRecs);
    if (state->lpPhiScHighGain == NULL || state->lpPhiScLowGain == NULL || state->lpPhiSc == NULL)
    {
        fprintf(state->processingLogFile, "%sUnable to allocate memory for LP interpolation.", infoHeader);
        return TIICT_MEMORY;
    }
    bzero(state->lpPhiScHighGain, sizeof(float) * state->nRecs);
    bzero(state->lpPhiScLowGain, sizeof(float) * state->nRecs);
    bzero(state->lpPhiSc, sizeof(float) * state->nRecs);

    double *tiiTime = (double*)state->dataBuffers[0];

    interpolate(lpTimes2Hz, lpVsHg, state->nLpRecs, tiiTime, state->nRecs, state->lpPhiScHighGain);
    interpolate(lpTimes2Hz, lpVsLg, state->nLpRecs, tiiTime, state->nRecs, state->lpPhiScLowGain);
    interpolate(lpTimes2Hz, lpVs, state->nLpRecs, tiiTime, state->nRecs, state->lpPhiSc);

    free(lpTimes2Hz);
    free(lpVsHg);
    free(lpVsLg);
    free(lpVs);

    return TIICT_OK;

}

// Adapted from SLIDEM
int loadLpInputs(const char *cdfFile, double **lpTime, double **lpPhiScHighGain, double **lpPhiScLowGain, double **lpPhiSc, long *numberOfRecords)
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

    void *newMem = NULL;

    status = CDFopenCDF(cdfFile, &cdfId);
    if (status != CDF_OK) 
    {
        return TIICT_CDF_READ;
    }

    status = CDFgetFormat(cdfId, &format);
    status = CDFgetDecoding(cdfId, &decoding);
    status = CDFinquireCDF(cdfId, &numDims, dimSizes, &encoding, &majority, &maxrRec, &numrVars, &maxzRec, &numzVars, &numAttrs);
    if (status != CDF_OK)
    {
        closeCdf(cdfId);
        return TIICT_CDF_READ;
    }
    int nVariables = 4;
    char * variables[4] = {"Timestamp", "Vs_hgn", "Vs_lgn", "U_SC"};

    for (uint8_t i = 0; i<nVariables; i++)
    {
        status = CDFconfirmzVarExistence(cdfId, variables[i]);
        if (status != CDF_OK)
        {
            closeCdf(cdfId);
            return TIICT_CDF_READ;
        }
    }
    
    for (uint8_t i = 0; i < nVariables; i++)
    {
        varNum = CDFgetVarNum(cdfId, variables[i]);
        status = CDFreadzVarAllByVarID(cdfId, varNum, &numRecs, &dataType, &numElems, &numDims, dimSizes, &recVary, dimVarys, &data);
        if (status != CDF_OK)
        {
            closeCdf(cdfId);
            CDFdataFree(data);
            return TIICT_CDF_READ;
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
                newMem = realloc(*lpTime, totalBytes);
                if (newMem == NULL)
                    return TIICT_MEMORY;
                *lpTime = (double*) newMem;
                memcpy(*lpTime + (*numberOfRecords), data, numBytesToAdd);
                break;
            case 1: // Vs_hgn
                newMem = (double*) realloc(*lpPhiScHighGain, totalBytes);
                if (newMem == NULL)
                    return TIICT_MEMORY;
                *lpPhiScHighGain = (double*) newMem;
                memcpy(*lpPhiScHighGain + (*numberOfRecords), data, numBytesToAdd);
                break;
            case 2: // Vs_lgn
                newMem = (double*) realloc(*lpPhiScLowGain, totalBytes);
                if (newMem == NULL)
                    return TIICT_MEMORY;
                *lpPhiScLowGain = (double*) newMem;
                memcpy(*lpPhiScLowGain + (*numberOfRecords), data, numBytesToAdd);
                break;
            case 3: // U_SC
                newMem = (double*) realloc(*lpPhiSc, totalBytes);
                if (newMem == NULL)
                    return TIICT_MEMORY;
                *lpPhiSc = (double*) newMem;
                memcpy(*lpPhiSc + (*numberOfRecords), data, numBytesToAdd);
                break;
            default:
                return TIICT_CDF_READ;
        }
        CDFdataFree(data);
    }
    // close CDF
    closeCdf(cdfId);

    // Update number of records found
    *numberOfRecords += numRecs;

    return TIICT_OK;

}

int loadTiiData(ProcessorState *state)
{
    // Store processing date
    int year = state->args.year;
    int month = state->args.month;
    int day = state->args.day;

    // Get data for prior, requested, and following days
    for (int8_t i = -1; i < 2; i++)
    {
        state->args.day = day + i;
        loadTiiCalDataFromDate(i, state);
        loadTracisDataFromDate(i, state);
    }
    // Reset processing date
    state->args.year = year;
    state->args.month = month;
    state->args.day = day;

    if (state->nRecs < 16*SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING)
    {
        fprintf(state->processingLogFile, "%sNot enough calibration records\n", infoHeader);
        return TIICT_NOT_ENOUGH_CALIBRATION_RECORDS;
    }

    // During commissioning there were science operations at image rates of 1/256.
    // Ensure we have at least two TRACIS records
    if (state->nTracisRecs < (long) ceil(1.0/256. * SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING + 2))
    {
        fprintf(state->processingLogFile, "%sNot enough TRACIS records\n", infoHeader);
        return TIICT_NOT_ENOUGH_TRACIS_RECORDS;
    }

    return TIICT_OK;
}

void loadTiiCalDataFromDate(const DayType dayType, ProcessorState *state)
{
    // c time manipulation: see  https://fresh2refresh.com/c-programming/c-time-related-functions/
    struct tm timestructure;
    time_t date;

    // day or month can be outside nominal ranges due to date arithmetic, so make them canonical
    timestructure.tm_year = state->args.year - 1900;
    timestructure.tm_mon = state->args.month - 1;
    timestructure.tm_mday = state->args.day;
    timestructure.tm_hour = 0;
    timestructure.tm_min = 0;
    timestructure.tm_sec = 0;
    timestructure.tm_isdst = 0;
    mktime(&timestructure);

    state->args.year = timestructure.tm_year + 1900;
    state->args.month = timestructure.tm_mon + 1;
    state->args.day = timestructure.tm_mday;
    setCalibrationFileName(state, state->args.year, state->args.month, state->args.day);

    // Variables
    const char* variables[] = {
        "epoch",                // index 0
        "1st Moment - H",       // 1
        "1st Moment - V",       // 2
        "Det_H__vX",            // 3
        "MLT",                  // 4
        "QDLat",                // 5
        "QDLon",                // 6
        "Latitude",             // 7
        "Long Tle",             // 8
        "Radius",               // 9
        "Det H CorVx",          // 10
        "Sat_Vel_N",            // 11
        "B-North",              // 12
        "Bias_Grid_H",          // 13
        "Faceplate_Volt_Mon_H"  // 14
        };

    loadCdf(dayType, state->calibrationFileName, variables, NUM_CAL_VARIABLES, state, 16, state->dataBuffers, &state->nRecs, &state->memoryAllocated);

    return;
}

void loadTracisDataFromDate(const DayType dayType, ProcessorState *state)
{
    // c time manipulation: see  https://fresh2refresh.com/c-programming/c-time-related-functions/
    struct tm timestructure;
    time_t date;

    // day or month can be outside nominal ranges due to date arithmetic, so make them canonical
    timestructure.tm_year = state->args.year - 1900;
    timestructure.tm_mon = state->args.month - 1;
    timestructure.tm_mday = state->args.day;
    timestructure.tm_hour = 0;
    timestructure.tm_min = 0;
    timestructure.tm_sec = 0;
    timestructure.tm_isdst = 0;
    mktime(&timestructure);

    state->args.year = timestructure.tm_year + 1900;
    state->args.month = timestructure.tm_mon + 1;
    state->args.day = timestructure.tm_mday;
    setTracisFileName(state, state->args.year, state->args.month, state->args.day);

    // Variables
    const char* variables[] = {
        "Timestamp",            // index 0
        "Image_anomaly_flags_H",// 1
        "Image_anomaly_flags_H" // 2
        };

    loadCdf(dayType, state->tracisFileName, variables, NUM_TRACIS_VARIABLES, state, 1.0/256.0, state->tracisDataBuffers, &state->nTracisRecs, &state->tracisMemoryAllocated);

    return;
}

void loadCdf(const DayType dayType, const char * filename, const char *variables[], int nVars, ProcessorState *state, float sampleRate, uint8_t **dataBuffers, size_t *numberOfRecords, size_t *memoryAllocated)
{
    // Open the CDF file with validation
    CDFsetValidate(VALIDATEFILEon);
    CDFid calCdfId = 0;
    CDFstatus status = 0;
    status = CDFopenCDF(filename, &calCdfId);
    if (status != CDF_OK) 
    {
        // Not necessarily an error. For example, some dates will have not calibration data, hence no file.
        return;
    }

    // Attributes
    long attrN;
    long entryN;
    char attrName[CDF_ATTR_NAME_LEN256+1];
    long attrScope; long maxEntry; long dataType; long numElems;

    // Check CDF info
    long numDims, decoding, encoding, majority, maxrRec, numrVars, maxzRec, numzVars, numAttrs, format;
    long dimSizes[CDF_MAX_DIMS];

    status = CDFgetFormat(calCdfId, &format);
    status = CDFgetDecoding(calCdfId, &decoding);

    status = CDFinquireCDF(calCdfId, &numDims, dimSizes, &encoding, &majority, &maxrRec, &numrVars, &maxzRec, &numzVars, &numAttrs);
    if (status != CDF_OK)
    {
        fprintf(state->processingLogFile, "%sProblem reading ^%s.\n", infoHeader, filename);
        closeCdf(calCdfId);
        return;
    }
    long nRecs, calibrationMemorySize = 0;
    // First variable should be the time
    status = CDFgetzVarAllocRecords(calCdfId, CDFgetVarNum(calCdfId, (char *)variables[0]), &nRecs);
    if (status != CDF_OK)
    {
        fprintf(state->processingLogFile, "%sProblem reading %s.\n", infoHeader, filename);
        closeCdf(calCdfId);
        return;
    }
    if (nRecs < (long)ceil(sampleRate * SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING))
    {
        // Not enough to do anything useful 
        // TODO: maybe increase this threshold to require a larger number of points each day?
        CDFcloseCDF(calCdfId);
        return;
    }

    // Determine start and stop record numbers to load
    long startRecord;
    long stopRecord;
    double timeReference, recordTime;
    // get epoch variable number
    long epochNum = CDFgetVarNum(calCdfId, (char *)variables[0]);
    if (epochNum < 0)
        return;
    switch (dayType)
    {
        case PREVIOUS_DAY:
            stopRecord = nRecs-1;
            // Work backward from end of file to find first record greater than or equal to 95 minutes from end of day
            timeReference = computeEPOCH(state->args.year, state->args.month, state->args.day + 1, 0, 0, -(SECONDS_OF_BOUNDARY_DATA_REQUIRED_FOR_PROCESSING), 0);
            for (startRecord = stopRecord; startRecord >= 0; startRecord--)
            {
                status = CDFgetzVarRecordData(calCdfId, epochNum, startRecord, &recordTime);
                if (status != CDF_OK)
                {
                    fprintf(state->processingLogFile, "%sCould not read epoch record from CDF file.\n", infoHeader);
                    closeCdf(calCdfId);
                    return;
                }
                if (recordTime < timeReference)
                {
                    break;
                }
            }
            break;
        case REQUESTED_DAY:
            // Get all records
            startRecord = 0;
            stopRecord = nRecs - 1;
            break;
        case NEXT_DAY: 
            // Work forward from start of file to find last record less than or equal to 95 minutes from start of day
            startRecord = 0;
            timeReference = computeEPOCH(state->args.year, state->args.month, state->args.day, 0, 0, (SECONDS_OF_BOUNDARY_DATA_REQUIRED_FOR_PROCESSING), 0);
            for (stopRecord = startRecord; stopRecord < nRecs; stopRecord++)
            {
                status = CDFgetzVarRecordData(calCdfId, epochNum, stopRecord, &recordTime);
                if (status != CDF_OK)
                {
                    fprintf(state->processingLogFile, "%sCould not read time record from CDF file.\n", infoHeader);
                    closeCdf(calCdfId);
                    return;
                }
                if (recordTime > timeReference)
                {
                    break;
                }
            }
            break;
        default:
            fprintf(state->processingLogFile, "%sError: Day type must be one of PREVIOUS_DAY, REQUESTED_DAY, or NEXT_DAY.\n", infoHeader);
            closeCdf(calCdfId);
            return;
    }

    // Update number of records
    nRecs = stopRecord - startRecord + 1;
    if ((dayType == REQUESTED_DAY && nRecs < (long)ceil(sampleRate*SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING)) || ((dayType == PREVIOUS_DAY || dayType == NEXT_DAY) && nRecs < (long)ceil(sampleRate*SECONDS_OF_BOUNDARY_DATA_REQUIRED_FOR_PROCESSING)))
    {
        // Not enough to do anything useful 
        CDFcloseCDF(calCdfId);
        return;
    }

    for (uint8_t i = 0; i<nVars; i++)
    {
        // fprintf(state->processingLogFile, "%s%20s ", infoHeader, variables[i]);
        status = CDFconfirmzVarExistence(calCdfId, variables[i]);
        if (status != CDF_OK)
        {

            fprintf(state->processingLogFile, "%sVariable %s not found in %s.\n", infoHeader, variables[i], filename);
            closeCdf(calCdfId);
            return;
        }
        else
        {
            // fprintf(state->processingLogFile, "%s OK\n", infoHeader);
        }
    }
    
    long varNum, numValues, numVarBytes;
    long numBytesPrev, numBytesToAdd, numBytesNew;
    void *newMem = NULL;

    for (uint8_t i = 0; i < nVars; i++)
    {
        varNum = CDFgetVarNum(calCdfId, (char *)variables[i]);
        if (varNum < CDF_OK)
        {
            printErrorMessage(varNum);
            fprintf(state->processingLogFile, "%sError reading %s variable ID from %s.\n", infoHeader, variables[i], filename);
            closeCdf(calCdfId);
            return;
        }
        status = CDFgetzVarNumDims(calCdfId, varNum, &numDims);
        status = CDFgetzVarDimSizes(calCdfId, varNum, dimSizes);
        status = CDFgetzVarDataType(calCdfId, varNum, &dataType);
        // Calculate new size of memory to allocate
        status = CDFgetDataTypeSize(dataType, &numVarBytes);
        numValues = 1;
        for (uint8_t j = 0; j < numDims; j++)
        {
            numValues *= dimSizes[j];
        }
        numBytesPrev = numValues * (state->nRecs) * numVarBytes;
        numBytesToAdd = numValues * nRecs * numVarBytes;
        numBytesNew = numBytesPrev + numBytesToAdd;
        calibrationMemorySize += numBytesNew;
        newMem = realloc(dataBuffers[i], (size_t) numBytesNew);
        if (newMem == NULL)
            return;
        dataBuffers[i] = (uint8_t*) newMem;
        status = CDFgetzVarRangeRecordsByVarID(calCdfId, varNum, startRecord, stopRecord, dataBuffers[i] + numBytesPrev);
        if (status != CDF_OK)
        {

            fprintf(state->processingLogFile, "%sError loading %s from %s.\n", infoHeader, variables[i], filename);
            closeCdf(calCdfId);
            return;
        }
    }
    // close CDF
    closeCdf(calCdfId);

    // Update number of records found and memory allocated
    if (numberOfRecords != NULL)
        *numberOfRecords += nRecs;
    if (memoryAllocated != NULL)
        *memoryAllocated = calibrationMemorySize;

    return;
}

void setCalibrationFileName(ProcessorState *state, int year, int month, int day)
{
    Arguments *a = &state->args;
    snprintf(state->calibrationFileName, CDF_PATHNAME_LEN, "%s/%s/%04d/Swarm_%s/%02d/TiiClbr%s_Swarm_%s_%04d_%02d_%02d.cdf", a->calDir, a->calVersion, year, a->satellite, month, a->calVersion, a->satellite, year, month, day);

    return;
}

void setTracisFileName(ProcessorState *state, int year, int month, int day)
{
    Arguments *a = &state->args;
    snprintf(state->tracisFileName, CDF_PATHNAME_LEN, "%s/%s/SW_OPER_EFI%cTISL1B_%04d%02d%02dT000000_%04d%02d%02dT235959_%s.cdf", a->tracisDir, a->tracisVersion, a->satellite[0], year, month, day, year, month, day, a->tracisVersion);

    return;
}

int checkCalDataAvailability(ProcessorState *state)
{
    setCalibrationFileName(state, state->args.year, state->args.month, state->args.day);
    if (access(state->calibrationFileName, F_OK) != 0)
    {
        fprintf(state->processingLogFile, "%sTII calibration file %s not found.\n", infoHeader, state->calibrationFileName);
        return TIICT_NO_CAL_FILE;
    }
    CDFid calCdfId;
    CDFstatus status;
    status = CDFopenCDF(state->calibrationFileName, &calCdfId);
    if (status != CDF_OK) 
    {
        fprintf(state->processingLogFile, "%sUnable to open %s.\n", infoHeader, state->calibrationFileName);
        return TIICT_CDF_READ;
    }

    // Get number of records for zVar "epoch"
    long nRecords = 0;
    status = CDFgetzVarAllocRecords(calCdfId, CDFgetVarNum(calCdfId, "epoch"), &nRecords);
    if (status != CDF_OK) 
    {
        printErrorMessageToFile(state->processingLogFile, status);
        return TIICT_CDF_READ;
    }
    closeCdf(calCdfId);

    if (nRecords < (16 * SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING))
    {
        fprintf(state->processingLogFile, "%sLess than %.0f s of data available.\n", infoHeader, (float)SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING);
        return TIICT_NOT_ENOUGH_CALIBRATION_RECORDS;
    }
    else
    {
        fprintf(state->processingLogFile, "%sProcessing %ld calibration records.\n", infoHeader, nRecords);
    }

    return TIICT_OK;

}

int loadData(ProcessorState *state)
{
    int status = TIICT_OK;

    if ((status = loadTiiData(state)) != TIICT_OK)
        return status;
 
    if ((status = loadLpCalData(state)) != TIICT_OK)
        return status;

    fprintf(state->processingLogFile, "%sNumber of records: TII Cal: %ld\tTRACIS: %ld\tLP: %ld\n", infoHeader, state->nRecs, state->nTracisRecs, state->nLpRecs);
    fflush(state->processingLogFile);

    return TIICT_OK;
}