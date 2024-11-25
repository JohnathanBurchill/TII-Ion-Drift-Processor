/*

    TII Cross-Track Ion Drift Processor: loadData.c

    Copyright (C) 2024  Johnathan K Burchill

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
#include "settings.h"
#include "errors.h"
#include "state.h"
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
    // Free memory, since this can be called multiple times in interactive mode

    int status = TIICT_OK;
    // Only used for version 0401 or greater, as earlier versions of calibration files did not include potential estimates
    state->usePotentials &= (strcmp(state->args.exportVersion, "0401") >= 0);
    // Read in the potentials in any case so they can be stored in the CDF file
    status = getLpData(state);
    if (status == TIICT_OK)
    {
        if (state->vars16hz.nLpRecs < LP_MIN_NUMBER_OF_POTENTIALS)
        {
            if (state->writeLogFiles) {
                fprintf(state->processingLogFile, "%sNot enough (%lu) LP potentials imported.\n", infoHeader, state->vars16hz.nLpRecs);
            }
            status =  TIICT_NO_LP_HM_DATA;
        }
        else
        {
            if (state->writeLogFiles) {
                fprintf(state->processingLogFile, "%sLoaded %lu LP potentials, and interpolated them to the TII times.\n", infoHeader, state->vars16hz.nLpRecs);
            }
        }
    }

    return status;

}

int getLpData(ProcessorState *state)
{
    freeVariable(state->vars16hz.lpPhiSc);
    freeVariable(state->vars16hz.lpPhiScLowGain);
    freeVariable(state->vars16hz.lpPhiScHighGain);
    state->vars16hz.nLpRecs = 0;

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
            if (state->writeLogFiles) {
                fprintf(state->processingLogFile, "%sNo LP data for %4d%02d%02d\n", infoHeader, date.tm_year+1900, date.tm_mon+1, date.tm_mday);
            }
            date.tm_mday = date.tm_mday + 1;
            continue;
        }
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sLoading LP data from %s\n", infoHeader, lpFile);
        }

        status = loadLpInputs(lpFile, &lpTimes2Hz, &lpVsHg, &lpVsLg, &lpVs, &state->vars16hz.nLpRecs);
        if (status == TIICT_MEMORY)
        {
            if (state->writeLogFiles) {
                fprintf(state->processingLogFile, "%sUnable to allocate memory for LP data. Skipping processing.", infoHeader);
            }
            return status;
        }

        date.tm_mday = date.tm_mday + 1;

    }
    if (state->vars16hz.nLpRecs == 0)
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sNo LP records found. Skipping processing.", infoHeader);
        }
        return TIICT_NO_LP_HM_DATA;
    }

    // interpolate LP data to TII times
    state->vars16hz.lpPhiScHighGain = malloc(sizeof(float) * state->vars16hz.nRecs);
    state->vars16hz.lpPhiScLowGain = malloc(sizeof(float) * state->vars16hz.nRecs);
    state->vars16hz.lpPhiSc = malloc(sizeof(float) * state->vars16hz.nRecs);
    if (state->vars16hz.lpPhiScHighGain == NULL || state->vars16hz.lpPhiScLowGain == NULL || state->vars16hz.lpPhiSc == NULL)
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sUnable to allocate memory for LP interpolation. Skipping processing.", infoHeader);
        }
        return TIICT_MEMORY;
    }
    bzero(state->vars16hz.lpPhiScHighGain, sizeof(float) * state->vars16hz.nRecs);
    bzero(state->vars16hz.lpPhiScLowGain, sizeof(float) * state->vars16hz.nRecs);
    bzero(state->vars16hz.lpPhiSc, sizeof(float) * state->vars16hz.nRecs);

    double *tiiTime = state->vars16hz.timestamp;
    interpolate(lpTimes2Hz, lpVsHg, state->vars16hz.nLpRecs, tiiTime, state->vars16hz.nRecs, state->vars16hz.lpPhiScHighGain);
    interpolate(lpTimes2Hz, lpVsLg, state->vars16hz.nLpRecs, tiiTime, state->vars16hz.nRecs, state->vars16hz.lpPhiScLowGain);
    interpolate(lpTimes2Hz, lpVs, state->vars16hz.nLpRecs, tiiTime, state->vars16hz.nRecs, state->vars16hz.lpPhiSc);

    free(lpTimes2Hz);
    free(lpVsHg);
    free(lpVsLg);
    free(lpVs);

    return TIICT_OK;

}

// Adapted from SLIDEM
int loadLpInputs(const char *cdfFile, double **lpTime, double **lpPhiScHighGain, double **lpPhiScLowGain, double **lpPhiSc, size_t *numberOfRecords)
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

    // Update number of records found and memory allocated
    *numberOfRecords += numRecs;
    // *totalMemoryAllocated = fpMemorySize;

    return TIICT_OK;

}

int loadCalData(ProcessorState *state)
{
    // Store processing date
    int year = state->args.year;
    int month = state->args.month;
    int day = state->args.day;

    ProcessorVariables_t *v = &state->vars16hz;
    CalibrationVariable_t variables[] = {
        {"epoch", (void**)&v->timestamp, 0},
        {"Latitude", (void**)&v->latitude, 0},
        {"Longitude", (void**)&v->longitude, 0},
        {"Radius", (void**)&v->radius, 0},
        {"Radius", (void**)&v->radius, 0},
        {"QDLat", (void**)&v->qdlat, 0},
        {"MLT", (void**)&v->mlt, 0},
        {"1st Moment - H", (void**)&v->mxh, 0},
        {"1st Moment - H", (void**)&v->myh, 1},
        {"1st Moment - V", (void**)&v->mxv, 0},
        {"1st Moment - V", (void**)&v->myv, 1},
        {"Det_H__vX", (void**)&v->vsatx, 0},
        {"Det_H__vX", (void**)&v->vsaty, 1},
        {"Det_H__vX", (void**)&v->vsatz, 2},
        {"Det H CorVx", (void**)&v->vicrx, 0},
        {"Det H CorVx", (void**)&v->vicry, 1},
        {"Det H CorVx", (void**)&v->vicrz, 2},
        {"Sat_Vel_N", (void**)&v->vsatn, 0},
        {"Sat_Vel_N", (void**)&v->vsate, 1},
        {"Sat_Vel_N", (void**)&v->vsatc, 2},
        {"B-North", (void**)&v->bn, 0},
        {"B-North", (void**)&v->be, 1},
        {"B-North", (void**)&v->bc, 2},
        {"MCP_Voltage_H", (void**)&v->vmcph, 0},
        {"MCP_Voltage_V", (void**)&v->vmcpv, 0},
        {"Bias_Grid_H", (void**)&v->vbiash, 0},
        {"Bias_Grid_V", (void**)&v->vbiasv, 0},
        {"Faceplate_Volt_Mon_H", (void**)&v->vfp, 0},
    };
    const int nVariables = sizeof(variables) / sizeof(CalibrationVariable_t);

    // Get data for prior, requested, and following days
    for (int8_t i = -1; i < 2; i++)
    {
        state->args.year = year;
        state->args.month= month;
        state->args.day = day + i;
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sLoading calibration data for %04d%02d%02d\n", infoHeader, state->args.year, state->args.month, state->args.day);
        }

        loadCalDataFromDate(i, state, variables, nVariables);
    }

    // Reset processing date
    state->args.year = year;
    state->args.month = month;
    state->args.day = day;
    if (state->writeLogFiles) {
        fprintf(state->processingLogFile, "%sNumber of records: %ld\n", infoHeader, state->vars16hz.nRecs);
        fprintf(state->processingLogFile, "%sLoaded %ld bytes (%ld MB) of calibration data.\n", infoHeader, state->vars16hz.memoryAllocated, state->vars16hz.memoryAllocated / 1024 / 1024);
        fflush(state->processingLogFile);
    }

    if (state->vars16hz.nRecs < 16*SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING) {
        return TIICT_NOT_ENOUGH_CALIBRATION_RECORDS;
    }

    // Ensure all 16 Hz variables have sufficient memory
    reallocVariables(&state->vars16hz, state->vars16hz.nRecs);

    return TIICT_OK;
}

void loadCalDataFromDate(const DayType dayType, ProcessorState *state, CalibrationVariable_t *variables, int nVariables)
{
    // c time manipulation: see  https://fresh2refresh.com/c-programming/c-time-related-functions/
    struct tm timestructure = {0};
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

    setCalibrationFileName(state, state->args.year, state->args.month, state->args.day, CAL_FILE_TII);
    char *calibrationFileName = state->tiiCalibrationFileName;
    if (state->writeLogFiles) {
        fprintf(state->processingLogFile, "%s from %s\n", infoHeader, calibrationFileName);
    }

    // Open the CDF file with validation
    CDFsetValidate(VALIDATEFILEon);
    CDFid calCdfId;
    CDFstatus status;
    status = CDFopenCDF(calibrationFileName, &calCdfId);
    if (status != CDF_OK)
    {
        // Not necessarily an error. For example, some dates will have not calibration data.
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sSkipping this date.\n", infoHeader);
        }
        return;
    }

    if (state->writeLogFiles) {
        fprintf(state->processingLogFile, "%sFound CDF file.\n", infoHeader);
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
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sProblem with calibration file. Skipping this date.\n", infoHeader);
        }
        closeCdf(calCdfId);
        return;
    }
    long nRecs = 0;
    long calibrationMemorySize = 0;
    status = CDFgetzVarAllocRecords(calCdfId, CDFgetVarNum(calCdfId, "epoch"), &nRecs);
    if (status != CDF_OK)
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sProblem with calibration file. Skipping this date.\n", infoHeader);
        }
        closeCdf(calCdfId);
        return;
    }
    if (nRecs < (16*SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING))
    {
        // Not enough to do anything useful
        // TODO: maybe increase this threshold to require a larger number of points each day?
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sFewer than %.0f s of data. Skipping this date.\n", infoHeader, (float)SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING);
        }
        CDFcloseCDF(calCdfId);
        return;
    }

    // Determine start and stop record numbers to load
    long startRecord;
    long stopRecord;
    double timeReference, recordTime;
    // get epoch variable number
    long epochNum = CDFgetVarNum(calCdfId, "epoch");
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
                    if (state->writeLogFiles) {
                        fprintf(state->processingLogFile, "%sCould not read epoch record from CDF file. Skipping this calibration date.\n", infoHeader);
                    }
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
                    if (state->writeLogFiles) {
                        fprintf(state->processingLogFile, "%sCould not read epoch record from CDF file. Skipping this calibration date.\n", infoHeader);
                    }
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
            if (state->writeLogFiles) {
                fprintf(state->processingLogFile, "%sError: Day type must be one of PREVIOUS_DAY, REQUESTED_DAY, or NEXT_DAY. Skipping this calibration data.\n", infoHeader);
            }
            closeCdf(calCdfId);
            return;
    }

    // Update number of records
    nRecs = stopRecord - startRecord + 1;
    if ((dayType == REQUESTED_DAY && nRecs < (16*SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING)) || ((dayType == PREVIOUS_DAY || dayType == NEXT_DAY) && nRecs < (16*SECONDS_OF_BOUNDARY_DATA_REQUIRED_FOR_PROCESSING)))
    {
        // Not enough to do anything useful
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sFewer than %.0f s of data meet constraints. Skipping this date.\n", infoHeader, (float)SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING);
        }
        CDFcloseCDF(calCdfId);
        return;
    }

    for (int i = 0; i < nVariables; i++) {
        status = loadCdfVariable(state, calCdfId, &variables[i], startRecord, stopRecord, &calibrationMemorySize);
        if (status != TIICT_OK) {
            if (state->writeLogFiles) {
                fprintf(state->processingLogFile, "%sError reading input variable %s.\n", infoHeader, variables[i].name);
            }
            closeCdf(calCdfId);
            return;
        }
    }

    // close CDF
    closeCdf(calCdfId);
    // Number of records obtained for this date
    if (state->writeLogFiles) {
        switch (dayType)
        {
            case PREVIOUS_DAY:
                fprintf(state->processingLogFile, "%sGot %ld s of data for previous day\n", infoHeader, nRecs / 16);
                break;
            case REQUESTED_DAY:
                fprintf(state->processingLogFile, "%sGot %ld s of data for requested day\n", infoHeader, nRecs / 16);
                break;
            case NEXT_DAY:
                fprintf(state->processingLogFile, "%sGot %ld s of data for next day\n", infoHeader, nRecs / 16);
                break;
            default:
                break;
        }
    }
    // Update number of records found and memory allocated
    state->vars16hz.nRecs += nRecs;
    state->vars16hz.memoryAllocated = calibrationMemorySize;

    return;
}

int loadCdfVariable(ProcessorState *state, CDFid calCdfId, CalibrationVariable_t *variable, long startRecord, long stopRecord, long *calibrationMemorySize)
{
    int status = TIICT_OK;
    if (variable == NULL) {
        return TIICT_ARGS_BAD;
    }

    // Check CDF info
    long varNum = CDFgetVarNum(calCdfId, variable->name);
    if (varNum < CDF_OK)
    {
        printErrorMessage(varNum);
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sVariable %s not found.\n", infoHeader, variable->name);
        }
        return TIICT_CDF_READ;
    }

    long dimensionsPerRecord = 0;
    status = CDFgetzVarNumDims(calCdfId, varNum, &dimensionsPerRecord);
    long dimensionSizes[CDF_MAX_DIMS];
    status = CDFgetzVarDimSizes(calCdfId, varNum, dimensionSizes);
    long valuesPerRecord = 1;
    for (uint8_t j = 0; j < dimensionsPerRecord; j++)
    {
        valuesPerRecord *= dimensionSizes[j];
    }

    // Calculate size of memory to allocate for loading
    long dataType = 0;
    status = CDFgetzVarDataType(calCdfId, varNum, &dataType);
    long bytesPerValue = 0;
    status = CDFgetDataTypeSize(dataType, &bytesPerValue);

    long nRecords = stopRecord - startRecord + 1;
    long nBytesToLoad = nRecords * valuesPerRecord * bytesPerValue;
    void *mem = malloc(nBytesToLoad);
    if (mem == NULL) {
        return TIICT_MEMORY;
    }
    status = CDFgetzVarRangeRecordsByVarID(calCdfId, varNum, startRecord, stopRecord, mem);
    if (status != CDF_OK)
    {
        return TIICT_CDF_READ;
    }

    long nRecordsAlreadyLoaded = state->vars16hz.nRecs;
    long nBytesAlreadyLoaded = nRecordsAlreadyLoaded * valuesPerRecord * bytesPerValue;
    long nBytesToAdd = nRecords * bytesPerValue;
    calibrationMemorySize += nBytesToAdd;
    long nBytesNew = nBytesAlreadyLoaded + nBytesToAdd;
    status = reallocVariable(variable->memoryPointer, (size_t) nBytesNew);
    if (status != TIICT_OK) {
        return status;
    }

    for (int i = 0; i < nRecords; i++) {
        switch(dataType) {
            case CDF_REAL8:
            case CDF_EPOCH:
                ((double*)*variable->memoryPointer)[nRecordsAlreadyLoaded + i] = ((double*)mem)[i*valuesPerRecord + variable->recordValueOffset];
                break;
            default:
                // EfiCalCdfs have only epoch and real4 types
                ((float*)*variable->memoryPointer)[nRecordsAlreadyLoaded + i] = ((float*)mem)[i*valuesPerRecord + variable->recordValueOffset];
                break;
        }
    }

    free(mem);

    return TIICT_OK;
}

void setCalibrationFileName(ProcessorState *state, int year, int month, int day, CalibrationFileType_enum fileType)
{
    Arguments *a = &state->args;
    switch (fileType) {
        case CAL_FILE_TII:
            snprintf(state->tiiCalibrationFileName, CDF_PATHNAME_LEN, "%s/%s/%04d/Swarm_%s/%02d/TiiClbr%s_Swarm_%s_%04d_%02d_%02d.cdf", a->calDir, a->calVersion, year, a->satellite, month, a->calVersion, a->satellite, year, month, day);
            break;
        case CAL_FILE_LP:
            // TODO Optimize this to not have to check for latest LP file version
            getLpInputFilename(state->args.satellite[0], year, month, day, state->args.lpDir, state->lpCalibrationFileName);
            break;
        default:
            break;
    }

    return;
}

int checkCalDataAvailability(ProcessorState *state)
{
    setCalibrationFileName(state, state->args.year, state->args.month, state->args.day, CAL_FILE_TII);
    if (access(state->tiiCalibrationFileName, F_OK) != 0)
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sCalibration file %s not found. Skipping this date.\n", infoHeader, state->tiiCalibrationFileName);
        }
        return TIICT_NO_CAL_FILE;
    }
    CDFid calCdfId;
    CDFstatus status;
    status = CDFopenCDF(state->tiiCalibrationFileName, &calCdfId);
    if (status != CDF_OK)
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sUnable to open %s. Skipping this date.\n", infoHeader, state->tiiCalibrationFileName);
        }
        return TIICT_CDF_READ;
    }

    // Get number of records for zVar "epoch"
    long nRecords = 0;
    status = CDFgetzVarAllocRecords(calCdfId, CDFgetVarNum(calCdfId, "epoch"), &nRecords);
    if (status != CDF_OK)
    {
        if (state->writeLogFiles) {
            printErrorMessageToFile(state->processingLogFile, status);
        }
        return TIICT_CDF_READ;
    }
    closeCdf(calCdfId);

    if (nRecords < (16 * SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING))
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sLess than %.0f s of data available. Skipping this date.\n", infoHeader, (float)SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING);
        }
        return TIICT_NOT_ENOUGH_CALIBRATION_RECORDS;
    }
    else
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sProcessing %ld calibration records for this date.\n", infoHeader, nRecords);
        }
    }

    return TIICT_OK;

}

void freeVariable(void *var)
{
    free(var);
    var = NULL;

    return;
}

void freeVariables(ProcessorVariables_t *vars)
{
    freeVariable(vars->timestamp);
    freeVariable(vars->latitude);
    freeVariable(vars->longitude);
    freeVariable(vars->radius);
    freeVariable(vars->qdlat);
    freeVariable(vars->mlt);
    freeVariable(vars->mxh);
    freeVariable(vars->myh);
    freeVariable(vars->mxv);
    freeVariable(vars->myv);
    freeVariable(vars->vsatx);
    freeVariable(vars->vsaty);
    freeVariable(vars->vsatz);
    freeVariable(vars->enhRaw);
    freeVariable(vars->envRaw);
    freeVariable(vars->enh);
    freeVariable(vars->env);
    freeVariable(vars->vixh);
    freeVariable(vars->vixherror);
    freeVariable(vars->vixv);
    freeVariable(vars->vixverror);
    freeVariable(vars->viy);
    freeVariable(vars->viyerror);
    freeVariable(vars->viz);
    freeVariable(vars->vizerror);
    freeVariable(vars->vsatn);
    freeVariable(vars->vsate);
    freeVariable(vars->vsatc);
    freeVariable(vars->ectxh);
    freeVariable(vars->ectyh);
    freeVariable(vars->ectzh);
    freeVariable(vars->ectxv);
    freeVariable(vars->ectyv);
    freeVariable(vars->ectzv);
    freeVariable(vars->bctx);
    freeVariable(vars->bcty);
    freeVariable(vars->bctz);
    freeVariable(vars->bn);
    freeVariable(vars->be);
    freeVariable(vars->bc);
    freeVariable(vars->vicrx);
    freeVariable(vars->vicry);
    freeVariable(vars->vicrz);
    freeVariable(vars->flags);
    freeVariable(vars->fitInfo);
    freeVariable(vars->geoelectricPotential);
    freeVariable(vars->geoelectricPotentialDifference);
    freeVariable(vars->maxAbsGeoelectricPotentialBaselineSlope);
    freeVariable(vars->ehxAdjusted);
    freeVariable(vars->ehxAdjustmentParameter);
    freeVariable(vars->geoelectricPotentialDetrended);
    freeVariable(vars->maxAbsGeoelectricPotentialDetrendedBaselineSlope);
    freeVariable(vars->orbitRegion);
    freeVariable(vars->lpTimes);
    freeVariable(vars->lpPhiScHighGain);
    freeVariable(vars->lpPhiScLowGain);
    freeVariable(vars->lpPhiSc);
    freeVariable(vars->potentials);
    freeVariable(vars->xhat);
    freeVariable(vars->yhat);
    freeVariable(vars->zhat);
    freeVariable(vars->geoelectricPotential);
    freeVariable(vars->maxAbsGeoelectricPotentialBaselineSlope);
    freeVariable(vars->geoelectricPotentialDifference);
    freeVariable(vars->geoelectricPotentialDetrended);
    freeVariable(vars->maxAbsGeoelectricPotentialDetrendedBaselineSlope);
    freeVariable(vars->ehxAdjusted);
    freeVariable(vars->ehxAdjustmentParameter);

    return;
}

int reallocVariable(void **var, size_t newSize)
{

    if (var == NULL) {
        return TIICT_ARGS_BAD;
    }
    void *mem = realloc(*var, newSize * sizeof **var);
    if (mem == NULL) {
        return TIICT_MEMORY;
    }
    *var = mem;

    return TIICT_OK;
}

int reallocVariables(ProcessorVariables_t *vars, size_t newSize)
{
    int status = TIICT_OK;
    status |= reallocVariable((void*)&vars->timestamp, newSize);
    status |= reallocVariable((void*)&vars->latitude, newSize);
    status |= reallocVariable((void*)&vars->longitude, newSize);
    status |= reallocVariable((void*)&vars->radius, newSize);
    status |= reallocVariable((void*)&vars->qdlat, newSize);
    status |= reallocVariable((void*)&vars->mlt, newSize);
    status |= reallocVariable((void*)&vars->mxh, newSize);
    status |= reallocVariable((void*)&vars->myh, newSize);
    status |= reallocVariable((void*)&vars->mxv, newSize);
    status |= reallocVariable((void*)&vars->myv, newSize);
    status |= reallocVariable((void*)&vars->vmcph, newSize);
    status |= reallocVariable((void*)&vars->vmcpv, newSize);
    status |= reallocVariable((void*)&vars->vbiash, newSize);
    status |= reallocVariable((void*)&vars->vbiasv, newSize);
    status |= reallocVariable((void*)&vars->vfp, newSize);
    status |= reallocVariable((void*)&vars->vsatx, newSize);
    status |= reallocVariable((void*)&vars->vsaty, newSize);
    status |= reallocVariable((void*)&vars->vsatz, newSize);
    status |= reallocVariable((void*)&vars->enhRaw, newSize);
    status |= reallocVariable((void*)&vars->envRaw, newSize);
    status |= reallocVariable((void*)&vars->enh, newSize);
    status |= reallocVariable((void*)&vars->env, newSize);
    status |= reallocVariable((void*)&vars->vixh, newSize);
    status |= reallocVariable((void*)&vars->vixherror, newSize);
    status |= reallocVariable((void*)&vars->vixv, newSize);
    status |= reallocVariable((void*)&vars->vixverror, newSize);
    status |= reallocVariable((void*)&vars->viy, newSize);
    status |= reallocVariable((void*)&vars->viyerror, newSize);
    status |= reallocVariable((void*)&vars->viz, newSize);
    status |= reallocVariable((void*)&vars->vizerror, newSize);
    status |= reallocVariable((void*)&vars->vsatn, newSize);
    status |= reallocVariable((void*)&vars->vsate, newSize);
    status |= reallocVariable((void*)&vars->vsatc, newSize);
    status |= reallocVariable((void*)&vars->bn, newSize);
    status |= reallocVariable((void*)&vars->be, newSize);
    status |= reallocVariable((void*)&vars->bc, newSize);
    status |= reallocVariable((void*)&vars->vicrx, newSize);
    status |= reallocVariable((void*)&vars->vicry, newSize);
    status |= reallocVariable((void*)&vars->vicrz, newSize);
    status |= reallocVariable((void*)&vars->ectxh, newSize);
    status |= reallocVariable((void*)&vars->ectyh, newSize);
    status |= reallocVariable((void*)&vars->ectzh, newSize);
    status |= reallocVariable((void*)&vars->ectxv, newSize);
    status |= reallocVariable((void*)&vars->ectyv, newSize);
    status |= reallocVariable((void*)&vars->ectzv, newSize);
    status |= reallocVariable((void*)&vars->bctx, newSize);
    status |= reallocVariable((void*)&vars->bcty, newSize);
    status |= reallocVariable((void*)&vars->bctz, newSize);
    status |= reallocVariable((void*)&vars->flags, newSize);
    status |= reallocVariable((void*)&vars->fitInfo, newSize);
    status |= reallocVariable((void*)&vars->geoelectricPotential, newSize);
    status |= reallocVariable((void*)&vars->geoelectricPotentialDifference, newSize);
    status |= reallocVariable((void*)&vars->maxAbsGeoelectricPotentialBaselineSlope, newSize);
    status |= reallocVariable((void*)&vars->ehxAdjusted, newSize);
    status |= reallocVariable((void*)&vars->ehxAdjustmentParameter, newSize);
    status |= reallocVariable((void*)&vars->geoelectricPotentialDetrended, newSize);
    status |= reallocVariable((void*)&vars->maxAbsGeoelectricPotentialDetrendedBaselineSlope, newSize);
    status |= reallocVariable((void*)&vars->orbitRegion, newSize);
    status |= reallocVariable((void*)&vars->lpTimes, newSize);
    status |= reallocVariable((void*)&vars->lpPhiScHighGain, newSize);
    status |= reallocVariable((void*)&vars->lpPhiScLowGain, newSize);
    status |= reallocVariable((void*)&vars->lpPhiSc, newSize);
    status |= reallocVariable((void*)&vars->potentials, newSize);
    status |= reallocVariable((void*)&vars->xhat, newSize * 3);
    status |= reallocVariable((void*)&vars->yhat, newSize * 3);
    status |= reallocVariable((void*)&vars->zhat, newSize * 3);

    return status;
}


