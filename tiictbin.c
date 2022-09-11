/*

    TII Cross-Track Ion Drift Processor: tiictbin.c

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
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <dirent.h>

#include "utilities.h"
#include "tiictbin.h"
#include <cdf.h>

// https://www.gnu.org/software/gsl/doc/html/
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>


char infoHeader[50];

int main(int argc, char* argv[])
{
    time_t currentTime;
    struct tm * timeParts;
    time(&currentTime);
    timeParts = localtime(&currentTime);
    char * dateString = asctime(timeParts);
    char date[255];
    snprintf(date, strlen(dateString), "%s", dateString);

    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--about") == 0)
        {
            fprintf(stdout, "tiictbin - prints specified parameter binned in QD latitude and MLT within specified limits. Version %s.\n", SOFTWARE_VERSION);
            fprintf(stdout, "Copyright (C) 2022  Johnathan K Burchill\n");
            fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
            fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
            fprintf(stdout, "under the terms of the GNU General Public License.\n");

            exit(0);
        }
    }


    if (argc != 10)
    {
        fprintf(stdout, "usage: %s directory satelliteLetter parameterName qdlatmin qdlatmax deltaqdlat mltmin mltmax deltamlt\n", argv[0]);
        exit(1);
    }
    const char *directory = argv[1];
    const char *satelliteLetter = argv[2];
    const char *parameterName = argv[3];
    double qdlatmin = atof(argv[4]);
    double qdlatmax = atof(argv[5]);
    double deltaqdlat = atof(argv[6]);
    double mltmin = atof(argv[7]);
    double mltmax = atof(argv[8]);
    double deltamlt = atof(argv[9]);
    int qdlatIndex = 0;
    int mltIndex = 0;

    if (deltaqdlat <= 0.0)
    {
        fprintf(stderr, "%s: deltaqdlat must be greater than 0.0\n", argv[0]);
        exit(0);
    }
    if (deltamlt <= 0.0)
    {
        fprintf(stderr, "%s: deltamlt must be greater than 0.0\n", argv[0]);
        exit(0);
    }

    int nQDLats = (int)floor((qdlatmax - qdlatmin) / deltaqdlat);
    int nMLTs = (int)floor((mltmax - mltmin) / deltamlt);
    if (nQDLats <= 0)
    {
        fprintf(stderr, "%s: invalid QD latitude bin specification.\n", argv[0]);
        exit(0);
    }
    if (nMLTs <= 0)
    {
        fprintf(stderr, "%s: invalid MLT bin specification.\n", argv[0]);
        exit(0);
    }

    // Allocate memory for binning
    double *bins = (double *)malloc((size_t)(nQDLats * nMLTs) * sizeof(double));
    long *binNumbers = (long *)malloc((size_t)(nQDLats * nMLTs) * sizeof(long));
    if (bins == NULL)
    {
        fprintf(stderr, "Unable to allocate memory.\n");
        exit(0);
    }
    // Access bins with bins[mltIndex * nQDLats + qdlatIndex];


    // Turn off GSL failsafe error handler. We typically check the GSL return codes.
    gsl_set_error_handler_off();

    // Count files to process
    DIR *dir = opendir(directory);
    struct dirent *entry;
    if (dir == NULL)
    {
        fprintf(stdout, "Could not open directory %s for reading.\n", directory);
        exit(1);
    }
    char *filename;
    long nFiles = 0;
    while ((entry = readdir(dir)) != NULL)
    {
        filename = entry->d_name;
        uint16_t length = strlen(filename);
        if (length != 59)
            continue;
        if (strcmp(filename+length-3, "cdf") != 0)
            continue;
        const char satellite = *(filename+strlen(filename)-48);
        if (satellite == satelliteLetter[0])
            nFiles++;
    }
    closedir(dir);
    // Process files
    dir = opendir(directory);
    if (dir == NULL)
    {
        fprintf(stdout, "Could not open directory %s for reading.\n", directory);
        exit(1);
    }
    long processedFiles = 0;
    int percentDone = 0;
    while ((entry = readdir(dir)) != NULL)
    {
        filename = entry->d_name;
        uint16_t length = strlen(filename);
        if (length != 59)
        {
            continue;
        }
        if (strcmp(filename+length-3, "cdf") != 0)
        {
            continue;
        }
        const char satellite = *(filename+strlen(filename)-48);
        if (satellite == satelliteLetter[0])
        {
            // Do the processing
            CDFstatus status;

            // The memory pointers
            uint8_t * dataBuffers[NUM_DATA_VARIABLES];
            for (uint8_t i = 0; i < NUM_DATA_VARIABLES; i++)
            {
                dataBuffers[i] = NULL;
            }
            long nRecs = 0;
            char fullPath[CDF_PATHNAME_LEN];
            sprintf(fullPath, "%s/%s", directory, filename);
            bool fourByteCalFlag = false;
            loadCrossTrackData(fullPath, dataBuffers, &nRecs, &fourByteCalFlag, parameterName);
            double *pEpoch = (double*) dataBuffers[0];
            long long timeIndex = 0;

            // The macros TIME(), QDLAT(), etc. give the value at index timeIndex.
            double epoch0 = TIME();

            uint8_t minorVersion = getMinorVersion(filename);
            // Number of records
            // Include all measurements for Swarm C, which are set to 0 flag always,
            // except those for which baseline calibration was not done or was problematic: CALFLAG() == 0
            for (timeIndex = 0; timeIndex < nRecs; timeIndex++)
            {
                if (FLAG() == 4 || (satellite == 'C' && minorVersion == 1) || (satellite == 'C' && fourByteCalFlag == true && (CALFLAG() == 0)))
                {
                    // Access bins with bins[mltIndex * nQDLats + qdlatIndex];
                    if (isfinite(PARAMETER()))
                    {
                        mltIndex = (int) floor((MLT() - mltmin) / deltamlt);
                        qdlatIndex = (int) floor((QDLAT() - qdlatmin) / deltaqdlat);
                        if (mltIndex >= 0 && mltIndex < nMLTs && qdlatIndex >=0 && qdlatIndex < nQDLats)
                        {
                            // TODO handle vector parameters
                            bins[mltIndex * nQDLats + qdlatIndex] += PARAMETER();
                            binNumbers[mltIndex * nQDLats + qdlatIndex] += 1;
                        }
                    }
                }
            }

            // free memory
            for (uint8_t i = 0; i < NUM_DATA_VARIABLES; i++)
            {
                free(dataBuffers[i]);
            }
            processedFiles++;
            percentDone = (int) floor((double)processedFiles / (double)nFiles * 100.0);
            if (percentDone % 5 == 0)
                fprintf(stderr, "Processed %ld of %ld files (%d%%)\n", processedFiles, nFiles, percentDone);
        }
    }

    double qdlat = 0.0;
    double mlt = 0.0;

    for (int q = 0; q < nQDLats; q++)
    {
        for (int m = 0; m < nMLTs; m++)
        {
            if (binNumbers[m * nQDLats + q] > 0)
                bins[m * nQDLats + q] /= (double) binNumbers[m * nQDLats + q];

            qdlat = qdlatmin + deltaqdlat * ((double)q + 0.5);
            mlt = mltmin + deltamlt* ((double)m + 0.5);

            fprintf(stdout, "QD: %lf MLT: %lf VAL: %lf N: %ld\n", qdlat, mlt, bins[m * nQDLats + q], binNumbers[m * nQDLats + q]);
        }
    }

    closedir(dir);

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

void loadCrossTrackData(const char *filename, uint8_t **dataBuffers, long *numberOfRecords, bool *fourByteCalFlag, const char *parameterName)
{
    char validationFileName[CDF_PATHNAME_LEN];
    snprintf(validationFileName, strlen(filename)-3, "%s", filename);

    uint8_t minorVersion = getMinorVersion(filename);

    // Open the CDF file with validation
    CDFsetValidate(VALIDATEFILEon);
    CDFid calCdfId;
    CDFstatus status;
    status = CDFopenCDF(validationFileName, &calCdfId);
    if (status != CDF_OK) 
    {
        printErrorMessage(status);
        // Not necessarily an error. For example, some dates will have not calibration data.
        fprintf(stdout, "%sSkipping this file.\n", infoHeader);
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
        printErrorMessage(status);
        fprintf(stdout, "%sProblem with data file. Skipping this file.\n", infoHeader);
        closeCdf(calCdfId);
        return;
    }
    long nRecs, memorySize = 0;
    status = CDFgetzVarAllocRecords(calCdfId, CDFgetVarNum(calCdfId, "Timestamp"), &nRecs);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        fprintf(stdout, "%sProblem with data file. Skipping this file.\n", infoHeader);
        closeCdf(calCdfId);
        return;
    }

    // Variables
    uint8_t nVars = NUM_DATA_VARIABLES;
    char* variables[NUM_DATA_VARIABLES] = {
        "Timestamp",
        "MLT",
        "QDLatitude",
        "Quality_flags",
        "Calibration_flags",
        ""
    };
    variables[NUM_DATA_VARIABLES-1] = (char *) parameterName;
    if (minorVersion == 1)
    {
        variables[2] = "flags";
    }
    for (uint8_t i = 0; (i<nVars-1) || (i == nVars-1 && minorVersion == 2); i++)
    {
        status = CDFconfirmzVarExistence(calCdfId, variables[i]);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            fprintf(stdout, "%sError reading variable %s. Skipping this file.\n", infoHeader, variables[i]);
            closeCdf(calCdfId);
            exit(1);
        }
    }
    
    long varNum, numValues, numVarBytes;
    long numBytesPrev, numBytesToAdd, numBytesNew;

    for (uint8_t i = 0; (i<nVars-1) || (i == nVars-1 && minorVersion == 2); i++)
    {
        varNum = CDFgetVarNum(calCdfId, variables[i]);
        if (varNum < CDF_OK)
        {
            printErrorMessage(varNum);
            fprintf(stdout, "%sError reading variable ID for %s. Skipping this file.\n", infoHeader, variables[i]);
            closeCdf(calCdfId);
            exit(1);
        }
        status = CDFgetzVarNumDims(calCdfId, varNum, &numDims);
        status = CDFgetzVarDimSizes(calCdfId, varNum, dimSizes);
        status = CDFgetzVarDataType(calCdfId, varNum, &dataType);
        // Calculate new size of memory to allocate
        status = CDFgetDataTypeSize(dataType, &numVarBytes);
        if (i == 3 && numVarBytes == 4)
        {
            *fourByteCalFlag = true;            
        }

        numValues = 1;
        for (uint8_t j = 0; j < numDims; j++)
        {
            numValues *= dimSizes[j];
        }
        memorySize = numValues * nRecs * numVarBytes;
        dataBuffers[i] = (uint8_t*) malloc((size_t) memorySize);
        status = CDFgetzVarAllRecordsByVarID(calCdfId, varNum, dataBuffers[i]);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            fprintf(stdout, "%sError loading data for %s. Skipping this file.\n", infoHeader, variables[i]);
            closeCdf(calCdfId);
            return;
        }
    }
    // close CDF
    closeCdf(calCdfId);
    // Update number of records found and memory allocated
    *numberOfRecords = nRecs;

}

uint8_t getMinorVersion(const char *filename)
{
    char minorVersionChar = *(filename + strlen(filename)-5);
    uint8_t minorVersion = (uint8_t) strtol(&minorVersionChar, NULL, 10);
    return minorVersion;
}
