/*

    TII Cross-Track Ion Drift Processor: tiictqualityflags.c

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
#include "tiictqualityflags.h"
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
            fprintf(stdout, "tiictqualityflags - prints info about quality flags produced by the TII Cross-track ion drift processor, version %s.\n", SOFTWARE_VERSION);
            fprintf(stdout, "Copyright (C) 2022  Johnathan K Burchill\n");
            fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
            fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
            fprintf(stdout, "under the terms of the GNU General Public License.\n");

            exit(0);
        }
    }


    if (argc != 3)
    {
        fprintf(stdout, "usage: %s directory satelliteLetter\n", argv[0]);
        exit(1);
    }
    const char *directory = argv[1];
    const char *satelliteLetter = argv[2];

    // Turn off GSL failsafe error handler. We typically check the GSL return codes.
    gsl_set_error_handler_off();


    DIR *dir = opendir(directory);
    struct dirent *entry;
    if (dir == NULL)
    {
        fprintf(stdout, "Could not open directory %s for reading.\n", directory);
        exit(1);
    }
    char *filename;
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
            loadCrossTrackData(fullPath, dataBuffers, &nRecs, &fourByteCalFlag);
            double *pEpoch = (double*) dataBuffers[0];
            long long timeIndex = 0;

            // The macros TIME(), QDLAT(), etc. give the value at index timeIndex.
            double epoch0 = TIME();

            uint8_t minorVersion = getMinorVersion(filename);
            // Number of records
            // Include all measurements for Swarm C, which are set to 0 flag always
            long nFlag4 = 0, nTotal = 0;
            for (timeIndex = 0; timeIndex < nRecs; timeIndex++)
            {
                nTotal++;

                if (FLAG() == 4 || (satellite == 'C' && minorVersion == 1) || (satellite == 'C' && fourByteCalFlag == true && (CALFLAG() == 0)))
                {
                    nFlag4++;
                }
            }

            // Allocate arrays
            double * time = malloc((size_t)((long)nFlag4 * sizeof(double)));
            double * viy4 = malloc((size_t)((long)nFlag4 * sizeof(double)));
            double * gslWorkspace = malloc((size_t)((long)nFlag4 * sizeof(double)));

            nFlag4 = 0;
            double meanTime = 0.;
            for (timeIndex = 0; timeIndex < nRecs; timeIndex++)
            {
                if (FLAG() == 4 || (satellite == 'C' && minorVersion == 1) || (satellite == 'C' && fourByteCalFlag == true && (CALFLAG() == 0)))
                {
                    time[nFlag4] = TIME()/1000.;
                    meanTime += time[nFlag4];
                    viy4[nFlag4] = (double)VIY();
                    nFlag4++;
                }
            }
            meanTime /= (double)nFlag4;
            double viyMad = gsl_stats_mad(viy4, 1, nFlag4, gslWorkspace);
            fprintf(stdout, "%f %ld %ld %f\n", meanTime, nTotal, nFlag4, viyMad);

            // free memory
            for (uint8_t i = 0; i < NUM_DATA_VARIABLES; i++)
            {
                free(dataBuffers[i]);
            }
            free(time);
            free(viy4);
            free(gslWorkspace);
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

void loadCrossTrackData(const char * filename, uint8_t **dataBuffers, long *numberOfRecords, bool *fourByteCalFlag)
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
        "QDLatitude"
        "Viy",
        "Viy_error",
        "Quality_flags",
        "Calibration_flags"
    };
    if (minorVersion == 1)
    {
        variables[4] = "flags";
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
