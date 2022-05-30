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

    if (satelliteLetter[0] == 'C')
    {
        fprintf(stdout, "Swarm C TCT version 0302 quality is not flagged.\n");
        exit(0);
    }


    DIR *dir = opendir(directory);
    struct dirent *entry;
    if (dir == NULL)
    {
        fprintf(stdout, "Could not open directory %s for reading.\n", directory);
        exit(1);
    }
    char *filename;

    long nFlag4 = 0, nTotal = 0, nHighLat = 0;

    // Count number of files in directory:
    long nFiles = 0;
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
            nFiles++;
    }

    if (nFiles == 0)
    {
        closedir(dir);
        fprintf(stdout, "Swarm %c: no TCT CDF files found.\n", satelliteLetter[0]);
        exit(0);
    }

    // Process the files
    rewinddir(dir);
    long statusInterval = (long) (STATUS_INTERVAL_FRACTION * (double) nFiles);
    long counter = 0;
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
            if (counter % statusInterval == 0)
            {
                fprintf(stdout, "\r%3.0lf%%", 100.0 * (double)counter / (double)nFiles);
		fflush(stdout);
	    }
	    counter++;

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
            long timeIndex = 0;

            // Number of records
            // Include all measurements for Swarm C, which are set to 0 flag always
            for (timeIndex = 0; timeIndex < nRecs; timeIndex++)
            {
                nTotal++;
                if (fabsf(QDLAT()) > 50.0)
                {
	            nHighLat++;
                    if (FLAG() == 4)
                    {
                        nFlag4++;
                    }
                }
            }

            // free memory
            for (uint8_t i = 0; i < NUM_DATA_VARIABLES; i++)
            {
                free(dataBuffers[i]);
            }
        }
    }
    closedir(dir);

    fprintf(stdout, "\rSwarm %c\n", satelliteLetter[0]);
    fprintf(stdout, "Total records:\t\t%ld\n", nTotal);
    fprintf(stdout, "High-lat records:\t%ld\n", nFlag4);
    fprintf(stdout, "Good high-lat fraction:\t%.3lf\n", (double)nFlag4 / (double)nHighLat);
    fprintf(stdout, "Good total fraction:\t%.3lf\n", (double)nFlag4 / (double)nTotal);

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
        CDFcloseCDF(calCdfId);
        return;
    }
    long nRecs, memorySize = 0;
    status = CDFgetzVarAllocRecords(calCdfId, CDFgetVarNum(calCdfId, "Timestamp"), &nRecs);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        fprintf(stdout, "%sProblem with data file. Skipping this file.\n", infoHeader);
        CDFcloseCDF(calCdfId);
        return;
    }

    // Variables
    uint8_t nVars = NUM_DATA_VARIABLES;
    char* variables[NUM_DATA_VARIABLES] = {
        "Timestamp",
        "QDLatitude",
        "Quality_flags"
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
            CDFcloseCDF(calCdfId);
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
            CDFcloseCDF(calCdfId);
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
            CDFcloseCDF(calCdfId);
            return;
        }
    }
    // close CDF
    CDFcloseCDF(calCdfId);
    // Update number of records found and memory allocated
    *numberOfRecords = *numberOfRecords + nRecs;

}

uint8_t getMinorVersion(const char *filename)
{
    char minorVersionChar = *(filename + strlen(filename)-5);
    uint8_t minorVersion = (uint8_t) strtol(&minorVersionChar, NULL, 10);
    return minorVersion;
}
