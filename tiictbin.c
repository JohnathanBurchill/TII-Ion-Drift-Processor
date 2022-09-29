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
#include "statistics.h"
#include <cdf.h>

// https://www.gnu.org/software/gsl/doc/html/
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <time.h>

char infoHeader[50];

int main(int argc, char* argv[])
{
    time_t currentTime;
    struct tm * timeParts;
    time(&currentTime);
    timeParts = gmtime(&currentTime);
    char * dateString = asctime(timeParts);
    char date[255];
    snprintf(date, strlen(dateString), "%s", dateString);

    bool viyToEastward = false;

    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--about") == 0)
        {
            fprintf(stdout, "tiictbin - Calculates and prints requested statistics of specified parameter per QD latitude and MLT bin. Version %s.\n", SOFTWARE_VERSION);
            fprintf(stdout, "Copyright (C) 2022  Johnathan K Burchill\n");
            fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
            fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
            fprintf(stdout, "under the terms of the GNU General Public License.\n");

            exit(0);
        }
        if (strcmp(argv[i], "--viy-to-eastward") == 0)
            viyToEastward = true;

        if (strcmp(argv[i], "--available-statistics") == 0)
        {
            fprintf(stderr, "Available statistics:\n");
            printAvailableStatistics(stderr);
            exit(1);
        }

        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-help") == 0)
        {
            usage(argv[0]);
            exit(1);
        }
    }


    if (argc < 11 || argc > 14)
    {
        usage(argv[0]);
        exit(1);
    }
    const char *directory = argv[1];
    const char *satelliteLetter = argv[2];
    const char *parameterName = argv[3];
    const char *statistic = argv[4];
    if (!validStatistic(statistic))
    {
        fprintf(stderr, "Invalid statistic '%s'\n", statistic);
        fprintf(stderr, "Must be one of:\n");
        printAvailableStatistics(stderr);
        exit(1);
    }

    float qdlatmin = atof(argv[5]);
    float qdlatmax = atof(argv[6]);
    float deltaqdlat = atof(argv[7]);
    float mltmin = atof(argv[8]);
    float mltmax = atof(argv[9]);
    float deltamlt = atof(argv[10]);
    char *firstDate = "20131208";
    char lastDate[255] = {0};
    time_t today = time(NULL);
    timeParts = gmtime(&today);
    sprintf(lastDate, "%4d%02d%02d", timeParts->tm_year + 1900, timeParts->tm_mon + 1, timeParts->tm_mday);
    if (argc > 11 && strncmp(argv[11], "--", 2) != 0)
        firstDate = argv[11];
    if (argc > 12 && strncmp(argv[12], "--", 2) != 0)
        sprintf(lastDate, "%s", argv[12]);

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
    float **binStorage = NULL;
    size_t *binSizes = NULL;
    size_t *binMaxSizes = NULL;
    if (allocateBinStorage(&binStorage, &binSizes, &binMaxSizes, nMLTs, nQDLats, BIN_STORAGE_BLOCK_SIZE))
    {
        fprintf(stderr, "Could not allocate bin storage memory.\n");
        exit(1);
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
        const char satellite = *(filename+length-48);
        if (satellite == satelliteLetter[0] && strncmp(filename + length - 40, firstDate, 8) >=0 && strncmp(filename + length - 40, lastDate, 8) <= 0)
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
    float percentDone = 0.0;
    int percentCheck = (int) ceil(0.1 * (float)nFiles);
    float value = 0.0;

    size_t index = 0;
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
        if (satellite == satelliteLetter[0] && strncmp(filename + length - 40, firstDate, 8) >=0 && strncmp(filename + length - 40, lastDate, 8) <= 0)
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
            status = loadCrossTrackData(fullPath, dataBuffers, &nRecs, &fourByteCalFlag, parameterName);
            if (status != CDF_OK)
            {
                continue;
            }
            double *pEpoch = (double*) dataBuffers[0];
            long long timeIndex = 0;

            // The macros TIME(), QDLAT(), etc. give the value at index timeIndex.
            double epoch0 = TIME()/1000.0;

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
                            value = PARAMETER();
                            if (viyToEastward)
                            {
                                // Flip sign of Viy to make positive viy always eastward as requested 
                                if (VSATN() < 0)
                                {
                                    value = -value;
                                }
                            }
                            index = mltIndex * nQDLats + qdlatIndex;
                            if (binSizes[index] >= binMaxSizes[index])
                            {
                                if(adjustBinStorage(binStorage, binMaxSizes, index, BIN_STORAGE_BLOCK_SIZE))
                                {
                                    fprintf(stderr, "Unable to allocate additional bin storage.\n");
                                    exit(1);
                                }
 
                            }
                            binStorage[index][binSizes[index]] = value;
                            binSizes[index]++;
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
            percentDone = (float)processedFiles / (float)nFiles * 100.0;
            if (processedFiles % percentCheck == 0)
                fprintf(stderr, "Processed %ld of %ld files (%3.0f%%)\n", processedFiles, nFiles, percentDone);
        }
    }
 
    float qdlat = 0.0;
    float mlt = 0.0;
    float result = 0.0;
    long nVals = 0;
    fprintf(stdout, "QDLat\tMLT\t%s(%s)\tCount\n", statistic, parameterName);

    for (size_t q = 0; q < nQDLats; q++)
    {
        for (size_t m = 0; m < nMLTs; m++)
        {
            index = m * nQDLats + q;
            qdlat = qdlatmin + deltaqdlat * ((float)q + 0.5);
            mlt = mltmin + deltamlt* ((float)m + 0.5);
            if (calculateStatistic(statistic, binStorage, binSizes, index, (void*) &result))
                result = GSL_NAN;
            nVals += binSizes[index];
            fprintf(stdout, "%.2f\t%.2f\t%.2f\t%ld\n", qdlat, mlt, result, binSizes[index]);
        }
    }

    fprintf(stderr, "%ld values binned.\n", nVals);

    freeBinStorage(binStorage, binSizes, binMaxSizes, nMLTs, nQDLats);

    closedir(dir);

}

CDFstatus loadCrossTrackData(const char *filename, uint8_t **dataBuffers, long *numberOfRecords, bool *fourByteCalFlag, const char *parameterName)
{
    CDFstatus status = CDF_OK;
    char validationFileName[CDF_PATHNAME_LEN];
    snprintf(validationFileName, strlen(filename)-3, "%s", filename);

    if (numberOfRecords != NULL)
        *numberOfRecords = 0;

    uint8_t minorVersion = getMinorVersion(filename);

    // Open the CDF file with validation
    CDFsetValidate(VALIDATEFILEon);
    CDFid calCdfId;
    status = CDFopenCDF(validationFileName, &calCdfId);
    if (status != CDF_OK) 
    {
        printErrorMessage(status);
        // Not necessarily an error. For example, some dates will have not calibration data.
        fprintf(stderr, "%sSkipping this file.\n", infoHeader);
        return status;
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
        return status;
    }
    long nRecs, memorySize = 0;
    status = CDFgetzVarAllocRecords(calCdfId, CDFgetVarNum(calCdfId, "Timestamp"), &nRecs);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        fprintf(stdout, "%sProblem with data file. Skipping this file.\n", infoHeader);
        closeCdf(calCdfId);
        return status;
    }

    // Variables
    uint8_t nVars = NUM_DATA_VARIABLES;
    char* variables[NUM_DATA_VARIABLES] = {
        "Timestamp",
        "MLT",
        "QDLatitude",
        "Quality_flags",
        "Calibration_flags",
        "VsatN"
        ""
    };
    variables[NUM_DATA_VARIABLES-1] = (char *) parameterName;
    if (minorVersion == 1)
    {
        variables[3] = "flags";
    }
    for (uint8_t i = 0; (i<nVars-1) || (i == nVars-1 && minorVersion == 2); i++)
    {
        status = CDFconfirmzVarExistence(calCdfId, variables[i]);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            fprintf(stdout, "%sError reading variable %s. Skipping this file.\n", infoHeader, variables[i]);
            closeCdf(calCdfId);
            return status;
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
            return status;
        }
        status = CDFgetzVarNumDims(calCdfId, varNum, &numDims);
        status = CDFgetzVarDimSizes(calCdfId, varNum, dimSizes);
        status = CDFgetzVarDataType(calCdfId, varNum, &dataType);
        // Calculate new size of memory to allocate
        status = CDFgetDataTypeSize(dataType, &numVarBytes);
        if (i == 4 && numVarBytes == 4)
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
            return status;
        }
    }
    // close CDF
    closeCdf(calCdfId);
    // Update number of records found and memory allocated
    if (numberOfRecords != NULL)
        *numberOfRecords = nRecs;

    return status;
}

uint8_t getMinorVersion(const char *filename)
{
    char minorVersionChar = *(filename + strlen(filename)-5);
    uint8_t minorVersion = (uint8_t) strtol(&minorVersionChar, NULL, 10);
    return minorVersion;
}

void usage(char *name)
{
    fprintf(stdout, "usage: %s directory satelliteLetter parameterName statistic qdlatmin qdlatmax deltaqdlat mltmin mltmax deltamlt [firstDate] [lastDate] [--viy-to-eastward] [--help] [--about]\n", name);
    fprintf(stdout, "Options:\n");
    fprintf(stdout, "\t--help or -h\t\tprints this message.\n");
    fprintf(stdout, "\t--about \t\tdescribes the program, declares license.\n");
    fprintf(stdout, "\t--available-statistics\tprints a list of statistics to calculate. Pass one statistic per call.\n");
    fprintf(stdout, "\t--viy-to-eastward\tflips sign of Viy for descending part of the orbit so that positive ion drift is always eastward.\n");

    return;    
}
