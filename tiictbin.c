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
    int32_t qualityFlagMask = 0;
    bool qualityMaskIsAnd = true;

    int nOptions = 0;

    bool showFileProgress = true;

    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--about") == 0)
        {
            fprintf(stdout, "tiictbin - Calculates and prints requested statistics of specified parameter per QD latitude and MLT bin. Version %s.\n", SOFTWARE_VERSION);
            fprintf(stdout, "Copyright (C) 2022  Johnathan K Burchill\n");
            fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
            fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
            fprintf(stdout, "under the terms of the GNU General Public License.\n");

            exit(EXIT_SUCCESS);
        }
        else if (strcmp(argv[i], "--available-statistics") == 0)
        {
            fprintf(stderr, "Available statistics:\n");
            printAvailableStatistics(stderr);
            exit(EXIT_FAILURE);
        }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-help") == 0)
        {
            usage(argv[0]);
            exit(EXIT_FAILURE);
        }
        else if (strcmp(argv[i], "--no-file-progress") == 0)
        {
            nOptions++;
            showFileProgress = false;
        }
        else if (strcmp(argv[i], "--viy-to-eastward") == 0)
        {
            nOptions++;
            viyToEastward = true;
        }
        else if (strncmp(argv[i], "--qualityflagmask=", 18) == 0)
        {
            nOptions++;
            int base = 10;
            int sign = 1;
            size_t offset = 18;
            size_t len = strlen(argv[i]);
            if (len < 19)
            {
                fprintf(stderr, "Invalid quality flag mask value.\n");
                exit(EXIT_FAILURE);
            }
            if ((argv[i] + offset)[0] == '-')
            {
                sign = -1;
                offset++;
            }
            if (len > 20)
            {
                if (strncmp(argv[i] + offset, "0b", 2) == 0)
                {
                    offset += 2;
                    base = 2;
                }
                else if (strncmp(argv[i] + offset, "0x", 2) == 0)
                {
                    offset += 2;
                    base = 16;
                }
            }
            qualityFlagMask = (int32_t) strtol(argv[i] + offset, (char **)NULL, base);
            qualityFlagMask *= sign;
        }
        else if (strncmp(argv[i], "--qualityflagmasktype=", 22) == 0)
        {
            nOptions++;
            if (strcmp(argv[i] + 22, "OR") == 0)
                qualityMaskIsAnd = false;
            else
                qualityMaskIsAnd = true;
        }
        else if (strncmp(argv[i], "--", 2) == 0)
        {
            fprintf(stderr, "Unknown option %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }

    }


    if (argc != 11 + nOptions)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
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
        exit(EXIT_FAILURE);
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
        exit(EXIT_FAILURE);
    }
    if (deltamlt <= 0.0)
    {
        fprintf(stderr, "%s: deltamlt must be greater than 0.0\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    int nQDLats = (int)floor((qdlatmax - qdlatmin) / deltaqdlat);
    int nMLTs = (int)floor((mltmax - mltmin) / deltamlt);
    if (nQDLats <= 0)
    {
        fprintf(stderr, "%s: invalid QD latitude bin specification.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (nMLTs <= 0)
    {
        fprintf(stderr, "%s: invalid MLT bin specification.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // Flag mask
    char *flagParams[4] = {"Vixh", "Vixv", "Viy", "Viz"};
    if (qualityFlagMask == 0)
    {
        fprintf(stderr, "Quality mask: selecting %s irrespective of data quality.\n", parameterName);
    }
    else
    {
        if (qualityFlagMask > 0)
            fprintf(stderr, "Quality mask: INCLUDING %s for good values of ", parameterName);
        else
            fprintf(stderr, "Quality mask: EXCLUDING %s for good values of ", parameterName);

        int maskBit = 0;
        int nFlaggedParams = 0;
        int unsignedMask = abs(qualityFlagMask);
        for (int b = 0; b < 4; b++)
        {
            maskBit = (unsignedMask >> b) & 0x1;
            if (maskBit)
            {
                nFlaggedParams++;
                if (nFlaggedParams > 1)
                    fprintf(stderr, "%s", qualityMaskIsAnd ? " AND " : " OR ");
                fprintf(stderr, "%s", flagParams[b]);
            }
        }
        fprintf(stderr, "\n");                

    }

    uint16_t positiveQualityFlagMask = (uint16_t) abs(qualityFlagMask);

    // Allocate memory for binning
    float **binStorage = NULL;
    size_t *binSizes = NULL;
    size_t *binMaxSizes = NULL;
    if (allocateBinStorage(&binStorage, &binSizes, &binMaxSizes, nMLTs, nQDLats, BIN_STORAGE_BLOCK_SIZE))
    {
        fprintf(stderr, "Could not allocate bin storage memory.\n");
        exit(EXIT_FAILURE);
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
        exit(EXIT_FAILURE);
    }
    char *filename;
    long nFiles = 0;

    long nValsRead = 0;
    long nValsWithinBinLimits = 0;
    long nValsBinned = 0;

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
        exit(EXIT_FAILURE);
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

            bool includeValue = false;
            uint16_t maskedValue = 0;

            for (timeIndex = 0; timeIndex < nRecs; timeIndex++)
            {
                nValsRead++;

                uint16_t flag = FLAG();
                maskedValue = (FLAG() & positiveQualityFlagMask);
                // Is flag matched by mask?
                if (positiveQualityFlagMask == 0)
                    includeValue = true;
                else if (qualityMaskIsAnd)
                    includeValue = (maskedValue == positiveQualityFlagMask);
                else
                    includeValue = (maskedValue > 0);
                
                // Exclude this match?
                if (qualityFlagMask < 0)
                    includeValue = !includeValue;

                // Access bins with bins[mltIndex * nQDLats + qdlatIndex];
                if (isfinite(PARAMETER()))
                {
                    mltIndex = (int) floor((MLT() - mltmin) / deltamlt);
                    qdlatIndex = (int) floor((QDLAT() - qdlatmin) / deltaqdlat);
                    if (mltIndex >= 0 && mltIndex < nMLTs && qdlatIndex >=0 && qdlatIndex < nQDLats)
                    {
                        nValsWithinBinLimits++;
                        if (includeValue)
                        {
                            nValsBinned++;
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
                                    exit(EXIT_FAILURE);
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
            if (showFileProgress)
            {
                percentDone = (float)processedFiles / (float)nFiles * 100.0;
                if (processedFiles % percentCheck == 0)
                    fprintf(stderr, "Processed %ld of %ld files (%3.0f%%)\n", processedFiles, nFiles, percentDone);
            }
        }
    }
 
    float qdlat = 0.0;
    float mlt = 0.0;
    float result = 0.0;
    char statExpression[50];
    snprintf(statExpression, 50, "%s(%s)", statistic, parameterName);
    fprintf(stdout, "  MLT\t QDLat\t%35s\t       Count\n", statExpression);

    for (size_t q = 0; q < nQDLats; q++)
    {
        for (size_t m = 0; m < nMLTs; m++)
        {
            index = m * nQDLats + q;
            qdlat = qdlatmin + deltaqdlat * ((float)q + 0.5);
            mlt = mltmin + deltamlt* ((float)m + 0.5);
            if (calculateStatistic(statistic, binStorage, binSizes, index, (void*) &result))
                result = GSL_NAN;
            fprintf(stdout, "%5.2f\t%6.2f\t%35.2f\t%12ld\n", mlt, qdlat, result, binSizes[index]);
        }
    }

    fprintf(stderr, "Values read: %12ld\tValues within bin limits: %12ld\tValues binned: %12ld (%6.2lf%% of those within bin limits)\n", nValsRead, nValsWithinBinLimits, nValsBinned, 100.0 * (double)nValsBinned / (double)nValsWithinBinLimits);

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
    for (uint8_t i = 0; i < NUM_DATA_VARIABLES; i++)
    {
        status = CDFconfirmzVarExistence(calCdfId, variables[i]);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            fprintf(stdout, "%sError reading variable %s. Skipping this file.\n", infoHeader, variables[i]);
            fflush(stdout);
            closeCdf(calCdfId);
            return status;
        }
    }
    
    long varNum, numValues, numVarBytes;
    long numBytesPrev, numBytesToAdd, numBytesNew;

    for (uint8_t i = 0; i<NUM_DATA_VARIABLES; i++)
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
    fprintf(stdout, "usage: %s directory satelliteLetter parameterName statistic qdlatmin qdlatmax deltaqdlat mltmin mltmax deltamlt [firstDate] [lastDate] [--viy-to-eastward] [--qualityflagmask=mask] [--qualityflagmasktype=type] [--no-file-progress] [--help] [--about]\n", name);
    fprintf(stdout, "Options:\n");
    fprintf(stdout, "\t--help or -h\t\tprints this message.\n");
    fprintf(stdout, "\t--about \t\tdescribes the program, declares license.\n");
    fprintf(stdout, "\t--available-statistics\tprints a list of statistics to calculate. Pass one statistic per call.\n");
    fprintf(stdout, "\t--viy-to-eastward\tflips sign of Viy for descending part of the orbit so that positive ion drift is always eastward.\n");
    fprintf(stdout, "\t--qualityflagmask=value\tselects (mask > 0) or rejects (mask < 0) measurements with quality flag bitwise-and-matching abs(mask) according to the mask type given by --qualityflagmasktype, e.g., --qualityflagmask=0b0110 or --qualityflagmask=-15\n");
    fprintf(stdout, "\t--qualityflagmasktype={AND|OR}\tinterpret --qualityflagmask values as bitwise AND or OR\n");
    fprintf(stdout, "\t--no-file-progress\tdo not print progress of files being processed\n");

    return;    
}
