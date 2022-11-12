/*

    TII Cross-Track Ion Drift Processor: cdfbin.c

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

#include "cdfbin.h"
#include "utilities.h"
#include "statistics.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <dirent.h>

#include <cdf.h>

// https://www.gnu.org/software/gsl/doc/html/
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <time.h>

char infoHeader[50];

static inline float solidAngle(float lat1, float lat2, float mlt1, float mlt2)
{
    // lat2 is more northward than lat1
    return (mlt2 - mlt1) * (cosf((90.0 - lat2) * M_PI / 180.0) - cosf((90.0 - lat1) * M_PI / 180.0));
}

static inline int nRingBins(float ringSolidAngle, float solidAngleUnit)
{
    return (int) fabsf(roundf(ringSolidAngle / solidAngleUnit));
} 

int main(int argc, char* argv[])
{
    time_t currentTime;
    struct tm * timeParts;
    time(&currentTime);
    timeParts = gmtime(&currentTime);
    char * dateString = asctime(timeParts);
    char date[255];
    snprintf(date, strlen(dateString), "%s", dateString);

    bool flipParamWhenDescending = false;
    int32_t qualityFlagMask = 0;
    bool qualityMaskIsAnd = true;
    bool goodFlagIsOne = true;

    int nOptions = 0;

    bool showFileProgress = true;

    char firstDate[255] = {0};
    snprintf(firstDate, 9, "%s", "20131208");
    char lastDate[255] = {0};
    time_t today = time(NULL);
    timeParts = gmtime(&today);
    sprintf(lastDate, "%4d%02d%02d", timeParts->tm_year + 1900, timeParts->tm_mon + 1, timeParts->tm_mday);

    bool useEqualArea = false;
    bool tctData = false;

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
            fprintf(stdout, "Available statistics:\n");
            printAvailableStatistics(stdout);
            exit(EXIT_FAILURE);
        }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-help") == 0)
        {
            usage(argv[0]);
            exit(EXIT_FAILURE);
        }
        else if (strncmp(argv[i], "--first-date=", 13) == 0)
        {
            if (strlen(argv[i]) != 21)
            {
                fprintf(stderr, "%s: invalid date\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            nOptions++;
            snprintf(firstDate, 9, "%s", argv[i] + 13);
        }
        else if (strncmp(argv[i], "--last-date=", 12) == 0)
        {
            if (strlen(argv[i]) != 20)
            {
                fprintf(stderr, "%s: invalid date\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            nOptions++;
            snprintf(lastDate, 9, "%s", argv[i] + 12);
        }
        else if (strcmp(argv[i], "--no-file-progress") == 0)
        {
            nOptions++;
            showFileProgress = false;
        }
        else if (strcmp(argv[i], "--equal-area-bins") == 0)
        {
            nOptions++;
            useEqualArea = true;
        }
        else if (strcmp(argv[i], "--flip-when-descending") == 0)
        {
            nOptions++;
            flipParamWhenDescending = true;
        }
        else if (strncmp(argv[i], "--quality-flag-mask=", 20) == 0)
        {
            nOptions++;
            int base = 10;
            int sign = 1;
            size_t offset = 20;
            size_t len = strlen(argv[i]);
            if (len < 21)
            {
                fprintf(stderr, "Invalid quality flag mask value.\n");
                exit(EXIT_FAILURE);
            }
            if ((argv[i] + offset)[0] == '-')
            {
                sign = -1;
                offset++;
            }
            if (len > 22)
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
        else if (strncmp(argv[i], "--quality-flag-mask-type=", 25) == 0)
        {
            nOptions++;
            if (strcmp(argv[i] + 25, "OR") == 0)
                qualityMaskIsAnd = false;
            else
                qualityMaskIsAnd = true;
        }
        else if (strcmp(argv[i], "--quality-flag-zero-is-good") == 0)
        {
            nOptions++;
            goodFlagIsOne = false;
        }
        else if (strcmp(argv[i], "--tct-data") == 0)
        {
            nOptions++;
            tctData = true;
        }
        else if (strncmp(argv[i], "--", 2) == 0)
        {
            fprintf(stderr, "Unknown or incomplete option %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    if (argc != 12 + nOptions)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    const char *directory = argv[1];
    const char *satelliteLetter = argv[2];
    const char *parameterName = argv[3];
    const char *flagName = argv[4];
    const char *statistic = argv[5];
    if (!validStatistic(statistic))
    {
        fprintf(stderr, "Invalid statistic '%s'\n", statistic);
        fprintf(stderr, "Must be one of:\n");
        printAvailableStatistics(stderr);
        exit(EXIT_FAILURE);
    }

    float qdlatmin = atof(argv[6]);
    float qdlatmax = atof(argv[7]);
    float deltaqdlat = atof(argv[8]);
    float mltmin = atof(argv[9]);
    float mltmax = atof(argv[10]);
    float deltamlt = atof(argv[11]);
    int qdlatIndex = 0;
    int mltIndex = 0;

    int nQDLats = (int)floor((qdlatmax - qdlatmin) / deltaqdlat);
    int nMLTs = (int)floor((mltmax - mltmin) / deltamlt);
    if (nQDLats <= 0)
    {
        fprintf(stderr, "%s: invalid QD latitude bin specification.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (nMLTs <= 0 && !useEqualArea)
    {
        fprintf(stderr, "%s: invalid MLT bin specification.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    fprintf(stdout, "Analyzing Swarm %s files between %s and %s\n", satelliteLetter, firstDate, lastDate);

    // Flag mask
    if (tctData)
    {
        char *flagParams[4] = {"Vixh", "Vixv", "Viy", "Viz"};
        if (qualityFlagMask == 0)
        {
            fprintf(stdout, "Quality mask: selecting %s irrespective of data quality.\n", parameterName);
        }
        else
        {
            if (qualityFlagMask > 0)
                fprintf(stdout, "Quality mask: INCLUDING %s for good values of ", parameterName);
            else
                fprintf(stdout, "Quality mask: EXCLUDING %s for good values of ", parameterName);

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
                        fprintf(stdout, "%s", qualityMaskIsAnd ? " AND " : " OR ");
                    fprintf(stdout, "%s", flagParams[b]);
                }
            }
            fprintf(stdout, "\n");                

        }
    }

    uint32_t positiveQualityFlagMask = (uint32_t) abs(qualityFlagMask);

    // Allocate memory for binning
    float **binStorage = NULL;
    size_t *binSizes = NULL;
    size_t *binValidSizes = NULL;
    size_t *binMaxSizes = NULL;

    size_t nBins = 0;

    int *nMltsVsLatitude = (int*) calloc(nQDLats, sizeof(int));
    int *cumulativeMltsVsLatitude = (int*) calloc(nQDLats, sizeof(int)); 
    if (nMltsVsLatitude == NULL || cumulativeMltsVsLatitude == NULL)
    {
        fprintf(stderr, "Unable to allocate memory for nMltsVsLatitude[].\n");
        exit(EXIT_FAILURE);
    }

    // Approximate bin area
    // Will adjust deltaMlt to get bin areas closest to this
    float solidAngleUnit = solidAngle(90.0 - deltaqdlat, 90.0, mltmin, mltmax) / (float)nMLTs;
    float ringSolidAngle = 0.0;

    float lat = 0.0;
    // q = 0 is middle of lowest latitude bin
    cumulativeMltsVsLatitude[0] = 0;
    for (int q = 0; q < nQDLats; q++)
    {
        if (!useEqualArea)
            nMltsVsLatitude[q] = nMLTs;
        else
        {
            ringSolidAngle = solidAngle(qdlatmin + q * deltaqdlat, qdlatmin + (q+1) * deltaqdlat, mltmin, mltmax);
            nMltsVsLatitude[q] = nRingBins(ringSolidAngle, solidAngleUnit);
        }        
        nBins += nMltsVsLatitude[q];
        if (q > 0)
            cumulativeMltsVsLatitude[q] = cumulativeMltsVsLatitude[q-1] + nMltsVsLatitude[q - 1];
    }

    if (allocateBinStorage(&binStorage, &binSizes, &binValidSizes, &binMaxSizes, nBins, BIN_STORAGE_BLOCK_SIZE, useEqualArea))
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
    int percentCheck = (int) ceil(0.01 * (float)nFiles);
    float mlt = 0.0;
    float qdlat = 0.0;
    float value = 0.0;
    bool includeValue = false;
    uint32_t maskedValue = 0;
    uint32_t flag;

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
            long variableTypes[NUM_DATA_VARIABLES] = {0};
            sprintf(fullPath, "%s/%s", directory, filename);
            status = loadCdfData(fullPath, dataBuffers, &nRecs, parameterName, flagName, variableTypes);
            if (status != CDF_OK)
                continue;

            double *pEpoch = (double*) dataBuffers[0];
            long long timeIndex = 0;

            // The macros TIME(), QDLAT(), etc. give the value at index timeIndex.
            double epoch0 = TIME()/1000.0;

            timeIndex = 0;
            qdlat =  variableTypes[QDLAT_INDEX] == CDF_REAL8 ? (float)QDLAT8() : QDLAT();
            float lastQDLat = qdlat;
            float dir = 0.0;
            if (nRecs > 2)
            {
                timeIndex = 1;
                qdlat =  variableTypes[QDLAT_INDEX] == CDF_REAL8 ? (float)QDLAT8() : QDLAT();
                dir = qdlat - lastQDLat;
            }
            for (timeIndex = 0; timeIndex < nRecs; timeIndex++)
            {
                nValsRead++;

                switch(variableTypes[FLAG_INDEX])
                {
                    case CDF_UINT1:
                    case CDF_BYTE:
                    case CDF_UCHAR:
                        flag = (uint32_t)*((uint8_t*)dataBuffers[3] + (size_t)timeIndex);
                        break;

                    case CDF_UINT2:
                        flag = (uint32_t)*((uint16_t*)dataBuffers[3] + (size_t)timeIndex);
                        break;

                    case CDF_CHAR:
                    case CDF_INT1:
                        flag = (uint32_t)*((int8_t*)dataBuffers[3] + (size_t)timeIndex);
                        break;
    
                    case CDF_INT2:
                        flag = (uint32_t)*((int16_t*)dataBuffers[3] + (size_t)timeIndex);
                        break;

                    case CDF_INT4:
                        flag = (uint32_t)*((int32_t*)dataBuffers[3] + (size_t)timeIndex);
                        break;

                    default:
                        flag = *((uint32_t*)dataBuffers[3] + (size_t)timeIndex);
                        break;
                }

                maskedValue = (flag & positiveQualityFlagMask);
                // Is flag matched by mask?
                if (positiveQualityFlagMask == 0)
                    includeValue = true;
                else if (qualityMaskIsAnd)
                {
                    if (goodFlagIsOne)
                        includeValue = (maskedValue == positiveQualityFlagMask);
                    else
                        includeValue = (maskedValue == 0);
                }
                else
                {
                    if (goodFlagIsOne)
                        includeValue = (maskedValue > 0);
                    else
                        includeValue = (maskedValue < positiveQualityFlagMask);
                }
                
                // Exclude this match?
                if (qualityFlagMask < 0)
                    includeValue = !includeValue;

                qdlat =  variableTypes[QDLAT_INDEX] == CDF_REAL8 ? (float)QDLAT8() : QDLAT();
                dir = qdlat - lastQDLat;
                lastQDLat = qdlat;
                // Access bins as bins[cumulativeMltsVsLatitude[qdlatIndex] + mltIndex];
                // TODO handle vector parameters
                value = variableTypes[PARAMETER_INDEX] == CDF_REAL8 ? (float)PARAMETER8() : PARAMETER();
                if (isfinite(value))
                {
                    qdlatIndex = (int) floor((qdlat - qdlatmin) / deltaqdlat);
                    if (qdlatIndex < 0 || qdlatIndex >= nQDLats)
                        continue;
                    deltamlt = (mltmax - mltmin) / (float)nMltsVsLatitude[qdlatIndex];
                    mlt = variableTypes[MLT_INDEX] == CDF_REAL8 ? (float)MLT8() : MLT();
                    mltIndex = (int) floor((mlt - mltmin) / deltamlt);
                    if (mltIndex >= 0 && mltIndex < nMltsVsLatitude[qdlatIndex])
                    {
                        index = cumulativeMltsVsLatitude[qdlatIndex] + mltIndex;
                        // Measurement lies within a QDLat and MLT bin
                        binValidSizes[index]++;
                        nValsWithinBinLimits++;

                        if (includeValue)
                        {
                            // Flip sign of parameter when moving southward, i.e. to make Viy eastward and Vixh or Vixv northward

                            if (flipParamWhenDescending && dir < 0.0)
                                value = -value;

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
                            nValsBinned++;
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
                    fprintf(stderr, "\rProcessed %ld of %ld files (%3.0f%%)", processedFiles, nFiles, percentDone);
            }
        }
    }

    if (showFileProgress)
        fprintf(stderr, "\r\n");

    float qdlat1 = 0.0;
    float qdlat2 = 0.0;
    float mlt1 = 0.0;
    float mlt2 = 0.0;
    float result = 0.0;

    fprintf(stdout, "Time range is inclusive. Bin specification for remaining quantities x and bin boundaries x1 and x2: x1 <= x < x2\n");
    fprintf(stdout, "Row legend:\n");
    fprintf(stdout, "firstDate lastDate MLT1 MLT2 QDLat1 QDLat2 %s(%s) binCount validRegionFraction totalReadFraction\n", statistic, parameterName);

    float denomBinValidSizes = 0.0;
    float denomNValsRead = nValsRead > 0 ? (float) nValsRead : 1.0;

    for (size_t q = 0; q < nQDLats; q++)
    {
        for (size_t m = 0; m < nMltsVsLatitude[q]; m++)
        {
            index = cumulativeMltsVsLatitude[q] + m;
            qdlat1 = qdlatmin + deltaqdlat * ((float)q);
            qdlat2 = qdlat1 + deltaqdlat;
            deltamlt = (mltmax - mltmin) / (float)nMltsVsLatitude[q];
            mlt1 = mltmin + deltamlt * (float)m;
            mlt2 = mlt1 + deltamlt;
            if (calculateStatistic(statistic, binStorage, binSizes, index, (void*) &result))
                result = GSL_NAN;

            denomBinValidSizes = binValidSizes[index] > 0 ? (float) binValidSizes[index] : 1.0;

            fprintf(stdout, "%8s %8s %5.2f %5.2f %6.2f %6.2f %f %ld %f %f\n", firstDate, lastDate, mlt1, mlt2, qdlat1, qdlat2, result, binSizes[index], (float)binSizes[index] / denomBinValidSizes, (float)binSizes[index] / denomNValsRead);
        }
    }

    fprintf(stdout, "Summary of counts\n");
    fprintf(stdout, "\tValues read: %ld; Values within bin limits: %ld; Values binned: %ld (%6.2lf%% of those within bin limits)\n", nValsRead, nValsWithinBinLimits, nValsBinned, 100.0 * (double)nValsBinned / (double)nValsWithinBinLimits);

    freeBinStorage(binStorage, binSizes, binValidSizes, binMaxSizes, nBins);
    if (nMltsVsLatitude != NULL)
        free(nMltsVsLatitude);
    if (cumulativeMltsVsLatitude != NULL)
        free(cumulativeMltsVsLatitude);

    closedir(dir);

}

CDFstatus loadCdfData(const char *filename, uint8_t **dataBuffers, long *numberOfRecords, const char *parameterName, const char *flagVarName, long *variableTypes)
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
    const char* variables[NUM_DATA_VARIABLES] = {
        "Timestamp",
        "MLT",
        "QDLatitude",
        flagVarName,
        parameterName
    };
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
        varNum = CDFgetVarNum(calCdfId, (char*)variables[i]);
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
        if (variableTypes != 0)
            variableTypes[i] = dataType;

        // Calculate new size of memory to allocate
        status = CDFgetDataTypeSize(dataType, &numVarBytes);

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
    fprintf(stdout, "usage: %s <directory> <satelliteLetter> <parameterName> <qualityFlagName> <statistic> <qdlatmin> <qdlatmax> <deltaqdlat> <mltmin> <mltmax> <deltamlt> [--first-date=yyyymmdd] [--last-date=yyyymmdd] [--equal-area-bins] [--flip-when-descending] [--quality-flag-mask=mask] [--quality-flag-mask-type=type] [--quality-flag-zero-is-good] [--tct-data] [--no-file-progress] [--help] [--about]\n", name);
    fprintf(stdout, "Options:\n");
    fprintf(stdout, "\t--help or -h\t\tprints this message.\n");
    fprintf(stdout, "\t--about \t\tdescribes the program, declares license.\n");
    fprintf(stdout, "\t--available-statistics\tprints a list of statistics to calculate. Pass one statistic per call.\n");
    fprintf(stdout, "\t--first-date=yyyymmdd\tFirst date to include in statistics.\n");
    fprintf(stdout, "\t--last-date=yyyymmdd\tLast date to include in statistics.\n");
    fprintf(stdout, "\t--flip-when-descending\tflips sign of Viy for descending part of the orbit so that positive ion drift is always eastward.\n");
    fprintf(stdout, "\t--quality-flag-mask=value\tselects (mask > 0) or rejects (mask < 0) measurements with quality flag bitwise-and-matching abs(mask) according to the mask type given by --quality-flag-mask-type, e.g., --quality-flag-mask=0b0110 or --quality-flag-mask=-15.\n");
    fprintf(stdout, "\t--quality-flag-mask-type={AND|OR}\tinterpret --qualityflagmask values as bitwise AND or OR.\n");
    fprintf(stdout, "\t--quality-flag-zero-is-good\tnon-zero quality flag value signifies an issue.\n");
    fprintf(stdout, "\t--no-file-progress\tdo not print progress of files being processed.\n");
    fprintf(stdout, "\t--equal-area-bins\tgenerate an equal-area grid centered on the magnetic pole. In this case deltaMlt determines the number of MLTs in a polar cap with half-angle deltaqdlat spanning mltmin to mltmax.\n");
    fprintf(stdout, "\t--tct-data\tTII cross-track ion drift data: print extra flag information.\n");
   
    return;    
}


