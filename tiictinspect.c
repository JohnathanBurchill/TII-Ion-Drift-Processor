/*

    TII Cross-Track Ion Drift Processor: tiictinspect.c

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
#include "tiictinspect.h"
#include <cdf.h>

// https://www.gnu.org/software/gsl/doc/html/
#include <gsl/gsl_errno.h>
//#include <gsl/gsl_fit.h>
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
    if (argc != 4)
    {
        fprintf(stdout, "usage: %s directory satelliteLetter type\n", argv[0]);
        exit(1);
    }
    const char *directory = argv[1];
    const char *satelliteLetter = argv[2];
    const char *analysisTypeArg = argv[3];

    // Turn off GSL failsafe error handler. We typically check the GSL return codes.
    gsl_set_error_handler_off();


    DIR *dir = opendir(directory);
    struct dirent *entry;
    if (dir == NULL)
    {
        fprintf(stdout, "Could not open directory %s for reading.\n", directory);
        exit(1);
    }

    AnalysisType analysisType = NO_ANALYSIS;
    if (strcmp(analysisTypeArg, "mad") == 0)
        analysisType = MAD_ANALYSIS;
    else if (strcmp(analysisTypeArg, "calibrationgraphic") == 0)
        analysisType = CALIBRATION_ANALYSIS_GRAPHIC;
    else if (strcmp(analysisTypeArg, "calibrationnumeric") == 0)
        analysisType = CALIBRATION_ANALYSIS_NUMERIC;
    else if (strcmp(analysisTypeArg, "timegap") == 0)
        analysisType = CALIBRATION_ANALYSIS_TIMEGAP;
    else if (strcmp(analysisTypeArg, "timegapsummary") == 0)
        analysisType = CALIBRATION_ANALYSIS_TIMEGAPSUMMARY;

    if (analysisType != MAD_ANALYSIS && analysisType != CALIBRATION_ANALYSIS_GRAPHIC && analysisType != CALIBRATION_ANALYSIS_NUMERIC && analysisType != CALIBRATION_ANALYSIS_TIMEGAP && analysisType != CALIBRATION_ANALYSIS_TIMEGAPSUMMARY)
    {
        fprintf(stdout, "Analysis type must be one of \"mad\", \"calibrationgraphic\", \"calibrationnumeric\", \"timegap\", or \"timgapsummary\".\n");
        exit(1);
    }

    char *filename;
    long totalCounts = 0;
    long gapCounts = 0;
    long largeGapCounts = 0;
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
        {
            nFiles++;
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

            uint8_t majorVersion = getMajorVersion(filename);
            uint8_t minorVersion = getMinorVersion(filename);
            switch (analysisType) {
                case MAD_ANALYSIS:
                    analyzeMAD(satellite, majorVersion, minorVersion, fourByteCalFlag, dataBuffers, nRecs);
                    break;
                case CALIBRATION_ANALYSIS_GRAPHIC:
                    analyzeCalibrationFlag(satellite, majorVersion, minorVersion, fourByteCalFlag, dataBuffers, nRecs, false);
                    break;
                case CALIBRATION_ANALYSIS_NUMERIC:
                    analyzeCalibrationFlag(satellite, majorVersion, minorVersion, fourByteCalFlag, dataBuffers, nRecs, true);
                    break;
                case CALIBRATION_ANALYSIS_TIMEGAP:
                    analyzeTimeGaps(filename, dataBuffers, nRecs, &totalCounts, &gapCounts);
                    break;
                case CALIBRATION_ANALYSIS_TIMEGAPSUMMARY:
                    summarizeTimeGaps(filename, dataBuffers, nRecs, &totalCounts, &gapCounts, &largeGapCounts);
                    break;
                case NO_ANALYSIS:
                    break;
            }

            // free memory
            for (uint8_t i = 0; i < NUM_DATA_VARIABLES; i++)
            {
                free(dataBuffers[i]);
            }
        }
    }
    if (analysisType == CALIBRATION_ANALYSIS_TIMEGAP)
    {
        fprintf(stdout, "Time gap analysis: \n  Number of files: %ld\n  Number of points: %ld\n  Number of gaps: %ld (%.3f%%)\n", nFiles, totalCounts, gapCounts, (float)gapCounts / (float)totalCounts * 100.0);
    }
    if (analysisType == CALIBRATION_ANALYSIS_TIMEGAPSUMMARY)
    {
        fprintf(stdout, "Time gap summary: \n  Number of files: %ld\n  Number of points: %ld\n  Number of gaps > 0.5 s: %ld (%.3f%%)\n  Number of gaps > 5.0 s: %ld\n", nFiles, totalCounts, gapCounts, (float)gapCounts / (float)totalCounts * 100.0, largeGapCounts);
    }
    closedir(dir);

}

void analyzeMAD(const char satellite, uint8_t majorVersion, uint8_t minorVersion, bool fourByteCalFlag, uint8_t * dataBuffers[], long nRecs)
{
    // Number of records
    // Include all measurements for Swarm C for which its viy calibration flag is 0
    long nFlag4 = 0, nTotal = 0;
    long long timeIndex;
    for (timeIndex = 0; timeIndex < nRecs; timeIndex++)
    {
        nTotal++;

        if (FLAG() == 4 || (satellite == 'C' && majorVersion == 2 &&  minorVersion == 1) || (satellite == 'C' && fourByteCalFlag == true && (((CALFLAG() >> 16) & 0xff) == 0)))
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
    double meanViyError = 0.;
    for (timeIndex = 0; timeIndex < nRecs; timeIndex++)
    {
        if (FLAG() == 4 || (satellite == 'C' && majorVersion == 2 && minorVersion == 1) || (satellite == 'C' && fourByteCalFlag == true && (((CALFLAG() >> 16) & 0xff) == 0)))
        {
            time[nFlag4] = TIME()/1000.;
            meanTime += time[nFlag4];
            meanViyError += (double)VIYERROR();
            viy4[nFlag4] = (double)VIY();
            nFlag4++;
        }
    }
    meanTime /= (double)nFlag4;
    meanViyError /= (double)nFlag4;
    double viyMad = gsl_stats_mad(viy4, 1, nFlag4, gslWorkspace);
    fprintf(stdout, "%18.6f %7ld %7ld %5.1f %8.1f %8.1f\n", meanTime, nTotal, nFlag4, (float)nFlag4/(float)nTotal*100.0, viyMad, meanViyError);
    free(time);
    free(viy4);
    free(gslWorkspace);

}

void analyzeTimeGaps(const char *filename, uint8_t * dataBuffers[], long nRecs, long *totalCount, long *gapCount)
{
    // Number of records
    // Include all measurements for Swarm C, which are set to 0 flag always

    long nGaps = 0, nTotal = 0;
    long long timeIndex = 0;
    double lastTime = TIME()/1000.0;
    double time;
    double timeDelta = 0.5;
    double firstTime = lastTime;
    char fileDate[9];
    getFileDate(filename, fileDate);

    if (nRecs < 1025)
	return;

    for (timeIndex = 0; timeIndex < nRecs; timeIndex++)
    {
        nTotal++;
        time = TIME()/1000.;	

	if ((time - lastTime) > timeDelta && firstTime > 0.0) 
        {
            fprintf(stdout, "%s %18.5f %lld %18.5f %18.5f %8.2f %8.2f\n", fileDate, firstTime, timeIndex, lastTime, time, lastTime-firstTime, time-firstTime);
            nGaps++;
        }
        lastTime = time;
    }

    fprintf(stdout, "  %s %7ld %7ld %7.3f\n", fileDate, nTotal, nGaps, (float)nGaps/(float)nTotal*100.0);
    *totalCount += nTotal;
    *gapCount += nGaps;
}

void summarizeTimeGaps(const char *filename, uint8_t * dataBuffers[], long nRecs, long *totalCount, long *gapCount, long *largeGapCount)
{
    // Number of records
    // Include all measurements for Swarm C, which are set to 0 flag always

    long nGaps = 0, nTotal = 0, nLargeGaps = 0;
    long long timeIndex = 0;
    double lastTime = TIME()/1000.0;
    double time;
    double timeDelta = 0.5;
    double largeTimeDelta = 5.0;
    double firstTime = lastTime;
    char fileDate[9];
    getFileDate(filename, fileDate);

    if (nRecs < 1025)
	return;

    for (timeIndex = 0; timeIndex < nRecs; timeIndex++)
    {
        nTotal++;
        time = TIME()/1000.;	

	if ((time - lastTime) > timeDelta && firstTime > 0.0) 
        {
            nGaps++;
	    if ((time - lastTime) > largeTimeDelta && firstTime > 0.0) 
	    {
            	nLargeGaps++;
            }
        }
        lastTime = time;
    }

    *totalCount += nTotal;
    *gapCount += nGaps;
    *largeGapCount += nLargeGaps;
}

void analyzeCalibrationFlag(const char satellite, uint8_t majorVersion, uint8_t minorVersion, bool fourByteCalFlag, uint8_t * dataBuffers[], long nRecs, bool numericOutput)
{
    // Number of records
    // Include all measurements for Swarm C, which are set to 0 flag always
    if (fourByteCalFlag == false)
    {
        return;
    }

    long long timeIndex;
    for (timeIndex = 0; timeIndex < nRecs; timeIndex++)
    {
        if (numericOutput == true)
        {
            fprintf(stdout, "%f %ld\n", TIME()/1000., (long) CALFLAG());
        }
        else
        {
            fprintf(stdout, "%f ", TIME()/1000.);
            uint8_t bit = 32;
            for (uint8_t k = 0; k < 4; k++)
            {
                fprintf(stdout, "(");
                for (uint8_t f = 0; f < 8; f++)
                {
                    uint8_t flag = ((CALFLAG() >> --bit) & 1);
                    if (flag == 1)
                        fprintf(stdout, "|");
                    else
                        fprintf(stdout, ".");

                }
                fprintf(stdout, ")");
            }
            fprintf(stdout, "\n");
        }
    }

}

void loadCrossTrackData(const char * filename, uint8_t **dataBuffers, long *numberOfRecords, bool *fourByteCalFlag)
{
    char validationFileName[CDF_PATHNAME_LEN];
    snprintf(validationFileName, strlen(filename)-3, "%s", filename);

    uint8_t majorVersion = getMajorVersion(filename);
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
        fprintf(stdout, "%sCDFinquireCDF: Problem with data file. Skipping this file.\n", infoHeader);
        closeCdf(calCdfId);
        return;
    }
    long nRecs, memorySize = 0;
    status = CDFgetzVarAllocRecords(calCdfId, CDFgetVarNum(calCdfId, "Timestamp"), &nRecs);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        fprintf(stdout, "%sCDFgetzVarAllocRecords: Problem with data file. Skipping this file.\n", infoHeader);
        closeCdf(calCdfId);
        return;
    }

    // Variables
    uint8_t nVars = NUM_DATA_VARIABLES;
    char* variables[NUM_DATA_VARIABLES] = {
        "Timestamp",
        "Viy",
        "Viy_error",
        "QDLatitude",
        "Quality_flags",
        "Calibration_flags"
    };
    if ((majorVersion == 2 && minorVersion == 1))
    {
        variables[3] = "flags";
    }
    for (uint8_t i = 0; (i<nVars-1) || (i == nVars-1 && (majorVersion > 2 || (majorVersion == 2 && minorVersion == 2))); i++)
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

    for (uint8_t i = 0; (i<nVars-1) || (i == nVars-1 && (majorVersion > 2 || (majorVersion == 2 && minorVersion == 2))); i++)
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
        if (i == 5 && numVarBytes == 4)
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

uint8_t getMajorVersion(const char *filename)
{
    char majorVersionChar = *(filename + strlen(filename)-7);
    uint8_t majorVersion = (uint8_t) strtol(&majorVersionChar, NULL, 10);
    return majorVersion;
}

uint8_t getMinorVersion(const char *filename)
{
    char minorVersionChar = *(filename + strlen(filename)-5);
    uint8_t minorVersion = (uint8_t) strtol(&minorVersionChar, NULL, 10);
    return minorVersion;
}

void getFileDate(const char *filename, char *startDate)
{
    strncpy(startDate, (filename + strlen(filename) - 40), 8);
    startDate[8] = '\0';
}


