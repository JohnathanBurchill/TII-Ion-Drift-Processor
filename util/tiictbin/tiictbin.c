/*

    TIICT Processor: util/tiictbin/tiictbin.c

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

// Equal-area bin statistics adapted from spaceseries
// https://github.com/JohnathanBurchill/spaceseries

#include "tiictbin.h"
#include "statistics.h"

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <fts.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <libgen.h>
#include <errno.h>

#include <cdf.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#define N_QUALITY_FLAG_BITS 4

static char* qualityFlagInfo[N_QUALITY_FLAG_BITS] = {
    "Vixh data quality OK",
    "Vixv data quality OK",
    "Viy data quality OK",
    "Viz data quality OK",
};

int main(int argc, char *argv[])
{

    int status = 0;

    ProcessingParameters params = {0};

    params.satelliteLetter = 'X';
    params.verbose = false;
    params.showFileProgress = true;
    params.cdfDirectory = ".";
    params.inputFile = NULL;

    // Defaults to equal-area binning
    params.binningState.equalArea = true;
    params.binningState.qdlatmin = 50.0;
    params.binningState.qdlatmax = 90.0;
    params.binningState.deltaqdlat = 5.0;
    params.binningState.mltmin = 0.0;
    params.binningState.mltmax = 24.0;
    params.binningState.deltamlt = 8.0;
    params.binningState.flipParamWhenDescending = false;
    params.binningState.flipParamWhenDawnward = false;
    params.binningState.flipParamWhenSunward = false;

    // check options and arguments
    parseCommandLine(&params, argc, argv);

    status = initBinningState(&params.binningState);
    if (status != BIN_OK)
        exit(EXIT_FAILURE);

    // Turn off GSL failsafe error handler. Check the GSL return codes.
    gsl_set_error_handler_off();

    long processedFiles = 0;
    float percentDone = 0.0;

    void *mem = NULL;

    // Analyze files in requested directory
    char *dir[2] = {params.cdfDirectory, NULL};
    long nFiles = 0;
    char *filename = NULL;
    char fullpath[FILENAME_MAX + 1] = {0};
    char *cwdresult = getcwd(fullpath, FILENAME_MAX);
    size_t pathlen = strlen(fullpath);
    if (cwdresult == NULL)
    {
        fprintf(stderr, "Unable to get current directory.\n");
        exit(EXIT_FAILURE);
    }

    // Count files
    FTS *f = fts_open(dir, FTS_PHYSICAL | FTS_NOSTAT, NULL);
    FTSENT *e = fts_read(f);
    while (e != NULL)
    {
        if (fileMatch(e, &params) == true)
            nFiles++;
        e = fts_read(f);
    }
    fts_close(f);

    int percentCheck = (int) ceil(0.01 * (float)nFiles);

    // Reopen directory listing to do the analysis
    f = fts_open(dir, FTS_PHYSICAL | FTS_NOSTAT, NULL);
    e = fts_read(f);
    while (e != NULL)
    {
        if (fileMatch(e, &params) == true)
        {
            if (params.verbose)
                fprintf(stderr, "\nAnalyzing %s\n", e->fts_name);

            // snprintf(fullpath + pathlen, FILENAME_MAX - pathlen, "/%s", e->fts_path);
            params.inputFile = e->fts_path;
            processFile(&params);
            processedFiles++;

            if (params.showFileProgress)
            {
                percentDone = (float)processedFiles / (float)nFiles * 100.0;
                if (processedFiles % percentCheck == 0)
                    fprintf(stderr, "\r%s: %ld of %ld files processed (%3.0f%%)", argv[0], processedFiles, nFiles, percentDone);
            }
        }

        e = fts_read(f);
    }
    fts_close(f);

    if (params.showFileProgress)
        fprintf(stderr, "\r\n");

    printBinningResults(&params.binningState, params.parameter, params.statistic);

    freeBinStorage(&params.binningState);

    return EXIT_SUCCESS;
}

void usage(char *name)
{
    fprintf(stdout, "usage: %s <satLetter> <parameter> <statistic> <startDate> <stopDate>\n", name);
    fprintf(stdout, "%35s - %s\n", "--help", "print this message");
    fprintf(stdout, "%35s - %s\n", "--about", "print program and license info");
    fprintf(stdout, "%35s - %s\n", "--verbose", "extra processing information");
    fprintf(stdout, "%35s - %s\n", "--available-statistics", "print list of supported statistics");
    fprintf(stdout, "%35s - %s\n", "--no-file-progress", "suppress printing file progress");
    fprintf(stdout, "%35s - %s\n", "--equal-length-bins", "use standard binning rather than equal-area");
    fprintf(stdout, "%35s - %s\n", "--qdlatmin=<value>", "minimum quasi-dipole magnetic latitude");
    fprintf(stdout, "%35s - %s\n", "--qdlatmax=<value>", "maximum quasi-dipole magnetic latitude");
    fprintf(stdout, "%35s - %s\n", "--deltaqdlat=<value>", "quasi-dipole magnetic latitude bin width");
    fprintf(stdout, "%35s - %s\n", "--mltmin=<value>", "minimum magnetic local time");
    fprintf(stdout, "%35s - %s\n", "--mltmax=<value>", "maximum magnetic local time");
    fprintf(stdout, "%35s - %s\n", "--deltamlt=<value>", "magnetic local time bin width (at the polar cap if for equal-area binning)");
    fprintf(stdout, "%35s - %s\n", "--flip-when-descending", "change sign of value when satellite is moving magnetically southward");
    fprintf(stdout, "%35s - %s\n", "--flip-when-sunward", "change sign of value when satellite is moving magnetically sunward");
    fprintf(stdout, "%35s - %s\n", "--flip-when-dawnward", "change sign of value when satellite is moving magnetically dawnward");
    fprintf(stdout, "%35s - %s\n", "--cdf-input-directory=<dir>", "path to directory containing binary input files");
    fprintf(stdout, "%35s - %s\n", "--flag-ignore-mask=<mask>", "ignores the given flag bits for determining data quality, e.g. --flag-ignore-mask=0b00000110 or --flag-ignore-mask=16");
    fprintf(stdout, "%35s - %s\n", "--flag-mask-type={AND|OR}", "interpret --flag-mask values as bitwise AND or as bitwise OR");
    fprintf(stdout, "%35s - %s\n", "--flag-raised-is-bad", "flag bit 1 signifies an issue. Default: bit equals 0 signifies an issue");
    fprintf(stdout, "%35s - %s\n", "--list-quality-flag-descriptions", "print a table of quality flags");

    return;
}

void about(void)
{
    fprintf(stdout, "tiictbin\n");
    fprintf(stdout, "Copyright (C) 2024 Johnathan K. Burchill.\n");
    fprintf(stdout, " Bins the requested TIICT parameter for statistical analysis.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
    fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
    fprintf(stdout, "under the terms of the GNU General Public License.\n");
    fprintf(stdout, "See the file LICENSE in the source repository for details.\n");
    return;
}

void parseCommandLine(ProcessingParameters *params, int argc, char *argv[])
{
    params->flagIgnoreMask = 0;
    params->flagMaskIsAnd = true;
    params->flagRaisedIsGood = true;

    for (int i = 0; i < argc; i++)
    {
        if (strcmp("--help", argv[i]) == 0)
        {
            usage(argv[0]);
            exit(EXIT_SUCCESS);
        }
        else if (strcmp("--about", argv[i])==0)
        {
            // call about function then exit
            about();
            exit(EXIT_SUCCESS);
        }
        else if (strcmp("--verbose", argv[i])==0)
        {
            params->nOptions++;
            params->verbose = true;
        }
        else if (strcmp(argv[i], "--available-statistics") == 0)
        {
            fprintf(stdout, "Available statistics:\n");
            printAvailableStatistics(stdout);
            exit(EXIT_SUCCESS);
        }
        else if (strcmp(argv[i], "--no-file-progress") == 0)
        {
            params->nOptions++;
            params->showFileProgress = false;
        }
        else if (strcmp(argv[i], "--equal-length-bins") == 0)
        {
            params->nOptions++;
            params->binningState.equalArea = false;
        }
        else if (strncmp(argv[i], "--qdlatmin=", 11) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 11)
                params->binningState.qdlatmin = atof(argv[i]+11);
            else
            {
                fprintf(stderr, "Could not parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp(argv[i], "--qdlatmax=", 11) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 11)
                params->binningState.qdlatmax = atof(argv[i]+11);
            else
            {
                fprintf(stderr, "Could not parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp(argv[i], "--deltaqdlat=", 13) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 13)
                params->binningState.deltaqdlat = atof(argv[i]+13);
            else
            {
                fprintf(stderr, "Could not parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp(argv[i], "--mltmin=", 9) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 9)
                params->binningState.mltmin = atof(argv[i]+9);
            else
            {
                fprintf(stderr, "Could not parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp(argv[i], "--mltmax=", 9) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 9)
                params->binningState.mltmax = atof(argv[i]+9);
            else
            {
                fprintf(stderr, "Could not parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp(argv[i], "--deltamlt=", 11) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 11)
                params->binningState.deltamlt = atof(argv[i]+11);
            else
            {
                fprintf(stderr, "Could not parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--flip-when-descending") == 0)
        {
            params->nOptions++;
            params->binningState.flipParamWhenDescending = true;
        }
        else if (strcmp(argv[i], "--flip-when-sunward") == 0)
        {
            params->nOptions++;
            params->binningState.flipParamWhenSunward = true;
        }
        else if (strcmp(argv[i], "--flip-when-dawnward") == 0)
        {
            params->nOptions++;
            params->binningState.flipParamWhenDawnward = true;
        }
        else if (strncmp(argv[i], "--cdf-input-directory=", 25) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 25)
                params->cdfDirectory = argv[i] + 25;
            else
            {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        // From cdfbin.c in TII-Ion-Drift-Processor on github
        else if (strncmp(argv[i], "--flag-ignore-mask=", 19) == 0)
        {
            params->nOptions++;
            int base = 10;
            int sign = 1;
            size_t offset = 19;
            size_t len = strlen(argv[i]);
            if (len < 20)
            {
                fprintf(stderr, "Invalid quality flag ignore-mask value.\n");
                exit(EXIT_FAILURE);
            }
            if ((argv[i] + offset)[0] == '-')
            {
                sign = -1;
                offset++;
            }
            if (len > 21)
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
            params->flagIgnoreMask = (int32_t) strtol(argv[i] + offset, (char **)NULL, base);
            params->flagIgnoreMask *= sign;
        }
        else if (strncmp(argv[i], "--flag-mask-type=", 17) == 0)
        {
            params->nOptions++;
            if (strcmp(argv[i] + 17, "OR") == 0)
                params->flagMaskIsAnd = false;
            else
                params->flagMaskIsAnd = true;
        }
        else if (strcmp(argv[i], "--flag-raised-is-bad") == 0)
        {
            params->nOptions++;
            params->flagRaisedIsGood = false;
        }
        else if (strcmp(argv[i], "--list-quality-flag-descriptions") == 0)
        {
            params->nOptions++;
            printQualityFlagTable();
            exit(EXIT_SUCCESS);
        }
        else if (strncmp("--", argv[i], 2) == 0)
        {
            fprintf(stderr, "Unrecognized option %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    if (argc - params->nOptions != 6)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    // Process files in directory based on time range
    params->satelliteLetter = argv[1][0];
    params->parameter = argv[2];
    params->statistic = argv[3];
    params->firstTimeString = argv[4];
    params->lastTimeString = argv[5];
    params->firstTime = parseEPOCH4(params->firstTimeString);
    params->lastTime = parseEPOCH4(params->lastTimeString);

    if (!validStatistic(params->statistic))
    {
        fprintf(stderr, "Invalid statistic '%s'\n", params->statistic);
        fprintf(stderr, "Must be one of:\n");
        printAvailableStatistics(stderr);
        exit(EXIT_FAILURE);
    }

    if (params->verbose)
    {
        if (params->flagIgnoreMask == 0)
            fprintf(stdout, "INCLUDING data with quality flag 0 (good data)\n");
        else
        {
            if (params->flagIgnoreMask > 0)
                fprintf(stdout, "INCLUDING good data while ignoring these quality flags\n");
            else
                fprintf(stdout, "INCLUDING good data while ignoring quality flags other than\n");

            int maskBit = 0;
            int nFlaggedParams = 0;
            int64_t unsignedMask = params->flagIgnoreMask;
            if (unsignedMask < 0) {
                unsignedMask = -unsignedMask;
            }
            for (int b = 0; b < N_QUALITY_FLAG_BITS; b++)
            {
                maskBit = (unsignedMask >> b) & 0x1;
                if (maskBit)
                {
                    nFlaggedParams++;
                    if (nFlaggedParams > 1)
                        fprintf(stdout, "%s", params->flagMaskIsAnd ? " AND " : "  OR ");
                    else
                        fprintf(stdout, "     ");
                    fprintf(stdout, "\"%s\"\n", qualityFlagInfo[b]);
                }
            }
            fprintf(stdout, "\n");
        }
    }

    int64_t posVal = params->flagIgnoreMask;
    if (posVal < 0) {
        posVal = -posVal;
    }
    params->positiveFlagMask = posVal;

    return;
}

bool fileMatch(FTSENT *e, ProcessingParameters *params)
{
    // Split filename by "_"
    // Adapted from the strsep manual page example
    char *inputstring = strdup(e->fts_name);
    if (inputstring == NULL) {
        return false;
    }

    char **ap = NULL;
    char *argv[10] = {0};
    int nparts = 0;

    for (ap = argv; (*ap = strsep(&inputstring, "_")) != NULL;nparts++) {
        if (**ap != '\0') {
            if (++ap >= &argv[10]) {
                break;
            }
        }
    }

    if (nparts != 7) {
        return false;
    }

    if (strcmp("SW", argv[0]) != 0) {
        return false;
    }

    if (strlen(argv[2]) != 4 || strlen(argv[3]) != 5) {
        return false;
    }

    if (strncmp("TCT", argv[3], 3) != 0) {
        return false;
    }

    if (argv[2][3] != params->satelliteLetter) {
        return false;
    }

    long y1, m1, d1, h1, min1, s1;
    long y2, m2, d2, h2, min2, s2;
    int nAssigned = 0;
    nAssigned = sscanf(argv[4], "%4ld%2ld%2ldT%2ld%2ld%2ld", &y1, &m1, &d1, &h1, &min1, &s1);
    if (nAssigned != 6) {
        return false;
    }
    nAssigned = sscanf(argv[5], "%4ld%2ld%2ldT%2ld%2ld%2ld", &y2, &m2, &d2, &h2, &min2, &s2);
    if (nAssigned != 6) {
        return false;
    }

    double fileFirstTime = computeEPOCH(y1, m1, d1, h1, min1, s1, 0);
    if (fileFirstTime == ILLEGAL_EPOCH_VALUE) {
        return false;
    }
    double fileLastTime = computeEPOCH(y2, m2, d2, h2, min2, s2, 0);
    if (fileLastTime == ILLEGAL_EPOCH_VALUE) {
        return false;
    }

    double firstTime = params->firstTime;
    double lastTime = params->lastTime;

    bool match = ((strcmp(".cdf", e->fts_name + e->fts_namelen - 4) == 0) && ((firstTime >= fileFirstTime && firstTime <= fileLastTime) || (lastTime >= fileFirstTime && lastTime <= fileLastTime) || (firstTime < fileFirstTime && lastTime > fileLastTime)));

    free(inputstring);

    return match;
}

int processFile(ProcessingParameters *params)
{
    if (params == NULL || params->inputFile == NULL)
        return STATISTICS_POINTER;

    int status = loadTiictData(params);

    if (status != CDF_OK)
        return status;

    float mlt = 0.0;
    float lastMlt = 0.0;
    float qdlat = 0.0;
    float lastQdLat = 0.0;
    float value = 0.0;
    uint32_t flag = 0;
    uint32_t flagMask = ~params->flagIgnoreMask;

    bool includeValue = false;
    float qdDisplacement = 0.0;
    float dawnwardDisplacement = 0.0;
    float sunwardDisplacement = 0.0;

    // Track direction in magnetic coordinates
    float xmag = 0.0;
    float lastXmag =0.0;
    float ymag = 0.0;
    float lastYmag = 0.0;
    float theta = 0.0;

    for (int i = 0; i < params->nRecords; i++)
    {
        qdlat = params->qdlat[i];
        if (i > 0) {
            qdDisplacement = qdlat - lastQdLat;
        }
        else {
            qdDisplacement = 0.0;
        }
        lastQdLat = qdlat;

        mlt = params->mlt[i];

        // Track direction of orbit in magnetic coordinates
        if (i > 0) {
            // xhat is antisunward, yhat is dawnward
            theta = 2.0 * M_PI * mlt / 24.0;
            xmag = fabsf(qdlat) * cosf(theta);
            ymag = fabsf(qdlat) * sinf(theta);
            sunwardDisplacement = -(xmag - lastXmag);
            dawnwardDisplacement = ymag - lastYmag;
        }
        else {
            sunwardDisplacement = 0.0;
            dawnwardDisplacement = 0.0;
        }
        lastQdLat = qdlat;
        lastMlt = mlt;
        lastXmag = xmag;
        lastYmag = ymag;

        value = params->values[i];

        if (params->flags != NULL)
            flag = params->flags[i];
        else flag = 0;

        if (params->binningState.flipParamWhenDescending && qdDisplacement < 0.0) {
            value = -value;
        }
        else if (params->binningState.flipParamWhenSunward && sunwardDisplacement > 0) {
            value = -value;
        }
        else if (params->binningState.flipParamWhenDawnward && dawnwardDisplacement > 0) {
            value = -value;
        }
        params->binningState.nValsRead++;

        // Filter based on flag ignore mask
        // Ignore the masked flag bits. All other bits must be 1.
        // TODO Fix this
        //includeValue = (flag & flagMask) == 0;
        includeValue = true;

        binData(&params->binningState, qdlat, mlt, value, includeValue);
    }

    free(params->time);
        params->time = NULL;
    free(params->mlt);
        params->mlt = NULL;
    free(params->qdlat);
        params->qdlat = NULL;
    free(params->values);
        params->values = NULL;
    free(params->flags);
        params->flags = NULL;

    return STATISTICS_OK;
}

int loadTiictData(ProcessingParameters *params)
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
    long decoding, encoding, majority, maxrRec, numrVars, maxzRec, numzVars, numAttrs, format, numDims, dimSizes[CDF_MAX_DIMS];

    status = CDFopenCDF(params->inputFile, &cdfId);
    if (status != CDF_OK)
        return status;

    status = CDFgetFormat(cdfId, &format);
    status = CDFgetDecoding(cdfId, &decoding);
    status = CDFinquireCDF(cdfId, &numDims, dimSizes, &encoding, &majority, &maxrRec, &numrVars, &maxzRec, &numzVars, &numAttrs);
    if (status != CDF_OK)
    {
    	CDFcloseCDF(cdfId);
        return status;
    }
    uint8_t nVars = numzVars;

    status = loadCdfVariable(cdfId, "Timestamp", (void**)&params->time, &params->nRecords);
    if (status != CDF_OK)
    {
        CDFcloseCDF(cdfId);
        return status;
    }

    status = loadCdfVariable(cdfId, "MLT", (void**)&params->mlt, NULL);
    if (status != CDF_OK)
    {
        CDFcloseCDF(cdfId);
        return status;
    }

    status = loadCdfVariable(cdfId, "QDLatitude", (void**)&params->qdlat, NULL);
    if (status != CDF_OK)
    {
        CDFcloseCDF(cdfId);
        return status;
    }

    status = loadCdfVariable(cdfId, params->parameter, (void**)&params->values, NULL);
    if (status != CDF_OK)
    {
        CDFcloseCDF(cdfId);
        return status;
    }

    status = loadCdfVariable(cdfId, "Quality_flags", (void**)&params->flags, NULL);

	CDFcloseCDF(cdfId);

	return status;

}

CDFstatus loadCdfVariable(CDFid cdfId, char *variable, void **mem, long *nRecords)
{
    if (mem == NULL)
        return -1;

    long numBytesToAdd, numVarBytes, numValuesPerRec;
    long varNum, dataType, numElems, numRecs, numDims, recVary;
    long dimSizes[CDF_MAX_DIMS], dimVarys[CDF_MAX_DIMS];
    CDFdata data;

	varNum = CDFgetVarNum(cdfId, variable);

	CDFstatus status = CDFreadzVarAllByVarID(cdfId, varNum, &numRecs, &dataType, &numElems, &numDims, dimSizes, &recVary, dimVarys, &data);
	if (status != CDF_OK)
	{
    	CDFcloseCDF(cdfId);
		CDFdataFree(data);
		return status;
	}

    if (nRecords != NULL)
        *nRecords = numRecs;

	// Calculate new size of memory to allocate
	status = CDFgetDataTypeSize(dataType, &numVarBytes);
	numValuesPerRec = 1;
	for (uint8_t j = 0; j < numDims; j++)
	{
		numValuesPerRec *= dimSizes[j];
	}
	numBytesToAdd = numValuesPerRec * numRecs * numVarBytes;
    *mem = malloc(numBytesToAdd);
    if (*mem == NULL)
    {
        printf("Could not allocate heap.\n");
        CDFdataFree(data);
        CDFcloseCDF(cdfId);
        exit(42);
    }
    memcpy(*mem, data, numBytesToAdd);
	CDFdataFree(data);

    return CDF_OK;
}

void printQualityFlagTable(void)
{
    fprintf(stdout, "Quality flag bit = 1 indicates nominal measurement.\n");
    fprintf(stdout, "%27s%s\n", "Flag value", " Description");
    int value = 0;
    for (int i = 0; i < N_QUALITY_FLAG_BITS; i++) {
        value = 1 << i;
        for (int b = N_QUALITY_FLAG_BITS - 1; b >= 0; b--) {
            fprintf(stdout, "%1d", (value >> b) & 0x01);
        }
        fprintf(stdout, " %7d %s\n", value, qualityFlagInfo[i]);
    }
}
