#include "lpData.h"

#include "indexing.h"
#include "utilities.h"

#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <cdf.h>
#include <fts.h>
#include <stdio.h>
#include <string.h>

extern char infoHeader[50];

int getLpData(const char *lpDir, const char *satellite, const int year, const int month, const int day, uint8_t **dataBuffers, double **lpPhiScHighGain, double **lpPhiScLowGain)
{

    int status = 0;

    // Get data from previous day, requested date, and next day    
    int y = year;
    int m = month;
    int d = day;
    
    struct tm date = {0};
    date.tm_year = y - 1900;
    date.tm_mon = m - 1;
    date.tm_mday = d - 1;
    char lpFile[FILENAME_MAX];
    int res = 0;

    double *lpTimes2Hz = NULL;
    double *lpVsHg = NULL;
    double *lpVsLg = NULL;
    long nLpRecs = 0;

    for (int i = 0; i < 3; i++)
    {
        timegm(&date);
        bzero(lpFile, FILENAME_MAX);
        res = getLpInputFilename(satellite[0], date.tm_year+1900, date.tm_mon+1, date.tm_mday, lpDir, lpFile);
        if (res != 0)
        {
            fprintf(stderr, "%sNo LP data for %4d%02d%02d\n", infoHeader, date.tm_year+1900, date.tm_mon+1, date.tm_mday);
            continue;
        }
        fprintf(stderr, "%sLoading LP data from %s\n", infoHeader, lpFile);
        
        loadLpInputs(lpFile, &lpTimes2Hz, &lpVsHg, &lpVsLg, &nLpRecs);

        date.tm_mday = date.tm_mday + 1;
        
    }

    fprintf(stderr, "%sLoaded %ld LP HM records\n", infoHeader, nLpRecs);


    long ym, mm, dm, hm, minm, sm, msm;

    for (long i = 0; i < nLpRecs; i+=1000)
    {
        EPOCHbreakdown(lpTimes2Hz[i], &ym, &mm, &dm, &hm, &minm, &sm, &msm);
        fprintf(stderr, "t=%4ld%02ld%02ldT%02ld%02ld%02ld.%03ld\t%f\t%f\n", ym, mm, dm, hm, minm, sm, msm, lpVsHg[i], lpVsLg[i]);
    }

    // interpolate LP data to TII times

    


    return status;

}

// Adapted from SLIDEM
int getLpInputFilename(const char satelliteLetter, long year, long month, long day, const char *path, char *filename)
{
	char *searchPath[2] = {NULL, NULL};
    searchPath[0] = (char *)path;

	FTS * fts = fts_open(searchPath, FTS_PHYSICAL | FTS_NOCHDIR, NULL);	
	if (fts == NULL)
	{
		printf("Could not open directory %s for reading.", path);
		return -1;
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
    {
        return 0;
    }
    else
    {
        return -1;
    }

}



// Adapted from SLIDEM
void loadLpInputs(const char *cdfFile, double **lpTime, double **lpPhiScHighGain, double **lpPhiScLowGain, long *numberOfRecords)
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

    status = CDFopenCDF(cdfFile, &cdfId);
    if (status != CDF_OK) 
    {
        printErrorMessage(status);
        fprintf(stdout, "%s Could not open CDF file. Skipping this date.\n", infoHeader);
        return;
    }

    status = CDFgetFormat(cdfId, &format);
    status = CDFgetDecoding(cdfId, &decoding);
    status = CDFinquireCDF(cdfId, &numDims, dimSizes, &encoding, &majority, &maxrRec, &numrVars, &maxzRec, &numzVars, &numAttrs);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        fprintf(stdout, "\n%s Problem with CDF file. Skipping this date.\n", infoHeader);
        closeCdf(cdfId);
        return;
    }
    int nVariables = 3;
    char * variables[3] = {"Timestamp", "Vs_hgn", "Vs_lgn"};

    for (uint8_t i = 0; i<nVariables; i++)
    {
        status = CDFconfirmzVarExistence(cdfId, variables[i]);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            fprintf(stdout, "\n%s Error reading variable %s from CDF file. Skipping this date.\n", infoHeader, variables[i]);
            closeCdf(cdfId);
                return;
        }
        else
        {
            // fprintf(stdout, "%s OK\n", infoHeader);
        }
    }
    
    for (uint8_t i = 0; i < nVariables; i++)
    {
        varNum = CDFgetVarNum(cdfId, variables[i]);
        status = CDFreadzVarAllByVarID(cdfId, varNum, &numRecs, &dataType, &numElems, &numDims, dimSizes, &recVary, dimVarys, &data);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            fprintf(stdout, "%s Error loading data for %s. Skipping this date.\n", infoHeader, variables[i]);
            closeCdf(cdfId);
            CDFdataFree(data);
            return;
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
                *lpTime = (double*) realloc(*lpTime, totalBytes);
                memcpy(*lpTime + (*numberOfRecords), data, numBytesToAdd);
                break;
            case 1: // Vs_hgn
                *lpPhiScHighGain = (double*) realloc(*lpPhiScHighGain, totalBytes);
                memcpy(*lpPhiScHighGain + (*numberOfRecords), data, numBytesToAdd);
                break;
            case 2: // Vs_lgn
                *lpPhiScLowGain = (double*) realloc(*lpPhiScLowGain, totalBytes);
                memcpy(*lpPhiScLowGain + (*numberOfRecords), data, numBytesToAdd);
                break;
            default:
                return;
        }
        CDFdataFree(data);
    }
    // close CDF
    closeCdf(cdfId);

    // Update number of records found and memory allocated
    *numberOfRecords += numRecs;
    // *totalMemoryAllocated = fpMemorySize;

}
