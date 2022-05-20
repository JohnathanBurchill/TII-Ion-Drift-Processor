/*

    TII Cross-Track Ion Drift Processor: tiict.c

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

#include "tiict.h"
#include "settings.h"
#include "utilities.h"
#include "processing.h"
#include "indexing.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

// https://cdf.gsfc.nasa.gov/
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

    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--about") == 0)
        {
            fprintf(stdout, "tiict - TII Cross-track ion drift processor, version %s.\n", SOFTWARE_VERSION);
            fprintf(stdout, "Copyright (C) 2022  Johnathan K Burchill\n");
            fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
            fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
            fprintf(stdout, "under the terms of the GNU General Public License.\n");

            exit(0);
        }
    }

    if (argc != 9)
    {
        fprintf(stdout, "usage: %s satLetter year month day calversionString exportVersionString calDir exportDir\n", argv[0]);
        exit(1);
    }

    const char* satellite = argv[1];
    const int year = atoi(argv[2]);
    const int month = atoi(argv[3]);
    const int day = atoi(argv[4]);
    const char* calVersion = argv[5];
    const char* exportVersion = argv[6];
    const char* calDir = argv[7];
    const char* exportDir = argv[8];

    char fitLogFileName[CDF_PATHNAME_LEN + 8];

    // Check satellite letter
    if (strlen(satellite) != 1 || (satellite[0] != 'A' && satellite[0] != 'B' && satellite[0] != 'C'))
    {
        fprintf(stdout, "Satellite must be one of 'A', 'B', or 'C' (no quotes).\n");
        exit(1);
    }

    // set up info header
    sprintf(infoHeader, "TIICT %c%s %04d-%02d-%02d: ", satellite[0], exportVersion, year, month, day);
    fprintf(stdout, "\n%s-------------------------------------------------\n", infoHeader);
    fprintf(stdout, "%sVersion 0302 20220519\n", infoHeader);
    fprintf(stdout, "%sProcessing date: %s\n", infoHeader, asctime(timeParts));

    // Confirm requested date has records. Abort otherwise.
    long numAvailableRecords = numberOfAvailableRecordsForDate(satellite, year, month, day, calDir, calVersion);
    if (numAvailableRecords < (16 * SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING))
    {
        fprintf(stdout, "%sLess than %.0f s of data available. Skipping this date.\n", infoHeader, (float)SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING);
        exit(1);
    }
    else
    {
        fprintf(stdout, "%sProcessing %ld records for this date.\n", infoHeader, numAvailableRecords);
    }

    sprintf(fitLogFileName, "%s/%s/TCT16/SW_EXPT_EFI%s_TCT16_%04d%02d%02dT000000_%04d%02d%02dT235959_%s_fits.txt", exportDir, exportVersion, satellite, year, month, day, year, month, day, exportVersion);
    fprintf(stdout, "%sFit log file: %s\n", infoHeader, fitLogFileName);

    CDFstatus status;

    // Turn off GSL failsafe error handler. We typically check the GSL return codes.
    gsl_set_error_handler_off();

    // Offset model parameters
    offset_model_fit_arguments fitargs[4] = {
        {0, "Northern ascending", 44.0, 50.0, 50.0, 44.0},
        {1, "Equatorial descending", 44.0, 38.0, -38.0, -44.0},
        {2, "Southern descending", -44.0, -50.0, -50.0, -44.0},
        {3, "Equatorial ascending", -44.0, -38.0, 38.0, 44.0},
    };

    FILE *fitFile = fopen(fitLogFileName, "w");
    if (fitFile == NULL)
    {
        fprintf(stdout, "%sCould not open fit log for writing:\n  %s\nAborting.\n", infoHeader, fitLogFileName);
        exit(1);
    }
    fprintf(fitFile, "EFI TII CrossTrackCalibration fit results by fit region.\n");
    fprintf(fitFile, "Each region consists of two mid-latitude segments denoted by CDF_EPOCH times T11, T12, T21, and T22.\n");
    fprintf(fitFile, "Linear models based on robust least squares (GNU Scientific Library) are subtracted from each region for which a fit can be obtained.\n");
    fprintf(fitFile, "Regions:\n");
    for (uint8_t ind = 0; ind < 4; ind++)
    {
        fprintf(fitFile, "%d %21s: (% 5.1f, % 5.1f) -> (% 5.1f, % 5.1f)\n", fitargs[ind].regionNumber, fitargs[ind].regionName, fitargs[ind].lat1, fitargs[ind].lat2, fitargs[ind].lat3, fitargs[ind].lat4);
    }
    fprintf(fitFile, "\n");
    fprintf(fitFile, "The columns are:\n");
    fprintf(fitFile, "regionNumber fitNumber numPoints1 numPoints2 T11 T12 T21 T22 offsetHX slopeHX adjRsqHX rmseHX medianHX1 medianHX2 madHX madHX1 madHX2 offsetHY slopeHY adjRsqHY rmseHY medianHY1 medianHY2 madHY madHY1 madHY2 offsetVX slopeVX adjRsqVX rmseVX medianVX1 medianVX2 madVX madVX1 madVX2 offsetVY slopeVY adjRsqVY rmseVY medianVY1 medianVY2 madVY madVY1 madVY2\n");
    fprintf(fitFile, "\n");
    fflush(fitFile);
    fflush(stdout);

    // load calibration data
    // The memory pointers
    uint8_t * dataBuffers[NUM_CAL_VARIABLES];
    for (uint8_t i = 0; i < NUM_CAL_VARIABLES; i++)
    {
        dataBuffers[i] = NULL;
    }
    long nRecs = 0;

    loadCalData(calDir, calVersion, satellite, year, month, day, dataBuffers, &nRecs);

    fflush(stdout);

    // Do the calibration
    double *pEpoch = (double*) dataBuffers[0];
    long long timeIndex = 0;

    // The macros TIME(), QDLAT(), etc. give the value at index timeIndex.
    double epoch0 = TIME();

    // Strategy is to perform calculations in place in memory without allocating space unecessarily.
    // Update 
    //  0) Adjust times by + 1./32 - 0.0875 s
    //  1) VSatXYZ from km/s to m/s
    //  2) bias voltage to nearest of -100 V, -62 V or left alone if too different from those values
    //  3) apply Level 1 calibration with bias dependence
    float value, innerDomeBias;

    if (nRecs < 16*SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING)
    {
        fprintf(stdout, "%sFewer than %d records found (got %ld). Skipping this date.\n", infoHeader, 16*5, nRecs);
        fclose(fitFile);
        exit(1);
    }

    for (long timeIndex = 0; timeIndex < nRecs; timeIndex++)
    {
        // Adjust times for sample lag
        ((double*)dataBuffers[0])[timeIndex] += ((1./32 - 0.0875)*1000.); // ms
        // Bias Voltage
        // TODO: need a more accurate replacement: i.e., early in mission the voltage was ~-60 V, not -62.
        if (VBIAS() < -95.0)
            *ADDR(13, 0, 1) = -100.;
        else if(VBIAS() > -65. && VBIAS() < -59.0)
            *ADDR(13, 0, 1) = -62.0;
//        fprintf(stdout, "%s; after: %f\n", infoHeader, VBIAS());

        // VSatXYZ to m/s
        *ADDR(3, 0, 3) *= 1000.0;
        *ADDR(3, 1, 3) *= 1000.0;
        *ADDR(3, 2, 3) *= 1000.0;

        // VCorot to m/s
        *ADDR(10, 0, 3) *= 1000.0;
        *ADDR(10, 1, 3) *= 1000.0;
        *ADDR(10, 2, 3) *= 1000.0;

        // Apply level1 calibration with bias dependence
        innerDomeBias = VBIAS() - VFP();

        // Get scaling parameter
        float shx, shy, svx, svy;
        switch(satellite[0])
        {
            case 'A': // 20190930 slew experiment, plus simulations
                shx = 574.0;
                shy = 574.0;
                svx = 712.0;
                svy = 712.0;
                if (innerDomeBias >= -63.0 && innerDomeBias < -59.0)
                {
                    shx *= .76;
                    shy *= .76;
                    svx *= .76;
                    svy *= .76;
                }
                break;
            case 'B': // 20210616 slew experiment
                shx = 553.0;
                shy = 553.0;
                svx = 548.0;
                svy = 548.0;
                if (innerDomeBias >= -63.0 && innerDomeBias < -59.0)
                {
                    shx *= 450.5;
                    shy = 450.5; // Calibration 20210624, inner dome bias at -62 V.
                    svx *= 451.0;
		    svy = 525.0; // Calibration 20210624, inner dome bias at -62 V.
                }
                break;
            case 'C': // 20191126 slew experiment, plus simulations
                shx = 679.0;
                shy = 679.0;
                svx = 2377.0;
                svy = 2377.0;
                if (innerDomeBias >= -63.0 && innerDomeBias < -59.0)
                {
                    shx *= .76;
                    shy *= .76;
                    svx *= .76;
                    svy *= .76;
                }
                break;
        }
        // Change sign to get flow directions correct, then apply scaling and remove known satellite signal
//        fprintf(stdout, "%sshx: %f\n", infoHeader, shx);
        *ADDR(1, 0, 2) = -1.0 * (MXH() - 32.5) * shx - VSATX(); 
        *ADDR(1, 1, 2) = -1.0 * (MYH() - 32.5) * shy - VSATY();
        *ADDR(2, 0, 2) = -1.0 * (MXV() - 32.5) * svx - VSATX();
        *ADDR(2, 1, 2) = -1.0 * (MYV() - 32.5) * svy - VSATZ();

    }

    fprintf(stdout, "%sPrepared calibration data.\n", infoHeader);
    fflush(stdout);

    // Set up new memory
    float *xhat = (float*) malloc((size_t) (nRecs * sizeof(float) * 3));
    float *yhat = (float*) malloc((size_t) (nRecs * sizeof(float) * 3));
    float *zhat = (float*) malloc((size_t) (nRecs * sizeof(float) * 3));
    float *ectFieldH = (float*) malloc((size_t) (nRecs * sizeof(float) * 3));
    float *ectFieldV = (float*) malloc((size_t) (nRecs * sizeof(float) * 3));
    float *bctField = (float*) malloc((size_t) (nRecs * sizeof(float) * 3));
    // Quality flag and fitInfo flag initialized to zero
    uint16_t *flags = (uint16_t*) malloc((size_t) (nRecs * sizeof(uint16_t)));
    uint32_t *fitInfo = (uint32_t*) malloc((size_t) (nRecs * sizeof(uint32_t)));
    // Error estimates from Mean Absolute Deviation (MAD): default is -42. :)
    float *viErrors = (float*) malloc((size_t) (nRecs * sizeof(float) * 4));
    for (long ind = 0; ind < nRecs; ind++)
    {
        viErrors[4*ind+0] = DEFAULT_VI_ERROR;
        viErrors[4*ind+1] = DEFAULT_VI_ERROR;
        viErrors[4*ind+2] = DEFAULT_VI_ERROR;
        viErrors[4*ind+3] = DEFAULT_VI_ERROR;
        flags[ind] = 0;
        // FITINFO_OFFSET_NOT_REMOVED = 1 and FITINFO_INCOMPLETE_REGION = 1 are the defaults for fitInfo
        // Set for each velocity component
        for (uint8_t k = 0; k < 4; k++)
        {
            fitInfo[ind] |= ((FITINFO_OFFSET_NOT_REMOVED | FITINFO_INCOMPLETE_REGION)) << (k * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT);
        }
    }

    // Remove offsets and calculate flags
    // 1. Ascending northern polar region to descending northern polar region
    // 2. Descending low-lat region
    // 3. Descending southern polar region to ascending southern polar region
    // 4. Ascending low-lat region
    // 5. Repeat. Does not have to be sequential
    for (uint8_t ind = 0; ind < 4; ind++)
    {
        removeOffsetsAndSetFlags(satellite, fitargs[ind], nRecs, dataBuffers, viErrors, flags, fitInfo, fitFile);
    }

    fprintf(stdout, "%sRemoved offsets and calculated flags.\n", infoHeader);
    fflush(fitFile);
    fflush(stdout);

    // Calculate fields
    long long ind;
    for (timeIndex = 0; timeIndex < nRecs; timeIndex++)
    {
        // Calculate xhat, yhat, zhat
        // xhat parallel to satellite velocity 
        float magVsat = sqrtf(VSATN()*VSATN() + VSATE()*VSATE() + VSATC() * VSATC());
        ind = 3*timeIndex;
        xhat[ind + 0] = VSATN() / magVsat;
        xhat[ind + 1] = VSATE() / magVsat;
        xhat[ind + 2] = VSATC() / magVsat;
        // yhat is xhat cross {0, 0, -1}
        yhat[ind + 0] = -xhat[ind + 1];
        yhat[ind + 1] = xhat[ind + 0];
        yhat[ind + 2] = 0.0;
        float magyhat = sqrtf(yhat[ind + 0] * yhat[ind + 0] + yhat[ind + 1] * yhat[ind + 1] + yhat[ind + 2] * yhat[ind + 2]);
        yhat[ind + 0] /= magyhat;
        yhat[ind + 1] /= magyhat;
        yhat[ind + 2] /= magyhat;
        // zhat is the cross-product of x into y
        zhat[ind+0] = xhat[ind+1] * yhat[ind+2] - xhat[ind+2] * yhat[ind+1];
        zhat[ind+1] = -1.0 * xhat[ind+0] * yhat[ind+2] + xhat[ind+2] * yhat[ind+0];
        zhat[ind+2] = xhat[ind+0] * yhat[ind+1] - xhat[ind+1] * yhat[ind+0];
        float magzhat = sqrtf(zhat[ind + 0] * zhat[ind + 0] + zhat[ind + 1] * zhat[ind + 1] + zhat[ind + 2] * zhat[ind + 2]);
        zhat[ind + 0] /= magzhat;
        zhat[ind + 1] /= magzhat;
        zhat[ind + 2] /= magzhat;

        // B field in cross-track frame, nT:
        bctField[ind + 0] = BN() * xhat[ind + 0] + BE() * xhat[ind + 1] + BC() * xhat[ind + 2];
        bctField[ind + 1] = BN() * yhat[ind + 0] + BE() * yhat[ind + 1] + BC() * yhat[ind + 2];
        bctField[ind + 2] = BN() * zhat[ind + 0] + BE() * zhat[ind + 1] + BC() * zhat[ind + 2];

        // E field from H sensor X, in cross-track frame, mV/m:
        ectFieldH[ind + 0] = -1.0 * (MYH() * bctField[ind + 2] - MYV() * bctField[ind + 1]) / 1000000000.0 * 1000.0; 
        ectFieldH[ind + 1] = -1.0 * (-1.0 * MXH() * bctField[ind + 2] + MYV() * bctField[ind + 0]) / 1000000000.0 * 1000.0; 
        ectFieldH[ind + 2] = -1.0 * (MXH() * bctField[ind + 1] - MYH() * bctField[ind + 0]) / 1000000000.0 * 1000.0; 

        // E field from V sensor X, in cross-track frame, mV/m:
        ectFieldV[ind + 0] = -1.0 * (MYH() * bctField[ind + 2] - MYV() * bctField[ind + 1]) / 1000000000.0 * 1000.0; 
        ectFieldV[ind + 1] = -1.0 * (-1.0 * MXV() * bctField[ind + 2] + MYV() * bctField[ind + 0]) / 1000000000.0 * 1000.0; 
        ectFieldV[ind + 2] = -1.0 * (MXV() * bctField[ind + 1] - MYH() * bctField[ind + 0]) / 1000000000.0 * 1000.0; 

    }

    fprintf(stdout, "%sCalculated fields.\n", infoHeader);
    fflush(stdout);

    // Export the cross-track flow CDF file
    long exportDimSizes[1] = {0};
    long exportVarNumbers[29];
    static long recVary = {VARY};
    static long dimVary = {VARY};
    static long dimNoVary = {NOVARY};

    // Keep only data from the requested day
    double minTime, maxTime, startTime, stopTime, prevTime;
    long startIndex = 0, stopIndex = nRecs;
    minTime = computeEPOCH(year, month, day, 0, 0, 0, 0);
    maxTime = computeEPOCH(year, month, day + 1, 0, 0, 0, 0);
    bool gotStartTime = false, gotStopTime = false;
    char errorMessage[CDF_STATUSTEXT_LEN + 1];
    bool exporting = true;
    long recordsExported = 0;
    double minutesExported = 0.;
    uint16_t filesExported = 0;

    // Find first index with time on requested day
    timeIndex = 0;
    while (TIME() < minTime) timeIndex++;

    startIndex = timeIndex;
    startTime = TIME();
    // Find last index not containing a gap of more than 10 minutes, or last index of the day
    stopIndex = startIndex;
    stopTime = startTime;
    
    // case where there are no records on the requested or following days
    if (startIndex == nRecs)
    {
        fprintf(stdout, "%sNo records found for requested date.\n", infoHeader);
        exporting = false;
    }
    while (exporting)
    {
        double gapTime = (TIME() - stopTime)/1000.;
        double duration = (stopTime - startTime)/1000.;
        if ((gapTime > MAX_ALLOWED_CDF_GAP_SECONDS) || (TIME() >= maxTime) || (timeIndex == (nRecs-1)))
        {
            // Export data and continue
            if (duration >= SECONDS_OF_DATA_REQUIRED_FOR_EXPORTING)
            {
                // 16 Hz data
                exportTCT16Cdfs(startTime, stopTime, exportDir, exportVersion, calVersion, satellite, startIndex, stopIndex, dataBuffers, ectFieldH, ectFieldV, bctField, viErrors, flags, fitInfo);
                // 2 Hz data
                exportTCT02Cdfs(startTime, stopTime, exportDir, exportVersion, calVersion, satellite, startIndex, stopIndex, dataBuffers, ectFieldH, ectFieldV, bctField, viErrors, flags, fitInfo);
                recordsExported += (stopIndex - startIndex + 1);
                minutesExported += (stopTime - startTime)/1000./60.;
                filesExported++; // Only count 16 Hz files
            }
            else
            {
                char startString[EPOCH_STRING_LEN+1], stopString[EPOCH_STRING_LEN+1];
                toEncodeEPOCH(startTime, 0, startString);
                toEncodeEPOCH(TIME(), 0, stopString);
                fprintf(stdout, "%sInterval spanning %s to %s has a duration of less than %d seconds: not exporting %ld records.\n", infoHeader, startString, stopString, SECONDS_OF_DATA_REQUIRED_FOR_EXPORTING, stopIndex - startIndex + 1);
            }

            // Try next interval
            startIndex = timeIndex;
            startTime = TIME();
            if (TIME() >= maxTime || timeIndex == nRecs)
            {
                exporting = false;
            }
        }
        stopIndex = timeIndex;
        stopTime = TIME();
        timeIndex++;
        if (timeIndex == nRecs)
        {
            exporting = false;
        }
    }
    // report
    fprintf(stdout, "%sExported %.0f orbits (%ld 16 Hz records) of science data in %d files. %.1f%% coverage.\n", infoHeader, minutesExported/94., recordsExported, filesExported, minutesExported/1440.0*100.0); 

    // Close fit log file
    fclose(fitFile);
    fflush(fitFile);
    fflush(stdout);

    // Free the memory
    for (uint8_t i = 0; i < NUM_CAL_VARIABLES; i++)
    {
        free(dataBuffers[i]);
    }

    free(xhat);
    free(yhat);
    free(zhat);
    free(ectFieldH);
    free(ectFieldV);
    free(bctField);
    free(flags);
 
    return 0;
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

CDFstatus createVarFrom1DVar(CDFid id, char *name, long dataType, long startIndex, long stopIndex, void *buffer)
{
    CDFstatus status;
    static long exportDimSizes[1] = {0};
    static long recVary = {VARY};
    static long dimNoVary = {NOVARY};
    long varNumber;
    status = CDFcreatezVar(id, name, dataType, 1, 0L, exportDimSizes, recVary, dimNoVary, &varNumber);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    status = CDFsetzVarSparseRecords(id, varNumber, NO_SPARSERECORDS);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    long dataTypeSize;
    status = CDFgetDataTypeSize(dataType, &dataTypeSize);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    status = CDFputVarRangeRecordsByVarName(id, name, 0, stopIndex-startIndex, (void*)((uint8_t*)buffer + (dataTypeSize*startIndex)));
    if (status != CDF_OK)
    {
        printErrorMessage(status);
    }
    return status;
}

CDFstatus createVarFrom2DVar(CDFid id, char *name, long dataType, long startIndex, long stopIndex, void *buffer, uint8_t index, uint8_t dimSize)
{
    CDFstatus status;
    long long nRecs = stopIndex - startIndex + 1;
    static long exportDimSizes[1] = {0};
    static long recVary = {VARY};
    static long dimNoVary = {NOVARY};
    long varNumber;
    status = CDFcreatezVar(id, name, dataType, 1, 0L, exportDimSizes, recVary, dimNoVary, &varNumber);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    status = CDFsetzVarSparseRecords(id, varNumber, NO_SPARSERECORDS);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    long dataSize;
    status = CDFgetDataTypeSize(dataType, &dataSize);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    float * exportBuffer = (float*)malloc((size_t)(nRecs * dataSize));
    long ind = 0;
    float * buf = (float *) buffer;
    for (long i = (dimSize * startIndex + index); i<(dimSize*(startIndex + nRecs)); i+=dimSize)
    {
        exportBuffer[ind++] = buf[i];
    }
    status = CDFputVarRangeRecordsByVarName(id, name, 0, stopIndex-startIndex, (void *)exportBuffer);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
    }

    free(exportBuffer);

    return status;
}

CDFstatus addgEntry(CDFid id, long attrNum, long entryNum, const char *entry)
{
    CDFstatus status = CDFputAttrgEntry(id, attrNum, entryNum, CDF_CHAR, strlen(entry), (void *)entry);
    return status;
}

CDFstatus addVariableAttributes(CDFid id, varAttr attr)
{
    CDFstatus status;
    char * variableName = attr.name;
    long varNum = CDFvarNum(id, variableName);
    status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "FIELDNAM"), varNum, CDF_CHAR, strlen(variableName), variableName);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "LABLAXIS"), varNum, CDF_CHAR, strlen(variableName), variableName);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "VAR_TYPE"), varNum, CDF_CHAR, 4, "data");
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    if (varNum != 0) // Everything but time
    {
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "DISPLAY_TYPE"), varNum, CDF_CHAR, 11, "time_series");
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            return status;
        }
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "DEPEND_0"), varNum, CDF_CHAR, 4, "Timestamp");
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            return status;
        }
    }
    else // Add the time base to Time
    {
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "TIME_BASE"), varNum, CDF_CHAR, 3, "AD0");
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            return status;
        }
    }
    status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "Type"), varNum, CDF_CHAR, strlen(attr.type), attr.type);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    if (attr.units[0] == '*')
    {
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "UNITS"), varNum, CDF_CHAR, 1, "");
    }
    else
    {
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "UNITS"), varNum, CDF_CHAR, strlen(attr.units), attr.units);
    }
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "CATDESC"), varNum, CDF_CHAR, strlen(attr.desc), attr.desc);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }

    // data type for valid min and max
    if (strcmp(attr.type, "CDF_EPOCH") == 0)
    {
        double val = attr.validMin;
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "VALIDMIN"), varNum, CDF_EPOCH, 1, &val);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            return status;
        }
        val = attr.validMax;
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "VALIDMAX"), varNum, CDF_EPOCH, 1, &val);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            return status;
        }
    }
    else if (strcmp(attr.type, "CDF_UINT2") == 0)
    {
        uint16_t val = (uint16_t) attr.validMin;
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "VALIDMIN"), varNum, CDF_UINT2, 1, &val);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            return status;
        }
        val = (uint16_t) attr.validMax;
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "VALIDMAX"), varNum, CDF_UINT2, 1, &val);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            return status;
        }
    }
    else if (strcmp(attr.type, "CDF_UINT4") == 0)
    {
        uint32_t val = (uint32_t) attr.validMin;
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "VALIDMIN"), varNum, CDF_UINT4, 1, &val);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            return status;
        }
        val = (uint32_t) attr.validMax;
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "VALIDMAX"), varNum, CDF_UINT4, 1, &val);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            return status;
        }
    }
    else if (strcmp(attr.type, "CDF_FLOAT") == 0)
    {
        float val = (float) attr.validMin;
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "VALIDMIN"), varNum, CDF_FLOAT, 1, &val);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            return status;
        }
        val = (float) attr.validMax;
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "VALIDMAX"), varNum, CDF_FLOAT, 1, &val);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            return status;
        }
    }

    return status;
}

void addAttributes(CDFid id, const char *dataset, const char *satellite, const char *version, double minTime, double maxTime)
{
    long attrNum;
    char buf[1000];

    bool highRate = (strcmp(dataset, "TCT16") == 0);

    // Global attributes
    CDFcreateAttr(id, "File_naming_convention", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "SW_EXPT_EFIX");
    CDFcreateAttr(id, "Logical_file_id", GLOBAL_SCOPE, &attrNum);
    sprintf(buf, "swarm%c_TII_H0__v%s", tolower(satellite[0]), version);
    addgEntry(id, attrNum, 0, buf);
    CDFcreateAttr(id, "Logical_source", GLOBAL_SCOPE, &attrNum);
    sprintf(buf, "Swarm%s_TII_H0", satellite);
    addgEntry(id, attrNum, 0, buf);
    CDFcreateAttr(id, "Logical_source_description", GLOBAL_SCOPE, &attrNum);
    sprintf(buf, "Swarm %s Thermal Ion Imager High resolution data", satellite);
    addgEntry(id, attrNum, 0, buf);
    CDFcreateAttr(id, "Mission_group", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "Swarm");
    CDFcreateAttr(id, "MODS", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "Second release of version 3, corrects data gaps associated with L0 overlaps.");
    CDFcreateAttr(id, "PI_name", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "David Knudsen");   
    CDFcreateAttr(id, "PI_affiliation", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "University of Calgary");
    CDFcreateAttr(id, "Acknowledgement", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "ESA Swarm EFI TII data are available from https://swarm-diss.eo.esa.int");
    CDFcreateAttr(id, "Source_name", GLOBAL_SCOPE, &attrNum);
    sprintf(buf, "Swarm%s>Swarm %s", satellite, satellite);
    addgEntry(id, attrNum, 0, buf);
    CDFcreateAttr(id, "Data_type", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "H0>High resolution data");
    CDFcreateAttr(id, "Data_version", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, version);
    CDFcreateAttr(id, "Descriptor", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "TII>Electric Field Instrument Thermal Ion Imager");
    CDFcreateAttr(id, "Discipline", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "Space Physics>Ionospheric Science");
    CDFcreateAttr(id, "Generated_by", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "University of Calgary EFI team");
    CDFcreateAttr(id, "Generation_date", GLOBAL_SCOPE, &attrNum);
    time_t created;
    time(&created);
    addgEntry(id, attrNum, 0, ctime(&created));
    CDFcreateAttr(id, "LINK_TEXT", GLOBAL_SCOPE, &attrNum);
    if (highRate)
    {
        addgEntry(id, attrNum, 0, "16 Hz EFI TII ion drift and electric field data available at");
    }
    else
    {
        addgEntry(id, attrNum, 0, "2 Hz EFI TII ion drift and electric field data available at");
    }
    CDFcreateAttr(id, "LINK_TITLE", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "ESA Swarm Data Access");
    CDFcreateAttr(id, "HTTP_LINK", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "https://swarm-diss.eo.esa.int");
    CDFcreateAttr(id, "Instrument_type", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "Electric Fields (space)");
    addgEntry(id, attrNum, 1, "Particles (space)");
    CDFcreateAttr(id, "Instrument_type", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "Plasma and Solar Wind");
    CDFcreateAttr(id, "TEXT", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "Swarm electric field instrument thermal ion imager ion drift and electric field data.");
    addgEntry(id, attrNum, 1, "Satellite-trajectory coordinates: X: parallel to the satellite velocity vector; Y: to the right looking forward; Z: nadir");
    addgEntry(id, attrNum, 2, "Knudsen, D.J., Burchill, J.K., Buchert, S.C., Eriksson, A.I., Gill, R., Wahlund, J.E., Ahlen, L., Smith, M. and Moffat, B., 2017. Thermal ion imagers and Langmuir probes in the Swarm electric field instruments. Journal of Geophysical Research: Space Physics, 122(2), pp.2655-2673.");
    addgEntry(id, attrNum, 3, "Burchill, J.K. and Knudsen, D.J, 2020. EFI TII Cross-Track Flow Data Release Notes. Swarm DISC document SW-RN-UoC-GS-004, Revision 6.");
    CDFcreateAttr(id, "Time_resolution", GLOBAL_SCOPE, &attrNum);
    if (highRate)
    {
        addgEntry(id, attrNum, 0, "0.0625 seconds");
    }
    else
    {
        addgEntry(id, attrNum, 0, "0.5 seconds");
    }
    CDFcreateAttr(id, "TITLE", GLOBAL_SCOPE, &attrNum);
    sprintf(buf, "Swarm %s EFI TII High resolution data.", satellite);
    addgEntry(id, attrNum, 0, buf);
    CDFcreateAttr(id, "Project", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, "ESA Living Planet Programme");
    CDFcreateAttr(id, "Software_version", GLOBAL_SCOPE, &attrNum);
    addgEntry(id, attrNum, 0, SOFTWARE_VERSION_STRING);
    CDFcreateAttr(id, "spase_DatasetResourceID", GLOBAL_SCOPE, &attrNum);
    if (highRate)
    {
        sprintf(buf, "spase://ESA/Instrument/Swarm%s/TII/0.0625s", satellite);
    }
    else
    {
        sprintf(buf, "spase://ESA/Instrument/Swarm%s/TII/0.5s", satellite);
    }
    addgEntry(id, attrNum, 0, buf);

    CDFcreateAttr(id, "FIELDNAM", VARIABLE_SCOPE, &attrNum);
    CDFcreateAttr(id, "CATDESC", VARIABLE_SCOPE, &attrNum);
    CDFcreateAttr(id, "Type", VARIABLE_SCOPE, &attrNum);
    CDFcreateAttr(id, "UNITS", VARIABLE_SCOPE, &attrNum);
    CDFcreateAttr(id, "VAR_TYPE", VARIABLE_SCOPE, &attrNum);
    CDFcreateAttr(id, "DEPEND_0", VARIABLE_SCOPE, &attrNum);
    CDFcreateAttr(id, "DISPLAY_TYPE", VARIABLE_SCOPE, &attrNum);
    CDFcreateAttr(id, "LABLAXIS", VARIABLE_SCOPE, &attrNum);
    CDFcreateAttr(id, "VALIDMIN", VARIABLE_SCOPE, &attrNum);
    CDFcreateAttr(id, "VALIDMAX", VARIABLE_SCOPE, &attrNum);
    CDFcreateAttr(id, "TIME_BASE", VARIABLE_SCOPE, &attrNum);

    const varAttr variableAttrs[NUM_EXPORT_VARIABLES] = {
        {"Timestamp", "CDF_EPOCH", "*", "UT", minTime, maxTime},
        {"Latitude", "CDF_FLOAT", "degrees", "Geocentric latitude.", -90., 90.},
        {"Longitude", "CDF_FLOAT", "degrees", "Geocentric longitude.", -180., 180.},
        {"Radius", "CDF_FLOAT", "m", "Geocentric radius.", 6400000., 7000000.},
        {"QDLatitude", "CDF_FLOAT", "degrees", "Quasi-dipole magnetic latitude.", -90., 90.},
        {"MLT", "CDF_FLOAT", "hour", "Magnetic local time.", 0., 24.},
        {"Vixh", "CDF_FLOAT", "m/s", "Along-track ion drift from horizontal TII sensor in satellite-track coordinates.", -8000., 8000.},
        {"Vixh_error", "CDF_FLOAT", "m/s", "Random error estimate for along-track ion drift from horizontal TII sensor in satellite-track coordinates. Negative value indicates no estimate available.", -42., 10000.},
        {"Vixv", "CDF_FLOAT", "m/s", "Along-track ion drift from vertical TII sensor in satellite-track coordinates.", -8000., 8000.},
        {"Vixv_error", "CDF_FLOAT", "m/s", "Random error estimate for along-track ion drift from vertical TII sensor in satellite-track coordinates. Negative value indicates no estimate available.", -42., 10000.},
        {"Viy", "CDF_FLOAT", "m/s", "Cross-track horizontal ion drift from horizontal TII sensor in satellite-track coordinates.", -8000., 8000.},
        {"Viy_error", "CDF_FLOAT", "m/s", "Random error estimate for cross-track horizontal ion drift from horizontal TII sensor in satellite-track coordinates. Negative value indicates no estimate available.", -42., 10000.},
        {"Viz", "CDF_FLOAT", "m/s", "Cross-track vertical ion drift from vertical TII sensor in satellite-track coordinates.", -8000., 8000.},
        {"Viz_error", "CDF_FLOAT", "m/s", "Random error estimate for cross-track vertical ion drift from vertical TII sensor in satellite-track coordinates. Negative value indicates no estimate available.", -42., 10000.},
        {"VsatN", "CDF_FLOAT", "m/s", "Satellite velocity N component in north-east-centre coordinates.", -8000., 8000.},
        {"VsatE", "CDF_FLOAT", "m/s", "Satellite velocity E component in north-east-centre coordinates.", -8000., 8000.},
        {"VsatC", "CDF_FLOAT", "m/s", "Satellite velocity C component in north-east-centre coordinates.", -200., 200.},
        {"Ehx", "CDF_FLOAT", "mV/m", "Electric field x component in satellite-track coordinates, derived from -VxB with along-track ion drift from horizontal sensor.", -400., 400.},
        {"Ehy", "CDF_FLOAT", "mV/m", "Electric field y component in satellite-track coordinates, derived from -VxB with along-track ion drift from horizontal sensor.", -400., 400.},
        {"Ehz", "CDF_FLOAT", "mV/m", "Electric field z component in satellite-track coordinates, derived from -VxB with along-track ion drift from horizontal sensor.", -400., 400.},
        {"Evx", "CDF_FLOAT", "mV/m", "Electric field x component in satellite-track coordinates, derived from -VxB with along-track ion drift from vertical sensor.", -400., 400.},
        {"Evy", "CDF_FLOAT", "mV/m", "Electric field y component in satellite-track coordinates, derived from -VxB with along-track ion drift from vertical sensor.", -400., 400.},
        {"Evz", "CDF_FLOAT", "mV/m", "Electric field z component in satellite-track coordinates, derived from -VxB with along-track ion drift from vertical sensor.", -400., 400.},
        {"Bx", "CDF_FLOAT", "nT", "Geomagnetic field x component in satellite-track coordinates, derived from the 1 Hz product.", -65000., 65000.},
        {"By", "CDF_FLOAT", "nT", "Geomagnetic field y component in satellite-track coordinates, derived from the 1 Hz product.", -65000., 65000.},
        {"Bz", "CDF_FLOAT", "nT", "Geomagnetic field z component in satellite-track coordinates, derived from the 1 Hz product.", -65000., 65000.},
        {"Vicrx", "CDF_FLOAT", "m/s", "Ion drift corotation signal x component in satellite-track coorinates. This has been removed from ion drift and electric field.", -1000., 1000.},
        {"Vicry", "CDF_FLOAT", "m/s", "Ion drift corotation signal y component in satellite-track coorinates. This has been removed from ion drift and electric field.", -1000., 1000.},
        {"Vicrz", "CDF_FLOAT", "m/s", "Ion drift corotation signal z component in satellite-track coorinates. This has been removed from ion drift and electric field.", -1000., 1000.},
        {"Quality_flags", "CDF_UINT2", "*", "Bitwise flag for each velocity component, where a value of 1 for a particular component signifies that calibration was successful, and that the baseline 1-sigma noise level is less than or equal to 100 m/s at 2 Hz. Electric field quality can be assessed from these flags according to -vxB. Bit0 (least significant) = Vixh, bit1 = Vixv, bit2 = Viy, bit3 = Viz. Refer to the release notes for details.", 0, 65535},
        {"Calibration_flags", "CDF_UINT4", "*", "Information about the calibration process. Refer to the release notes for details.", 0, 4294967295}
    };

    for (uint8_t i = 0; i < NUM_EXPORT_VARIABLES; i++)
    {
        addVariableAttributes(id, variableAttrs[i]);
    }

}


void loadCalData(const char *calDir, const char *calVersion, const char *satellite, const int year, const int month, const int day, uint8_t **dataBuffers, long *nRecs)
{
    // c time manipulation: see  https://fresh2refresh.com/c-programming/c-time-related-functions/
    struct tm timestructure;
    time_t date;
 
    timestructure.tm_year = year - 1900;
    timestructure.tm_mon = month - 1;
    timestructure.tm_mday = day - 1;
    timestructure.tm_hour = 0;
    timestructure.tm_min = 0;
    timestructure.tm_sec = 0;
    timestructure.tm_isdst = 0;

    long calibrationMemorySize = 0;
    for (int8_t i = 0; i < 3; i++)
    {
        mktime(&timestructure);
        fprintf(stdout, "%sLoading calibration data for %04d%02d%02d\n", infoHeader, timestructure.tm_year + 1900, timestructure.tm_mon + 1, timestructure.tm_mday);
        loadCalDataFromDate((DayType)(i-1), calDir, calVersion, satellite, timestructure, dataBuffers, nRecs, &calibrationMemorySize);
        timestructure.tm_mday++;
    }
    fprintf(stdout, "%sNumber of records: %ld\n", infoHeader, *nRecs);
    fprintf(stdout, "%sLoaded %ld bytes (%ld MB) of calibration data.\n", infoHeader, calibrationMemorySize, calibrationMemorySize / 1024 / 1024);

}

void setCalibrationFileName(const char *satellite, const int year, const int month, const int day, const char *calDir, const char *calVersion, char *calibrationFileName)
{
    sprintf(calibrationFileName, "%s/%s/%04d/Swarm_%s/%02d/TiiClbr%s_Swarm_%s_%04d_%02d_%02d.cdf", calDir, calVersion, year, satellite, month, calVersion, satellite, year, month, day);
}

void loadCalDataFromDate(const DayType dayType, const char *calDir, const char *calVersion, const char *satellite, struct tm timestructure,  uint8_t **dataBuffers, long *numberOfRecords, long *totalMemoryAllocated)
{
    char calibrationFileName[CDF_PATHNAME_LEN];
    int year = timestructure.tm_year + 1900;
    int month = timestructure.tm_mon + 1;
    int day = timestructure.tm_mday;
    setCalibrationFileName(satellite, year, month, day, calDir, calVersion, calibrationFileName);
    fprintf(stdout, "%s from %s\n", infoHeader, calibrationFileName);

    // Open the CDF file with validation
    CDFsetValidate(VALIDATEFILEon);
    CDFid calCdfId;
    CDFstatus status;
    status = CDFopenCDF(calibrationFileName, &calCdfId);
    if (status != CDF_OK) 
    {
        printErrorMessage(status);
        // Not necessarily an error. For example, some dates will have not calibration data.
        fprintf(stdout, "%sSkipping this date.\n", infoHeader);
        return;
    }

    fprintf(stdout, "%sFound CDF file.\n", infoHeader);

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
        fprintf(stdout, "%sProblem with calibration file. Skipping this date.\n", infoHeader);
        closeCdf(calCdfId);
        return;
    }
    long nRecs, calibrationMemorySize = 0;
    status = CDFgetzVarAllocRecords(calCdfId, CDFgetVarNum(calCdfId, "epoch"), &nRecs);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        fprintf(stdout, "%sProblem with calibration file. Skipping this date.\n", infoHeader);
        closeCdf(calCdfId);
        return;
    }
    if (nRecs < (16*SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING))
    {
        // Not enough to do anything useful 
        // TODO: maybe increase this threshold to require a larger number of points each day?
        fprintf(stdout, "%sFewer than %.0f s of data. Skipping this date.\n", infoHeader, (float)SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING);
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
            timeReference = computeEPOCH(year, month, day + 1, 0, 0, -(SECONDS_OF_BOUNDARY_DATA_REQUIRED_FOR_PROCESSING), 0);
            for (startRecord = stopRecord; startRecord >= 0; startRecord--)
            {
                status = CDFgetzVarRecordData(calCdfId, epochNum, startRecord, &recordTime);
                if (status != CDF_OK)
                {
                    fprintf(stdout, "%sCould not read epoch record from CDF file. Skipping this calibration date.\n", infoHeader);
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
            timeReference = computeEPOCH(year, month, day, 0, 0, (SECONDS_OF_BOUNDARY_DATA_REQUIRED_FOR_PROCESSING), 0);
            for (stopRecord = startRecord; stopRecord < nRecs; stopRecord++)
            {
                status = CDFgetzVarRecordData(calCdfId, epochNum, stopRecord, &recordTime);
                if (status != CDF_OK)
                {
                    fprintf(stdout, "%sCould not read epoch record from CDF file. Skipping this calibration date.\n", infoHeader);
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
            fprintf(stdout, "%sError: Day type must be one of PREVIOUS_DAY, REQUESTED_DAY, or NEXT_DAY. Skipping this calibration data.\n", infoHeader);
            closeCdf(calCdfId);
            return;
    }

    // Update number of records
    nRecs = stopRecord - startRecord + 1;
    if ((dayType == REQUESTED_DAY && nRecs < (16*SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING)) || ((dayType == PREVIOUS_DAY || dayType == NEXT_DAY) && nRecs < (16*SECONDS_OF_BOUNDARY_DATA_REQUIRED_FOR_PROCESSING)))
    {
        // Not enough to do anything useful 
        // TODO: maybe increase this threshold to require a larger number of points each day?
        fprintf(stdout, "%sFewer than %.0f s of data meet constraints. Skipping this date.\n", infoHeader, (float)SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING);
        CDFcloseCDF(calCdfId);
        return;
    }

    // Variables
    char* variables[] = {
        "epoch",                // index 0
        "1st Moment - H",       // 1
        "1st Moment - V",       // 2
        "Det_H__vX",            // 3
        "MLT",                  // 4
        "QDLat",                // 5
        "QDLon",                // 6
        "Latitude",             // 7
        "Long Tle",             // 8
        "Radius",               // 9
        "Det H CorVx",          // 10
        "Sat_Vel_N",            // 11
        "B-North",              // 12
        "Bias_Grid_H",          // 13
        "Faceplate_Volt_Mon_H"  // 14
        };

    uint8_t nVars = numzVars;
    if (nVars != NUM_CAL_VARIABLES)
    {
        fprintf(stdout, "%sError: number of calibration variables should be %d. Got %ld. Skipping this date.\n", infoHeader, (uint8_t) NUM_CAL_VARIABLES, numzVars);
        closeCdf(calCdfId);
        return;
    }
    fprintf(stdout, "%sChecking calibration file variables...", infoHeader);
    for (uint8_t i = 0; i<nVars; i++)
    {
        // fprintf(stdout, "%s%20s ", infoHeader, variables[i]);
        status = CDFconfirmzVarExistence(calCdfId, variables[i]);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            fprintf(stdout, "%sError reading variable %s. Skipping this date.\n", infoHeader, variables[i]);
            closeCdf(calCdfId);
            return;
        }
        else
        {
            // fprintf(stdout, "%s OK\n", infoHeader);
        }
    }
    fprintf(stdout, "%sOK\n", infoHeader);
    
    long varNum, numValues, numVarBytes;
    long numBytesPrev, numBytesToAdd, numBytesNew;

    for (uint8_t i = 0; i < nVars; i++)
    {
        varNum = CDFgetVarNum(calCdfId, variables[i]);
        if (varNum < CDF_OK)
        {
            printErrorMessage(varNum);
            fprintf(stdout, "%sError reading variable ID for %s. Skipping this date.\n", infoHeader, variables[i]);
            closeCdf(calCdfId);
            return;
        }
        status = CDFgetzVarNumDims(calCdfId, varNum, &numDims);
        status = CDFgetzVarDimSizes(calCdfId, varNum, dimSizes);
        status = CDFgetzVarDataType(calCdfId, varNum, &dataType);
        // Calculate new size of memory to allocate
        status = CDFgetDataTypeSize(dataType, &numVarBytes);
        numValues = 1;
        for (uint8_t j = 0; j < numDims; j++)
        {
            numValues *= dimSizes[j];
        }
        numBytesPrev = numValues * (*numberOfRecords) * numVarBytes;
        numBytesToAdd = numValues * nRecs * numVarBytes;
        numBytesNew = numBytesPrev + numBytesToAdd;
        calibrationMemorySize += numBytesNew;
        dataBuffers[i] = (uint8_t*) realloc(dataBuffers[i], (size_t) numBytesNew);
        status = CDFgetzVarRangeRecordsByVarID(calCdfId, varNum, startRecord, stopRecord, dataBuffers[i] + numBytesPrev);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            fprintf(stdout, "%sError loading data for %s. Skipping this date.\n", infoHeader, variables[i]);
            closeCdf(calCdfId);
            return;
        }
    }
    // close CDF
    closeCdf(calCdfId);
    // Number of records obtained for this date
    switch (dayType)
    {
        case PREVIOUS_DAY:
            fprintf(stdout, "%sGot %ld s of data for previous day\n", infoHeader, nRecs / 16);
            break;
        case REQUESTED_DAY:
            fprintf(stdout, "%sGot %ld s of data for requested day\n", infoHeader, nRecs / 16);
            break;
        case NEXT_DAY:
            fprintf(stdout, "%sGot %ld s of data for next day\n", infoHeader, nRecs / 16);
            break;
        default:
            break;
    }
    // Update number of records found and memory allocated
    *numberOfRecords += nRecs;
    *totalMemoryAllocated = calibrationMemorySize;

}
long numberOfAvailableRecordsForDate(const char *satellite, const int year, const int month, const int day, const char *calDir, const char *calVersion)
{
    char calibrationFileName[CDF_PATHNAME_LEN];
    setCalibrationFileName(satellite, year, month, day, calDir, calVersion, calibrationFileName);
    CDFid calCdfId;
    CDFstatus status;
    status = CDFopenCDF(calibrationFileName, &calCdfId);
    if (status != CDF_OK) 
    {
        // Not necessarily an error. For example, some dates will have not calibration data.
        printErrorMessage(status);
        return 0;
    }

    // Get number of records for zVar "epoch"
    long nRecords = 0;
    status = CDFgetzVarAllocRecords(calCdfId, CDFgetVarNum(calCdfId, "epoch"), &nRecords);
    if (status != CDF_OK) 
    {
        printErrorMessage(status);
        nRecords = 0;
    }
    closeCdf(calCdfId);
    return nRecords;

}

void exportTCT16Cdfs(double startTime, double stopTime, const char *exportDir, const char *exportVersion, const char *calVersion, const char *satellite, long startIndex, long stopIndex, uint8_t **dataBuffers, float *ectFieldH, float *ectFieldV, float* bctField, float *viErrors, uint16_t *flags, uint32_t *fitInfo)
{
    fprintf(stdout, "%sExporting 16 Hz data.\n",infoHeader);
    char epochString[EPOCH_STRING_LEN+1];
    toEncodeEPOCH(startTime, 0, epochString);
    fprintf(stdout, "%sStartTime: %s\n", infoHeader, epochString);
    toEncodeEPOCH(stopTime, 0, epochString);
    fprintf(stdout, "%sStoptime: %s\n", infoHeader, epochString);

    char cdfFileName[CDF_PATHNAME_LEN];
    constructExportFileName("TCT16", startTime, stopTime, exportDir, exportVersion, satellite, cdfFileName);
    
    CDFid exportCdfId;
    CDFstatus status;
    status = CDFcreateCDF(cdfFileName, &exportCdfId);
    if (status != CDF_OK) 
    {
        printErrorMessage(status);
        // Close export file
        closeCdf(exportCdfId);
    }
    else
    {
        // export variables
        createVarFrom1DVar(exportCdfId, "Timestamp", CDF_EPOCH, startIndex, stopIndex, dataBuffers[0]);
        createVarFrom1DVar(exportCdfId, "Latitude", CDF_REAL4, startIndex, stopIndex, dataBuffers[7]);
        createVarFrom1DVar(exportCdfId, "Longitude", CDF_REAL4, startIndex, stopIndex, dataBuffers[8]);
        createVarFrom1DVar(exportCdfId, "Radius", CDF_REAL4, startIndex, stopIndex, dataBuffers[9]);
        createVarFrom1DVar(exportCdfId, "QDLatitude", CDF_REAL4, startIndex, stopIndex, dataBuffers[5]);
        createVarFrom1DVar(exportCdfId, "MLT", CDF_REAL4, startIndex, stopIndex, dataBuffers[4]);
        createVarFrom2DVar(exportCdfId, "Vixh", CDF_REAL4, startIndex, stopIndex, dataBuffers[1], 0, 2);
        createVarFrom2DVar(exportCdfId, "Vixh_error", CDF_REAL4, startIndex, stopIndex, viErrors, 0, 4);
        createVarFrom2DVar(exportCdfId, "Vixv", CDF_REAL4, startIndex, stopIndex, dataBuffers[2], 0, 2);
        createVarFrom2DVar(exportCdfId, "Vixv_error", CDF_REAL4, startIndex, stopIndex, viErrors, 2, 4);
        createVarFrom2DVar(exportCdfId, "Viy", CDF_REAL4, startIndex, stopIndex, dataBuffers[1], 1, 2);
        createVarFrom2DVar(exportCdfId, "Viy_error", CDF_REAL4, startIndex, stopIndex, viErrors, 1, 4);
        createVarFrom2DVar(exportCdfId, "Viz", CDF_REAL4, startIndex, stopIndex, dataBuffers[2], 1, 2);
        createVarFrom2DVar(exportCdfId, "Viz_error", CDF_REAL4, startIndex, stopIndex, viErrors, 3, 4);
        createVarFrom2DVar(exportCdfId, "VsatN", CDF_REAL4, startIndex, stopIndex, dataBuffers[11], 0, 3);
        createVarFrom2DVar(exportCdfId, "VsatE", CDF_REAL4, startIndex, stopIndex, dataBuffers[11], 1, 3);
        createVarFrom2DVar(exportCdfId, "VsatC", CDF_REAL4, startIndex, stopIndex, dataBuffers[11], 2, 3);
        createVarFrom2DVar(exportCdfId, "Ehx", CDF_REAL4, startIndex, stopIndex, ectFieldH, 0, 3);
        createVarFrom2DVar(exportCdfId, "Ehy", CDF_REAL4, startIndex, stopIndex, ectFieldH, 1, 3);
        createVarFrom2DVar(exportCdfId, "Ehz", CDF_REAL4, startIndex, stopIndex, ectFieldH, 2, 3);
        createVarFrom2DVar(exportCdfId, "Evx", CDF_REAL4, startIndex, stopIndex, ectFieldV, 0, 3);
        createVarFrom2DVar(exportCdfId, "Evy", CDF_REAL4, startIndex, stopIndex, ectFieldV, 1, 3);
        createVarFrom2DVar(exportCdfId, "Evz", CDF_REAL4, startIndex, stopIndex, ectFieldV, 2, 3);
        createVarFrom2DVar(exportCdfId, "Bx", CDF_REAL4, startIndex, stopIndex, bctField, 0, 3);
        createVarFrom2DVar(exportCdfId, "By", CDF_REAL4, startIndex, stopIndex, bctField, 1, 3);
        createVarFrom2DVar(exportCdfId, "Bz", CDF_REAL4, startIndex, stopIndex, bctField, 2, 3);
        createVarFrom2DVar(exportCdfId, "Vicrx", CDF_REAL4, startIndex, stopIndex, dataBuffers[10], 0, 3);
        createVarFrom2DVar(exportCdfId, "Vicry", CDF_REAL4, startIndex, stopIndex, dataBuffers[10], 1, 3);
        createVarFrom2DVar(exportCdfId, "Vicrz", CDF_REAL4, startIndex, stopIndex, dataBuffers[10], 2, 3);
        createVarFrom1DVar(exportCdfId, "Quality_flags", CDF_UINT2, startIndex, stopIndex, flags);
        createVarFrom1DVar(exportCdfId, "Calibration_flags", CDF_UINT4, startIndex, stopIndex, fitInfo);

        // add attributes
        addAttributes(exportCdfId, "TCT16", satellite, exportVersion, startTime, stopTime);

        // Close export file
        closeCdf(exportCdfId);
        fprintf(stdout, "%sExported %ld records to %s.cdf\n", infoHeader, (stopIndex - startIndex + 1), cdfFileName);
        fflush(stdout);

    }
}

void exportTCT02Cdfs(double startTime, double stopTime, const char *exportDir, const char *exportVersion, const char *calVersion, const char *satellite, long startIndex, long stopIndex, uint8_t **dataBuffers, float *ectFieldH, float *ectFieldV, float* bctField, float *viErrors, uint16_t *flags, uint32_t *fitInfo)
{
    fprintf(stdout, "%sExporting 2 Hz data.\n",infoHeader);

    // Average the data to 2 Hz from 16 Hz
    // In principle all 16 Hz data come from a single instrument source packet (ISP) 
    // so we can safely assume the number of samples is a multiple of 8
    // This means we could simply average each 8 measurements together as assign the 
    // time for each measurement as the average time for the 8 samples
    // But this does not take acount day boundaries, where the first measurement of the day 
    // does not correspond to the beginning of a chunk of 8 samples because 
    // the times have been adjusted for sampling lag.

    // The averaging strategy is therefore to keep only those 2 Hz samples for which we have
    // 8 16 Hz samples in the interval [0, 0.5) or [0.5, 1.0) for each UT second.
    // Any samples for which 8 were not available are discarded, but can be obtained 
    // for science from the 16 Hz data.

    // No attempt is made in this version of the software to remove outliers.
    fprintf(stdout, "%sDown-sampling 16 Hz to 2 Hz.\n",infoHeader);

    long timeIndex;
    long n2HzSamples = 0;
    long storageIndex = startIndex;
    double t0;
    bool downSampled = false;
    for (timeIndex = startIndex; timeIndex <= stopIndex;) // Time index advanced below
    {
        t0 = floor(TIME()/1000.0); // UT second reference
        for (uint8_t halfSecond = 0; halfSecond < 2; halfSecond ++)
        {
            downSampled = downSampleHalfSecond(&timeIndex, storageIndex, t0 + 0.5 * halfSecond, stopIndex, dataBuffers, ectFieldH, ectFieldV, bctField, viErrors, flags, fitInfo);
            if (downSampled)
            {
                storageIndex++;
                n2HzSamples++;
            }
        }

    }
    // Update start and stop indexes and times
    stopIndex = startIndex + n2HzSamples - 1;
    startTime = *((double*)dataBuffers[0] + (startIndex));
    stopTime = *((double*)dataBuffers[0] + (stopIndex));

    char epochString[EPOCH_STRING_LEN+1];
    toEncodeEPOCH(startTime, 0, epochString);
    fprintf(stdout, "%sStartTime: %s\n", infoHeader, epochString);
    toEncodeEPOCH(stopTime, 0, epochString);
    fprintf(stdout, "%sStoptime: %s\n", infoHeader, epochString);

    CDFid exportCdfId;
    CDFstatus status;
    char cdfFileName[CDF_PATHNAME_LEN];
    constructExportFileName("TCT02", startTime, stopTime, exportDir, exportVersion, satellite, cdfFileName);
    status = CDFcreateCDF(cdfFileName, &exportCdfId);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        // Close export file
        closeCdf(exportCdfId);
    }
    else
    {
        // export variables
        createVarFrom1DVar(exportCdfId, "Timestamp", CDF_EPOCH, startIndex, stopIndex, dataBuffers[0]);
        createVarFrom1DVar(exportCdfId, "Latitude", CDF_REAL4, startIndex, stopIndex, dataBuffers[7]);
        createVarFrom1DVar(exportCdfId, "Longitude", CDF_REAL4, startIndex, stopIndex, dataBuffers[8]);
        createVarFrom1DVar(exportCdfId, "Radius", CDF_REAL4, startIndex, stopIndex, dataBuffers[9]);
        createVarFrom1DVar(exportCdfId, "QDLatitude", CDF_REAL4, startIndex, stopIndex, dataBuffers[5]);
        createVarFrom1DVar(exportCdfId, "MLT", CDF_REAL4, startIndex, stopIndex, dataBuffers[4]);
        createVarFrom2DVar(exportCdfId, "Vixh", CDF_REAL4, startIndex, stopIndex, dataBuffers[1], 0, 2);
        createVarFrom2DVar(exportCdfId, "Vixh_error", CDF_REAL4, startIndex, stopIndex, viErrors, 0, 4);
        createVarFrom2DVar(exportCdfId, "Vixv", CDF_REAL4, startIndex, stopIndex, dataBuffers[2], 0, 2);
        createVarFrom2DVar(exportCdfId, "Vixv_error", CDF_REAL4, startIndex, stopIndex, viErrors, 2, 4);
        createVarFrom2DVar(exportCdfId, "Viy", CDF_REAL4, startIndex, stopIndex, dataBuffers[1], 1, 2);
        createVarFrom2DVar(exportCdfId, "Viy_error", CDF_REAL4, startIndex, stopIndex, viErrors, 1, 4);
        createVarFrom2DVar(exportCdfId, "Viz", CDF_REAL4, startIndex, stopIndex, dataBuffers[2], 1, 2);
        createVarFrom2DVar(exportCdfId, "Viz_error", CDF_REAL4, startIndex, stopIndex, viErrors, 3, 4);
        createVarFrom2DVar(exportCdfId, "VsatN", CDF_REAL4, startIndex, stopIndex, dataBuffers[11], 0, 3);
        createVarFrom2DVar(exportCdfId, "VsatE", CDF_REAL4, startIndex, stopIndex, dataBuffers[11], 1, 3);
        createVarFrom2DVar(exportCdfId, "VsatC", CDF_REAL4, startIndex, stopIndex, dataBuffers[11], 2, 3);
        createVarFrom2DVar(exportCdfId, "Ehx", CDF_REAL4, startIndex, stopIndex, ectFieldH, 0, 3);
        createVarFrom2DVar(exportCdfId, "Ehy", CDF_REAL4, startIndex, stopIndex, ectFieldH, 1, 3);
        createVarFrom2DVar(exportCdfId, "Ehz", CDF_REAL4, startIndex, stopIndex, ectFieldH, 2, 3);
        createVarFrom2DVar(exportCdfId, "Evx", CDF_REAL4, startIndex, stopIndex, ectFieldV, 0, 3);
        createVarFrom2DVar(exportCdfId, "Evy", CDF_REAL4, startIndex, stopIndex, ectFieldV, 1, 3);
        createVarFrom2DVar(exportCdfId, "Evz", CDF_REAL4, startIndex, stopIndex, ectFieldV, 2, 3);
        createVarFrom2DVar(exportCdfId, "Bx", CDF_REAL4, startIndex, stopIndex, bctField, 0, 3);
        createVarFrom2DVar(exportCdfId, "By", CDF_REAL4, startIndex, stopIndex, bctField, 1, 3);
        createVarFrom2DVar(exportCdfId, "Bz", CDF_REAL4, startIndex, stopIndex, bctField, 2, 3);
        createVarFrom2DVar(exportCdfId, "Vicrx", CDF_REAL4, startIndex, stopIndex, dataBuffers[10], 0, 3);
        createVarFrom2DVar(exportCdfId, "Vicry", CDF_REAL4, startIndex, stopIndex, dataBuffers[10], 1, 3);
        createVarFrom2DVar(exportCdfId, "Vicrz", CDF_REAL4, startIndex, stopIndex, dataBuffers[10], 2, 3);
        createVarFrom1DVar(exportCdfId, "Quality_flags", CDF_UINT2, startIndex, stopIndex, flags);
        createVarFrom1DVar(exportCdfId, "Calibration_flags", CDF_UINT4, startIndex, stopIndex, fitInfo);

        // add attributes
        // update start and stop times to the averaged ones
        addAttributes(exportCdfId, "TCT02", satellite, exportVersion, startTime, stopTime);

        // Close export file
        closeCdf(exportCdfId);
        fprintf(stdout, "%sExported %ld records to %s.cdf\n", infoHeader, (stopIndex - startIndex + 1), cdfFileName);
        fflush(stdout);

    }
}

bool downSampleHalfSecond(long *index, long storageIndex, double t0, long maxIndex, uint8_t **dataBuffers, float *ectFieldH, float *ectFieldV, float *bctField, float *viErrors, uint16_t *flags, uint32_t *fitInfo)
{
    long timeIndex = *index;
    uint8_t nSamples = 0;
    double timeBuf = 0.0;
    uint16_t flagBuf = 65535;
    uint32_t fitInfoBuf = 0;
    float floatBuf[NUM_BUFFER_VARIABLES]; // Time, flags, and fitInfo are handled separately; need one extra working buffer each for lattitude and longitude averages
    float theta, phi, x, y, z;
    bool downSampled = false;

    for (uint8_t b = 0; b < NUM_BUFFER_VARIABLES; b++)
    {
        floatBuf[b] = 0.0;
    }

    while (((TIME()/1000. - t0) < 0.5) && (timeIndex <= maxIndex))
    {
        timeBuf += TIME();
        // Handle lat and lon and mlt in cartesian coordinates
        // For spherical coordinates
        theta = M_PI / 2.0 - LAT() * M_PI / 180.0;
        phi = LON() * M_PI / 180.0;
        floatBuf[0] += cosf(phi) * sinf(theta);
        floatBuf[1] += sinf(phi) * sinf(theta);
        floatBuf[2] += cosf(theta);
        floatBuf[3] += RADIUS();
        // For magnetic coordinates
        theta = M_PI / 2.0 - QDLAT() * M_PI / 180.0;
        phi = MLT() / 24.0 * 2.0 * M_PI;
        floatBuf[4] += cosf(phi) * sinf(theta);
        floatBuf[5] += sinf(phi) * sinf(theta);
        floatBuf[6] += cosf(theta);
        floatBuf[7] += MXH();
        floatBuf[8] += viErrors[timeIndex * 4 + 0];
        floatBuf[9] += MXV();
        floatBuf[10] += viErrors[timeIndex * 4 + 2];
        floatBuf[11] += MYH();
        floatBuf[12] += viErrors[timeIndex * 4 + 1];
        floatBuf[13] += MYV();
        floatBuf[14] += viErrors[timeIndex * 4 + 3];
        floatBuf[15] += VSATN();
        floatBuf[16] += VSATE();
        floatBuf[17] += VSATC();
        floatBuf[18] += ectFieldH[timeIndex * 3 + 0]; // EH xyz
        floatBuf[19] += ectFieldH[timeIndex * 3 + 1];
        floatBuf[20] += ectFieldH[timeIndex * 3 + 2];
        floatBuf[21] += ectFieldV[timeIndex * 3 + 0]; // EV xyz
        floatBuf[22] += ectFieldV[timeIndex * 3 + 1];
        floatBuf[23] += ectFieldV[timeIndex * 3 + 2];
        floatBuf[24] += bctField[timeIndex * 3 + 0]; // Bxyz
        floatBuf[25] += bctField[timeIndex * 3 + 1];
        floatBuf[26] += bctField[timeIndex * 3 + 2];
        floatBuf[27] += VCORX();
        floatBuf[28] += VCORY();
        floatBuf[29] += VCORZ();
        flagBuf &= flags[timeIndex];
        fitInfoBuf |= fitInfo[timeIndex];
        nSamples++;
        timeIndex++;
    }
    if (nSamples == 8)
    {
        // do the averaging and store result in original array
        *((double*)dataBuffers[0] + (storageIndex)) = timeBuf / 8.0; // Average time
        x = floatBuf[0] / 8.0;
        y = floatBuf[1] / 8.0;
        z = floatBuf[2] / 8.0;
        *((float*)dataBuffers[7] + (storageIndex)) = atan2f(z, sqrtf(x*x + y*y)) * 180.0 / M_PI; // Lat
        *((float*)dataBuffers[8] + (storageIndex)) = fmodf(atan2f(y, x) * 180.0 / M_PI, 360.0); // Lon
        *((float*)dataBuffers[9] + (storageIndex)) = floatBuf[3] / 8.0; // radius
        x = floatBuf[4] / 8.0;
        y = floatBuf[5] / 8.0;
        z = floatBuf[6] / 8.0;
        *((float*)dataBuffers[5] + (storageIndex)) = atan2f(z, sqrtf(x*x + y*y)) * 180.0 / M_PI; // QD Lat
        *((float*)dataBuffers[4] + (storageIndex)) = fmodf(atan2f(y, x) * 180.0 / M_PI + 360.0, 360.0) / 360. * 24.0; // MLT
        *((float*)dataBuffers[1] + (2*storageIndex) + 0) = floatBuf[7] / 8.0; // VxH
        viErrors[storageIndex * 4 + 0] = floatBuf[8] / 8.0 / sqrtf(8.0); // VxH error scaled for 2 Hz
        *((float*)dataBuffers[2] + (2*storageIndex) + 0) = floatBuf[9] / 8.0; // VxV
        viErrors[storageIndex * 4 + 2] = floatBuf[10] / 8.0 / sqrtf(8.0); // VxV error scaled for 2 Hz
        *((float*)dataBuffers[1] + (2*storageIndex) + 1) = floatBuf[11] / 8.0; // VyH
        viErrors[storageIndex * 4 + 1] = floatBuf[12] / 8.0 / sqrtf(8.0); // VyH error scaled for 2 Hz
        *((float*)dataBuffers[2] + (2*storageIndex) + 1) = floatBuf[13] / 8.0; // VyV
        viErrors[storageIndex * 4 + 3] = floatBuf[14] / 8.0 / sqrtf(8.0); // VyV error scaled for 2 Hz
        *((float*)dataBuffers[11] + (3*storageIndex) + 0) = floatBuf[15] / 8.0; // VSatN
        *((float*)dataBuffers[11] + (3*storageIndex) + 1) = floatBuf[16] / 8.0; // VSatE
        *((float*)dataBuffers[11] + (3*storageIndex) + 2) = floatBuf[17] / 8.0; // VSatC
        ectFieldH[storageIndex * 3 + 0] = floatBuf[18] / 8.0; // EHxyz
        ectFieldH[storageIndex * 3 + 1] = floatBuf[19] / 8.0;
        ectFieldH[storageIndex * 3 + 2] = floatBuf[20] / 8.0;
        ectFieldV[storageIndex * 3 + 0] = floatBuf[21] / 8.0; // EVxyz
        ectFieldV[storageIndex * 3 + 1] = floatBuf[22] / 8.0;
        ectFieldV[storageIndex * 3 + 2] = floatBuf[23] / 8.0;
        bctField[storageIndex * 3 + 0] = floatBuf[24] / 8.0; // Bxyz
        bctField[storageIndex * 3 + 1] = floatBuf[25] / 8.0;
        bctField[storageIndex * 3 + 2] = floatBuf[26] / 8.0;
        *((float*)dataBuffers[10] + (3*storageIndex) + 0) = floatBuf[27] / 8.0; // Vicrxyz
        *((float*)dataBuffers[10] + (3*storageIndex) + 1) = floatBuf[28] / 8.0;
        *((float*)dataBuffers[10] + (3*storageIndex) + 2) = floatBuf[29] / 8.0;
        // Flags set to 0 at 16 Hz based on magnitude of flow,
        // are not reset at 2 Hz, to ensure integrity of 2 Hz measurements
        // One can review 16 Hz measurements to examine details of flow where even a
        // single sample of the eight has a magnitude greater than 8 km/s
        flags[storageIndex] = flagBuf;
        fitInfo[storageIndex] = fitInfoBuf;

        downSampled = true;

    }

    *index = timeIndex;
    return downSampled;

}

