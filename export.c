/*

    TII Cross-Track Ion Drift Processor: export.c

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

#include "export.h"

#include "loadData.h"
#include "settings.h"
#include "errors.h"
#include "processing.h"
#include "state.h"
#include "utilities.h"

#include <cdf.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <ctype.h>

extern char infoHeader[50];

int exportCdfs(ProcessorState *state)
{

    int status = TIICT_OK;

    ProcessorVariables_t *v = &state->vars16hz;

    // Keep only data from the requested day
    double minTime, maxTime, startTime, stopTime, prevTime;
    long startIndex = 0, stopIndex = v->nRecs;
    minTime = computeEPOCH(state->args.year, state->args.month, state->args.day, 0, 0, 0, 0);
    maxTime = computeEPOCH(state->args.year, state->args.month, state->args.day + 1, 0, 0, 0, 0);
    bool gotStartTime = false, gotStopTime = false;
    char errorMessage[CDF_STATUSTEXT_LEN + 1];
    bool exporting = true;
    long recordsExported = 0;
    double minutesExported = 0.;
    uint16_t filesExported = 0;

    // Find first index with time on requested day
    int i = 0;
    while (i < v->nRecs && v->timestamp[i] < minTime) i++;
    if (i == v->nRecs)
        return TIICT_NO_RECORDS_TO_EXPORT;

    startIndex = i;
    startTime = v->timestamp[i];
    // Find last index not containing a gap of more than 10 minutes, or last index of the day
    stopIndex = startIndex;
    stopTime = startTime;

    // case where there are no records on the requested or following days
    if (startIndex == v->nRecs)
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sNo records found for requested date.\n", infoHeader);
        }
        exporting = false;
    }
    while (exporting)
    {
        double gapTime = (v->timestamp[i] - stopTime)/1000.;
        double duration = (stopTime - startTime)/1000.;
        if ((gapTime > MAX_ALLOWED_CDF_GAP_SECONDS) || (v->timestamp[i] >= maxTime) || (i == (v->nRecs-1)))
        {
            // Export data and continue even on error (for this interval)
            if (duration >= SECONDS_OF_DATA_REQUIRED_FOR_EXPORTING)
            {
                // 16 Hz data
                if (state->export16Hz) {
                    status |= exportTCT16Cdfs(state, startTime, stopTime, startIndex, stopIndex);
                }
                // 2 Hz data
                if (state->export2Hz) {
                    status |= exportTCT02Cdfs(state, startTime, stopTime, startIndex, stopIndex);
                }

                recordsExported += (stopIndex - startIndex + 1);
                minutesExported += (stopTime - startTime)/1000./60.;
                filesExported++; // Only count 16 Hz files
            }
            else
            {
                if (state->writeLogFiles) {
                    char startString[EPOCH_STRING_LEN+1], stopString[EPOCH_STRING_LEN+1];
                    toEncodeEPOCH(startTime, 0, startString);
                    toEncodeEPOCH(v->timestamp[i], 0, stopString);
                    fprintf(state->processingLogFile, "%sInterval spanning %s to %s has a duration of less than %d seconds: not exporting %ld records.\n", infoHeader, startString, stopString, SECONDS_OF_DATA_REQUIRED_FOR_EXPORTING, stopIndex - startIndex + 1);
                }
            }

            // Try next interval
            startIndex = i;
            startTime = v->timestamp[i];
            if (v->timestamp[i] >= maxTime || i == v->nRecs)
            {
                exporting = false;
            }
        }
        stopIndex = i;
        stopTime = v->timestamp[i];
        i++;
        if (i == v->nRecs)
        {
            exporting = false;
        }
    }
    // report
    if (state->writeLogFiles) {
        fprintf(state->processingLogFile, "%sExported %.0f orbits (%ld 16 Hz records) of science data in %d files. %.1f%% coverage.\n", infoHeader, minutesExported/94., recordsExported, filesExported, minutesExported/1440.0*100.0);
    }

    return status;

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
    else if (strcmp(attr.type, "CDF_UINT1") == 0)
    {
        uint8_t val = (uint8_t) attr.validMin;
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "VALIDMIN"), varNum, CDF_UINT1, 1, &val);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            return status;
        }
        val = (uint8_t) attr.validMax;
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "VALIDMAX"), varNum, CDF_UINT1, 1, &val);
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
    else if (strcmp(attr.type, "CDF_DOUBLE") == 0)
    {
        double val = (double) attr.validMin;
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "VALIDMIN"), varNum, CDF_DOUBLE, 1, &val);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            return status;
        }
        val = (double) attr.validMax;
        status = CDFputAttrzEntry(id, CDFgetAttrNum(id, "VALIDMAX"), varNum, CDF_DOUBLE, 1, &val);
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
        {"U_SC", "CDF_FLOAT", "V", "Satellite floating potential estimate from EXTD LP_HM dataset.", -50., 5.},
        {"Quality_flags", "CDF_UINT2", "*", "Bitwise flag for each velocity component, where a value of 1 for a particular component signifies that calibration was successful, and that the baseline 1-sigma noise level is less than or equal to 100 m/s at 2 Hz. Electric field quality can be assessed from these flags according to -vxB. Bit0 (least significant) = Vixh, bit1 = Vixv, bit2 = Viy, bit3 = Viz. Refer to the release notes for details.", 0, 65535},
        {"Calibration_flags", "CDF_UINT4", "*", "Information about the calibration process. Refer to the release notes for details.", 0, 4294967295},
        {"GeoelectricPotential", "CDF_FLOAT", "V", "Geoelectric potential estimate from Ehx.", -1000000., 1000000.},
        {"GeoelectricPotentialDifference", "CDF_FLOAT", "V", "Geoelectric potential at end of orbit region minus geoelectric potential at beginning of orbit region.", -1000000., 1000000.},
        {"MaxAbsGeoelectricPotentialBaselineSlope", "CDF_FLOAT", "mV/m", "Maximum absolute slope of the mid-latitude geoelectric potential estimates.", -1000000., 1000000.},
        {"EhxAdjusted", "CDF_FLOAT", "mV/m", "Along-track electric field component adjusted according to Zhu et al. (2020).", -400., 400.},
        {"EhxAdjustmentParameter", "CDF_FLOAT", "*", "Parameter used to adjust the along-track electric field component to remove mid-latitude potential offset (Zhu et al., 2020).", -400., 400.},
        {"GeoelectricPotentialDetrended", "CDF_FLOAT", "V", "Geoelectric potential linearly detrended within each orbit region.", -1000000., 1000000.},
        {"MaxAbsGeoelectricPotentialDetrendedBaselineSlope", "CDF_FLOAT", "mV/m", "Maximum absolute slope of the mid-latitude detrended geoelectric potential estimates.", -1000000., 1000000.},
        {"OrbitRegion", "CDF_UINT1", "*", "Orbit region poleward or equatorward of plus-or-minus 44.0 degrees quasi-dipole latitude. 0: northern polar ascending; 1: equatorial descending; 2: southern polar descending; 3: equatorial ascending; 255: incomplete or invalid.", 0, 3},
    };

    for (uint8_t i = 0; i < NUM_EXPORT_VARIABLES; i++)
    {
        addVariableAttributes(id, variableAttrs[i]);
    }

    return;

}

int exportTCT16Cdfs(ProcessorState *state, double startTime, double stopTime, long startIndex, long stopIndex)
{
    ProcessorVariables_t *v = &state->vars16hz;

    if (state->writeLogFiles) {
        fprintf(state->processingLogFile, "%sExporting 16 Hz data.\n",infoHeader);
        char epochString[EPOCH_STRING_LEN+1];
        toEncodeEPOCH(startTime, 0, epochString);
        fprintf(state->processingLogFile, "%sStartTime: %s\n", infoHeader, epochString);
        toEncodeEPOCH(stopTime, 0, epochString);
        fprintf(state->processingLogFile, "%sStoptime: %s\n", infoHeader, epochString);
    }

    char cdfFileName[CDF_PATHNAME_LEN];
    constructExportFileName("TCT16", startTime, stopTime, state->args.exportDir, state->args.exportVersion, state->args.satellite, cdfFileName);

    char zipFileName[FILENAME_MAX];
    sprintf(zipFileName, "%s.ZIP", cdfFileName);
    if (access(zipFileName, F_OK) == 0)
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sTIICT ZIP file exists. Not exporting.\n", infoHeader);
        }
        return TIICT_ZIP_EXISTS;
    }

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
        createVarFrom1DVar(exportCdfId, "Timestamp", CDF_EPOCH, startIndex, stopIndex, v->timestamp);
        createVarFrom1DVar(exportCdfId, "Latitude", CDF_REAL4, startIndex, stopIndex, v->latitude);
        createVarFrom1DVar(exportCdfId, "Longitude", CDF_REAL4, startIndex, stopIndex, v->longitude);
        createVarFrom1DVar(exportCdfId, "Radius", CDF_REAL4, startIndex, stopIndex, v->radius);
        createVarFrom1DVar(exportCdfId, "QDLatitude", CDF_REAL4, startIndex, stopIndex, v->qdlat);
        createVarFrom1DVar(exportCdfId, "MLT", CDF_REAL4, startIndex, stopIndex, v->mlt);
        createVarFrom1DVar(exportCdfId, "Vixh", CDF_REAL4, startIndex, stopIndex, v->vixh);
        createVarFrom1DVar(exportCdfId, "Vixh_error", CDF_REAL4, startIndex, stopIndex, v->vixherror);
        createVarFrom1DVar(exportCdfId, "Vixv", CDF_REAL4, startIndex, stopIndex, v->vixv);
        createVarFrom1DVar(exportCdfId, "Vixv_error", CDF_REAL4, startIndex, stopIndex, v->vixverror);
        createVarFrom1DVar(exportCdfId, "Viy", CDF_REAL4, startIndex, stopIndex, v->viy);
        createVarFrom1DVar(exportCdfId, "Viy_error", CDF_REAL4, startIndex, stopIndex, v->viyerror);
        createVarFrom1DVar(exportCdfId, "Viz", CDF_REAL4, startIndex, stopIndex, v->viz);
        createVarFrom1DVar(exportCdfId, "Viz_error", CDF_REAL4, startIndex, stopIndex, v->vizerror);
        createVarFrom1DVar(exportCdfId, "VsatN", CDF_REAL4, startIndex, stopIndex, v->vsatn);
        createVarFrom1DVar(exportCdfId, "VsatE", CDF_REAL4, startIndex, stopIndex, v->vsate);
        createVarFrom1DVar(exportCdfId, "VsatC", CDF_REAL4, startIndex, stopIndex, v->vsatc);
        createVarFrom1DVar(exportCdfId, "Ehx", CDF_REAL4, startIndex, stopIndex, v->ectxh);
        createVarFrom1DVar(exportCdfId, "Ehy", CDF_REAL4, startIndex, stopIndex, v->ectyh);
        createVarFrom1DVar(exportCdfId, "Ehz", CDF_REAL4, startIndex, stopIndex, v->ectzh);
        createVarFrom1DVar(exportCdfId, "Evx", CDF_REAL4, startIndex, stopIndex, v->ectxv);
        createVarFrom1DVar(exportCdfId, "Evy", CDF_REAL4, startIndex, stopIndex, v->ectyv);
        createVarFrom1DVar(exportCdfId, "Evz", CDF_REAL4, startIndex, stopIndex, v->ectzv);
        createVarFrom1DVar(exportCdfId, "Bx", CDF_REAL4, startIndex, stopIndex, v->bctx);
        createVarFrom1DVar(exportCdfId, "By", CDF_REAL4, startIndex, stopIndex, v->bcty);
        createVarFrom1DVar(exportCdfId, "Bz", CDF_REAL4, startIndex, stopIndex, v->bctz);
        createVarFrom1DVar(exportCdfId, "Vicrx", CDF_REAL4, startIndex, stopIndex, v->vicrx);
        createVarFrom1DVar(exportCdfId, "Vicry", CDF_REAL4, startIndex, stopIndex, v->vicry);
        createVarFrom1DVar(exportCdfId, "Vicrz", CDF_REAL4, startIndex, stopIndex, v->vicrz);
        createVarFrom1DVar(exportCdfId, "U_SC", CDF_REAL4, startIndex, stopIndex, v->potentials);
        createVarFrom1DVar(exportCdfId, "Quality_flags", CDF_UINT2, startIndex, stopIndex, v->flags);
        createVarFrom1DVar(exportCdfId, "Calibration_flags", CDF_UINT4, startIndex, stopIndex, v->fitInfo);
        createVarFrom1DVar(exportCdfId, "GeoelectricPotential", CDF_REAL4, startIndex, stopIndex, v->geoelectricPotential);
        createVarFrom1DVar(exportCdfId, "GeoelectricPotentialDifference", CDF_REAL4, startIndex, stopIndex, v->geoelectricPotentialDifference);
        createVarFrom1DVar(exportCdfId, "MaxAbsGeoelectricPotentialBaselineSlope", CDF_REAL4, startIndex, stopIndex, v->maxAbsGeoelectricPotentialBaselineSlope);
        createVarFrom1DVar(exportCdfId, "EhxAdjusted", CDF_REAL4, startIndex, stopIndex, v->ehxAdjusted);
        createVarFrom1DVar(exportCdfId, "EhxAdjustmentParameter", CDF_REAL4, startIndex, stopIndex, v->ehxAdjustmentParameter);
        createVarFrom1DVar(exportCdfId, "GeoelectricPotentialDetrended", CDF_REAL4, startIndex, stopIndex, v->geoelectricPotentialDetrended);
        createVarFrom1DVar(exportCdfId, "MaxAbsGeoelectricPotentialDetrendedBaselineSlope", CDF_REAL4, startIndex, stopIndex, v->maxAbsGeoelectricPotentialDetrendedBaselineSlope);
        createVarFrom1DVar(exportCdfId, "OrbitRegion", CDF_UINT1, startIndex, stopIndex, v->orbitRegion);

        // add attributes
        addAttributes(exportCdfId, "TCT16", state->args.satellite, state->args.exportVersion, startTime, stopTime);

        // Close export file
        closeCdf(exportCdfId);
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sExported %ld records to %s.cdf\n", infoHeader, (stopIndex - startIndex + 1), cdfFileName);
            fflush(state->processingLogFile);
        }

        if (state->exportZip) {
            status = zipCdfFile(state, cdfFileName);
        }
    }

    return status;

}

int exportTCT02Cdfs(ProcessorState *state, double startTime, double stopTime, long startIndex, long stopIndex)
{
    if (state->writeLogFiles) {
        fprintf(state->processingLogFile, "%sExporting 2 Hz data.\n",infoHeader);
    }

    ProcessorVariables_t *v = &state->vars16hz;

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
    if (state->writeLogFiles) {
        fprintf(state->processingLogFile, "%sDown-sampling 16 Hz to 2 Hz.\n",infoHeader);
    }

    long n2HzSamples = 0;
    long storageIndex = startIndex;
    double t0;
    bool downSampled = false;

    // request enough memory for 2 Hz data
    state->vars2hz.nRecs = state->vars16hz.nRecs / 8 + 1;
    // TODO set to actual number of records
    int reallocstatus = reallocVariables(&state->vars2hz, state->vars2hz.nRecs);
    if (reallocstatus != TIICT_OK) {
        fprintf(state->processingLogFile, "%sUnable to allocate memory for 2Hz data.\n",infoHeader);
        return reallocstatus;
    }

    for (long i = startIndex; i <= stopIndex;) // Time index advanced below
    {
        t0 = floor(v->timestamp[i] / 1000.0); // UT second reference
        for (uint8_t halfSecond = 0; halfSecond < 2; halfSecond ++)
        {
            downSampled = downSampleHalfSecond(state, &i, storageIndex, t0 + 0.5 * halfSecond, stopIndex);
            if (downSampled)
            {
                storageIndex++;
                n2HzSamples++;
            }
        }

    }
    // Update start and stop indexes and times
    v = &state->vars2hz;
    stopIndex = startIndex + n2HzSamples - 1;
    startTime = v->timestamp[startIndex];
    stopTime = v->timestamp[stopIndex];
    if (state->writeLogFiles) {
        char epochString[EPOCH_STRING_LEN+1];
        toEncodeEPOCH(startTime, 0, epochString);
        fprintf(state->processingLogFile, "%sStartTime: %s\n", infoHeader, epochString);
        toEncodeEPOCH(stopTime, 0, epochString);
        fprintf(state->processingLogFile, "%sStoptime: %s\n", infoHeader, epochString);
    }

    CDFid exportCdfId;
    CDFstatus status;
    char cdfFileName[CDF_PATHNAME_LEN];
    constructExportFileName("TCT02", startTime, stopTime, state->args.exportDir, state->args.exportVersion, state->args.satellite, cdfFileName);

    char zipFileName[FILENAME_MAX];
    sprintf(zipFileName, "%s.ZIP", cdfFileName);
    if (access(zipFileName, F_OK) == 0)
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sTIICT ZIP file exists. Not exporting.\n", infoHeader);
        }
        return TIICT_ZIP_EXISTS;
    }

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
        createVarFrom1DVar(exportCdfId, "Timestamp", CDF_EPOCH, startIndex, stopIndex, v->timestamp);
        createVarFrom1DVar(exportCdfId, "Latitude", CDF_REAL4, startIndex, stopIndex, v->latitude);
        createVarFrom1DVar(exportCdfId, "Longitude", CDF_REAL4, startIndex, stopIndex, v->longitude);
        createVarFrom1DVar(exportCdfId, "Radius", CDF_REAL4, startIndex, stopIndex, v->radius);
        createVarFrom1DVar(exportCdfId, "QDLatitude", CDF_REAL4, startIndex, stopIndex, v->qdlat);
        createVarFrom1DVar(exportCdfId, "MLT", CDF_REAL4, startIndex, stopIndex, v->mlt);
        createVarFrom1DVar(exportCdfId, "Vixh", CDF_REAL4, startIndex, stopIndex, v->vixh);
        createVarFrom1DVar(exportCdfId, "Vixh_error", CDF_REAL4, startIndex, stopIndex, v->vixherror);
        createVarFrom1DVar(exportCdfId, "Vixv", CDF_REAL4, startIndex, stopIndex, v->vixv);
        createVarFrom1DVar(exportCdfId, "Vixv_error", CDF_REAL4, startIndex, stopIndex, v->vixverror);
        createVarFrom1DVar(exportCdfId, "Viy", CDF_REAL4, startIndex, stopIndex, v->viy);
        createVarFrom1DVar(exportCdfId, "Viy_error", CDF_REAL4, startIndex, stopIndex, v->viyerror);
        createVarFrom1DVar(exportCdfId, "Viz", CDF_REAL4, startIndex, stopIndex, v->viz);
        createVarFrom1DVar(exportCdfId, "Viz_error", CDF_REAL4, startIndex, stopIndex, v->vizerror);
        createVarFrom1DVar(exportCdfId, "VsatN", CDF_REAL4, startIndex, stopIndex, v->vsatn);
        createVarFrom1DVar(exportCdfId, "VsatE", CDF_REAL4, startIndex, stopIndex, v->vsate);
        createVarFrom1DVar(exportCdfId, "VsatC", CDF_REAL4, startIndex, stopIndex, v->vsatc);
        createVarFrom1DVar(exportCdfId, "Ehx", CDF_REAL4, startIndex, stopIndex, v->ectxh);
        createVarFrom1DVar(exportCdfId, "Ehy", CDF_REAL4, startIndex, stopIndex, v->ectyh);
        createVarFrom1DVar(exportCdfId, "Ehz", CDF_REAL4, startIndex, stopIndex, v->ectzh);
        createVarFrom1DVar(exportCdfId, "Evx", CDF_REAL4, startIndex, stopIndex, v->ectxv);
        createVarFrom1DVar(exportCdfId, "Evy", CDF_REAL4, startIndex, stopIndex, v->ectyv);
        createVarFrom1DVar(exportCdfId, "Evz", CDF_REAL4, startIndex, stopIndex, v->ectzv);
        createVarFrom1DVar(exportCdfId, "Bx", CDF_REAL4, startIndex, stopIndex, v->bctx);
        createVarFrom1DVar(exportCdfId, "By", CDF_REAL4, startIndex, stopIndex, v->bcty);
        createVarFrom1DVar(exportCdfId, "Bz", CDF_REAL4, startIndex, stopIndex, v->bctz);
        createVarFrom1DVar(exportCdfId, "Vicrx", CDF_REAL4, startIndex, stopIndex, v->vicrx);
        createVarFrom1DVar(exportCdfId, "Vicry", CDF_REAL4, startIndex, stopIndex, v->vicry);
        createVarFrom1DVar(exportCdfId, "Vicrz", CDF_REAL4, startIndex, stopIndex, v->vicrz);
        createVarFrom1DVar(exportCdfId, "U_SC", CDF_REAL4, startIndex, stopIndex, v->potentials);
        createVarFrom1DVar(exportCdfId, "Quality_flags", CDF_UINT2, startIndex, stopIndex, v->flags);
        createVarFrom1DVar(exportCdfId, "Calibration_flags", CDF_UINT4, startIndex, stopIndex, v->fitInfo);
        createVarFrom1DVar(exportCdfId, "GeoelectricPotential", CDF_REAL4, startIndex, stopIndex, v->geoelectricPotential);
        createVarFrom1DVar(exportCdfId, "GeoelectricPotentialDifference", CDF_REAL4, startIndex, stopIndex, v->geoelectricPotentialDifference);
        createVarFrom1DVar(exportCdfId, "MaxAbsGeoelectricPotentialBaselineSlope", CDF_REAL4, startIndex, stopIndex, v->maxAbsGeoelectricPotentialBaselineSlope);
        createVarFrom1DVar(exportCdfId, "EhxAdjusted", CDF_REAL4, startIndex, stopIndex, v->ehxAdjusted);
        createVarFrom1DVar(exportCdfId, "EhxAdjustmentParameter", CDF_REAL4, startIndex, stopIndex, v->ehxAdjustmentParameter);
        createVarFrom1DVar(exportCdfId, "GeoelectricPotentialDetrended", CDF_REAL4, startIndex, stopIndex, v->geoelectricPotentialDetrended);
        createVarFrom1DVar(exportCdfId, "MaxAbsGeoelectricPotentialDetrendedBaselineSlope", CDF_REAL4, startIndex, stopIndex, v->maxAbsGeoelectricPotentialDetrendedBaselineSlope);
        createVarFrom1DVar(exportCdfId, "OrbitRegion", CDF_UINT1, startIndex, stopIndex, v->orbitRegion);

        // add attributes
        // update start and stop times to the averaged ones
        addAttributes(exportCdfId, "TCT02", state->args.satellite, state->args.exportVersion, startTime, stopTime);

        // close export file
        closeCdf(exportCdfId);
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sExported %ld records to %s.cdf\n", infoHeader, (stopIndex - startIndex + 1), cdfFileName);
            fflush(state->processingLogFile);
        }

        if (state->exportZip) {
            status = zipCdfFile(state, cdfFileName);
        }

    }

    return status;
}


int zipCdfFile(ProcessorState *state, char *cdfFileName)
{

    // Archive the CDF and HDR files in a ZIP file
    int status = system(NULL);
    if (status == 0)
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sSystem shell call not available. Not archiving CDF.\n", infoHeader);
        }
        return TIICT_SHELL;
    }
    status = system("zip -q 1 > /dev/null");
    if (WIFEXITED(status) && WEXITSTATUS(status) == 12)
    {
        char command[6*FILENAME_MAX + 100];
        sprintf(command, "zip -Z store -q -r -j %s.ZIP %s.cdf && rm %s.cdf", cdfFileName, cdfFileName, cdfFileName);
        status = system(command);
        if (WIFEXITED(status) && (WEXITSTATUS(status) == 0))
        {
            if (state->writeLogFiles) {
                fprintf(state->processingLogFile, "%sStored CDF file in %s.ZIP\n", infoHeader, cdfFileName);
            }
            status = TIICT_OK;
        }
        else
        {
            if (state->writeLogFiles) {
                fprintf(state->processingLogFile, "%sFailed to archive CDF file.\n", infoHeader);
            }
            status = TIICT_ZIP;
        }
    }
    else
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "zip is unusable. Not archiving CDF.\n");
        }
        status = TIICT_ZIP;
    }

    return status;
}

int initDirectories(ProcessorState *state)
{
    Arguments *args = &state->args;
    // Create output directory structure if it does not exist
    int dirStat = makeSureDirExists(args->exportDir, args->exportVersion, "TCT16");
    if (dirStat != 0)
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sTCT16 export directory is unavailable (or could not be created).\n", infoHeader);
        }
        return TIICT_EXPORT_DIRECTORY_TCT16;
    }
    dirStat = makeSureDirExists(args->exportDir, args->exportVersion, "TCT02");
    if (dirStat != 0)
    {
        if (state->writeLogFiles) {
            fprintf(state->processingLogFile, "%sTCT02 export directory is unavailable (or could not be created).\n", infoHeader);
        }
        return TIICT_EXPORT_DIRECTORY_TCT02;
    }

    return TIICT_OK;
}

int initLogFiles(ProcessorState *state)
{
    Arguments *a = &state->args;
    int status = makeSureDirExists(a->exportDir, a->exportVersion, "logs");
    if (status != TIICT_OK) {
        return status;
    }
    sprintf(state->fitLogFilename, "%s/%s/logs/%s%04d%02d%02d.fit", a->exportDir, a->exportVersion, a->satellite, a->year, a->month, a->day);
    sprintf(state->processingLogFilename, "%s/%s/logs/%s%04d%02d%02d.log", a->exportDir, a->exportVersion, a->satellite, a->year, a->month, a->day);

    state->fitFile = fopen(state->fitLogFilename, "a");
    if (state->fitFile == NULL) {
        return TIICT_LOG_WRITE;
    }
    fprintf(state->fitFile, "EFI TII CrossTrackCalibration fit results by fit region.\n");
    fprintf(state->fitFile, "Each region consists of two mid-latitude segments denoted by CDF_EPOCH times T11, T12, T21, and T22.\n");
    fprintf(state->fitFile, "Linear models based on robust least squares (GNU Scientific Library) are subtracted from each region for which a fit can be obtained.\n");
    fprintf(state->fitFile, "Regions:\n");
    for (uint8_t ind = 0; ind < 4; ind++)
    {
        fprintf(state->fitFile, "%d %21s: (% 5.1f, % 5.1f) -> (% 5.1f, % 5.1f)\n", state->fitargs[ind].regionNumber, state->fitargs[ind].regionName, state->fitargs[ind].lat1, state->fitargs[ind].lat2, state->fitargs[ind].lat3, state->fitargs[ind].lat4);
    }
    fprintf(state->fitFile, "\n");
    fprintf(state->fitFile, "The columns are:\n");
    fprintf(state->fitFile, "regionNumber fitNumber numPoints1 numPoints2 T11 T12 T21 T22 offsetHX slopeHX adjRsqHX rmseHX medianHX1 medianHX2 madHX madHX1 madHX2 offsetHY slopeHY adjRsqHY rmseHY medianHY1 medianHY2 madHY madHY1 madHY2 offsetVX slopeVX adjRsqVX rmseVX medianVX1 medianVX2 madVX madVX1 madVX2 offsetVY slopeVY adjRsqVY rmseVY medianVY1 medianVY2 madVY madVY1 madVY2\n");
    fprintf(state->fitFile, "\n");
    fflush(state->fitFile);

    if (state->writeLogFiles) {
        state->processingLogFile = fopen(state->processingLogFilename, "a");
    }
    else {
        state->processingLogFile = stdout;
    }
    if (state->processingLogFile == NULL)
        return TIICT_LOG_WRITE;
    fflush(state->processingLogFile);

    return TIICT_OK;
}

