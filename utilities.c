#include <stdio.h>
#include <cdf.h>

#include "utilities.h"

// Prefix for all fprintf messages
extern char infoHeader[50];


// Generates the filename for exported CDF file, with full path
void constructExportFileName(const char *dataset, double startTime, double stopTime, const char *exportDir, const char *exportVersion, const char *satellite, char *cdfFileName) 
{
 
    long startyear, startmonth, startday, starthour, startminute, startsecond, startmillisecond;
    long stopyear, stopmonth, stopday, stophour, stopminute, stopsecond, stopmillisecond;

    EPOCHbreakdown(startTime, &startyear, &startmonth, &startday, &starthour, &startminute, &startsecond, &startmillisecond);
    EPOCHbreakdown(stopTime, &stopyear, &stopmonth, &stopday, &stophour, &stopminute, &stopsecond, &stopmillisecond);

    sprintf(cdfFileName, "%s/%s/%s/SW_EXPT_EFI%s_%s_%04ld%02ld%02ldT%02ld%02ld%02ld_%04ld%02ld%02ldT%02ld%02ld%02ld_%s", exportDir, exportVersion, dataset, satellite, dataset, startyear, startmonth, startday, starthour, startminute, startsecond, stopyear, stopmonth, stopday, stophour, stopminute, stopsecond, exportVersion);

}

void printErrorMessage(CDFstatus status)
{
    char errorMessage[CDF_STATUSTEXT_LEN + 1];
    CDFgetStatusText(status, errorMessage);
    fprintf(stdout, "%s%s\n", infoHeader, errorMessage);
}
