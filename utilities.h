#ifndef UTILITIES_H
#define UTILITIES_H

#include <cdf.h>

// Constructs a full path to the export CDF file in the argument constructExportFileName.
void constructExportFileName(const char *dataset, double startTime, double stopTime, const char *exportDir, const char *exportVersion, const char *satellite, char *cdfFileName);

// Prints an error message from the CDFstatus
void printErrorMessage(CDFstatus status);

#endif // UTILITIES_H
