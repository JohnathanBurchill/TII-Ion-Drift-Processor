#ifndef _LPDATA_H
#define _LPDATA_H

#include <stdint.h>
#include <stddef.h>

int getLpData(const char *lpDir, const char *satellite, const int year, const int month, const int day, uint8_t **dataBuffers, long nRecs, float **lpPhiScHighGain, float **lpPhiScLowGain, float **lpPhiSc);

int getLpInputFilename(const char satelliteLetter, long year, long month, long day, const char *path, char *filename);

void loadLpInputs(const char *cdfFile, double **lpTime, double **lpPhiScHighGain, double **lpPhiScLowGain, double **lpPhiSc, long *numberOfRecords);

void interpolate(double *times, double *values, size_t nVals, double *requestedTimes, long nRequestedValues, float *newValues);

#endif // _LPDATA_H
