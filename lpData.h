#ifndef _LPDATA_H
#define _LPDATA_H

#include <stdint.h>

int getLpData(const char *lpDir, const char *satellite, const int year, const int month, const int day, uint8_t **dataBuffers, double **lpPhiScHighGain, double **lpPhiScLowGain);

int getLpInputFilename(const char satelliteLetter, long year, long month, long day, const char *path, char *filename);

void loadLpInputs(const char *cdfFile, double **lpTime, double **lpPhiSc1, double **lpPhiSc2, long *numberOfRecords);

#endif // _LPDATA_H
