#ifndef _TIICT_SAMPLE_H
#define _TIICT_SAMPLE_H

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

void zeroOrderInterpolateImageFlags(double *times, uint8_t *values, size_t nVals, double *requestedTimes, long nRequestedValues, uint8_t *newValues);

void interpolate(double *times, double *values, size_t nVals, double *requestedTimes, long nRequestedValues, float *newValues);

bool downSampleHalfSecond(long *index, long storageIndex, double t0, long maxIndex, uint8_t **dataBuffers, float *ectFieldH, float *ectFieldV, float *bctField, float *viErrors, float *potentials, uint16_t *flags, uint32_t *fitInfo, bool usePotentials);


#endif // _TIICT_SAMPLE_H
