#include "sample.h"
#include "indexing.h"
#include "settings.h"

#include <math.h>
#include <stdio.h>

void zeroOrderInterpolateU8(double *times, uint8_t *values, size_t nVals, double *requestedTimes, long nRequestedValues, uint8_t *newValues)
{
    
    size_t lastIndex = 0;
    double thisTime = 0;

    size_t interpIndex = 0;

    for (size_t i = 0; i < nRequestedValues; i++)
    {
        thisTime = requestedTimes[i];

        while (lastIndex < nVals && times[lastIndex] <= thisTime)
            lastIndex++;
        
        // Copy flag
        if (lastIndex > 0)
            interpIndex = lastIndex - 1;
        else
            interpIndex = 0;

        newValues[i] = values[interpIndex];

    }

    return;

}

void zeroOrderInterpolateU32(double *times, uint32_t *values, size_t nVals, double *requestedTimes, long nRequestedValues, uint32_t *newValues)
{
    
    size_t lastIndex = 0;
    double thisTime = 0;

    size_t interpIndex = 0;

    for (size_t i = 0; i < nRequestedValues; i++)
    {
        thisTime = requestedTimes[i];

        while (lastIndex < nVals && times[lastIndex] <= thisTime)
            lastIndex++;
        
        // Copy flag
        if (lastIndex > 0)
            interpIndex = lastIndex - 1;
        else
            interpIndex = 0;

        newValues[i] = values[interpIndex];

    }

    return;

}

// Copied and modified from TRACIS interpolate.c
void interpolate(double *times, double *values, size_t nVals, double *requestedTimes, long nRequestedValues, float *newValues)
{
    
    size_t lastIndex = 0;
    double thisTime = 0;
    double t1 = 0, t2 = 0, dt = 0;
    double fraction = 0;
    double v1 = 0, v2 = 0, dv = 0;
    double x = 0, y = 0, z = 0;
    double rad = 0;
    double lat = 0;
    double lon = 0;
    for (size_t i = 0; i < nRequestedValues; i++)
    {
        thisTime = requestedTimes[i];
        // 
        while (times[lastIndex] <= thisTime && lastIndex < nVals) 
        {
            lastIndex++;
        }
        // Extrapolate earlier or later, or the times are the same
        if (times[lastIndex] >= thisTime || lastIndex == nVals - 1)
        {
            newValues[i] = (float)values[lastIndex];
        }
        // Interpolate
        else
        {
            t1 = times[lastIndex];
            t2 = times[lastIndex+1];
            // Assumes t2 > t1
            dt = thisTime - t1;
            fraction = dt / (t2 - t1);
            newValues[i] = (float)(values[lastIndex] + (values[lastIndex+1] - values[lastIndex]) * fraction);
        }
    }

    return;

}

bool downSampleHalfSecond(long *index, long storageIndex, double t0, long maxIndex, uint8_t **dataBuffers, float *ectFieldH, float *ectFieldV, float *bctField, float *viErrors, float *potentials, uint16_t *flags, uint32_t *fitInfo, bool usePotentials)
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
        if (usePotentials)
            floatBuf[30] += potentials[timeIndex];
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
        if (usePotentials)
            potentials[storageIndex] = floatBuf[30] / 8.0; // Floating potential U_SC
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

