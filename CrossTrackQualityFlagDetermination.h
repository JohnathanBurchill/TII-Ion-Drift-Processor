#ifndef CROSSTRACKVALIDATION
#define CROSSTRACKVALIDATION

#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include <cdf.h>

#define SOFTWARE_VERSION_STRING "CrossTrackValidation 2020-05-11"

#define NUM_DATA_VARIABLES 6

void closeCdf(CDFid id);

void loadCrossTrackData(const char * filename, uint8_t **dataBuffers, long *numberOfRecords, bool *fourByteCalFlag);


#define ADDR(n, m) (((float*)dataBuffers[(n)]+(timeIndex + m)))
#define MEAS(n, m) ((float)(*(ADDR(n, m))))
#define TIME() ((double)*((double*)dataBuffers[0]+(timeIndex)))
#define QDLAT() (MEAS(1, 0))
#define VIY() (MEAS(2, 0))
#define VIYERROR() (MEAS(3, 0))
#define FLAG() ((uint16_t)*((uint16_t*)dataBuffers[4]+(timeIndex)))
#define CALFLAG() ((uint32_t)*((uint32_t*)dataBuffers[5]+(timeIndex)))


uint8_t getMinorVersion(const char *filename);


#endif // CROSSTRACKVALIDATION
