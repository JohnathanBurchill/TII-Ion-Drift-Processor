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
#define VIY() (MEAS(1, 0))
#define VIYERROR() (MEAS(2, 0))
#define QDLAT() (MEAS(3, 0))
#define FLAG() ((uint16_t)*((uint16_t*)dataBuffers[4]+(timeIndex)))
#define CALFLAG() ((uint32_t)*((uint32_t*)dataBuffers[5]+(timeIndex)))

void analyzeMAD(const char satellite, uint8_t majorVersion, uint8_t minorVersion, bool fourByteCalFlag, uint8_t * dataBuffers[], long nRecs);

void analyzeTimeGaps(const char *filename, uint8_t * dataBuffers[], long nRecs, long *totalCount, long *gapCount);
void summarizeTimeGaps(const char *filename, uint8_t * dataBuffers[], long nRecs, long *totalCount, long *gapCount, long *largeGapCount);

void analyzeCalibrationFlag(const char satellite, uint8_t majorVersion, uint8_t minorVersion, bool fourByteCalFlag, uint8_t * dataBuffers[], long nRecs, bool numericOutput);

uint8_t getMajorVersion(const char *filename);
uint8_t getMinorVersion(const char *filename);

void getFileDate(const char *filename, char *startDate);

typedef enum AnalysisType {
    NO_ANALYSIS = 0,
    MAD_ANALYSIS = 1,
    CALIBRATION_ANALYSIS_NUMERIC = 2,
    CALIBRATION_ANALYSIS_GRAPHIC = 3,
    CALIBRATION_ANALYSIS_TIMEGAP = 4,
    CALIBRATION_ANALYSIS_TIMEGAPSUMMARY = 5
} AnalysisType;

#endif // CROSSTRACKVALIDATION
