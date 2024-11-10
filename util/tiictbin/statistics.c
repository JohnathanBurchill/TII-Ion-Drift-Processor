/*

    spaceseries: statistics.c

    Copyright (C) 2023  Johnathan K Burchill

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

// Adapted from TII-Ion-Drift-Processor
// https://github.com/JohnathanBurchill/TII-Ion-Drift-Processor/tree/tracis_flagging

#include "statistics.h"

#include <string.h>
#include <strings.h>
#include <stdbool.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics.h>

static char *availableStatistics[NSTATISTICS] = {
    "Mean",
    "Median",
    "StandardDeviation",
    "MedianAbsoluteDeviation",
    "Min",
    "Max",
    "Count"
};

// Equal-area binnary functions
static inline double solidAngle(double lat1, double lat2, double mlt1, double mlt2)
{
    // lat2 is more northward than lat1
    return (mlt2 - mlt1) * (cos((90.0 - lat2) * M_PI / 180.0) - cos((90.0 - lat1) * M_PI / 180.0));
}
static inline int nRingBins(double ringSolidAngle, double solidAngleUnit)
{
    return (int) fabs(round(ringSolidAngle / solidAngleUnit));
} 

// Preliminary memory management for binning
int initBinningState(BinningState *binningState)
{
    binningState->nQDLats = (int)floor((binningState->qdlatmax - binningState->qdlatmin) / binningState->deltaqdlat);
    binningState->nMLTs = (int)floor((binningState->mltmax - binningState->mltmin) / binningState->deltamlt);
    if (binningState->nQDLats <= 0)
    {
        fprintf(stderr, "Invalid QDLat bin dimension specification.\n");
        return BIN_SPECIFICATION;
    }
    if (binningState->nMLTs <= 0 && !binningState->equalArea)
    {
        fprintf(stderr, "Invalid MLT bin dimension specification.\n");
        return BIN_SPECIFICATION;
    }

    // Allocate memory for binning
    binningState->nMltsVsLatitude = calloc(binningState->nQDLats, sizeof *binningState->nMltsVsLatitude);
    binningState->cumulativeMltsVsLatitude = calloc(binningState->nQDLats, sizeof *binningState->cumulativeMltsVsLatitude); 
    if (binningState->nMltsVsLatitude == NULL || binningState->cumulativeMltsVsLatitude == NULL)
    {
        fprintf(stderr, "Unable to allocate memory.\n");
        return BIN_MEMORY;
    }

    // Approximate bin area
    // Will adjust deltaMlt to get bin areas closest to this
    binningState->solidAngleUnit = solidAngle(90.0 - binningState->deltaqdlat, 90.0, binningState->mltmin, binningState->mltmax) / (double)binningState->nMLTs;
    binningState->ringSolidAngle = 0.0;

    double lat = 0.0;
    // q = 0 is middle of lowest latitude bin
    binningState->cumulativeMltsVsLatitude[0] = 0;
    for (int q = 0; q < binningState->nQDLats; q++)
    {
        if (!binningState->equalArea)
            binningState->nMltsVsLatitude[q] = binningState->nMLTs;
        else
        {
            binningState->ringSolidAngle = solidAngle(binningState->qdlatmin + q * binningState->deltaqdlat, binningState->qdlatmin + (q+1) * binningState->deltaqdlat, binningState->mltmin, binningState->mltmax);
            binningState->nMltsVsLatitude[q] = nRingBins(binningState->ringSolidAngle, binningState->solidAngleUnit);
        }        
        binningState->nBins += binningState->nMltsVsLatitude[q];
        if (q > 0)
            binningState->cumulativeMltsVsLatitude[q] = binningState->cumulativeMltsVsLatitude[q-1] + binningState->nMltsVsLatitude[q - 1];
    }

    if (allocateBinStorage(&binningState->binStorage, &binningState->binSizes, &binningState->binValidSizes, &binningState->binMaxSizes, binningState->nBins, BIN_STORAGE_BLOCK_SIZE, binningState->equalArea))
    {
        fprintf(stderr, "Could not allocate bin storage memory.\n");
        return BIN_MEMORY;
    }
    // Access bins with bins[mltIndex * nQDLats + qdlatIndex];

    return BIN_OK;
}

// Main bin storage
int allocateBinStorage(double ***bins, size_t **binSizes, size_t **binValidSizes, size_t **binMaxSizes, size_t nBins, size_t sizePerBin, bool useEqualArea)
{
    *bins = malloc(nBins * sizeof *bins);
    *binSizes = malloc(nBins * sizeof *binSizes);
    *binValidSizes = malloc(nBins * sizeof *binValidSizes);
    *binMaxSizes = malloc(nBins * sizeof *binMaxSizes);
    if (*bins == NULL || *binSizes == NULL || *binValidSizes == NULL || *binMaxSizes == NULL)
        return STATISTICS_MEM;
    for (size_t i = 0; i < nBins; i++)
    {
        (*bins)[i] = NULL;
        (*binSizes)[i] = 0;
        (*binValidSizes)[i] = 0;
        (*binMaxSizes)[i] = 0;
    }
    for (size_t i = 0; i < nBins; i++)
    {
        (*bins)[i] = malloc((sizePerBin * sizeof **bins));
        if ((*bins)[i] == NULL)
            return STATISTICS_MEM;
        bzero((*bins)[i], sizePerBin * sizeof **bins);
        (*binMaxSizes)[i] = sizePerBin;
    }

    return STATISTICS_OK;
}

// Grow the allocated bin storage
int adjustBinStorage(double **bins, size_t *binMaxSizes, int mltQdLatIndex, long numberOfElementsToAdd)
{
    void *mem = realloc(bins[mltQdLatIndex], (binMaxSizes[mltQdLatIndex] + numberOfElementsToAdd) * sizeof *bins[mltQdLatIndex]);
    if (mem == NULL)
        return STATISTICS_MEM;
    bins[mltQdLatIndex] = mem;

    binMaxSizes[mltQdLatIndex] += numberOfElementsToAdd;

    return STATISTICS_OK;
}

void freeBinStorage(BinningState *binningState)
{
    for (int i = 0; i < binningState->nBins; i++)
    {
        free(binningState->binStorage[i]);
    }

    free(binningState->binStorage);
    free(binningState->binSizes);
    free(binningState->binValidSizes);
    free(binningState->binMaxSizes);

    free(binningState->nMltsVsLatitude);
    free(binningState->cumulativeMltsVsLatitude);

    return;
}

// Perform binning
int binData(BinningState *binningState, double qdlat, double mlt, double value, bool includeValue)
{
    int qdlatIndex = 0;
    int mltIndex = 0;

    if (!isfinite(value))
        return BIN_VALUE_NOT_FINITE;

    size_t index = 0;

    qdlatIndex = (int) floor((qdlat - binningState->qdlatmin) / binningState->deltaqdlat);
    if (qdlatIndex < 0 || qdlatIndex >= binningState->nQDLats)
        return BIN_QDLAT_OUTOFRANGE;

    double deltamlt = (binningState->mltmax - binningState->mltmin) / (double)binningState->nMltsVsLatitude[qdlatIndex];
    mltIndex = (int) floor((mlt - binningState->mltmin) / deltamlt);
    if (mltIndex >= 0 && mltIndex < binningState->nMltsVsLatitude[qdlatIndex])
    {
        index = binningState->cumulativeMltsVsLatitude[qdlatIndex] + mltIndex;
        // Measurement lies within a QDLat and MLT bin
        binningState->binValidSizes[index]++;
        binningState->nValsWithinBinLimits++;

        if (includeValue)
        {
            // Flip sign of parameter when moving southward, i.e. to make Viy eastward and Vixh or Vixv northward

            if (binningState->binSizes[index] >= binningState->binMaxSizes[index])
            {
                if(adjustBinStorage(binningState->binStorage, binningState->binMaxSizes, index, BIN_STORAGE_BLOCK_SIZE))
                {
                    fprintf(stderr, "Unable to allocate additional bin storage.\n");
                    exit(EXIT_FAILURE);
                }

            }
            binningState->binStorage[index][binningState->binSizes[index]] = value;
            binningState->binSizes[index]++;
            binningState->nValsBinned++;
        }
    }

    return BIN_OK;
}


void printBinningResults(BinningState *binningState, char *parameter, char *statistic)
{
    double qdlat1 = 0.0;
    double qdlat2 = 0.0;
    double mlt1 = 0.0;
    double mlt2 = 0.0;
    double result = 0.0;

    fprintf(stdout, "Time range is inclusive. Bin specification for remaining quantities x and bin boundaries x1 and x2: x1 <= x < x2\n");
    fprintf(stdout, "Row legend:\n");
    fprintf(stdout, "MLT1 MLT2 QDLat1 QDLat2 %s(%s) binCount validRegionFraction totalReadFraction\n", statistic, parameter);

    double denomBinValidSizes = 0.0;
    double denomNValsRead = binningState->nValsRead > 0 ? (double) binningState->nValsRead : 1.0;

    size_t index = 0;

    for (size_t q = 0; q < binningState->nQDLats; q++)
    {
        for (size_t m = 0; m < binningState->nMltsVsLatitude[q]; m++)
        {
            index = binningState->cumulativeMltsVsLatitude[q] + m;
            qdlat1 = binningState->qdlatmin + binningState->deltaqdlat * ((double)q);
            qdlat2 = qdlat1 + binningState->deltaqdlat;
            binningState->deltamlt = (binningState->mltmax - binningState->mltmin) / (double)binningState->nMltsVsLatitude[q];
            mlt1 = binningState->mltmin + binningState->deltamlt * (double)m;
            mlt2 = mlt1 + binningState->deltamlt;
            if (calculateStatistic(statistic, binningState->binStorage, binningState->binSizes, index, (void*) &result))
                result = GSL_NAN;

            denomBinValidSizes = binningState->binValidSizes[index] > 0 ? (double) binningState->binValidSizes[index] : 1.0;

            fprintf(stdout, "%5.2lf %5.2lf %6.2lf %6.2lf %lf %ld %lf %lf\n", mlt1, mlt2, qdlat1, qdlat2, result, binningState->binSizes[index], (double)binningState->binSizes[index] / denomBinValidSizes, (double)binningState->binSizes[index] / denomNValsRead);
        }
    }

    fprintf(stdout, "Summary of counts\n");
    fprintf(stdout, "\tValues read: %ld; Values within bin limits: %ld; Values binned: %ld (%6.2lf%% of those within bin limits)\n", binningState->nValsRead, binningState->nValsWithinBinLimits, binningState->nValsBinned, 100.0 * (double)binningState->nValsBinned / (double)binningState->nValsWithinBinLimits);

}

void printAvailableStatistics(FILE *dest)
{
    for (int i = 0; i < NSTATISTICS; i++)
        fprintf(dest, "\t%s\n", availableStatistics[i]);

    return;
}

bool validStatistic(const char *statistic)
{
    bool valid = false;
    for (int i = 0; i < NSTATISTICS; i++)
    {
        if (strcmp(statistic, availableStatistics[i])==0)
        {
            valid = true;
            break;
        }

    }

    return valid;
}

// Calculate requested statistic for each bin
int calculateStatistic(const char *statistic, double **bins, size_t *binSizes, size_t mltQdIndex, void *returnValue)
{
    int status = STATISTICS_OK;
    double result = 0.0;
    if (binSizes[mltQdIndex] == 0)
        return STATISTICS_NO_DATA;
    if (returnValue == NULL)
        return STATISTICS_POINTER;

    if (strcmp(statistic, "Mean")==0)
    {
        *(double*)returnValue = gsl_stats_mean(bins[mltQdIndex], 1, binSizes[mltQdIndex]);
    }
    else if (strcmp(statistic, "Median")==0)
    {
        *(double*)returnValue = gsl_stats_median(bins[mltQdIndex], 1, binSizes[mltQdIndex]);
    }
    else if (strcmp(statistic, "StandardDeviation")==0)
    {
        *(double*)returnValue = gsl_stats_sd(bins[mltQdIndex], 1, binSizes[mltQdIndex]);
    }
    else if (strcmp(statistic, "MedianAbsoluteDeviation")==0)
    {
        double *data = (double*)malloc(binSizes[mltQdIndex] * sizeof *data);
        if (data == NULL)
        {
            status = STATISTICS_MEM;
        }
        else
        {
            *(double*)returnValue = gsl_stats_mad(bins[mltQdIndex], 1, binSizes[mltQdIndex], data);
            free(data);
        }
    }
    else if (strcmp(statistic, "Min")==0)
    {
        *(double*)returnValue = gsl_stats_min(bins[mltQdIndex], 1, binSizes[mltQdIndex]);
    }
    else if (strcmp(statistic, "Max")==0)
    {
        *(double*)returnValue = gsl_stats_max(bins[mltQdIndex], 1, binSizes[mltQdIndex]);
    }
    else if (strcmp(statistic, "Count")==0)
    {
        *(double*)returnValue = (double) binSizes[mltQdIndex];
    }
    else
    {
        status = STATISTICS_UNSUPPORTED_STATISTIC;
    }

    return status;
}

