/*

    TII Cross-Track Ion Drift Processor: statistics.c

    Copyright (C) 2022  Johnathan K Burchill

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

#include "statistics.h"

#include <string.h>
#include <strings.h>
#include <stdbool.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics_float.h>

static char *availableStatistics[NSTATISTICS] = {
    "Mean",
    "Median",
    "StandardDeviation",
    "MedianAbsoluteDeviation",
    "Min",
    "Max",
    "Count"
};

int allocateBinStorage(float ***bins, size_t **binSizes, size_t **binValidSizes, size_t **binMaxSizes, size_t nBins, size_t sizePerBin, bool useEqualArea)
{
    *bins = (float **) malloc(nBins * sizeof(float **));
    *binSizes = (size_t *) malloc(nBins * sizeof(size_t));
    *binValidSizes = (size_t *) malloc(nBins * sizeof(size_t));
    *binMaxSizes = (size_t *) malloc(nBins * sizeof(size_t));
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
        (*bins)[i] = (float*)malloc((size_t)(sizePerBin * sizeof(float)));
        if ((*bins)[i] == NULL)
            return STATISTICS_MEM;
        bzero((*bins)[i], sizePerBin * sizeof(float));
        (*binMaxSizes)[i] = sizePerBin;
    }

    return STATISTICS_OK;
}

int adjustBinStorage(float **bins, size_t *binMaxSizes, int mltQdLatIndex, long numberOfElementsToAdd)
{
    bins[mltQdLatIndex] = realloc(bins[mltQdLatIndex], (binMaxSizes[mltQdLatIndex] + numberOfElementsToAdd) * sizeof(float));
    if (bins[mltQdLatIndex] == NULL)
        return STATISTICS_MEM;

    binMaxSizes[mltQdLatIndex] += numberOfElementsToAdd;

    return STATISTICS_OK;
}

void freeBinStorage(float **bins, size_t *binSizes, size_t *binValidSizes, size_t *binMaxSizes, size_t nBins)
{
    for (int i = 0; i < nBins; i++)
    {
        free(bins[i]);
    }

    free(bins);
    free(binSizes);
    free(binValidSizes);
    free(binMaxSizes);

    return;
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

int calculateStatistic(const char *statistic, float **bins, size_t *binSizes, size_t mltQdIndex, void *returnValue)
{
    int status = STATISTICS_OK;
    float result = 0.0;
    if (binSizes[mltQdIndex] == 0)
        return STATISTICS_NO_DATA;
    if (returnValue == NULL)
        return STATISTICS_POINTER;

    if (strcmp(statistic, "Mean")==0)
    {
        *(float*)returnValue = gsl_stats_float_mean(bins[mltQdIndex], 1, binSizes[mltQdIndex]);
    }
    else if (strcmp(statistic, "Median")==0)
    {
        *(float*)returnValue = gsl_stats_float_median(bins[mltQdIndex], 1, binSizes[mltQdIndex]);
    }
    else if (strcmp(statistic, "StandardDeviation")==0)
    {
        *(float*)returnValue = gsl_stats_float_sd(bins[mltQdIndex], 1, binSizes[mltQdIndex]);
    }
    else if (strcmp(statistic, "MedianAbsoluteDeviation")==0)
    {
        double *data = (double*)malloc(binSizes[mltQdIndex] * sizeof(double));
        if (data == NULL)
        {
            status = STATISTICS_MEM;
        }
        else
        {
            *(float*)returnValue = gsl_stats_float_mad(bins[mltQdIndex], 1, binSizes[mltQdIndex], data);
            free(data);
        }
    }
    else if (strcmp(statistic, "Min")==0)
    {
        *(float*)returnValue = gsl_stats_float_min(bins[mltQdIndex], 1, binSizes[mltQdIndex]);
    }
    else if (strcmp(statistic, "Max")==0)
    {
        *(float*)returnValue = gsl_stats_float_max(bins[mltQdIndex], 1, binSizes[mltQdIndex]);
    }
    else if (strcmp(statistic, "Count")==0)
    {
        *(float*)returnValue = (float) binSizes[mltQdIndex];
    }
    else
    {
        status = STATISTICS_UNSUPPORTED_STATISTIC;
    }

    return status;
}
