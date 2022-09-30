/*

    TII Cross-Track Ion Drift Processor: processing.c

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

#include "processing.h"

#include "settings.h"
#include "indexing.h"

#include <stdio.h>
#include <stdbool.h>
#include <cdf.h>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>

extern char infoHeader[50];

void removeOffsetsAndSetFlags(const char* satellite, offset_model_fit_arguments fitargs, long nRecs, uint8_t **dataBuffers, float *viErrors, uint16_t *flags, uint32_t *fitInfo, FILE* fitFile, bool setFlags)
{
    // Robust linear least squares from GSL: https://www.gnu.org/software/gsl/doc/html/lls.html#examples
    long long timeIndex = 0;
    double epoch0 = TIME();
    uint16_t minDataPointsNeeded = 80;
    float previousQDLat = QDLAT();
    float driftValue;
    bool regionBegin = false;
    long long beginIndex0 = 0, beginIndex1 = 0, endIndex0 = 0, endIndex1 = 0, modelDataIndex = 0, modelDataMidPoint = 0;

    double c0, c1, cov00, cov01, cov11, sumsq;
    double seconds, fitTime;
    int gslStatus;
    long long numModelPoints, numModel1Points, numModel2Points; // data for fit, and for median calculations at each end
    bool gotFirstModelData = false;
    bool gotStartOfSecondModelData = false;
    bool gotSecondModelData = false;
    uint16_t numFits = 0;

    char startString[EPOCH_STRING_LEN+1], stopString[EPOCH_STRING_LEN+1];

    float dlatfirst = fitargs.lat2 - fitargs.lat1;
    float dlatsecond = fitargs.lat4 - fitargs.lat3;

    int8_t firstDirection = 0;
    if (dlatfirst > 0)
        firstDirection = 1;
    else if (dlatfirst < 0)
        firstDirection = -1;

    int8_t secondDirection = 0;
    if (dlatsecond > 0)
        secondDirection = 1;
    else if (dlatsecond < 0)
        secondDirection = -1;

    float location = (fitargs.lat2 + fitargs.lat3) / 2.;
    int8_t region; // For flagging
    if (location >= 44. )
        region = 1; // North
    else if (location <= -44.)
        region = -1; // South
    else
        region = 0; // Equator

    // Buffers for mid-latitude linear fit data
    const gsl_multifit_robust_type * fitType = gsl_multifit_robust_bisquare;
    gsl_matrix *modelTimesMatrix, *cov;
    gsl_vector *modelValues, *model1Values, *model2Values, *fitCoefficients, *work1, *work2;
    gsl_multifit_robust_workspace * gslFitWorkspace;
    gsl_multifit_robust_stats stats;
    const size_t p = 2; // linear fit

    double tregion11, tregion12, tregion21, tregion22;

    for (timeIndex = 0; timeIndex < nRecs; timeIndex++)
    {
        c0 = 0.0;
        c1 = 0.0;
        seconds = (TIME() - epoch0) / 1000.;
        if ((firstDirection == 1 && QDLAT() >= fitargs.lat1 && previousQDLat < fitargs.lat1) || (firstDirection == -1 && QDLAT() <= fitargs.lat1 && previousQDLat > fitargs.lat1))
        {
            // Start a new region search
            regionBegin = true; // Found start of region to remove offset from
            gotFirstModelData = false;
            gotStartOfSecondModelData = false;
            gotSecondModelData = false;
            beginIndex0 = timeIndex;
            tregion11 = TIME();

        }
        else if (regionBegin && ((firstDirection == 1 && QDLAT() >= fitargs.lat2 && previousQDLat < fitargs.lat2) || (firstDirection == -1 && QDLAT() <= fitargs.lat2 && previousQDLat > fitargs.lat2)))
        {
            if ((TIME() - tregion11)/1000. < (5400. / 2.)) // Should be within 1/2 an orbit of start of segment
            {
                gotFirstModelData = true;
                beginIndex1 = timeIndex;
                tregion12 = TIME();
            }
            else
            {
                // reset search
                gotFirstModelData = false;
                gotStartOfSecondModelData = false;
                gotSecondModelData = false;
                regionBegin = false;
            }
        }
        else if (gotFirstModelData && ((secondDirection == -1 && QDLAT() <= fitargs.lat3 && previousQDLat > fitargs.lat3) || (secondDirection == 1 && QDLAT() >= fitargs.lat3 && previousQDLat < fitargs.lat3)))
        {
            if ((TIME() - tregion12)/1000. < (5400. / 2.)) // Should be within 1/2 an orbit of start of segment
            {
                gotStartOfSecondModelData = true;
                endIndex0 = timeIndex;
                tregion21 = TIME();
            }
            else
            {
                // reset search
                gotFirstModelData = false;
                gotStartOfSecondModelData = false;
                gotSecondModelData = false;
                regionBegin = false;
            }
        }
        else if (gotStartOfSecondModelData && ((secondDirection == -1 && QDLAT() <= fitargs.lat4 && previousQDLat > fitargs.lat4) || (secondDirection == 1 && QDLAT() >= fitargs.lat4 && previousQDLat < fitargs.lat4)))
        {
            if ((TIME() - tregion21)/1000. < (5400. / 2.)) // Should be within 1/2 an orbit of start of segment
            {
                // We have a complete region - remove linear offset model
                endIndex1 = timeIndex;
                tregion22 = TIME();
                gotSecondModelData = true;
            }
            else
            {
                // reset search
                gotFirstModelData = false;
                gotStartOfSecondModelData = false;
                gotSecondModelData = false;
                regionBegin = false;
            }
            if (gotFirstModelData && gotSecondModelData)
            {
                numFits++;
                numModel1Points = beginIndex1 - beginIndex0;
                numModel2Points = endIndex1 - endIndex0;
                numModelPoints = numModel1Points + numModel2Points;
                fprintf(fitFile, "%d %d %lld %lld %f %f %f %f", fitargs.regionNumber, numFits, numModel1Points, numModel2Points, tregion11, tregion12, tregion21, tregion22);
                // Allocate fit buffers
                modelTimesMatrix = gsl_matrix_alloc(numModelPoints, p);
                model1Values = gsl_vector_alloc(numModel1Points);
                model2Values = gsl_vector_alloc(numModel2Points);
                work1 = gsl_vector_alloc(numModel1Points);
                work2 = gsl_vector_alloc(numModel2Points);
                modelValues = gsl_vector_alloc(numModelPoints);
                fitCoefficients = gsl_vector_alloc(p);
                cov = gsl_matrix_alloc(p, p);

                // Load times into the model data buffer
                modelDataIndex = 0;
                for (timeIndex = beginIndex0; timeIndex < beginIndex1; timeIndex++)
                {
                    fitTime = (TIME() - epoch0)/1000.;
                    gsl_matrix_set(modelTimesMatrix, modelDataIndex, 0, 1.0);
                    gsl_matrix_set(modelTimesMatrix, modelDataIndex++, 1, fitTime); // seconds from start of file
                }
                for (timeIndex = endIndex0; timeIndex < endIndex1; timeIndex++)
                {
                    fitTime = (TIME() - epoch0)/1000.;
                    gsl_matrix_set(modelTimesMatrix, modelDataIndex, 0, 1.0);
                    gsl_matrix_set(modelTimesMatrix, modelDataIndex++, 1, fitTime); // seconds from start of file
                }
                // Load values into model data buffer once each for HX, HY, VX, VY
                for (uint8_t k = 0; k < 4; k++)
                {
                    uint8_t flagIndex = k; // Defined flags as bit 0 -> HX, bit 1 -> VX, bit 2-> HY and bit 3-> VY. Data are stored in memory differently.
                    if (flagIndex == 1)
                        flagIndex = 2;
                    else if (flagIndex == 2)
                        flagIndex = 1;

                    modelDataIndex = 0;
                    for (timeIndex = beginIndex0; timeIndex < beginIndex1; timeIndex++)
                    {
                        gsl_vector_set(model1Values, modelDataIndex, *ADDR(1 + k / 2, k % 2, 2)); 
                        gsl_vector_set(modelValues, modelDataIndex++, *ADDR(1 + k / 2, k % 2, 2)); 
                    }
                    modelDataMidPoint = modelDataIndex;
                    for (timeIndex = endIndex0; timeIndex < endIndex1; timeIndex++)
                    {
                        gsl_vector_set(model2Values, modelDataIndex - modelDataMidPoint, *ADDR(1 + k / 2, k % 2, 2)); 
                        gsl_vector_set(modelValues, modelDataIndex++, *ADDR(1 + k / 2, k % 2, 2)); 
                    }
                    // Robust linear model fit and removal
                    gslFitWorkspace = gsl_multifit_robust_alloc(fitType, numModelPoints, p);
                    gslStatus = gsl_multifit_robust_maxiter(GSL_FIT_MAXIMUM_ITERATIONS, gslFitWorkspace);
                    if (gslStatus)
                    {
                        fprintf(stdout, "%sCould not set maximum GSL iterations.\n", infoHeader);
                    }
                    gslStatus = gsl_multifit_robust(modelTimesMatrix, modelValues, fitCoefficients, cov, gslFitWorkspace);
                    if (gslStatus)
                    {
                        toEncodeEPOCH(tregion11, 0, startString);
                        toEncodeEPOCH(tregion22, 0, stopString);
                        fprintf(stdout, "%s<GSL Fit Error: %s> for fit region from %s to %s spanning latitudes %.0f to %.0f.\n", infoHeader, gsl_strerror(gslStatus), startString, stopString, fitargs.lat1, fitargs.lat4);
                        // Print "-9999999999.GSLERRORNUMBER" for each of the nine fit parameters
                        fprintf(fitFile, " -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d", gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus);
                        if (setFlags)
                        {
                            for (timeIndex = beginIndex0; timeIndex < endIndex1; timeIndex++)
                            {
                                // Got a complete region, but had a fit error
                                fitInfo[timeIndex] &= ~(FITINFO_INCOMPLETE_REGION << (flagIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
                                fitInfo[timeIndex] |= (FITINFO_GSL_FIT_ERROR << (flagIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
                            }
                        }
                    }
                    else
                    {
                        c0 = gsl_vector_get(fitCoefficients, 0);
                        c1 = gsl_vector_get(fitCoefficients, 1);
                        stats = gsl_multifit_robust_statistics(gslFitWorkspace);
                        gsl_multifit_robust_free(gslFitWorkspace);
                        // check median absolute deviation and median of signal
                        // Note that median calculation sorts the array, so do this last
                        double mad = stats.sigma_mad; // For full data fitted
                        double mad1 = gsl_stats_mad(model1Values->data, 1, numModel1Points, work1->data); // For first segment
                        double mad2 = gsl_stats_mad(model2Values->data, 1, numModel2Points, work2->data); // For last segment
                        double median1 = gsl_stats_median(model1Values->data, 1, numModel1Points);
                        double median2 = gsl_stats_median(model2Values->data, 1, numModel2Points);
                        fprintf(fitFile, " %f %f %f %f %f %f %f %f %f", c0, c1, stats.adj_Rsq, stats.rmse, median1, median2, mad, mad1, mad2);
                        // Remove the offsets and assign flags for this region
                        for (timeIndex = beginIndex0; timeIndex < endIndex1; timeIndex++)
                        {
                            // remove offset
                            *ADDR(1 + k / 2, k % 2, 2) -= (((TIME() - epoch0)/1000.0) * c1 + c0);
                            driftValue = *ADDR(1 + k / 2, k % 2, 2);
                            // Assign error estimate
                            // TODO: use interpolated MADs derived start and end of pass?
                            viErrors[4*timeIndex + k] = mad;
                            // Set offset-removed flag (most significant bit), and complete region found
                            // Having a complete region can be determined from offset removed and gsl error flags, but
                            // the extra bit will make it logically straightforward to find incomplete regions.
                            // Toggle off the "offset not removed" and "incomplete region" flag bits
                            if (setFlags)
                            {
                                fitInfo[timeIndex] &= ~((FITINFO_OFFSET_NOT_REMOVED | FITINFO_INCOMPLETE_REGION) << (flagIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
                                // Update quality and calibration flags based on thresholds
                                updateDataQualityFlags(satellite, flagIndex, fitargs.regionNumber, driftValue, mad, timeIndex, flags, fitInfo);
                            }
                        }
                    }
                }
                fprintf(fitFile, "\n");
                gsl_matrix_free(modelTimesMatrix);
                gsl_vector_free(model1Values);
                gsl_vector_free(model2Values);
                gsl_vector_free(modelValues);
                gsl_vector_free(work1);
                gsl_vector_free(work2);
                gsl_vector_free(fitCoefficients);
                gsl_matrix_free(cov);
            }
            else
            {
                fprintf(stdout, "%s Fit error: did not get both endpoints of region defined for CDF_EPOCHS %f, %f, %f, %f: not fitting and not removing offsets.\n", infoHeader, tregion11, tregion12, tregion21, tregion22);
                // Fit region flag for incomplete region is already accounted for as complete_region bit is 0
            }
            
            regionBegin = false;
            gotFirstModelData = false;
            gotStartOfSecondModelData = false;
            gotSecondModelData = false;
            if (timeIndex == nRecs)
            {
                break;
            }
            else
            {
                timeIndex++;
            }
            
        }
        
        previousQDLat = QDLAT();
    
    }

}

void updateDataQualityFlags(const char *satellite, uint8_t sensorIndex, uint8_t regionNumber, float driftValue, float mad, long timeIndex, uint16_t *flags, uint32_t *fitInfo)
{
    // Swarm C flags all zero for now
    // Flag is zero if drift magnitude is greater than FLAGS_MAXIMUM_DRIFT_VALUE
    uint16_t flagMask = (1<<sensorIndex);
    bool madOK = (mad < madThreshold(satellite[0], sensorIndex));
    bool magOK = fabs(driftValue) <= FLAGS_MAXIMUM_DRIFT_VALUE;
    // Currently quality flag is set to 1 only for Swarm A and B viy at middle-to-high latitudes (region 0 or region 2)
    if (satellite[0] != 'C' && sensorIndex == 2 && (regionNumber == 0 || regionNumber == 2) && madOK && magOK)
    {
        flags[timeIndex] |= flagMask;
    }
    // Set calibration info flag. Only set to one if baseline offset was subtracted
    if (!madOK)
    {
        fitInfo[timeIndex] |= (FITINFO_MAD_EXCEEDED << (sensorIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
    }
    if (!magOK)
    {
        fitInfo[timeIndex] |= (FITINFO_DRIFT_MAGNITUDE_EXCEEDED << (sensorIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
    }

}

float madThreshold(char satellite, int sensorIndex)
{
    // sensorIndex, per flagging and CDF file content, not as stored in dataBuffers:
    // 0 -> Hx
    // 1 -> Vx
    // 2 -> Hy
    // 3 -> Vy

    // Absolute threshold for noise levels defined as 100 m/s 1-sigma equivalant at 2 Hz
    return 100.0 * sqrtf(8.0);
}
