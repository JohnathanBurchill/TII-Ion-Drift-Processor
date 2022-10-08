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

#include "state.h"
#include "settings.h"
#include "indexing.h"
#include "errors.h"
#include "loadData.h"
#include "export.h"
#include "sample.h"

#include <tii/detector.h>
#include <tii/isp.h>

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <cdf.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>

char infoHeader[50] = {0};

int initQualityData(ProcessorState *state)
{
    // Quality flag and fitInfo flag initialized to zero
    state->flags = (uint16_t*) malloc((size_t) (state->nRecs * sizeof(uint16_t)));
    state->fitInfo = (uint32_t*) malloc((size_t) (state->nRecs * sizeof(uint32_t)));
    // Error estimates from Mean Absolute Deviation (MAD): default is -42. :)
    state->viErrors = (float*) malloc((size_t) (state->nRecs * sizeof(float) * 4));

    if (state->flags == NULL || state->fitInfo == NULL || state->viErrors == NULL)
        return TIICT_MEMORY;

    for (long ind = 0; ind < state->nRecs; ind++)
    {
        state->viErrors[4*ind+0] = DEFAULT_VI_ERROR;
        state->viErrors[4*ind+1] = DEFAULT_VI_ERROR;
        state->viErrors[4*ind+2] = DEFAULT_VI_ERROR;
        state->viErrors[4*ind+3] = DEFAULT_VI_ERROR;
        state->flags[ind] = 0;
        // FITINFO_OFFSET_NOT_REMOVED = 1 and FITINFO_INCOMPLETE_REGION = 1 are the defaults for fitInfo
        // Set for each velocity component
        for (uint8_t k = 0; k < 4; k++)
        {
            state->fitInfo[ind] |= ((FITINFO_OFFSET_NOT_REMOVED | FITINFO_INCOMPLETE_REGION)) << (k * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT);
        }
    }

    return TIICT_OK;
}

int calibrateFlows(ProcessorState *state)
{
    // Allocate memory for quality parameters
    int status = initQualityData(state);
    if (status != TIICT_OK)
        return status;

    // Do the calibration
    uint8_t **dataBuffers = state->dataBuffers;
    double *pEpoch = (double*) dataBuffers[0];
    long long timeIndex = 0;

    // The macros TIME(), QDLAT(), etc. give the value at index timeIndex.
    double epoch0 = TIME();

    // Strategy is to perform calculations in place in memory without allocating space unecessarily.
    // Update 
    //  0) Adjust times by + 1./32 - 0.0875 s
    //  1) VSatXYZ from km/s to m/s
    //  2) bias voltage to nearest of -100 V, -62 V or left alone if too different from those values
    //  3) apply Level 1 calibration with bias dependence
    float value, innerDomeBias;

    // Adjust times for sample lag
    for (long timeIndex = 0; timeIndex < state->nRecs; timeIndex++)
    {
        ((double*)dataBuffers[0])[timeIndex] += ((1./32 - 0.0875)*1000.); // ms
    }

    // Choose a potential estimate (lpPhiSc, lpPhiScHighGain, or lpPhiScLowGain)
    state->potentials = state->lpPhiSc;

    float enxh = 0.0;
    float enxv = 0.0;
    float vixh = 0.0;
    float vixv = 0.0;
    float q = 1.602e-19;
    float mass = 2.67e-26; // O+

    float shx = 0.0, shy = 0.0, svx = 0.0, svy = 0.0;

    // Detector centres
    float xch = 0.0;
    float ych = 0.0;
    float xcv = 0.0;
    float ycv = 0.0;
    detectorCoordinates(state->args.satellite[0], H_SENSOR, &xch, &ych);
    detectorCoordinates(state->args.satellite[0], V_SENSOR, &xcv, &ycv);

    // Estimate raw velocities
    for (long timeIndex = 0; timeIndex < state->nRecs; timeIndex++)
    {
        // Bias Voltage
        // TODO: need a more accurate replacement: i.e., early in mission the voltage was ~-60 V, not -62.
        if (VBIAS() < -95.0)
            *ADDR(13, 0, 1) = -100.;
        else if(VBIAS() > -65. && VBIAS() < -59.0)
            *ADDR(13, 0, 1) = -62.0;

        // VSatXYZ to m/s
        *ADDR(3, 0, 3) *= 1000.0;
        *ADDR(3, 1, 3) *= 1000.0;
        *ADDR(3, 2, 3) *= 1000.0;

        // VCorot to m/s
        *ADDR(10, 0, 3) *= 1000.0;
        *ADDR(10, 1, 3) *= 1000.0;
        *ADDR(10, 2, 3) *= 1000.0;

        // Apply level1 calibration with bias dependence
        innerDomeBias = VBIAS() - VFP();

        // Get scaling parameter
        switch(state->args.satellite[0])
        {
            case 'A': // 20190930 slew experiment, plus simulations
                shx = 574.0;
                shy = 574.0;
                svx = 712.0;
                svy = 712.0;
                if (innerDomeBias >= -63.0 && innerDomeBias < -59.0)
                {
                    shx *= .76;
                    shy *= .76;
                    svx *= .76;
                    svy *= .76;
                }
                break;
            case 'B': // 20210616 slew experiment
                shx = 553.0;
                shy = 553.0;
                svx = 548.0;
                svy = 548.0;
                if (innerDomeBias >= -63.0 && innerDomeBias < -59.0)
                {
                    shx *= 450.5;
                    shy = 450.5; // Calibration 20210624, inner dome bias at -62 V.
                    svx *= 451.0;
        		    svy = 525.0; // Calibration 20210624, inner dome bias at -62 V.
                }
                break;
            case 'C': // 20191126 slew experiment, plus simulations
                shx = 679.0;
                shy = 679.0;
                // Using H sensor sensitivity for Vixv
                svx = 679.0;
                svy = 2377.0;
                if (innerDomeBias >= -63.0 && innerDomeBias < -59.0)
                {
                    shx *= .76;
                    shy *= .76;
                    svx *= .76;
                    svy *= .76;
                }
                break;
        }
        // Cross-track flows do not take into satellite potential
        // Change sign to get flow directions correct, then apply scaling
        // and subtract satellite velocity

        *ADDR(1, 1, 2) = -1.0 * (MYH() - ych) * shy - VSATY();
        *ADDR(2, 1, 2) = -1.0 * (MYV() - ycv) * svy - VSATZ();

        // Along-track flows take into account satellite potential 
        // by first converting flow to energy (eofr)
        // adding the satellite potential estimate from LP
        // then converting to velocity
        // We no longer use the cross-track empirical sensitivity formulas,
        // which are approximations anyway.
        // TODO: revise calibration files to get MCP voltage 
        // NOTE that inner dome bias is not measured but the setting for the H sensor.
        // Remove effect of satellite potential
        // TODO: take into account v-cross-B-dot-dl

        if (state->usePotentials)
        {
            // Calculate ion energy for each sensor (from only the x moment for now)
            // Add in the satellite potential
            // Then remove offsets from this
            // Then convert to flow velocity, adding ram energy of O+ before taking sqare root. 
            *ADDR(1, 0, 2) = eofr(MXH() - xch, innerDomeBias, 0.0) + state->potentials[timeIndex];
            *ADDR(2, 0, 2) = eofr(MXV() - xcv, innerDomeBias, 0.0) + state->potentials[timeIndex];

        }
        else
        {
            // Old way, estimates a proxy based on image moments
            // no potential correction, and this is our offset-biased flow estimate
            *ADDR(1, 0, 2) = -1.0 * (MXH() - 32.5) * shx - VSATX();
            *ADDR(2, 0, 2) = -1.0 * (MXV() - 32.5) * svx - VSATX();
        }

    }

    fprintf(state->processingLogFile, "%sPrepared calibration data.\n", infoHeader);
    fflush(state->processingLogFile);

    // Remove offsets and set calibration flags
    status = removeOffsetsAndSetFlags(state, true);
    if (status != TIICT_OK)
        return status;

    // If using satellite potential, calculate ion along-track drift from offset-corrected energies
    // We have effectively removed 4.8 eV from each energy by setting the energy to 0 at mid-latitude
    // Add it back in before calculating velocity, then remove satellite velocity
    float backgroundRamEnergyeV = 0.0;
    float factor = 0.5 * mass / q;
    if (state->usePotentials)
    {
        for (long timeIndex = 0; timeIndex < state->nRecs; timeIndex++)
        {
            backgroundRamEnergyeV = factor * VSATX() * VSATX();
            // Calculate vix assuming pure O+
            // Positive, is flow towards satellite, in direction of sensor x axis.
            vixh = sqrtf((MXH() + backgroundRamEnergyeV) / factor);
            vixv = sqrtf((MXV() + backgroundRamEnergyeV) / factor);
            // then calculate vi. Note that VSATX is positive toward direction of motion
            // HX
            *ADDR(1, 0, 2) = VSATX() - vixh;
            // HV
            *ADDR(2, 0, 2) = VSATX() - vixv;
        }
    }

    return TIICT_OK;

}

int removeOffsetsAndSetFlags(ProcessorState *state, bool setFlags)
{
    int status = TIICT_OK;

    // Remove offsets and calculate flags
    // 1. Ascending northern polar region to descending northern polar region
    // 2. Descending low-lat region
    // 3. Descending southern polar region to ascending southern polar region
    // 4. Ascending low-lat region
    // 5. Repeat. Does not have to be sequential
    for (uint8_t ind = 0; ind < 4; ind++)
    {
        // Remove offsets again and set flags
        status = removeOffsetsAndSetFlagsForInterval(state, ind, setFlags);
        if (status != TIICT_OK)
            return status;
    }

    fprintf(state->processingLogFile, "%sRemoved offsets and calculated flags.\n", infoHeader);
    fflush(state->fitFile);
    fflush(state->processingLogFile);

    return status;

}

int removeOffsetsAndSetFlagsForInterval(ProcessorState *state, uint8_t interval, bool setFlags)
{
    // Robust linear least squares from GSL: https://www.gnu.org/software/gsl/doc/html/lls.html#examples
    long long timeIndex = 0;
    uint8_t **dataBuffers = state->dataBuffers;
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

    offset_model_fit_arguments *fitargs = &state->fitargs[interval];

    float dlatfirst = fitargs->lat2 - fitargs->lat1;
    float dlatsecond = fitargs->lat4 - fitargs->lat3;

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

    float location = (fitargs->lat2 + fitargs->lat3) / 2.;
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

    for (timeIndex = 0; timeIndex < state->nRecs; timeIndex++)
    {
        c0 = 0.0;
        c1 = 0.0;
        seconds = (TIME() - epoch0) / 1000.;
        if ((firstDirection == 1 && QDLAT() >= fitargs->lat1 && previousQDLat < fitargs->lat1) || (firstDirection == -1 && QDLAT() <= fitargs->lat1 && previousQDLat > fitargs->lat1))
        {
            // Start a new region search
            regionBegin = true; // Found start of region to remove offset from
            gotFirstModelData = false;
            gotStartOfSecondModelData = false;
            gotSecondModelData = false;
            beginIndex0 = timeIndex;
            tregion11 = TIME();

        }
        else if (regionBegin && ((firstDirection == 1 && QDLAT() >= fitargs->lat2 && previousQDLat < fitargs->lat2) || (firstDirection == -1 && QDLAT() <= fitargs->lat2 && previousQDLat > fitargs->lat2)))
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
        else if (gotFirstModelData && ((secondDirection == -1 && QDLAT() <= fitargs->lat3 && previousQDLat > fitargs->lat3) || (secondDirection == 1 && QDLAT() >= fitargs->lat3 && previousQDLat < fitargs->lat3)))
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
        else if (gotStartOfSecondModelData && ((secondDirection == -1 && QDLAT() <= fitargs->lat4 && previousQDLat > fitargs->lat4) || (secondDirection == 1 && QDLAT() >= fitargs->lat4 && previousQDLat < fitargs->lat4)))
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
                fprintf(state->fitFile, "%d %d %lld %lld %f %f %f %f", fitargs->regionNumber, numFits, numModel1Points, numModel2Points, tregion11, tregion12, tregion21, tregion22);
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
                        fprintf(state->processingLogFile, "%sCould not set maximum GSL iterations.\n", infoHeader);
                    }
                    gslStatus = gsl_multifit_robust(modelTimesMatrix, modelValues, fitCoefficients, cov, gslFitWorkspace);
                    if (gslStatus)
                    {
                        toEncodeEPOCH(tregion11, 0, startString);
                        toEncodeEPOCH(tregion22, 0, stopString);
                        fprintf(state->processingLogFile, "%s<GSL Fit Error: %s> for fit region from %s to %s spanning latitudes %.0f to %.0f.\n", infoHeader, gsl_strerror(gslStatus), startString, stopString, fitargs->lat1, fitargs->lat4);
                        // Print "-9999999999.GSLERRORNUMBER" for each of the nine fit parameters
                        fprintf(state->fitFile, " -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d", gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus);
                        if (setFlags)
                        {
                            for (timeIndex = beginIndex0; timeIndex < endIndex1; timeIndex++)
                            {
                                // Got a complete region, but had a fit error
                                state->fitInfo[timeIndex] &= ~(FITINFO_INCOMPLETE_REGION << (flagIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
                                state->fitInfo[timeIndex] |= (FITINFO_GSL_FIT_ERROR << (flagIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
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
                        fprintf(state->fitFile, " %f %f %f %f %f %f %f %f %f", c0, c1, stats.adj_Rsq, stats.rmse, median1, median2, mad, mad1, mad2);
                        // Remove the offsets and assign flags for this region
                        for (timeIndex = beginIndex0; timeIndex < endIndex1; timeIndex++)
                        {
                            // remove offset
                            *ADDR(1 + k / 2, k % 2, 2) -= (((TIME() - epoch0)/1000.0) * c1 + c0);
                            driftValue = *ADDR(1 + k / 2, k % 2, 2);
                            // Assign error estimate
                            // TODO: use interpolated MADs derived start and end of pass?
                            state->viErrors[4*timeIndex + k] = mad;
                            // Set offset-removed flag (most significant bit), and complete region found
                            // Having a complete region can be determined from offset removed and gsl error flags, but
                            // the extra bit will make it logically straightforward to find incomplete regions.
                            // Toggle off the "offset not removed" and "incomplete region" flag bits
                            if (setFlags)
                            {
                                state->fitInfo[timeIndex] &= ~((FITINFO_OFFSET_NOT_REMOVED | FITINFO_INCOMPLETE_REGION) << (flagIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
                                // Update quality and calibration flags based on thresholds
                                updateDataQualityFlags(state, flagIndex, fitargs->regionNumber, driftValue, mad, timeIndex);
                            }
                        }
                    }
                }
                fprintf(state->fitFile, "\n");
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
                fprintf(state->processingLogFile, "%s Fit error: did not get both endpoints of region defined for CDF_EPOCHS %f, %f, %f, %f: not fitting and not removing offsets.\n", infoHeader, tregion11, tregion12, tregion21, tregion22);
                // Fit region flag for incomplete region is already accounted for as complete_region bit is 0
            }
            
            regionBegin = false;
            gotFirstModelData = false;
            gotStartOfSecondModelData = false;
            gotSecondModelData = false;
            if (timeIndex == state->nRecs)
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

    return TIICT_OK;

}

void updateDataQualityFlags(ProcessorState *state, uint8_t sensorIndex, uint8_t regionNumber, float driftValue, float mad, long timeIndex)
{
    uint16_t *flags = state->flags;
    uint32_t *fitInfo = state->fitInfo;
    uint8_t *imageFlags = state->tracisImageFlagsH;
    if (sensorIndex > 1)
        imageFlags = state->tracisImageFlagsV;

    // We lower the quality flag even for peripheral anomalies
    // This flag is obtained from zero-order interpolation of the
    // nearest previous TII full image, every 10 s to 30 s or so
    // for most of the Swarm mission.
    bool imageOK = imageFlags[timeIndex] == 0;

    // Flag is zero if drift magnitude is greater than FLAGS_MAXIMUM_DRIFT_VALUE
    // TODO incorporate ion density from LP_HM into flagging
    uint16_t flagMask = (1<<sensorIndex);
    bool madOK = (mad < madThreshold(state->args.satellite[0], sensorIndex));
    bool magOK = fabs(driftValue) <= FLAGS_MAXIMUM_DRIFT_VALUE;

    // Currently quality flag is set to 1 only for Swarm A and B viy at middle-to-high latitudes (region 0 or region 2)
    // Swarm C flags all zero for now

    if (state->args.satellite[0] != 'C' && sensorIndex == 2 && (regionNumber == 0 || regionNumber == 2) && madOK && magOK && imageOK)
    {
        flags[timeIndex] |= flagMask;
    }
    // Set calibration info flag. Only set to one if baseline offset was subtracted
    if (!imageOK)
    {
        fitInfo[timeIndex] |= (FITINFO_TII_IMAGE_QUALITY << (sensorIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
    }
    if (!madOK)
    {
        fitInfo[timeIndex] |= (FITINFO_MAD_EXCEEDED << (sensorIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
    }
    if (!magOK)
    {
        fitInfo[timeIndex] |= (FITINFO_DRIFT_MAGNITUDE_EXCEEDED << (sensorIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
    }

    return;

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

int initFields(ProcessorState *state)
{
    // Set up new memory
    state->xhat = (float*) malloc((size_t) (state->nRecs * sizeof(float) * 3));
    state->yhat = (float*) malloc((size_t) (state->nRecs * sizeof(float) * 3));
    state->zhat = (float*) malloc((size_t) (state->nRecs * sizeof(float) * 3));
    state->ectFieldH = (float*) malloc((size_t) (state->nRecs * sizeof(float) * 3));
    state->ectFieldV = (float*) malloc((size_t) (state->nRecs * sizeof(float) * 3));
    state->bctField = (float*) malloc((size_t) (state->nRecs * sizeof(float) * 3));

    if (state->xhat == NULL || state->yhat == NULL || state->zhat == NULL || state->ectFieldH == NULL || state->ectFieldV == NULL || state->bctField == NULL)
        return TIICT_MEMORY;

    return TIICT_OK;
}

int calculateFields(ProcessorState *state)
{
    // Init memory for fields
    int status = initFields(state);
    if (status != TIICT_OK)
        return status;

    // Calculate fields
    long long ind;
    uint8_t **dataBuffers = state->dataBuffers;
    float *xhat = state->xhat;
    float *yhat = state->yhat;
    float *zhat = state->zhat;
    float *ectFieldH = state->ectFieldH;
    float *ectFieldV = state->ectFieldV;
    float *bctField = state->bctField;

    for (long timeIndex = 0; timeIndex < state->nRecs; timeIndex++)
    {
        // Calculate xhat, yhat, zhat
        // xhat parallel to satellite velocity 
        float magVsat = sqrtf(VSATN()*VSATN() + VSATE()*VSATE() + VSATC() * VSATC());
        ind = 3*timeIndex;
        xhat[ind + 0] = VSATN() / magVsat;
        xhat[ind + 1] = VSATE() / magVsat;
        xhat[ind + 2] = VSATC() / magVsat;
        // yhat is xhat cross {0, 0, -1}
        yhat[ind + 0] = -xhat[ind + 1];
        yhat[ind + 1] = xhat[ind + 0];
        yhat[ind + 2] = 0.0;
        float magyhat = sqrtf(yhat[ind + 0] * yhat[ind + 0] + yhat[ind + 1] * yhat[ind + 1] + yhat[ind + 2] * yhat[ind + 2]);
        yhat[ind + 0] /= magyhat;
        yhat[ind + 1] /= magyhat;
        yhat[ind + 2] /= magyhat;
        // zhat is the cross-product of x into y
        zhat[ind+0] = xhat[ind+1] * yhat[ind+2] - xhat[ind+2] * yhat[ind+1];
        zhat[ind+1] = -1.0 * xhat[ind+0] * yhat[ind+2] + xhat[ind+2] * yhat[ind+0];
        zhat[ind+2] = xhat[ind+0] * yhat[ind+1] - xhat[ind+1] * yhat[ind+0];
        float magzhat = sqrtf(zhat[ind + 0] * zhat[ind + 0] + zhat[ind + 1] * zhat[ind + 1] + zhat[ind + 2] * zhat[ind + 2]);
        zhat[ind + 0] /= magzhat;
        zhat[ind + 1] /= magzhat;
        zhat[ind + 2] /= magzhat;

        // B field in cross-track frame, nT:
        bctField[ind + 0] = BN() * xhat[ind + 0] + BE() * xhat[ind + 1] + BC() * xhat[ind + 2];
        bctField[ind + 1] = BN() * yhat[ind + 0] + BE() * yhat[ind + 1] + BC() * yhat[ind + 2];
        bctField[ind + 2] = BN() * zhat[ind + 0] + BE() * zhat[ind + 1] + BC() * zhat[ind + 2];

        // E field from H sensor X, in cross-track frame, mV/m:
        ectFieldH[ind + 0] = -1.0 * (MYH() * bctField[ind + 2] - MYV() * bctField[ind + 1]) / 1000000000.0 * 1000.0; 
        ectFieldH[ind + 1] = -1.0 * (-1.0 * MXH() * bctField[ind + 2] + MYV() * bctField[ind + 0]) / 1000000000.0 * 1000.0; 
        ectFieldH[ind + 2] = -1.0 * (MXH() * bctField[ind + 1] - MYH() * bctField[ind + 0]) / 1000000000.0 * 1000.0; 

        // E field from V sensor X, in cross-track frame, mV/m:
        ectFieldV[ind + 0] = -1.0 * (MYH() * bctField[ind + 2] - MYV() * bctField[ind + 1]) / 1000000000.0 * 1000.0; 
        ectFieldV[ind + 1] = -1.0 * (-1.0 * MXV() * bctField[ind + 2] + MYV() * bctField[ind + 0]) / 1000000000.0 * 1000.0; 
        ectFieldV[ind + 2] = -1.0 * (MXV() * bctField[ind + 1] - MYH() * bctField[ind + 0]) / 1000000000.0 * 1000.0; 

    }

    fprintf(state->processingLogFile, "%sCalculated fields.\n", infoHeader);
    fflush(state->processingLogFile);

    return TIICT_OK;

}

int runProcessor(int argc, char *argv[])
{
    int status = TIICT_OK;
    ProcessorState s = {0};
    ProcessorState *state = &s;

    if ((status = initProcessor(argc, argv, state)) != TIICT_OK)
        return shutdown(status, state);

    if ((status = loadData(state)) != TIICT_OK)
        return shutdown(status, state);
 
    if ((status = calibrateFlows(state)) != TIICT_OK)
        return shutdown(status, state);

    if ((status = calculateFields(state)) != TIICT_OK)
        return shutdown(status, state);

    status = exportCdfs(state);

    status = shutdown(status, state);

    return status;


}

int initProcessor(int argc, char *argv[], ProcessorState *state)
{
    int status = TIICT_OK;
    void *args = &state->args;

    // Check arguments and abort if not right
    status = parseArguments(argc, argv, state);
    if (status != TIICT_OK)
        return status;

    status = initLogFiles(state);
    if (status != TIICT_OK)
        return status;

    // Prefix for messages
    initHeader(state);

    // Confirm requested date has records. Abort otherwise.
    status = checkCalDataAvailability(state);
    if (status != TIICT_OK)
        return status;

    // Ensure export directories exist, or abort.
    status = initDirectories(state);
    if (status != TIICT_OK)
        return status;

    // The calibration data memory pointers
    for (uint8_t i = 0; i < NUM_CAL_VARIABLES; i++)
    {
        state->dataBuffers[i] = NULL;
    }
    state->nRecs = 0;

    for (uint8_t i = 0; i < NUM_TRACIS_VARIABLES; i++)
    {
        state->tracisDataBuffers[i] = NULL;
    }
    state->nTracisRecs = 0;

    // Turn off GSL failsafe error handler. We typically check the GSL return codes.
    gsl_set_error_handler_off();

    return TIICT_OK;
}

int parseArguments(int argc, char **argv, ProcessorState *state)
{
    Arguments *args = &state->args;
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--about") == 0)
        {
            fprintf(stdout, "tiict - TII Cross-track ion drift processor, version %s.\n", SOFTWARE_VERSION);
            fprintf(stdout, "Copyright (C) 2022  Johnathan K Burchill\n");
            fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
            fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
            fprintf(stdout, "under the terms of the GNU General Public License.\n");
            return TIICT_ARGS_ABOUT;
        }
    }

    if (argc != 12)
    {
        fprintf(stdout, "usage: %s satLetter year month day calversionString tracisVersionString exportVersionString calDir tracisDir lpDir exportDir\n", argv[0]);
        return TIICT_ARGS_BAD;
    }

    args->satellite = argv[1];
    args->year = atoi(argv[2]);
    args->month = atoi(argv[3]);
    args->day = atoi(argv[4]);
    args->calVersion = argv[5];
    args->tracisVersion = argv[6];
    args->exportVersion = argv[7];
    args->calDir = argv[8];
    args->tracisDir = argv[9];
    args->lpDir = argv[10];
    args->exportDir = argv[11];

    // Check satellite letter
    if (strlen(args->satellite) != 1 || (args->satellite[0] != 'A' && args->satellite[0] != 'B' && args->satellite[0] != 'C'))
    {
        fprintf(stdout, "Satellite must be one of 'A', 'B', or 'C' (no quotes).\n");
        return TIICT_ARGS_SATELLITE;
    }

    return TIICT_OK;

}

void initHeader(ProcessorState *state)
{
    Arguments *args = &state->args;
    time_t currentTime;
    struct tm * timeParts;
    time(&currentTime);
    timeParts = localtime(&currentTime);

    // set up info header
    sprintf(infoHeader, "TIICT %c%s %04d-%02d-%02d: ", args->satellite[0], args->exportVersion, args->year, args->month, args->day);
    fprintf(state->processingLogFile, "\n%s-------------------------------------------------\n", infoHeader);
    fprintf(state->processingLogFile, "%sVersion %s %s\n", infoHeader, args->exportVersion, SOFTWARE_VERSION_STRING);
    fprintf(state->processingLogFile, "%sProcessing date: %s", infoHeader, asctime(timeParts));

    return;
}

int checkResult(int status, ProcessorState *state)
{
    if (status != TIICT_OK)
    {
        if (state->processingLogFile != NULL)
            fprintf(state->processingLogFile, "Error processing file. status = %d\n", status);
        state->returnStatus = status;
        status = shutdown(status, state);
    }

    return status;

}

int shutdown(int status, ProcessorState *state)
{
    if (state == NULL)
        return(EXIT_FAILURE);

    // Close fit log file
    if (state->fitFile != NULL)
    {
        fclose(state->fitFile);
        state->fitFile = NULL;
    }
    if (state->processingLogFile != NULL)
    {
        fprintf(state->processingLogFile, "%sProcessing stopped with status %d\n", infoHeader, status);
        fclose(state->processingLogFile);
        state->processingLogFile = NULL;
    }
    fflush(stdout);

    // Free memory
    for (uint8_t i = 0; i < NUM_CAL_VARIABLES; i++)
        if (state->dataBuffers[i] != NULL)
        {
            free(state->dataBuffers[i]);
            state->dataBuffers[i] = NULL;
        }

    for (uint8_t i = 0; i < NUM_TRACIS_VARIABLES; i++)
        if (state->tracisDataBuffers[i] != NULL)
        {
            free(state->tracisDataBuffers[i]);
            state->tracisDataBuffers[i] = NULL;
        }

    if (state->tracisImageFlagsH != NULL)
    {
        free(state->tracisImageFlagsH);
        state->tracisImageFlagsH = NULL;
    }

    if (state->tracisImageFlagsV != NULL)
    {
        free(state->tracisImageFlagsV);
        state->tracisImageFlagsV = NULL;
    }

    if (state->lpPhiScHighGain != NULL)
    {
        free(state->lpPhiScHighGain);
        state->lpPhiScHighGain = NULL;
    }
    if (state->lpPhiScLowGain != NULL)
    {
        free(state->lpPhiScLowGain);
        state->lpPhiScLowGain = NULL;
    }
    if (state->lpPhiSc != NULL)
    {
        free(state->lpPhiSc);
        state->lpPhiSc = NULL;
    }
    // state->potentials is just a pointer to one of the above potentials

    if (state->xhat != NULL)
    {
        free(state->xhat);
        state->xhat = NULL;
    }
    if (state->yhat != NULL)
    {
        free(state->yhat);
        state->yhat = NULL;
    }
    if (state->zhat != NULL)
    {
        free(state->zhat);
        state->zhat = NULL;
    }
    if (state->ectFieldH != NULL)
    {
        free(state->ectFieldH);
        state->ectFieldH = NULL;
    }
    if (state->ectFieldV != NULL)
    {
        free(state->ectFieldV);
        state->ectFieldV = NULL;
    }
    if (state->bctField != NULL)
    {
        free(state->bctField);
        state->bctField = NULL;
    }
    if (state->viErrors != NULL)
    {
        free(state->viErrors);
        state->viErrors = NULL;
    }
    if (state->flags != NULL)
    {
        free(state->flags);
        state->flags = NULL;
    }
    if (state->fitInfo != NULL)
    {
        free(state->fitInfo);
        state->fitInfo = NULL;
    }

    return status;
}
