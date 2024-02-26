/*

    TII Cross-Track Ion Drift Processor: processing.c

    Copyright (C) 2024  Johnathan K Burchill

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
    state->flags = malloc(state->nRecs * sizeof state->flags);
    state->fitInfo = malloc(state->nRecs * sizeof state->fitInfo);
    state->region = malloc(state->nRecs * sizeof state->region);
    // Error estimates from Mean Absolute Deviation (MAD): default is -42. :)
    state->viErrors = malloc(state->nRecs * sizeof state->viErrors * 4);
    if (state->flags == NULL || state->fitInfo == NULL || state->viErrors == NULL)
    {
        return TIICT_MEMORY;
    }
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

    // Strategy is to perform calculations in place in memory without allocating space unecessarily.
    // Update 
    //  0) Adjust times by + 1./32 - 0.0875 s
    //  1) VSatXYZ from km/s to m/s
    //  2) bias voltage to nearest of -100 V, -62 V or left alone if too different from those values
    //  3) apply Level 1 calibration with bias dependence
    float value, innerDomeBias;
    float vMcpH = -2000.0;
    float vMcpV = -2000.0;

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
                svx = 2377.0;
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
            *ADDR(1, 0, 2) = eofr(MXH() - xch, innerDomeBias, vMcpH) + state->potentials[timeIndex];
            *ADDR(2, 0, 2) = eofr(MXV() - xcv, innerDomeBias, vMcpV) + state->potentials[timeIndex];

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
    state->setFlags = true;
    // Linear offset model
    state->bgws.fitDegree = 2;
    status = removeOffsetsAndSetFlags(state, velocityBackgroundRemoval);
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

int removeOffsetsAndSetFlags(ProcessorState *state, void (*doInterestingStuff)(ProcessorState *))
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
        state->interval = ind;
        status = removeOffsetsAndSetFlagsForInterval(state, doInterestingStuff);
        if (status != TIICT_OK)
            return status;
    }

    fprintf(state->processingLogFile, "%sRemoved offsets and calculated flags.\n", infoHeader);
    fflush(state->fitFile);
    fflush(state->processingLogFile);

    return status;

}

int removeOffsetsAndSetFlagsForInterval(ProcessorState *state, void (*doInterestingStuff)(ProcessorState *))
{
    // Robust linear least squares from GSL: https://www.gnu.org/software/gsl/doc/html/lls.html#examples
    long long timeIndex = 0;
    uint8_t **dataBuffers = state->dataBuffers;
    uint16_t minDataPointsNeeded = 80;
    float previousQDLat = QDLAT();
    bool regionBegin = false;
    state->bgws.epoch0 = TIME();
    double seconds, fitTime;
    int gslStatus;
    bool gotFirstModelData = false;
    bool gotStartOfSecondModelData = false;
    bool gotSecondModelData = false;
    uint16_t numFits = 0;

    offset_model_fit_arguments *fitargs = &state->fitargs[state->interval];

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
    int8_t currentregion; // For flagging
    if (location >= 44. )
        currentregion = 1; // North
    else if (location <= -44.)
        currentregion = -1; // South
    else
        currentregion = 0; // Equator

    for (timeIndex = 0; timeIndex < state->nRecs; timeIndex++)
    {
        seconds = (TIME() - state->bgws.epoch0) / 1000.;
        if ((firstDirection == 1 && QDLAT() >= fitargs->lat1 && previousQDLat < fitargs->lat1) || (firstDirection == -1 && QDLAT() <= fitargs->lat1 && previousQDLat > fitargs->lat1))
        {
            // Start a new region search
            regionBegin = true; // Found start of region to remove offset from
            gotFirstModelData = false;
            gotStartOfSecondModelData = false;
            gotSecondModelData = false;
            state->bgws.beginIndex0 = timeIndex;
            state->bgws.tregion11 = TIME();

        }
        else if (regionBegin && ((firstDirection == 1 && QDLAT() >= fitargs->lat2 && previousQDLat < fitargs->lat2) || (firstDirection == -1 && QDLAT() <= fitargs->lat2 && previousQDLat > fitargs->lat2)))
        {
            if ((TIME() - state->bgws.tregion11)/1000. < (5400. / 2.)) // Should be within 1/2 an orbit of start of segment
            {
                gotFirstModelData = true;
                state->bgws.beginIndex1 = timeIndex;
                state->bgws.tregion12 = TIME();
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
            if ((TIME() - state->bgws.tregion12)/1000. < (5400. / 2.)) // Should be within 1/2 an orbit of start of segment
            {
                gotStartOfSecondModelData = true;
                state->bgws.endIndex0 = timeIndex;
                state->bgws.tregion21 = TIME();
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
            if ((TIME() - state->bgws.tregion21)/1000. < (5400. / 2.)) // Should be within 1/2 an orbit of start of segment
            {
                // We have a complete region - remove linear offset model
                state->bgws.endIndex1 = timeIndex;
                state->bgws.tregion22 = TIME();
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
                state->bgws.numModel1Points = state->bgws.beginIndex1 - state->bgws.beginIndex0;
                state->bgws.numModel2Points = state->bgws.endIndex1 - state->bgws.endIndex0;
                state->bgws.numModelPoints = state->bgws.numModel1Points + state->bgws.numModel2Points;
                fprintf(state->fitFile, "%d %d %lld %lld %f %f %f %f", fitargs->regionNumber, numFits, state->bgws.numModel1Points, state->bgws.numModel2Points, state->bgws.tregion11, state->bgws.tregion12, state->bgws.tregion21, state->bgws.tregion22);
                // Allocate fit buffers
                state->bgws.modelTimesMatrix = gsl_matrix_alloc(state->bgws.numModelPoints, state->bgws.fitDegree);
                state->bgws.modelTimes1Matrix = gsl_matrix_alloc(state->bgws.numModel1Points, state->bgws.fitDegree);
                state->bgws.modelTimes2Matrix = gsl_matrix_alloc(state->bgws.numModel2Points, state->bgws.fitDegree);
                state->bgws.model1Values = gsl_vector_alloc(state->bgws.numModel1Points);
                state->bgws.model2Values = gsl_vector_alloc(state->bgws.numModel2Points);
                state->bgws.work1 = gsl_vector_alloc(state->bgws.numModel1Points);
                state->bgws.work2 = gsl_vector_alloc(state->bgws.numModel2Points);
                state->bgws.modelValues = gsl_vector_alloc(state->bgws.numModelPoints);
                state->bgws.fitCoefficients = gsl_vector_alloc(state->bgws.fitDegree);

                // Load times into the model data buffer
                state->bgws.modelDataIndex = 0;
                for (timeIndex = state->bgws.beginIndex0; timeIndex < state->bgws.beginIndex1; timeIndex++)
                {
                    fitTime = (TIME() - state->bgws.epoch0)/1000.;
                    gsl_matrix_set(state->bgws.modelTimesMatrix, state->bgws.modelDataIndex, 0, 1.0);
                    gsl_matrix_set(state->bgws.modelTimes1Matrix, state->bgws.modelDataIndex, 0, 1.0);
                    gsl_matrix_set(state->bgws.modelTimesMatrix, state->bgws.modelDataIndex, 1, fitTime); // seconds from start of file
                    gsl_matrix_set(state->bgws.modelTimes1Matrix, state->bgws.modelDataIndex++, 1, fitTime); // seconds from start of file
                }
                state->bgws.modelDataMidPoint = state->bgws.modelDataIndex;
                for (timeIndex = state->bgws.endIndex0; timeIndex < state->bgws.endIndex1; timeIndex++)
                {
                    fitTime = (TIME() - state->bgws.epoch0)/1000.;
                    gsl_matrix_set(state->bgws.modelTimes2Matrix, state->bgws.modelDataIndex - state->bgws.modelDataMidPoint, 0, 1.0); 
                    gsl_matrix_set(state->bgws.modelTimesMatrix, state->bgws.modelDataIndex, 0, 1.0);
                    gsl_matrix_set(state->bgws.modelTimes2Matrix, state->bgws.modelDataIndex - state->bgws.modelDataMidPoint, 1, fitTime); 
                    gsl_matrix_set(state->bgws.modelTimesMatrix, state->bgws.modelDataIndex++, 1, fitTime); // seconds from start of file
                }
                // Perform background subtraction
                doInterestingStuff(state);
                fprintf(state->fitFile, "\n");

                gsl_matrix_free(state->bgws.modelTimes1Matrix);
                gsl_matrix_free(state->bgws.modelTimes2Matrix);
                gsl_matrix_free(state->bgws.modelTimesMatrix);
                gsl_vector_free(state->bgws.model1Values);
                gsl_vector_free(state->bgws.model2Values);
                gsl_vector_free(state->bgws.modelValues);
                gsl_vector_free(state->bgws.work1);
                gsl_vector_free(state->bgws.work2);
                gsl_vector_free(state->bgws.fitCoefficients);
            }
            else
            {
                fprintf(state->processingLogFile, "%s Fit error: did not get both endpoints of region defined for CDF_EPOCHS %f, %f, %f, %f: not fitting and not removing offsets.\n", infoHeader, state->bgws.tregion11, state->bgws.tregion12, state->bgws.tregion21, state->bgws.tregion22);
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
    state->xhat = malloc(state->nRecs * sizeof state->xhat * 3);
    state->yhat = malloc(state->nRecs * sizeof state->yhat * 3);
    state->zhat = malloc(state->nRecs * sizeof state->zhat * 3);
    state->ectFieldH = malloc(state->nRecs * sizeof state->ectFieldH* 3);
    state->ectFieldV = malloc(state->nRecs * sizeof state->ectFieldV * 3);
    state->bctField = malloc(state->nRecs * sizeof state->bctField * 3);
    state->geoPotential = malloc(state->nRecs * sizeof state->geoPotential);
    state->maxAbsGeopotentialSlope= malloc(state->nRecs * sizeof state->maxAbsGeopotentialSlope);

    if (state->xhat == NULL || state->yhat == NULL || state->zhat == NULL || state->ectFieldH == NULL || state->ectFieldV == NULL || state->bctField == NULL || state->geoPotential == NULL || state->maxAbsGeopotentialSlope == NULL)
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

    long timeIndex = 0;
    double previousTime = TIME();
    double currentTime = TIME();
    double deltaTime = 0.0;
    state->geoPotential[0] = 0.0;
    float previousExH = 0.0;

    for (timeIndex = 0; timeIndex < state->nRecs; timeIndex++)
    {
        // For geoelectric potential estimation
        currentTime = TIME();
        deltaTime = (currentTime - previousTime) / 1000.0; // seconds
        previousTime = currentTime;

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

        if (timeIndex > 0)
        {
            state->geoPotential[timeIndex] = state->geoPotential[timeIndex-1];
            previousExH = ectFieldH[ind + 0] / 1e3;
        }
        state->geoPotential[timeIndex] += previousExH * magVsat * deltaTime;
    }

    // Remove offsets and set calibration flags
    state->setFlags = true;
    // Linear offset model
    state->bgws.fitDegree = 2;
    status = removeOffsetsAndSetFlags(state, geoelectricPotentialBackgroundRemoval);
    if (status != TIICT_OK)
        return status;

    fprintf(state->processingLogFile, "%sCalculated fields.\n", infoHeader);
    fflush(state->processingLogFile);

    return TIICT_OK;

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


bool downSampleHalfSecond(long *index, long storageIndex, double t0, long maxIndex, uint8_t **dataBuffers, float *ectFieldH, float *ectFieldV, float *geoPotential, float *maxAbsGeopotentialSlope, float *bctField, float *viErrors, float *potentials, uint16_t *flags, uint32_t *fitInfo, uint8_t *region, bool usePotentials)
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
        floatBuf[31] += geoPotential[timeIndex];
        // Take the largest of each 8-sample interval
        if (floatBuf[32] < maxAbsGeopotentialSlope[timeIndex])
            floatBuf[32] = maxAbsGeopotentialSlope[timeIndex];
        // Take latest region in the half-second interval
        floatBuf[33] = region[timeIndex];
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
        geoPotential[storageIndex] = floatBuf[31] / 8.0; // Geoelectric potential H sensor
        maxAbsGeopotentialSlope[storageIndex] = floatBuf[32]; // Take the maximum value                                                        
        region[storageIndex] = floatBuf[33]; // Latest region in the sample                                                            
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

int runProcessor(int argc, char *argv[])
{
    int status = TIICT_OK;
    ProcessorState s = {0};
    ProcessorState *state = &s;

    if ((status = initProcessor(argc, argv, state)) != TIICT_OK)
        return shutdown(status, state);

    if ((status = loadTiiCalData(state)) != TIICT_OK)
        return shutdown(status, state);
 
    if ((status = loadLpCalData(state)) != TIICT_OK)
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

    // Turn off GSL failsafe error handler. We typically check the GSL return codes.
    gsl_set_error_handler_off();

    return TIICT_OK;
}

int parseArguments(int argc, char **argv, ProcessorState *state)
{
    Arguments *args = &state->args;
    // LP estimates of satellite potential
    state->usePotentials = true;

    state->nOptions = 0;
    for (int i = 1; i < argc; i++)
    {
        if (strcmp("--do-not-use-satellite-potential", argv[i]) == 0)
        {
            state->nOptions++;
            state->usePotentials = false;
        }
        else if (strcmp(argv[i], "--about") == 0)
        {
            fprintf(stdout, "tiict - TII Cross-track ion drift processor, version %s.\n", SOFTWARE_VERSION);
            fprintf(stdout, "Copyright (C) 2024  Johnathan K Burchill\n");
            fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
            fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
            fprintf(stdout, "under the terms of the GNU General Public License.\n");
            return TIICT_ARGS_ABOUT;
        }
    }

    if (argc - state->nOptions != 10)
    {
        fprintf(stdout, "usage: %s satLetter year month day calversionString exportVersionString calDir lpDir exportDir\n", argv[0]);
        return TIICT_ARGS_BAD;
    }

    args->satellite = argv[1];
    args->year = atoi(argv[2]);
    args->month = atoi(argv[3]);
    args->day = atoi(argv[4]);
    args->calVersion = argv[5];
    args->exportVersion = argv[6];
    args->calDir = argv[7];
    args->lpDir = argv[8];
    args->exportDir = argv[9];

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
    fprintf(state->processingLogFile, "%sVersion 0302 20220519\n", infoHeader);
    fprintf(state->processingLogFile, "%sProcessing date: %s\n", infoHeader, asctime(timeParts));

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
        fclose(state->processingLogFile);
        state->processingLogFile = NULL;
    }
    fflush(stdout);

    // Free the memory
    for (uint8_t i = 0; i < NUM_CAL_VARIABLES; i++)
        if (state->dataBuffers[i] != NULL)
        {
            free(state->dataBuffers[i]);
            state->dataBuffers[i] = NULL;
        }

    if (state->lpTimes != NULL)
    {
        free(state->lpTimes);
        state->lpTimes = NULL;
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
    if (state->geoPotential != NULL)
    {
        free(state->geoPotential);
        state->geoPotential = NULL;
    }
    if (state->maxAbsGeopotentialSlope != NULL)
    {
        free(state->maxAbsGeopotentialSlope);
        state->maxAbsGeopotentialSlope = NULL;
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
    if (state->region != NULL)
    {
        free(state->region);
        state->region = NULL;
    }

    return status;
}

void velocityBackgroundRemoval(ProcessorState *state)
{
    long long timeIndex = 0;
    uint8_t **dataBuffers = state->dataBuffers;

    int gslStatus = GSL_SUCCESS;

    // Buffers for mid-latitude linear fit data
    const gsl_multifit_robust_type * fitType = gsl_multifit_robust_bisquare;
    gsl_multifit_robust_workspace * gslFitWorkspace;
    gsl_multifit_robust_stats stats;
    BackgroundRemovalWorkspace_t *ws = &state->bgws;
    gsl_matrix *cov = gsl_matrix_alloc(ws->fitDegree, ws->fitDegree);
    double c0, c1, cov00, cov01, cov11, sumsq;
    float driftValue = 0.0;

    offset_model_fit_arguments *fitargs = &state->fitargs[state->interval];

    // Load values into model data buffer once each for HX, HY, VX, VY
    for (uint8_t k = 0; k < 4; k++)
    {
        uint8_t flagIndex = k; // Defined flags as bit 0 -> HX, bit 1 -> VX, bit 2-> HY and bit 3-> VY. Data are stored in memory differently.
        if (flagIndex == 1)
            flagIndex = 2;
        else if (flagIndex == 2)
            flagIndex = 1;

        ws->modelDataIndex = 0;
        for (timeIndex = ws->beginIndex0; timeIndex < ws->beginIndex1; timeIndex++)
        {
            gsl_vector_set(ws->model1Values, ws->modelDataIndex, *ADDR(1 + k / 2, k % 2, 2)); 
            gsl_vector_set(ws->modelValues, ws->modelDataIndex++, *ADDR(1 + k / 2, k % 2, 2)); 
        }
        ws->modelDataMidPoint = ws->modelDataIndex;
        for (timeIndex = ws->endIndex0; timeIndex < ws->endIndex1; timeIndex++)
        {
            gsl_vector_set(ws->model2Values, ws->modelDataIndex - ws->modelDataMidPoint, *ADDR(1 + k / 2, k % 2, 2)); 
            gsl_vector_set(ws->modelValues, ws->modelDataIndex++, *ADDR(1 + k / 2, k % 2, 2)); 
        }
        // Robust linear model fit and removal
        gslFitWorkspace = gsl_multifit_robust_alloc(fitType, ws->numModelPoints, ws->fitDegree);
        gslStatus = gsl_multifit_robust_maxiter(GSL_FIT_MAXIMUM_ITERATIONS, gslFitWorkspace);
        if (gslStatus)
        {
            fprintf(state->processingLogFile, "%sCould not set maximum GSL iterations.\n", infoHeader);
        }
        gslStatus = gsl_multifit_robust(ws->modelTimesMatrix, ws->modelValues, ws->fitCoefficients, cov, gslFitWorkspace);
        if (gslStatus)
        {
            toEncodeEPOCH(ws->tregion11, 0, ws->startString);
            toEncodeEPOCH(ws->tregion22, 0, ws->stopString);
            fprintf(state->processingLogFile, "%s<GSL Fit Error: %s> for fit region from %s to %s spanning latitudes %.0f to %.0f.\n", infoHeader, gsl_strerror(gslStatus), ws->startString, ws->stopString, fitargs->lat1, fitargs->lat4);
            // Print "-9999999999.GSLERRORNUMBER" for each of the nine fit parameters
            fprintf(state->fitFile, " -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d", gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus);
            if (state->setFlags)
            {
                for (timeIndex = ws->beginIndex0; timeIndex < ws->endIndex1; timeIndex++)
                {
                    // Got a complete region, but had a fit error
                    state->fitInfo[timeIndex] &= ~(FITINFO_INCOMPLETE_REGION << (flagIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
                    state->fitInfo[timeIndex] |= (FITINFO_GSL_FIT_ERROR << (flagIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
                }
            }
        }
        else
        {
            c0 = gsl_vector_get(ws->fitCoefficients, 0);
            c1 = gsl_vector_get(ws->fitCoefficients, 1);
            stats = gsl_multifit_robust_statistics(gslFitWorkspace);
            gsl_multifit_robust_free(gslFitWorkspace);
            // check median absolute deviation and median of signal
            // Note that median calculation sorts the array, so do this last
            double mad = stats.sigma_mad; // For full data fitted
            double mad1 = gsl_stats_mad(ws->model1Values->data, 1, ws->numModel1Points, ws->work1->data); // For first segment
            double mad2 = gsl_stats_mad(ws->model2Values->data, 1, ws->numModel2Points, ws->work2->data); // For last segment
            double median1 = gsl_stats_median(ws->model1Values->data, 1, ws->numModel1Points);
            double median2 = gsl_stats_median(ws->model2Values->data, 1, ws->numModel2Points);
            fprintf(state->fitFile, " %f %f %f %f %f %f %f %f %f", c0, c1, stats.adj_Rsq, stats.rmse, median1, median2, mad, mad1, mad2);
            // Remove the offsets and assign flags for this region
            for (timeIndex = ws->beginIndex0; timeIndex < ws->endIndex1; timeIndex++)
            {
                // remove offset
                *ADDR(1 + k / 2, k % 2, 2) -= (((TIME() - ws->epoch0)/1000.0) * c1 + c0);
                driftValue = *ADDR(1 + k / 2, k % 2, 2);
                // Assign error estimate
                // TODO: use interpolated MADs derived start and end of pass?
                state->viErrors[4*timeIndex + k] = mad;
                // Set offset-removed flag (most significant bit), and complete region found
                // Having a complete region can be determined from offset removed and gsl error flags, but
                // the extra bit will make it logically straightforward to find incomplete regions.
                // Toggle off the "offset not removed" and "incomplete region" flag bits
                if (state->setFlags)
                {
                    state->fitInfo[timeIndex] &= ~((FITINFO_OFFSET_NOT_REMOVED | FITINFO_INCOMPLETE_REGION) << (flagIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
                    // Update quality and calibration flags based on thresholds
                    updateDataQualityFlags(state->args.satellite, flagIndex, fitargs->regionNumber, driftValue, mad, timeIndex, state->flags, state->fitInfo);
                }
            }
        }
    }
    gsl_matrix_free(cov);

    return;

}

void geoelectricPotentialBackgroundRemoval(ProcessorState *state)
{
    long long timeIndex = 0;
    uint8_t **dataBuffers = state->dataBuffers;

    int gslStatus = GSL_SUCCESS;

    // Buffers for mid-latitude linear fit data
    const gsl_multifit_robust_type * fitType = gsl_multifit_robust_bisquare;
    gsl_multifit_robust_workspace * gslFitWorkspace;
    gsl_multifit_robust_workspace * gslFitWorkspace1;
    gsl_multifit_robust_workspace * gslFitWorkspace2;
    gsl_multifit_robust_stats stats;
    BackgroundRemovalWorkspace_t *ws = &state->bgws;
    gsl_matrix *cov = gsl_matrix_alloc(ws->fitDegree, ws->fitDegree);
    double c0, c1, cov00, cov01, cov11, sumsq;
    double slope1 = 0.0;
    double slope2 = 0.0;

    offset_model_fit_arguments *fitargs = &state->fitargs[state->interval];

    float vmag1 = 0.0;
    float vmag2 = 0.0;

    // Set the region number
    // 0 northern ascending, 1 equatorial descending
    // 2 southern descending, 3 equatorial ascending
    for (timeIndex = ws->beginIndex0; timeIndex < ws->endIndex1; timeIndex++)
    {
        state->region[timeIndex] = fitargs->regionNumber; 
    }

    // Load values into model data buffer once each for HX, HY, VX, VY
    ws->modelDataIndex = 0;
    vmag1 = 0.0;
    vmag2 = 0.0;
    for (timeIndex = ws->beginIndex0; timeIndex < ws->beginIndex1; timeIndex++)
    {
        gsl_vector_set(ws->model1Values, ws->modelDataIndex, state->geoPotential[timeIndex]); 
        gsl_vector_set(ws->modelValues, ws->modelDataIndex++, state->geoPotential[timeIndex]); 
        vmag1 += sqrtf(VSATN()*VSATN() + VSATE()*VSATE() + VSATC() * VSATC());
    }
    if (ws->numModel1Points > 0)
        vmag1 /= (float)ws->numModel1Points;
    else
        vmag1 = 0.0;
    ws->modelDataMidPoint = ws->modelDataIndex;
    for (timeIndex = ws->endIndex0; timeIndex < ws->endIndex1; timeIndex++)
    {
        gsl_vector_set(ws->model2Values, ws->modelDataIndex - ws->modelDataMidPoint, state->geoPotential[timeIndex]); 
        gsl_vector_set(ws->modelValues, ws->modelDataIndex++, state->geoPotential[timeIndex]); 
        vmag2 += sqrtf(VSATN()*VSATN() + VSATE()*VSATE() + VSATC() * VSATC());
    }
    if (ws->numModel2Points > 0)
        vmag2 /= (float)ws->numModel2Points;
    else
        vmag2 = 0.0;

    // Robust linear model fit and removal
    gslFitWorkspace = gsl_multifit_robust_alloc(fitType, ws->numModelPoints, ws->fitDegree);
    gslFitWorkspace1 = gsl_multifit_robust_alloc(fitType, ws->numModel1Points, ws->fitDegree);
    gslFitWorkspace2 = gsl_multifit_robust_alloc(fitType, ws->numModel2Points, ws->fitDegree);
    gslStatus = gsl_multifit_robust_maxiter(GSL_FIT_MAXIMUM_ITERATIONS, gslFitWorkspace);
    gslStatus = gsl_multifit_robust_maxiter(GSL_FIT_MAXIMUM_ITERATIONS, gslFitWorkspace1);
    gslStatus = gsl_multifit_robust_maxiter(GSL_FIT_MAXIMUM_ITERATIONS, gslFitWorkspace2);
    gslStatus = gsl_multifit_robust(ws->modelTimesMatrix, ws->modelValues, ws->fitCoefficients, cov, gslFitWorkspace);
    if (gslStatus)
    {
        toEncodeEPOCH(ws->tregion11, 0, ws->startString);
        toEncodeEPOCH(ws->tregion22, 0, ws->stopString);
        fprintf(state->processingLogFile, "%s<GSL Fit Error: %s> for fit region from %s to %s spanning latitudes %.0f to %.0f.\n", infoHeader, gsl_strerror(gslStatus), ws->startString, ws->stopString, fitargs->lat1, fitargs->lat4);
        // Print "-9999999999.GSLERRORNUMBER" for each of the nine fit parameters
        fprintf(state->fitFile, " -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d", gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus);
    }
    else
    {
        c0 = gsl_vector_get(ws->fitCoefficients, 0);
        c1 = gsl_vector_get(ws->fitCoefficients, 1);
        stats = gsl_multifit_robust_statistics(gslFitWorkspace);
        // check median absolute deviation and median of signal
        // Note that median calculation sorts the array, so do this last
        double mad = stats.sigma_mad; // For full data fitted
        double mad1 = gsl_stats_mad(ws->model1Values->data, 1, ws->numModel1Points, ws->work1->data); // For first segment
        double mad2 = gsl_stats_mad(ws->model2Values->data, 1, ws->numModel2Points, ws->work2->data); // For last segment
        double median1 = gsl_stats_median(ws->model1Values->data, 1, ws->numModel1Points);
        double median2 = gsl_stats_median(ws->model2Values->data, 1, ws->numModel2Points);
        gslStatus = gsl_multifit_robust(ws->modelTimes1Matrix, ws->model1Values, ws->fitCoefficients, cov, gslFitWorkspace1);
        if (gslStatus == GSL_SUCCESS && vmag1 > 0)
            slope1 = fabs(gsl_vector_get(ws->fitCoefficients, 1))/vmag1 * 1000.0; // mV / m
        else
            slope1 = 999999999.0;
        gslStatus = gsl_multifit_robust(ws->modelTimes2Matrix, ws->model2Values, ws->fitCoefficients, cov, gslFitWorkspace2);
        if (gslStatus == GSL_SUCCESS && vmag2 > 0)
            slope2 = fabs(gsl_vector_get(ws->fitCoefficients, 1))/vmag2 * 1000.0; // mV / m
        else
            slope2 = 999999999.0;
        // Remove the offsets and estimate max abs mean ex at mid-latitude 
        for (timeIndex = ws->beginIndex0; timeIndex < ws->endIndex1; timeIndex++)
        {
            // remove offset
            state->geoPotential[timeIndex] -= (((TIME() - ws->epoch0)/1000.0) * c1 + c0);
            state->maxAbsGeopotentialSlope[timeIndex] = slope1 > slope2 ? slope1 : slope2;
            // TODO convert to V / m?
        }
    }
    gsl_multifit_robust_free(gslFitWorkspace);
    gsl_multifit_robust_free(gslFitWorkspace1);
    gsl_multifit_robust_free(gslFitWorkspace2);
    gsl_matrix_free(cov);

    return;

}

