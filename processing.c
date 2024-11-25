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
#include "errors.h"
#include "loadData.h"
#include "export.h"
#include "visualize.h"

#include <time.h>
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

int initQualityData(ProcessorVariables_t *var)
{
    // Quality flag and fitInfo flag initialized to zero
    var->flags = malloc(var->nRecs * sizeof *var->flags);
    var->fitInfo = malloc(var->nRecs * sizeof *var->fitInfo);
    var->orbitRegion = malloc(var->nRecs * sizeof *var->orbitRegion);
    // Error estimates from Mean Absolute Deviation (MAD): default is -42. :)
    var->vixherror = malloc(var->nRecs * sizeof *var->vixherror);
    var->vixverror = malloc(var->nRecs * sizeof *var->vixverror);
    var->viyerror = malloc(var->nRecs * sizeof *var->viyerror);
    var->vizerror = malloc(var->nRecs * sizeof *var->vizerror);
    if (var->flags == NULL || var->fitInfo == NULL || var->vixherror == NULL || var->vixverror == NULL || var->viyerror == NULL || var->vizerror == NULL)
    {
        return TIICT_MEMORY;
    }
    for (long ind = 0; ind < var->nRecs; ind++)
    {
        var->vixherror[ind] = DEFAULT_VI_ERROR;
        var->vixverror[ind] = DEFAULT_VI_ERROR;
        var->viyerror[ind] = DEFAULT_VI_ERROR;
        var->vizerror[ind] = DEFAULT_VI_ERROR;
        var->flags[ind] = 0;
        var->orbitRegion[ind] = 255; // Invalid or incomplete region
        // FITINFO_OFFSET_NOT_REMOVED = 1 and FITINFO_INCOMPLETE_REGION = 1 are the defaults for fitInfo
        // Set for each velocity component
        for (uint8_t k = 0; k < 4; k++)
        {
            var->fitInfo[ind] |= ((FITINFO_OFFSET_NOT_REMOVED | FITINFO_INCOMPLETE_REGION)) << (k * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT);
        }
    }

    return TIICT_OK;
}

int calibrateFlows(ProcessorState *state)
{
    // Allocate memory for quality parameters
    int status = initQualityData(&state->vars16hz);
    if (status != TIICT_OK)
        return status;

    // Do the calibration
    double *pEpoch = state->vars16hz.timestamp;
    long long i = 0;

    // Strategy is to perform calculations in place in memory without allocating space unecessarily.
    // Update
    //  0) Adjust times by + 1./32 - 0.0875 s
    //  1) VSatXYZ from km/s to m/s
    //  2) bias voltage to nearest of -100 V, -62 V or left alone if too different from those values
    //  3) apply Level 1 calibration with bias dependence
    float value, innerDomeBias;

    // TII calibration file times have already been adjusted for sample lag

    // Choose a potential estimate (lpPhiSc, lpPhiScHighGain, or lpPhiScLowGain)
    switch (state->vars16hz.lpPotentialSource) {
        case LP_POTENTIAL_HIGHGAIN:
            state->vars16hz.potentials = state->vars16hz.lpPhiScHighGain;
            break;
        case LP_POTENTIAL_LOWGAIN:
            state->vars16hz.potentials = state->vars16hz.lpPhiScLowGain;
            break;
        case LP_POTENTIAL_U_SC:
            state->vars16hz.potentials = state->vars16hz.lpPhiSc;
            break;
        default:
            state->vars16hz.potentials = NULL;
            break;
    }

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
    ProcessorVariables_t *v = &state->vars16hz;
    for (long i = 0; i < v->nRecs; i++)
    {
        // Bias Voltage
        // TODO: need a more accurate replacement: i.e., early in mission the voltage was ~-60 V, not -62.
        //
        if (v->vbiash[i] < -95.0) {
            v->vbiash[i] = -100.0;
        }
        else if(v->vbiash[i] > -65. && v->vbiash[i] < -55.0) {
            v->vbiash[i] = -62.0;
        }
        if (v->vbiasv[i] < -95.0) {
            v->vbiasv[i] = -100.0;
        }
        else if(v->vbiasv[i] > -65. && v->vbiasv[i] < -55.0) {
            v->vbiasv[i] = -62.0;
        }

        // VSatXYZ to m/s
        v->vsatx[i] *= 1000.0;
        v->vsaty[i] *= 1000.0;
        v->vsatz[i] *= 1000.0;

        // VCorot to m/s
        v->vicrx[i] *= 1000.0;
        v->vicry[i] *= 1000.0;
        v->vicrz[i] *= 1000.0;

        // Apply level1 calibration with bias dependence
        // Assumes H and V sensors have same bias, which is okay
        // since we set it to a set number due to known noise levels in monitors.
        // Just use H sensor here
        innerDomeBias = v->vbiash[i] - v->vfp[i];

        // Get scaling parameter
        switch(state->args.satellite[0])
        {
            case 'A': // 20190930 slew experiment, plus simulations
                shx = 574.0;
                shy = 574.0;
                svx = 712.0;
                svy = 712.0;
                if (innerDomeBias >= -63.0 && innerDomeBias < -55.0)
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
                if (innerDomeBias >= -63.0 && innerDomeBias < -55.0)
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
                if (innerDomeBias >= -63.0 && innerDomeBias < -55.0)
                {
                    shx *= .76;
                    shy *= .76;
                    svx *= .76;
                    svy *= .76;
                }
                break;
        }
        // Cross-track flows do not take into account variations in satellite potential
        // Change sign to get flow directions correct, then apply scaling
        // and subtract satellite velocity

        v->viy[i] = -1.0 * (v->myh[i] - ych) * shy - v->vsaty[i];
        v->viz[i] = -1.0 * (v->myv[i] - ycv) * svy - v->vsatz[i];

        // Along-track flows should take into account satellite potential
        // by first converting flow to energy (eofr)
        // adding the satellite potential estimate from LP then converting to velocity
        // We no longer use the cross-track empirical sensitivity formulas for along-track drift

        if (state->useEofR)
        {
            // Calculate ion energy for each sensor (from only the x moment for now)
            // Add in the satellite potential
            // Then remove offsets from this
            // Then convert to flow velocity, adding ram energy of O+ before taking sqare root.
            // Use eofr estimate, no correction for variations in satellite potential
            v->enh[i] = v->enhRaw[i] = eofr(v->mxh[i] - xch, innerDomeBias, v->vmcph[i]);
            v->env[i] = v->envRaw[i] = eofr(v->mxv[i] - xcv, innerDomeBias, v->vmcpv[i]);
            if (state->usePotentials)
            {
                // TODO include emf?
                v->enh[i] += v->potentials[i];
                v->env[i] += v->potentials[i];
            }
        }
        else
        {
            // Old way, estimates a proxy based on image moments
            // no potential correction, and this is our offset-biased flow estimate
            v->vixh[i] = -1.0 * (v->mxh[i] - 32.5) * shx - v->vsatx[i];
            v->vixv[i] = -1.0 * (v->mxv[i] - 32.5) * svx - v->vsatx[i];
            v->enhRaw[i] = v->envRaw[i] = v->enh[i] = v->env[i] = 0.0;
        }

    }

    if (state->writeLogFiles) {
        fprintf(state->processingLogFile, "%sPrepared calibration data.\n", infoHeader);
        fflush(state->processingLogFile);
    }

    // Remove offsets and set calibration flags
    state->setFlags = true;
    // Linear offset model
    state->bgws.fitDegree = 2;
    status = removeOffsetsAndSetFlags(state, backgroundRemoval);
    if (status != TIICT_OK)
        return status;

    // If using EofR method, calculate ion along-track drift from offset-corrected energies
    // We have effectively removed 4.8 eV from each energy by setting the energy to 0 at mid-latitude
    // Add it back in before calculating velocity, then remove satellite velocity
    float backgroundRamEnergyeV = 0.0;
    float factor = 0.5 * mass / q;
    if (state->useEofR)
    {
        for (long i = 0; i < v->nRecs; i++)
        {
            backgroundRamEnergyeV = factor * v->vsatx[i] * v->vsatx[i];
            // Calculate vix assuming pure O+
            // Positive, is flow towards satellite, in direction of sensor x axis.
            // Then calculate vi. Note that VSATX is positive toward direction of motion
            v->vixh[i] = v->vsatx[i] - sqrtf((v->enh[i] + backgroundRamEnergyeV) / factor);
            v->vixv[i] = v->vsatx[i] - sqrtf((v->env[i] + backgroundRamEnergyeV) / factor);
        }
    }

    return TIICT_OK;

}

int removeOffsetsAndSetFlags(ProcessorState *state, int (*doInterestingStuff)(ProcessorState *))
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

    if (state->writeLogFiles) {
        fprintf(state->processingLogFile, "%sRemoved offsets and calculated flags.\n", infoHeader);
        fflush(state->fitFile);
        fflush(state->processingLogFile);
    }

    return status;

}

int removeOffsetsAndSetFlagsForInterval(ProcessorState *state, int (*processRegion)(ProcessorState *))
{
    // Robust linear least squares from GSL: https://www.gnu.org/software/gsl/doc/html/lls.html#examples
    int status = TIICT_OK;
    uint16_t minDataPointsNeeded = 80;
    ProcessorVariables_t *v = &state->vars16hz;
    float previousQDLat = v->qdlat[0];
    bool regionBegin = false;
    state->bgws.epoch0 = v->timestamp[0];
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

    for (int i = 0; i < v->nRecs; i++)
    {
        seconds = (v->timestamp[i] - state->bgws.epoch0) / 1000.;
        if ((firstDirection == 1 && v->qdlat[i] >= fitargs->lat1 && previousQDLat < fitargs->lat1) || (firstDirection == -1 && v->qdlat[i] <= fitargs->lat1 && previousQDLat > fitargs->lat1))
        {
            // Start a new region search
            regionBegin = true; // Found start of region to remove offset from
            gotFirstModelData = false;
            gotStartOfSecondModelData = false;
            gotSecondModelData = false;
            state->bgws.beginIndex0 = i;
            state->bgws.tregion11 = v->timestamp[i];

        }
        else if (regionBegin && ((firstDirection == 1 && v->qdlat[i] >= fitargs->lat2 && previousQDLat < fitargs->lat2) || (firstDirection == -1 && v->qdlat[i] <= fitargs->lat2 && previousQDLat > fitargs->lat2)))
        {
            if ((v->timestamp[i] - state->bgws.tregion11)/1000. < (5400. / 2.)) // Should be within 1/2 an orbit of start of segment
            {
                gotFirstModelData = true;
                state->bgws.beginIndex1 = i;
                state->bgws.tregion12 = v->timestamp[i];
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
        else if (gotFirstModelData && ((secondDirection == -1 && v->qdlat[i] <= fitargs->lat3 && previousQDLat > fitargs->lat3) || (secondDirection == 1 && v->qdlat[i] >= fitargs->lat3 && previousQDLat < fitargs->lat3)))
        {
            if ((v->timestamp[i] - state->bgws.tregion12)/1000. < (5400. / 2.)) // Should be within 1/2 an orbit of start of segment
            {
                gotStartOfSecondModelData = true;
                state->bgws.endIndex0 = i;
                state->bgws.tregion21 = v->timestamp[i];
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
        else if (gotStartOfSecondModelData && ((secondDirection == -1 && v->qdlat[i] <= fitargs->lat4 && previousQDLat > fitargs->lat4) || (secondDirection == 1 && v->qdlat[i] >= fitargs->lat4 && previousQDLat < fitargs->lat4)))
        {
            if ((v->timestamp[i] - state->bgws.tregion21)/1000. < (5400. / 2.)) // Should be within 1/2 an orbit of start of segment
            {
                // We have a complete region - remove linear offset model
                state->bgws.endIndex1 = i;
                state->bgws.tregion22 = v->timestamp[i];
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
                if (state->writeLogFiles) {
                    fprintf(state->fitFile, "%d %d %lld %lld %f %f %f %f", fitargs->regionNumber, numFits, state->bgws.numModel1Points, state->bgws.numModel2Points, state->bgws.tregion11, state->bgws.tregion12, state->bgws.tregion21, state->bgws.tregion22);
                }
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
                for (i = state->bgws.beginIndex0; i < state->bgws.beginIndex1; i++)
                {
                    fitTime = (v->timestamp[i] - state->bgws.epoch0)/1000.;
                    gsl_matrix_set(state->bgws.modelTimesMatrix, state->bgws.modelDataIndex, 0, 1.0);
                    gsl_matrix_set(state->bgws.modelTimes1Matrix, state->bgws.modelDataIndex, 0, 1.0);
                    gsl_matrix_set(state->bgws.modelTimesMatrix, state->bgws.modelDataIndex, 1, fitTime); // seconds from start of file
                    gsl_matrix_set(state->bgws.modelTimes1Matrix, state->bgws.modelDataIndex++, 1, fitTime); // seconds from start of file
                }
                state->bgws.modelDataMidPoint = state->bgws.modelDataIndex;
                for (i = state->bgws.endIndex0; i < state->bgws.endIndex1; i++)
                {
                    fitTime = (v->timestamp[i] - state->bgws.epoch0)/1000.;
                    gsl_matrix_set(state->bgws.modelTimes2Matrix, state->bgws.modelDataIndex - state->bgws.modelDataMidPoint, 0, 1.0);
                    gsl_matrix_set(state->bgws.modelTimesMatrix, state->bgws.modelDataIndex, 0, 1.0);
                    gsl_matrix_set(state->bgws.modelTimes2Matrix, state->bgws.modelDataIndex - state->bgws.modelDataMidPoint, 1, fitTime);
                    gsl_matrix_set(state->bgws.modelTimesMatrix, state->bgws.modelDataIndex++, 1, fitTime); // seconds from start of file
                }
                // Perform regional analysis
                status = processRegion(state);
                if (status != TIICT_OK) {
                    return status;
                }
                if (state->writeLogFiles) {
                    fprintf(state->fitFile, "\n");
                }

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
                if (state->writeLogFiles) {
                    fprintf(state->processingLogFile, "%s Fit error: did not get both endpoints of region defined for CDF_EPOCHS %f, %f, %f, %f: not fitting and not removing offsets.\n", infoHeader, state->bgws.tregion11, state->bgws.tregion12, state->bgws.tregion21, state->bgws.tregion22);
                }
                // Fit region flag for incomplete region is already accounted for as complete_region bit is 0
            }

            regionBegin = false;
            gotFirstModelData = false;
            gotStartOfSecondModelData = false;
            gotSecondModelData = false;
            if (i == v->nRecs)
            {
                break;
            }
            else
            {
                i++;
            }

        }

        previousQDLat = v->qdlat[i];

    }

    return TIICT_OK;

}

void updateDataQualityFlags(const char *satellite, uint8_t sensorIndex, uint8_t regionNumber, float driftValue, float mad, long i, uint16_t *flags, uint32_t *fitInfo)
{
    // Swarm C flags all zero for now
    // Flag is zero if drift magnitude is greater than FLAGS_MAXIMUM_DRIFT_VALUE
    uint16_t flagMask = (1<<sensorIndex);
    bool madOK = (mad < madThreshold(satellite[0], sensorIndex));
    bool magOK = fabs(driftValue) <= FLAGS_MAXIMUM_DRIFT_VALUE;
    // Currently quality flag is set to 1 only for Swarm A and B viy at middle-to-high latitudes (region 0 or region 2)
    if (satellite[0] != 'C' && sensorIndex == 2 && (regionNumber == 0 || regionNumber == 2) && madOK && magOK)
    {
        flags[i] |= flagMask;
    }
    // Set calibration info flag. Only set to one if baseline offset was subtracted
    if (!madOK)
    {
        fitInfo[i] |= (FITINFO_MAD_EXCEEDED << (sensorIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
    }
    if (!magOK)
    {
        fitInfo[i] |= (FITINFO_DRIFT_MAGNITUDE_EXCEEDED << (sensorIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
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

int initFields(ProcessorVariables_t *var)
{
    // Set up new memory
    var->xhat = malloc(var->nRecs * (sizeof *var->xhat) * 3);
    var->yhat = malloc(var->nRecs * (sizeof *var->yhat) * 3);
    var->zhat = malloc(var->nRecs * (sizeof *var->zhat) * 3);
    var->ectxh = malloc(var->nRecs * sizeof *var->ectxh);
    var->ectyh = malloc(var->nRecs * sizeof *var->ectyh);
    var->ectzh = malloc(var->nRecs * sizeof *var->ectzh);
    var->ectxv = malloc(var->nRecs * sizeof *var->ectxv);
    var->ectyv = malloc(var->nRecs * sizeof *var->ectyv);
    var->ectzv = malloc(var->nRecs * sizeof *var->ectzv);
    var->bctx = malloc(var->nRecs * sizeof *var->bctx);
    var->bcty = malloc(var->nRecs * sizeof *var->bcty);
    var->bctz = malloc(var->nRecs * sizeof *var->bctz);
    var->geoelectricPotential = malloc(var->nRecs * sizeof *var->geoelectricPotential);
    var->geoelectricPotentialDifference = malloc(var->nRecs * sizeof *var->geoelectricPotentialDifference);
    var->ehxAdjusted = malloc(var->nRecs * sizeof *var->ehxAdjusted);
    var->ehxAdjustmentParameter= malloc(var->nRecs * sizeof *var->ehxAdjustmentParameter);
    var->maxAbsGeoelectricPotentialBaselineSlope = malloc(var->nRecs * sizeof *var->maxAbsGeoelectricPotentialBaselineSlope);
    var->geoelectricPotentialDetrended = malloc(var->nRecs * sizeof *var->geoelectricPotentialDetrended);
    var->maxAbsGeoelectricPotentialDetrendedBaselineSlope = malloc(var->nRecs * sizeof *var->maxAbsGeoelectricPotentialDetrendedBaselineSlope);

    if (var->xhat == NULL || var->yhat == NULL || var->zhat == NULL || var->ectxh == NULL || var->ectyh == NULL || var->ectzh == NULL
        || var->ectxv == NULL || var->ectyv == NULL || var->ectzv == NULL || var->bctx == NULL || var->bctx == NULL || var->bctx == NULL
        || var->geoelectricPotential == NULL || var->maxAbsGeoelectricPotentialBaselineSlope == NULL
        || var->maxAbsGeoelectricPotentialDetrendedBaselineSlope == NULL) {
        return TIICT_MEMORY;
    }

    return TIICT_OK;
}

int calculateFields(ProcessorState *state)
{
    ProcessorVariables_t *v = &state->vars16hz;

    // Init memory for fields
    int status = initFields(v);
    if (status != TIICT_OK)
        return status;

    // Calculate fields
    long long ind;

    float magVsat = sqrtf(v->vsatn[0]*v->vsatn[0] + v->vsate[0]*v->vsate[0] + v->vsatc[0] * v->vsatc[0]);

    for (int i = 0; i < v->nRecs; i++)
    {
        // Calculate xhat, yhat, zhat
        // xhat parallel to satellite velocity
        magVsat = sqrtf(v->vsatn[i]*v->vsatn[i] + v->vsate[i]*v->vsate[i] + v->vsatc[i] * v->vsatc[i]);
        ind = 3*i;
        v->xhat[ind + 0] = v->vsatn[i] / magVsat;
        v->xhat[ind + 1] = v->vsate[i] / magVsat;
        v->xhat[ind + 2] = v->vsatc[i] / magVsat;
        // yhat is xhat cross {0, 0, -1}
        v->yhat[ind + 0] = -v->xhat[ind + 1];
        v->yhat[ind + 1] = v->xhat[ind + 0];
        v->yhat[ind + 2] = 0.0;
        float magyhat = sqrtf(v->yhat[ind + 0] * v->yhat[ind + 0] + v->yhat[ind + 1] * v->yhat[ind + 1] + v->yhat[ind + 2] * v->yhat[ind + 2]);
        v->yhat[ind + 0] /= magyhat;
        v->yhat[ind + 1] /= magyhat;
        v->yhat[ind + 2] /= magyhat;
        // zhat is the cross-product of x into y
        v->zhat[ind+0] = v->xhat[ind+1] * v->yhat[ind+2] - v->xhat[ind+2] * v->yhat[ind+1];
        v->zhat[ind+1] = -1.0 * v->xhat[ind+0] * v->yhat[ind+2] + v->xhat[ind+2] * v->yhat[ind+0];
        v->zhat[ind+2] = v->xhat[ind+0] * v->yhat[ind+1] - v->xhat[ind+1] * v->yhat[ind+0];
        float magzhat = sqrtf(v->zhat[ind + 0] * v->zhat[ind + 0] + v->zhat[ind + 1] * v->zhat[ind + 1] + v->zhat[ind + 2] * v->zhat[ind + 2]);
        v->zhat[ind + 0] /= magzhat;
        v->zhat[ind + 1] /= magzhat;
        v->zhat[ind + 2] /= magzhat;

        // B field in cross-track frame, nT:
        v->bctx[ind] = v->bn[i] * v->xhat[ind + 0] + v->be[i] * v->xhat[ind + 1] + v->bc[i] * v->xhat[ind + 2];
        v->bcty[ind] = v->bn[i] * v->yhat[ind + 0] + v->be[i] * v->yhat[ind + 1] + v->bc[i] * v->yhat[ind + 2];
        v->bctz[ind] = v->bn[i] * v->zhat[ind + 0] + v->be[i] * v->zhat[ind + 1] + v->bc[i] * v->zhat[ind + 2];

        // E field from H sensor X, in cross-track frame, mV/m:
        v->ectxh[ind] = -1.0 * (v->myh[i] * v->bctz[ind] - v->myv[i] * v->bcty[ind]) / 1000000000.0 * 1000.0;
        v->ectyh[ind] = -1.0 * (-1.0 * v->mxh[i] * v->bctz[ind] + v->myv[i] * v->bctx[ind]) / 1000000000.0 * 1000.0;
        v->ectzh[ind] = -1.0 * (v->mxh[i] * v->bcty[ind] - v->myh[i] * v->bctx[ind]) / 1000000000.0 * 1000.0;

        // E field from V sensor X, in cross-track frame, mV/m:
        v->ectxv[ind] = -1.0 * (v->myh[i] * v->bctz[ind] - v->myv[i] * v->bcty[ind]) / 1000000000.0 * 1000.0;
        v->ectyv[ind] = -1.0 * (-1.0 * v->mxv[i] * v->bctz[ind] + v->myv[i] * v->bctx[ind]) / 1000000000.0 * 1000.0;
        v->ectzv[ind] = -1.0 * (v->mxv[i] * v->bcty[ind] - v->myh[i] * v->bctx[ind]) / 1000000000.0 * 1000.0;

    }

    state->setFlags = false;
    status = removeOffsetsAndSetFlags(state, geoelectricPotentialEstimator);
    if (status != TIICT_OK)
        return status;

    if (state->writeLogFiles) {
        fprintf(state->processingLogFile, "%sCalculated fields.\n", infoHeader);
        fflush(state->processingLogFile);
    }

    return TIICT_OK;

}

int integrateField(ProcessorState *state, float *sourceField, float scaleFactor, float *targetPotential, long startInd, long stopInd, bool vsDistance, bool positive, bool absoluteValue, bool removeMedianFromStart, float *medianDifference, float *firstSlope, float *lastSlope)
{
    int status = TIICT_OK;

    ProcessorVariables_t *v = &state->vars16hz;

    float previousField = 0.0;
    double previousTime = v->timestamp[startInd];
    double currentTime = v->timestamp[startInd];
    double deltaTime = 0.0;
    targetPotential[startInd] = 0.0;
    float magVsat = sqrtf(v->vsatn[startInd]*v->vsatn[startInd] + v->vsate[startInd]*v->vsate[startInd] + v->vsatc[startInd] * v->vsatc[startInd]);
    float previousMagVsat = magVsat;

    float signfactor = positive ? 1.0 : -1.0;
    float distanceFactor = 0.0;
    float nextContribution = 0.0;
    for (int i = startInd + 1; i < stopInd; i++)
    {
        currentTime = v->timestamp[i];
        deltaTime = (currentTime - previousTime) / 1000.0; // seconds
        magVsat = sqrtf(v->vsatn[i]*v->vsatn[i] + v->vsate[i]*v->vsate[i] + v->vsatc[i] * v->vsatc[i]);
        previousField = sourceField[(i-1)] * scaleFactor;
        nextContribution = signfactor * (absoluteValue ? fabsf(previousField) : previousField) * deltaTime;
        if (vsDistance)
            nextContribution *= magVsat;
        targetPotential[i] = targetPotential[i-1] + nextContribution;
        previousTime = currentTime;
        previousMagVsat = magVsat;
    }

    // Estimate median values at endpoints
    float firstMedian = 0.0;
    float lastMedian = 0.0;
    status = regionMetrics(state, state->bgws.beginIndex0, state->bgws.beginIndex1, targetPotential, &firstMedian, firstSlope);
    if (status != TIICT_OK)
        return status;
    status = regionMetrics(state, state->bgws.endIndex0, state->bgws.endIndex1, targetPotential, &lastMedian, lastSlope);
    if (status != TIICT_OK)
        return status;

    if (removeMedianFromStart)
    {
        for (int i = startInd; i < stopInd; i++)
        {
            targetPotential[i] -= firstMedian;
        }
    }
    if (medianDifference != NULL)
    {
        *medianDifference = lastMedian - firstMedian;
    }

    return TIICT_OK;
}

int regionMetrics(ProcessorState *state, long startInd, long stopInd, float *parameter, float *median, float *slope)
{
    ProcessorVariables_t *v = &state->vars16hz;

    double t0 = v->timestamp[startInd];
    size_t n = stopInd - startInd + 1;

    // GSL median sorts buffer, so use a working buffer
    gsl_matrix *distanceMatrix = gsl_matrix_alloc(n, 2);
    gsl_vector *values = gsl_vector_alloc(n);
    const gsl_multifit_robust_type * fitType = gsl_multifit_robust_bisquare;
    gsl_multifit_robust_workspace * gslFitWorkspace = gsl_multifit_robust_alloc(fitType, n, 2);
    gsl_vector *fitCoefficients = gsl_vector_alloc(2);
    gsl_matrix *cov = gsl_matrix_alloc(2, 2);
    gsl_multifit_robust_maxiter(GSL_FIT_MAXIMUM_ITERATIONS, gslFitWorkspace);

    // Load times into the model data buffer
    double fitDistance = 0.0;
    double vMag = 0.0;
    long ind = 0;
    for (int i = startInd; i < stopInd; i++)
    {
        ind = i - startInd;
        vMag = sqrtf(v->vsatn[i]*v->vsatn[i] + v->vsate[i]*v->vsate[i] + v->vsatc[i] * v->vsatc[i]);
        fitDistance = (v->timestamp[i] - t0)/1000. * (double)vMag;
        gsl_matrix_set(distanceMatrix, ind, 0, 1.0);
        gsl_matrix_set(distanceMatrix, ind, 1, fitDistance);
        gsl_vector_set(values, ind, (double)parameter[i]);
    }

    if (slope != NULL)
    {
        int gslStatus = gsl_multifit_robust(distanceMatrix, values, fitCoefficients, cov, gslFitWorkspace);
        if (gslStatus != GSL_SUCCESS)
        {
            double t1 = v->timestamp[stopInd];
            fprintf(stderr, "%sregionMetrics: unable to estimate linear fit between epochs %lf and %lf\n", infoHeader, t0, t1);
        }
        *slope = (float)gsl_vector_get(fitCoefficients, 1);
    }

    if (median != NULL)
    {
        *median = (float)gsl_stats_median(values->data, 1, n);
    }

    gsl_vector_free(fitCoefficients);
    gsl_matrix_free(distanceMatrix);
    gsl_vector_free(values);
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


bool downSampleHalfSecond(ProcessorState *state, long *index, long storageIndex, double t0, long maxIndex)
{
    long i = *index;
    uint8_t nSamples = 0;
    double timeBuf = 0.0;
    uint8_t orbitRegionBuf = 0;
    uint16_t flagBuf = 65535;
    uint32_t fitInfoBuf = 0;
    float theta, phi, x, y, z;
    bool downSampled = false;

    float floatBuf[100];
    memset(floatBuf, 0, 100 * sizeof *floatBuf);

    ProcessorVariables_t *v16 = &state->vars16hz;
    ProcessorVariables_t *v2 = &state->vars2hz;
    // TODO allocate memory for vars2hz

    int k = 0;

    while (((v16->timestamp[i]/1000. - t0) < 0.5) && (i <= maxIndex))
    {
        k = 0;
        timeBuf += v16->timestamp[i];
        // Handle lat and lon and mlt in cartesian coordinates
        // For spherical coordinates
        theta = M_PI / 2.0 - v16->latitude[i] * M_PI / 180.0;
        phi = v16->longitude[i] * M_PI / 180.0;
        floatBuf[k++] += cosf(phi) * sinf(theta);
        floatBuf[k++] += sinf(phi) * sinf(theta);
        floatBuf[k++] += cosf(theta);
        floatBuf[k++] += v16->radius[i];
        // For magnetic coordinates
        theta = M_PI / 2.0 - v16->qdlat[i] * M_PI / 180.0;
        phi = v16->mlt[i] / 24.0 * 2.0 * M_PI;
        floatBuf[k++] += cosf(phi) * sinf(theta);
        floatBuf[k++] += sinf(phi) * sinf(theta);
        floatBuf[k++] += cosf(theta);
        floatBuf[k++] += v16->mxh[i];
        floatBuf[k++] += v16->myh[i];
        floatBuf[k++] += v16->mxv[i];
        floatBuf[k++] += v16->myv[i];
        floatBuf[k++] += v16->vmcph[i];
        floatBuf[k++] += v16->vmcpv[i];
        floatBuf[k++] += v16->vbiash[i];
        floatBuf[k++] += v16->vbiasv[i];
        floatBuf[k++] += v16->vfp[i];
        floatBuf[k++] += v16->vsatx[i];
        floatBuf[k++] += v16->vsaty[i];
        floatBuf[k++] += v16->vsatz[i];
        floatBuf[k++] += v16->enhRaw[i];
        floatBuf[k++] += v16->envRaw[i];
        floatBuf[k++] += v16->enh[i];
        floatBuf[k++] += v16->env[i];
        floatBuf[k++] += v16->vixh[i];
        floatBuf[k++] += v16->vixherror[i];
        floatBuf[k++] += v16->vixv[i];
        floatBuf[k++] += v16->vixverror[i];
        floatBuf[k++] += v16->viy[i];
        floatBuf[k++] += v16->viyerror[i];
        floatBuf[k++] += v16->viz[i];
        floatBuf[k++] += v16->vizerror[i];
        floatBuf[k++] += v16->vsatn[i];
        floatBuf[k++] += v16->vsate[i];
        floatBuf[k++] += v16->vsatc[i];
        floatBuf[k++] += v16->ectxh[i]; // EH xyz
        floatBuf[k++] += v16->ectyh[i];
        floatBuf[k++] += v16->ectzh[i];
        floatBuf[k++] += v16->ectxv[i]; // EV xyz
        floatBuf[k++] += v16->ectyv[i];
        floatBuf[k++] += v16->ectzv[i];
        floatBuf[k++] += v16->bctx[i]; // Bxyz
        floatBuf[k++] += v16->bcty[i];
        floatBuf[k++] += v16->bctz[i];
        floatBuf[k++] += v16->bn[i];
        floatBuf[k++] += v16->be[i];
        floatBuf[k++] += v16->bc[i];
        floatBuf[k++] += v16->vicrx[i];
        floatBuf[k++] += v16->vicry[i];
        floatBuf[k++] += v16->vicrz[i];
        floatBuf[k++] += v16->geoelectricPotential[i];
        floatBuf[k++] += v16->geoelectricPotentialDifference[i];
        // Take the largest of each 8-sample interval
        if (floatBuf[k] < v16->maxAbsGeoelectricPotentialBaselineSlope[i]) {
            floatBuf[k++] = v16->maxAbsGeoelectricPotentialBaselineSlope[i];
        }
        floatBuf[k++] += v16->ehxAdjusted[i];
        // Latest value since this is constant over a pass
        floatBuf[k++] = v16->ehxAdjustmentParameter[i];
        // Take the largest of each 8-sample interval
        floatBuf[k++] += v16->geoelectricPotentialDetrended[i];
        if (floatBuf[k] < v16->maxAbsGeoelectricPotentialDetrendedBaselineSlope[i]) {
            floatBuf[k++] = v16->maxAbsGeoelectricPotentialDetrendedBaselineSlope[i];
        }
        // Take latest region in the half-second interval
        orbitRegionBuf = v16->orbitRegion[i];

        flagBuf &= v16->flags[i];
        fitInfoBuf |= v16->fitInfo[i];

        if (state->usePotentials) {
            floatBuf[k++] += v16->lpPhiScHighGain[i];
            floatBuf[k++] += v16->lpPhiScLowGain[i];
            floatBuf[k++] += v16->lpPhiSc[i];
            floatBuf[k++] += v16->potentials[i];
        }

        // Ignoring xhat, yhat, zhat


        nSamples++;
        i++;
    }
    if (nSamples == 8)
    {
        k = 0;
        // do the averaging and store result in original array
        v2->timestamp[storageIndex] = timeBuf / 8.0; // Average time
        x = floatBuf[k++] / 8.0;
        y = floatBuf[k++] / 8.0;
        z = floatBuf[k++] / 8.0;
        v2->latitude[storageIndex] = atan2f(z, sqrtf(x*x + y*y)) * 180.0 / M_PI; // Lat
        v2->longitude[storageIndex] = fmodf(atan2f(y, x) * 180.0 / M_PI, 360.0); // Lon
        v2->radius[storageIndex] = floatBuf[k++] / 8.0; // radius
        x = floatBuf[k++] / 8.0;
        y = floatBuf[k++] / 8.0;
        z = floatBuf[k++] / 8.0;
        v2->qdlat[storageIndex] = atan2f(z, sqrtf(x*x + y*y)) * 180.0 / M_PI; // QD Lat
        v2->mlt[storageIndex] = fmodf(atan2f(y, x) * 180.0 / M_PI + 360.0, 360.0) / 360. * 24.0; // MLT
        v2->mxh[storageIndex] = floatBuf[k++] / 8.0;
        v2->myh[storageIndex] = floatBuf[k++] / 8.0;
        v2->mxv[storageIndex] = floatBuf[k++] / 8.0;
        v2->myv[storageIndex] = floatBuf[k++] / 8.0;
        v2->vmcph[storageIndex] = floatBuf[k++] / 8.0;
        v2->vmcpv[storageIndex] = floatBuf[k++] / 8.0;
        v2->vbiash[storageIndex] = floatBuf[k++] / 8.0;
        v2->vbiasv[storageIndex] = floatBuf[k++] / 8.0;
        v2->vfp[storageIndex] = floatBuf[k++] / 8.0;
        v2->vsatx[storageIndex] = floatBuf[k++] / 8.0;
        v2->vsaty[storageIndex] = floatBuf[k++] / 8.0;
        v2->vsatz[storageIndex] = floatBuf[k++] / 8.0;
        v2->enhRaw[storageIndex] = floatBuf[k++] / 8.0;
        v2->envRaw[storageIndex] = floatBuf[k++] / 8.0;
        v2->enh[storageIndex] = floatBuf[k++] / 8.0;
        v2->env[storageIndex] = floatBuf[k++] / 8.0;
        v2->vixh[storageIndex] = floatBuf[k++] / 8.0;
        v2->vixherror[storageIndex] = floatBuf[k++] / 8.0 / sqrtf(8.0);
        v2->vixv[storageIndex] = floatBuf[k++] / 8.0;
        v2->vixverror[storageIndex] = floatBuf[k++] / 8.0 / sqrtf(8.0);
        v2->viy[storageIndex] = floatBuf[k++] / 8.0;
        v2->viyerror[storageIndex] = floatBuf[k++] / 8.0 / sqrtf(8.0);
        v2->viz[storageIndex] = floatBuf[k++] / 8.0;
        v2->vizerror[storageIndex] = floatBuf[k++] / 8.0 / sqrtf(8.0);
        v2->vsatn[storageIndex] = floatBuf[k++] / 8.0;
        v2->vsate[storageIndex] = floatBuf[k++] / 8.0;
        v2->vsatc[storageIndex] = floatBuf[k++] / 8.0;
        v2->ectxh[storageIndex] = floatBuf[k++] / 8.0;
        v2->ectyh[storageIndex] = floatBuf[k++] / 8.0;
        v2->ectzh[storageIndex] = floatBuf[k++] / 8.0;
        v2->ectxv[storageIndex] = floatBuf[k++] / 8.0;
        v2->ectyv[storageIndex] = floatBuf[k++] / 8.0;
        v2->ectzv[storageIndex] = floatBuf[k++] / 8.0;
        v2->bctx[storageIndex] = floatBuf[k++] / 8.0;
        v2->bcty[storageIndex] = floatBuf[k++] / 8.0;
        v2->bctz[storageIndex] = floatBuf[k++] / 8.0;
        v2->bn[storageIndex] = floatBuf[k++] / 8.0;
        v2->be[storageIndex] = floatBuf[k++] / 8.0;
        v2->bc[storageIndex] = floatBuf[k++] / 8.0;
        v2->vicrx[storageIndex] = floatBuf[k++] / 8.0;
        v2->vicry[storageIndex] = floatBuf[k++] / 8.0;
        v2->vicrz[storageIndex] = floatBuf[k++] / 8.0;
        v2->geoelectricPotential[storageIndex] = floatBuf[k++] / 8.0;
        v2->geoelectricPotentialDifference[storageIndex] = floatBuf[k++] / 8.0;
        v2->maxAbsGeoelectricPotentialBaselineSlope[storageIndex] = floatBuf[k++];
        v2->ehxAdjusted[storageIndex] = floatBuf[k++] / 8.0;
        v2->ehxAdjustmentParameter[storageIndex] = floatBuf[k++];
        v2->geoelectricPotentialDetrended[storageIndex] = floatBuf[k++] / 8.0;
        v2->maxAbsGeoelectricPotentialDetrendedBaselineSlope[storageIndex] = floatBuf[k++];
        v2->orbitRegion[storageIndex] = orbitRegionBuf;
        // Flags set to 0 at 16 Hz based on magnitude of flow,
        // are not reset at 2 Hz, to ensure integrity of 2 Hz measurements
        // One can review 16 Hz measurements to examine details of flow where even a
        // single sample of the eight has a magnitude greater than 8 km/s
        v2->flags[storageIndex] = flagBuf;
        v2->fitInfo[storageIndex] = fitInfoBuf;

        if (state->usePotentials) {
            v2->lpPhiScHighGain[storageIndex] = floatBuf[k++] / 8.0;
            v2->lpPhiScLowGain[storageIndex] = floatBuf[k++] / 8.0;
            v2->lpPhiSc[storageIndex] = floatBuf[k++] / 8.0;
            v2->potentials[storageIndex] = floatBuf[k++] / 8.0;
        }


        downSampled = true;

    }

    *index = i;
    return downSampled;

}

int runProcessor(int argc, char *argv[], ProcessorState **result)
{
    int status = TIICT_OK;
    ProcessorState *state = NULL;
    if (result != NULL && *result != NULL) {
        // re-use the provided state
        state = *result;
    }
    else {
        state = initState(argc, argv);
    }

    status = initProcessor(state);
    if (status != TIICT_OK) {
        goto cleanup;
    }

    status = loadCalData(state);
    if (status != TIICT_OK) {
        goto cleanup;
    }

    status = loadLpCalData(state);
    if (status != TIICT_OK) {
        goto cleanup;
    }

    status = calibrateFlows(state);
    if (status != TIICT_OK) {
        goto cleanup;
    }

    status = calculateFields(state);
    if (status != TIICT_OK) {
        goto cleanup;
    }

    if (state->export2Hz || state->export16Hz) {
        status = exportCdfs(state);
        if (status != TIICT_OK) {
            goto cleanup;
        }
    }

    if (state->visualizeResults) {
        state->vars = &state->vars16hz;
        status = visualizeResults(state);
        if (status != TIICT_OK) {
            goto cleanup;
        }
    }

cleanup:

    // Close files
    if (state->writeLogFiles) {
        closeFiles(state);
    }

    if (result != NULL) {
        // Return with results if requested
        *result = state;
        return status;
    }
    else {
        // free memory
        status = shutdown(state);
        free(state);
    }
    printf("9\n");

    return status;

}

ProcessorState *initState(int argc, char **argv)
{
    // Allocate memory for the processor state
    ProcessorState *state = malloc(sizeof *state);
    if (state == NULL) {
        fprintf(stderr, "Unable to allocate memory for processor state.\n");
        return NULL;
    }
    memset(state, 0, sizeof *state);
    Arguments *args = &state->args;
    args->argc = argc;
    args->argv = argv;

    // Check arguments and abort if not right
    int status = parseArguments(state);
    if (status != TIICT_OK) {
        return NULL;
    }

    return state;
}

int initProcessor(ProcessorState *state)
{
    int status = TIICT_OK;
    Arguments *args = &state->args;

    // Prefix for messages
    initHeader(state);

    // Offset model parameters
    offset_model_fit_arguments f[4] = {
        {0, "Northern ascending", 44.0, 50.0, 50.0, 44.0},
        {1, "Equatorial descending", 44.0, 38.0, -38.0, -44.0},
        {2, "Southern descending", -44.0, -50.0, -50.0, -44.0},
        {3, "Equatorial ascending", -44.0, -38.0, 38.0, 44.0},
    };
    memcpy(state->fitargs, f, 4 * sizeof(offset_model_fit_arguments));

    // Normalize date (where 32 Jan is 1 Feb)
    double epoch = computeEPOCH(args->year, args->month, args->day, 0, 0, 0, 0.0);
    long year, month, day, h, m, s, ms;
    EPOCHbreakdown(epoch, &year, &month, &day, &h, &m, &s, &ms);
    args->year = year;
    args->month = month;
    args->day = day;

    if (state->writeLogFiles) {
        status = initLogFiles(state);
        if (status != TIICT_OK) {
            return status;
        }

        fprintf(state->processingLogFile, "\n%s-------------------------------------------------\n", infoHeader);
        fprintf(state->processingLogFile, "%sVersion 0401 20241110\n", infoHeader);
        fprintf(state->processingLogFile, "%sProcessing date: %s", infoHeader, state->processingDateString);

        // Print command line
        fprintf(state->processingLogFile, "%sCalled as '", infoHeader);
        for (int i = 0; i < args->argc; i++) {
            fprintf(state->processingLogFile, "%s", args->argv[i]);
            if (i < args->argc - 1) {
                fprintf(state->processingLogFile, " ");
            }
        }
        fprintf(state->processingLogFile, "'\n");
    }

    // Confirm requested date has records. Abort otherwise.
    status = checkCalDataAvailability(state);
    if (status != TIICT_OK)
        return status;

    // Ensure export directories exist, or abort.
    status = initDirectories(state);
    if (status != TIICT_OK)
        return status;

    // Turn off GSL failsafe error handler. We typically check the GSL return codes.
    gsl_set_error_handler_off();
    return TIICT_OK;
}

int parseArguments(ProcessorState *state)
{
    Arguments *args = &state->args;
    int argc = args->argc;
    char **argv = args->argv;

    // E-of-R method for along-track drift
    state->useEofR = true;
    // LP estimates of satellite potential are disabled by default
    state->usePotentials = false;
    state->vars16hz.lpPotentialSource = LP_POTENTIAL_NONE;
    state->vars2hz.lpPotentialSource = LP_POTENTIAL_NONE;
    state->frames = NULL;

    state->export2Hz = true;
    state->export16Hz = true;
    state->exportZip = true;
    state->exportVideo = false;

    state->videoOutputDir = ".";
    state->videoFilename[0] = '\0';
    state->printVideoFilename = false;
    state->visualizeResults = false;
    state->frameWidth = 960;
    state->frameHeight = 540;
    state->framesPerSecond = 15;

    state->plotCommand = "QDLat,-90,90,1" ";PhiSc,-5,0,1" ";Vixh,-4,4,0.001" ";Vixv,-4,4,0.001" ";Viy,-2,2,0.001" ";Viz,-2,2,0.001";
    state->defaultPlotHeight = DEFAULT_PLOT_HEIGHT;
    state->maxPlotsPerScreen = 0;
    state->writeLogFiles = true;

    // Default automatically to first and last times for video export
    double epoch = computeEPOCH(state->args.year, state->args.month, state->args.day, 0, 0, 0, 0);
    state->plotT0 = epoch;
    state->plotT1 = epoch + 86400 * 1000.0;

    state->nOptions = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp("--no-log-files", argv[i]) == 0) {
            state->nOptions++;
            state->writeLogFiles = false;
        }
        else if (strcmp("--no-16hz-export", argv[i]) == 0) {
            state->nOptions++;
            state->export16Hz = false;
        }
        else if (strcmp("--no-2hz-export", argv[i]) == 0) {
            state->nOptions++;
            state->export2Hz = false;
        }
        else if (strcmp("--no-zip-export", argv[i]) == 0) {
            state->nOptions++;
            state->exportZip = false;
        }
        else if (strcmp("--export-video", argv[i]) == 0) {
            state->nOptions++;
            state->exportVideo = true;
        }
        else if (strcmp("--no-eofr", argv[i]) == 0) {
            state->nOptions++;
            state->useEofR = false;
        }
        else if (strcmp("--visualize", argv[i]) == 0) {
            state->nOptions++;
            state->visualizeResults = true;
        }
        else if (strncmp("--frame-width=", argv[i], 14) == 0) {
            state->nOptions++;
            if (strlen(argv[i]) < 15) {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                return TIICT_ARGS_BAD;
            }
            state->frameWidth = atoi(argv[i] + 14);
        }
        else if (strncmp("--frame-height=", argv[i], 15) == 0) {
            state->nOptions++;
            if (strlen(argv[i]) < 16) {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                return TIICT_ARGS_BAD;
            }
            state->frameHeight = atoi(argv[i] + 15);
        }
        else if (strncmp("--fps=", argv[i], 6) == 0) {
            state->nOptions++;
            if (strlen(argv[i]) < 7) {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                return TIICT_ARGS_BAD;
            }
            state->framesPerSecond = atoi(argv[i] + 6);
        }
        else if (strncmp("--plots=", argv[i], 8) == 0) {
            state->nOptions++;
            if (strlen(argv[i]) < 9) {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                return TIICT_ARGS_BAD;
            }
            state->plotCommand = argv[i] + 8;
        }
        else if (strncmp("--plot-height=", argv[i], 14) == 0) {
            state->nOptions++;
            if (strlen(argv[i]) < 15) {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                return TIICT_ARGS_BAD;
            }
            state->defaultPlotHeight = atoi(argv[i] + 14);
        }
        else if (strncmp("--max-screen-plots=", argv[i], 19) == 0) {
            state->nOptions++;
            if (strlen(argv[i]) < 20) {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                return TIICT_ARGS_BAD;
            }
            state->maxPlotsPerScreen = atoi(argv[i] + 19);
        }
        else if (strncmp("--video-dir=", argv[i], 12) == 0) {
            state->nOptions++;
            if (strlen(argv[i]) < 13) {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                return TIICT_ARGS_BAD;
            }
            state->videoOutputDir = argv[i] + 12;
        }
        else if (strncmp("--t0=", argv[i], 5) == 0) {
            state->nOptions++;
            if (strlen(argv[i]) < 6) {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                return TIICT_ARGS_BAD;
            }

            state->plotT0 = parseEPOCH4(argv[i] + 5);
        }
        else if (strncmp("--t1=", argv[i], 5) == 0) {
            state->nOptions++;
            if (strlen(argv[i]) < 6) {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                return TIICT_ARGS_BAD;
            }

            state->plotT1 = parseEPOCH4(argv[i] + 5);
        }
        else if (strncmp("--video-filename=", argv[i], 17) == 0) {
            state->nOptions++;
            if (strlen(argv[i]) < 18) {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                return TIICT_ARGS_BAD;
            }
            snprintf(state->videoFilename, FILENAME_MAX, "%s", argv[i] + 17);
        }
        else if (strcmp("--print-video-filename", argv[i]) == 0) {
            state->nOptions++;
            state->printVideoFilename = true;
        }
        else if (strncmp("--floating-potential-source=", argv[i], 28) == 0) {
            state->nOptions++;
            if (strlen(argv[i]) != 29) {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                return TIICT_ARGS_BAD;
            }
            char source = (argv[i] + 28)[0];
            switch (source) {
                case 'H':
                case 'h':
                    state->vars16hz.lpPotentialSource = LP_POTENTIAL_HIGHGAIN;
                    state->usePotentials = true;
                    break;
                case 'L':
                case 'l':
                    state->vars16hz.lpPotentialSource = LP_POTENTIAL_LOWGAIN;
                    state->usePotentials = true;
                    break;
                case 'U':
                case 'u':
                    state->vars16hz.lpPotentialSource = LP_POTENTIAL_U_SC;
                    state->usePotentials = true;
                    break;
                case 'N':
                case 'n':
                    state->vars16hz.lpPotentialSource = LP_POTENTIAL_NONE;
                    state->usePotentials = false;
                    break;
                default:
                    state->vars16hz.lpPotentialSource = LP_POTENTIAL_UNKNOWN;
                    state->usePotentials = false;
                    break;
            }
            state->vars2hz.lpPotentialSource = state->vars16hz.lpPotentialSource;
        }
        else if (strcmp("--help", argv[i]) == 0) {
            cmdUsage(argv[0]);
            return TIICT_ARGS_HELP;
        }
        else if (strcmp(argv[i], "--about") == 0) {
            fprintf(stdout, "tiict - TII Cross-track ion drift processor, version %s.\n", SOFTWARE_VERSION);
            fprintf(stdout, "Copyright (C) 2024  Johnathan K Burchill\n");
            fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
            fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
            fprintf(stdout, "under the terms of the GNU General Public License.\n");
            return TIICT_ARGS_ABOUT;
        }
        else if (strncmp(argv[i], "--", 2) == 0) {
            cmdUsage(argv[0]);
            return TIICT_ARGS_BAD;
        }
    }

    if (argc - state->nOptions != 10) {
        cmdUsage(argv[0]);
        return TIICT_ARGS_BAD;
    }

    if (!state->export16Hz && !state->export2Hz && !state->visualizeResults) {
        fprintf(stdout, "No record export or visualization has been requested. Not processing.\n");
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

void cmdUsage(char *name)
{
    fprintf(stdout, "usage: %s satLetter year month day calversionString exportVersionString calDir lpDir exportDir\n", name);
    fprintf(stdout, "options:\n");
    fprintf(stdout, "%40s - %s\n", "--no-log-files", "do not write log files");
    fprintf(stdout, "%40s - %s\n", "--no-16hz-export", "do not export 16 Hz dataset");
    fprintf(stdout, "%40s - %s\n", "--no-2hz-export", "do not export 2 Hz dataset");
    fprintf(stdout, "%40s - %s\n", "--no-zip-export", "do not export zip archive");
    fprintf(stdout, "%40s - %s\n", "--export-video", "export video of results");
    fprintf(stdout, "%40s - %s\n", "--visualize", "generate visualizations of results");
    fprintf(stdout, "%40s - %s\n", "--floating-potential-source=<source>", "satellite floating potential source");
    fprintf(stdout, "%40s %s\n", "", "'N' or 'n': none (default)");
    fprintf(stdout, "%40s %s\n", "", "'U' or 'u': EXTD U_SC (blended)");
    fprintf(stdout, "%40s %s\n", "", "'H' or 'h': EXTD high-gain probe");
    fprintf(stdout, "%40s %s\n", "", "'L' or 'l': EXTD low-gain probe");
    fprintf(stdout, "%40s - %s\n", "--plot-command=<cmd>", "plot instructions");
    fprintf(stdout, "%40s %s\n", "", "Semicolon-separated commands of the form <param>,<ymin>,<ymax>,<yscale>[,<plotheight>]");
    fprintf(stdout, "%40s %s\n", "", "params:");
    fprintf(stdout, "%40s   %s\n", "", "QDLat, MLT, PhiSc, Vixh, Vixv, Viy, Viz");
    fprintf(stdout, "%40s - %s\n", "--plot-height=<value>", "set the default plot height in pixels");
    fprintf(stdout, "%40s - %s\n", "--max-screen-plots=<value>", "limit to <value> plots per screen. Defaule: 0 (automatic)");
    fprintf(stdout, "%40s - %s\n", "--video-dir", "video output directory; default: '.'");
    fprintf(stdout, "%40s - %s\n", "--video-filename", "video filename; default: 'results.mp4'");
    fprintf(stdout, "%40s - %s\n", "--print-video-filename", "print the video filename");
    fprintf(stdout, "%40s - %s\n", "--no-eofr", "use legacy method for estimating along-track drift");
    fprintf(stdout, "%40s - %s\n", "--about", "print copyright and license information");
    fprintf(stdout, "%40s - %s\n", "--help", "print this message");

    return;
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

    // Store the processing start time string
    asctime_r(timeParts, state->processingDateString);

    return;
}

int checkResult(int status, ProcessorState *state)
{
    if (status != TIICT_OK)
    {
        if (state->processingLogFile != NULL) {
            if (state->writeLogFiles) {
                fprintf(state->processingLogFile, "Error processing file. status = %d\n", status);
            }
        }
        state->returnStatus = status;
        status = shutdown(state);
    }

    return status;

}

void closeFiles(ProcessorState *state)
{
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

    return;
}

int shutdown(ProcessorState *state)
{
    if (state == NULL)
        return(TIICT_ARGS_BAD);

    // Free the memory
    if (state->vars16hz.memoryAllocated) {
        freeVariables(&state->vars16hz);
        state->vars16hz.memoryAllocated = 0;
        state->vars16hz.nRecs = 0;
        state->vars16hz.nLpRecs = 0;
    }
    if (state->vars2hz.memoryAllocated) {
        freeVariables(&state->vars2hz);
        state->vars2hz.memoryAllocated = 0;
        state->vars2hz.nRecs = 0;
        state->vars2hz.nLpRecs = 0;
    }

    for (int i = 0; i < state->nVideoFrames; i++) {
        free(state->frames[i].pixels);
    }
    free(state->frames);
    state->frames = NULL;
    state->nVideoFrames = 0;

    return TIICT_OK;
}

int backgroundRemoval(ProcessorState *state)
{
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

    ProcessorVariables_t *v = &state->vars16hz;

    float *param = NULL;
    float *paramErrors = NULL;
    // Load values into model data buffer once each for HX, HY, VX, VY

    for (uint8_t k = 0; k < 4; k++)
    {
        if (state->measurementType == VELOCITY_MEASUREMENT) {
            switch (k) {
                case 0:
                    param = v->vixh;
                    paramErrors = v->vixherror;
                    break;
                case 1:
                    param = v->viy;
                    paramErrors = v->viyerror;
                    break;
                case 2:
                    param = v->vixv;
                    paramErrors = v->vixverror;
                    break;
                case 3:
                    param = v->viz;
                    paramErrors = v->vizerror;
                    break;
                default:
                    break;
            }
        }
        else {
            // TODO estimate energy errors
            switch (k) {
                case 0:
                    param = v->enhRaw;
                    paramErrors = NULL;
                    break;
                case 1:
                    param = v->envRaw;
                    paramErrors = NULL;
                    break;
                case 2:
                    param = v->enh;
                    paramErrors = NULL;
                    break;
                case 3:
                    param = v->env;
                    paramErrors = NULL;
                    break;
                default:
                    break;
            }
        }

        uint8_t flagIndex = k; // Defined flags as bit 0 -> HX, bit 1 -> VX, bit 2-> HY and bit 3-> VY. Data are stored in memory differently.
        if (flagIndex == 1)
            flagIndex = 2;
        else if (flagIndex == 2)
            flagIndex = 1;

        ws->modelDataIndex = 0;
        for (int i = ws->beginIndex0; i < ws->beginIndex1; i++)
        {
            gsl_vector_set(ws->model1Values, ws->modelDataIndex, param[i]);
            gsl_vector_set(ws->modelValues, ws->modelDataIndex++, param[i]);
        }
        ws->modelDataMidPoint = ws->modelDataIndex;
        for (int i = ws->endIndex0; i < ws->endIndex1; i++)
        {
            gsl_vector_set(ws->model2Values, ws->modelDataIndex - ws->modelDataMidPoint, param[i]);
            gsl_vector_set(ws->modelValues, ws->modelDataIndex++, param[i]);
        }
        // Robust linear model fit and removal
        gslFitWorkspace = gsl_multifit_robust_alloc(fitType, ws->numModelPoints, ws->fitDegree);
        gslStatus = gsl_multifit_robust_maxiter(GSL_FIT_MAXIMUM_ITERATIONS, gslFitWorkspace);
        if (gslStatus && state->writeLogFiles)
        {
            fprintf(state->processingLogFile, "%sCould not set maximum GSL iterations.\n", infoHeader);
        }
        gslStatus = gsl_multifit_robust(ws->modelTimesMatrix, ws->modelValues, ws->fitCoefficients, cov, gslFitWorkspace);
        if (gslStatus)
        {
            if (state->writeLogFiles) {
                toEncodeEPOCH(ws->tregion11, 0, ws->startString);
                toEncodeEPOCH(ws->tregion22, 0, ws->stopString);
                fprintf(state->processingLogFile, "%s<GSL Fit Error: %s> for fit region from %s to %s spanning latitudes %.0f to %.0f.\n", infoHeader, gsl_strerror(gslStatus), ws->startString, ws->stopString, fitargs->lat1, fitargs->lat4);
                // Print "-9999999999.GSLERRORNUMBER" for each of the nine fit parameters
                fprintf(state->fitFile, " -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d", gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus);
            }
            if (state->setFlags)
            {
                for (int i = ws->beginIndex0; i < ws->endIndex1; i++)
                {
                    // Got a complete region, but had a fit error
                    v->fitInfo[i] &= ~(FITINFO_INCOMPLETE_REGION << (flagIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
                    v->fitInfo[i] |= (FITINFO_GSL_FIT_ERROR << (flagIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
                }
            }
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
            if (state->writeLogFiles) {
                fprintf(state->fitFile, " %f %f %f %f %f %f %f %f %f", c0, c1, stats.adj_Rsq, stats.rmse, median1, median2, mad, mad1, mad2);
            }
            // Remove the offsets and assign flags for this region
            for (int i = ws->beginIndex0; i < ws->endIndex1; i++)
            {
                // remove offset
                param[i] -= (((v->timestamp[i] - ws->epoch0)/1000.0) * c1 + c0);
                driftValue = param[i];
                // Assign error estimate
                // TODO: use interpolated MADs derived start and end of pass?
                if (paramErrors != NULL) {
                    paramErrors[i] = mad;
                }
                // Set offset-removed flag (most significant bit), and complete region found
                // Having a complete region can be determined from offset removed and gsl error flags, but
                // the extra bit will make it logically straightforward to find incomplete regions.
                // Toggle off the "offset not removed" and "incomplete region" flag bits
                if (state->setFlags)
                {
                    v->fitInfo[i] &= ~((FITINFO_OFFSET_NOT_REMOVED | FITINFO_INCOMPLETE_REGION) << (flagIndex * MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT));
                    // Update quality and calibration flags based on thresholds
                    updateDataQualityFlags(state->args.satellite, flagIndex, fitargs->regionNumber, driftValue, mad, i, v->flags, v->fitInfo);
                }
            }
        }
        gsl_multifit_robust_free(gslFitWorkspace);
    }
    gsl_matrix_free(cov);

    return TIICT_OK;

}

void geoelectricPotentialBackgroundRemoval(ProcessorState *state)
{

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

    ws->modelDataIndex = 0;
    vmag1 = 0.0;
    vmag2 = 0.0;

    ProcessorVariables_t *v = &state->vars16hz;
    for (int i = ws->beginIndex0; i < ws->beginIndex1; i++)
    {
        gsl_vector_set(ws->model1Values, ws->modelDataIndex, v->geoelectricPotential[i]);
        gsl_vector_set(ws->modelValues, ws->modelDataIndex++, v->geoelectricPotential[i]);
        vmag1 += sqrtf(v->vsatn[i]*v->vsatn[i] + v->vsate[i]*v->vsate[i] + v->vsatc[i] * v->vsatc[i]);
    }
    if (ws->numModel1Points > 0)
        vmag1 /= (float)ws->numModel1Points;
    else
        vmag1 = 0.0;
    ws->modelDataMidPoint = ws->modelDataIndex;
    for (int i = ws->endIndex0; i < ws->endIndex1; i++)
    {
        gsl_vector_set(ws->model2Values, ws->modelDataIndex - ws->modelDataMidPoint, v->geoelectricPotential[i]);
        gsl_vector_set(ws->modelValues, ws->modelDataIndex++, v->geoelectricPotential[i]);
        vmag2 += sqrtf(v->vsatn[i]*v->vsatn[i] + v->vsate[i]*v->vsate[i] + v->vsatc[i] * v->vsatc[i]);
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
    if (gslStatus )
    {
        if (state->writeLogFiles) {
            toEncodeEPOCH(ws->tregion11, 0, ws->startString);
            toEncodeEPOCH(ws->tregion22, 0, ws->stopString);
            fprintf(state->processingLogFile, "%s<GSL Fit Error: %s> for fit region from %s to %s spanning latitudes %.0f to %.0f.\n", infoHeader, gsl_strerror(gslStatus), ws->startString, ws->stopString, fitargs->lat1, fitargs->lat4);
            // Print "-9999999999.GSLERRORNUMBER" for each of the nine fit parameters
            fprintf(state->fitFile, " -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d -9999999999.%d", gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus, gslStatus);
        }
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
        for (int i = ws->beginIndex0; i < ws->endIndex1; i++)
        {
            // Define median1 of first region as zero potential
            v->maxAbsGeoelectricPotentialBaselineSlope[i] = slope1 > slope2 ? slope1 : slope2;
            // linear detrend
            v->geoelectricPotentialDetrended[i] = v->geoelectricPotential[i] - (((v->timestamp[i] - ws->epoch0)/1000.0) * c1 + c0);
            v->geoelectricPotential[i] -= median1;
            // Remove median of start region
            v->geoelectricPotentialDifference[i] = median2 - median1;
        }

        // Now estimate detrended slopes
        for (int i = ws->beginIndex0; i < ws->endIndex0; i++)
        {
            gsl_vector_set(ws->model1Values, i-ws->beginIndex0, v->geoelectricPotentialDetrended[i]);
        }
        for (int i = ws->beginIndex1; i < ws->endIndex1; i++)
        {
            gsl_vector_set(ws->model2Values, i-ws->beginIndex1, v->geoelectricPotentialDetrended[i]);
        }
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
        for (int i = ws->beginIndex0; i < ws->endIndex1; i++)
        {
            v->maxAbsGeoelectricPotentialDetrendedBaselineSlope[i] = slope1 > slope2 ? slope1 : slope2;
        }
    }
    gsl_multifit_robust_free(gslFitWorkspace);
    gsl_multifit_robust_free(gslFitWorkspace1);
    gsl_multifit_robust_free(gslFitWorkspace2);
    gsl_matrix_free(cov);

    return;

}

int geoelectricPotentialEstimator(ProcessorState *state)
{
    int status = TIICT_OK;

    ProcessorVariables_t *v = &state->vars16hz;

    long i = state->bgws.beginIndex0;
    double previousTime = v->timestamp[i];
    double currentTime = v->timestamp[i];
    double deltaTime = 0.0;
    v->geoelectricPotential[i] = 0.0;

    float deltaPhi1 =0.0;
    float deltaPhi2 =0.0;
    float deltaPhi3 =0.0;
    size_t bInd0 = state->bgws.beginIndex0;
    size_t eInd1 = state->bgws.endIndex1;

    // Set the region number
    // 0 northern ascending, 1 equatorial descending
    // 2 southern descending, 3 equatorial ascending
    offset_model_fit_arguments *fitargs = &state->fitargs[state->interval];
    for (i = bInd0; i < eInd1; i++)
    {
        v->orbitRegion[i] = fitargs->regionNumber;
    }

    float s1 = 0.0;
    float s2 = 0.0;
    float s1p = 0.0;
    float s2p = 0.0;
    status = integrateField(state, v->ectxh, 1/1000.0, v->geoelectricPotential, bInd0, eInd1, true, false, false, true, &deltaPhi1, &s1, &s2);
    if (status != TIICT_OK)
        return status;
    status = integrateField(state, v->ectxh, 1/1000.0, v->geoelectricPotentialDetrended, bInd0, eInd1, true, true, false, true, &deltaPhi2, NULL, NULL);

    if (status != TIICT_OK)
        return status;
    status = integrateField(state, v->ectxh, 1/1000.0, v->ehxAdjusted, bInd0, eInd1, true, true, true, true, &deltaPhi3, NULL, NULL);
    if (status != TIICT_OK)
        return status;

    // Zhu et al. (2020) DMSP rescaling factor
    float c = - deltaPhi2 / deltaPhi3;
    for (i = bInd0; i < eInd1; i++)
    {
        v->ehxAdjusted[i] = v->ectxh[i] + c * fabsf(v->ectxh[i]);
        v->ehxAdjustmentParameter[i] = c;
        v->geoelectricPotentialDifference[i] = deltaPhi1;
    }

    // Calculated adjusted geopotential
    status = integrateField(state, v->ehxAdjusted, 1/1000.0, v->geoelectricPotentialDetrended, bInd0, eInd1, true, false, false, true, NULL, &s1p, &s2p);
    if (status != TIICT_OK)
        return status;

    s1 = fabsf(s1);
    s2 = fabsf(s2);
    s1p = fabsf(s1p);
    s2p = fabsf(s2p);
    for (i = bInd0; i < eInd1; i++)
    {
        v->maxAbsGeoelectricPotentialBaselineSlope[i] = 1000.0 * (s1 > s2 ? s1 : s2);
        v->maxAbsGeoelectricPotentialDetrendedBaselineSlope[i] = 1000.0 * (s1p > s2p ? s1p : s2p);
    }

    return status;
}

