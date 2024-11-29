#include "visualize.h"
#include "loadData.h"
#include "state.h"
#include "errors.h"

#include <stdio.h>
#include <tii/utility.h>
#include <tiigraphics/draw.h>
#include <tiigraphics/tiigraphics.h>
#include <tiigraphics/video.h>
#include <tiigraphics/colors.h>
#include <tiigraphics/fonts.h>

#include <stdlib.h>
#include <string.h>

int visualizeResults(ProcessorState *state)
{
    if (state == NULL) {
        return TIICT_ARGS_BAD;
    }
    int frameCounter = 0;
    static int timesReached = 0;

    int plotdy = 25;
    int plotX0 = 100;
    int topMargin = 10;

    Image image = {0};
    if (allocImage(&image, state->frameWidth, state->frameHeight, 1) != DRAW_OK)
    {
        printf("Could not allocate memory for image.\n");
        goto cleanup;
    }
    // Set background
    memset(image.pixels, BACKGROUND_COLOR, image.numberOfBytes);

    // Set by caller for 16 Hz or 2 Hz
    ProcessorVariables_t *v = state->vars;

    if (v->timestamp == NULL) {
        addAnnotations(state, &image, topMargin);
        annotate("No data", 24, state->frameWidth/2 - fontwidth(24)*strlen("No data")/2, state->frameHeight/2 + fontheight(24)/2, &image);
        storeImage(state, &image);
        goto cleanup;
    }
    int firstIndex = 0;
    int lastIndex = 0;

    // Find first index
    while (firstIndex < v->nRecs - 1 && v->timestamp[firstIndex] < state->plotT0) {
        firstIndex++;
    }
    lastIndex = firstIndex;
    while (lastIndex < v->nRecs - 1 && v->timestamp[lastIndex] < state->plotT1) {
        lastIndex++;
    }

    if (strlen(state->videoFilename) == 0) {
        char format[EPOCHx_FORMAT_MAX], t0String[EPOCHx_STRING_MAX], t1String[EPOCHx_STRING_MAX];
        snprintf(format, EPOCHx_FORMAT_MAX, "<year><mm.02><dom.02>T<hour><min><sec>");
        encodeEPOCHx(state->plotT0, format, t0String);
        encodeEPOCHx(state->plotT1, format, t1String);
        snprintf(state->videoFilename, FILENAME_MAX, "%s/Swarm%s_TIICT_%s_%s_%s.mp4", state->videoOutputDir, state->args.satellite, t0String, t1String, state->args.exportVersion);
    }

    // Draw PA and measles time series
    int plotWidth = 700;
    int plotHeight0 = state->defaultPlotHeight;
    int plotHeight = plotHeight0;
    int ox = 100;
    int oy = 140;
    int dotSize = 2;
    char xlabel[255];
    sprintf(xlabel, "%s", "Hours from start of file");

    dotSize = 2; // full day
    sprintf(xlabel, "%s", "UT hours");

    // Plots

    // parse plot command
    // from strsep docs
    char *plotOptions, *string, *tofree;
    tofree = string = strdup(state->plotCommand);
    int maxPlots = state->maxPlotsPerScreen;
    int nPlots = 0;
    int nParams = 0;
    float yr0 = 0.0;
    float yr1 = 0.0;
    float valueScale = 1.0;
    float valueOffset = 0.0;
    char yr0str[255];
    char yr1str[255];
    int plotY0 = plotHeight0 + plotdy + topMargin;
    int plotYOffset = plotY0;
    int yText = topMargin;
    int plotsMade = 0;
    char *xLabel = "";
    float *parameter = NULL;
    char *parameterLabel = "";
    bool gotParameter = false;
    int stride = 1;
    int tupleLength = 0;
    int tupleIndex = 0;
    // Count number of requested plots
    while ((plotOptions = strsep(&string, ";")) != NULL) {
        nPlots++;
    }
    free(tofree);

    int fontSize = 12;

    // Plot requested plots

    if (state->exportVideo) {
        int status = initVideo(state->videoFilename, state->frameWidth, state->frameHeight, state->framesPerSecond);
        if (status < 0)
        {
            fprintf(stderr, "Problem intializing video: got status %d.\n", status);
            goto cleanup;
        }
    }

    tofree = string = strdup(state->plotCommand);
    while ((plotOptions = strsep(&string, ";")) != NULL && plotsMade < nPlots) {
        gotParameter = true;
        tupleLength = 1;
        tupleIndex = 0;
        // Parse plot parameters
        char **ap, *params[10];
        nParams = 0;
        for (ap = params; (*ap = strsep(&plotOptions, ",")) != NULL;) {
            if (**ap != '\0') {
                nParams++;
                if (++ap >= &params[10]) {
                    break;
                }
            }
        }
        if (nParams < 4) {
            // Expected format: ...;<param>,<ymin>,<ymax>,<yscale>[,<plotheight>];...
            // Skip empty plots, otherwise warn if not enough parameters
            if (nParams > 0) {
                fprintf(stderr, "Warning: plot command was %s: expected at least %d plot parameters, skipping this plot\n", plotOptions, nParams);
            }
            gotParameter = false;
            continue;
        }

        yr0 = atof(params[1]) * state->yScaleFactor;
        yr1 = atof(params[2]) * state->yScaleFactor;
        valueScale = atof(params[3]);
        valueOffset = 0.0;
        int oldPlotHeight = plotHeight;
        if (nParams == 5) {
            plotHeight = atoi(params[4]);
        }
        else {
            plotHeight = plotHeight0;
        }
        plotYOffset += (plotHeight - oldPlotHeight);

        if (strcmp("QDLat", params[0]) == 0) {
            parameter = v->qdlat;
            parameterLabel = "QD Lat";
        } else if (strcmp("MLT", params[0]) == 0) {
            parameter = v->mlt;
            parameterLabel = "MLT";
        } else if (strcmp("Lat", params[0]) == 0) {
            parameter = v->latitude;
            parameterLabel = "Lat";
        } else if (strcmp("Lon", params[0]) == 0) {
            parameter = v->longitude;
            parameterLabel = "Lon";
        } else if (strcmp("Radius", params[0]) == 0) {
            parameter = v->radius;
            parameterLabel = "Radius";
        } else if (strcmp("Alt", params[0]) == 0) {
            parameter = v->radius;
            parameterLabel = "Alt";
            valueOffset = 6372000.0;
        } else if (strcmp("mxh", params[0]) == 0) {
            parameter = v->mxh;
            parameterLabel = "mxh";
        } else if (strcmp("myh", params[0]) == 0) {
            parameter = v->myh;
            parameterLabel = "myh";
        } else if (strcmp("mxv", params[0]) == 0) {
            parameter = v->mxv;
            parameterLabel = "mxv";
        } else if (strcmp("myv", params[0]) == 0) {
            parameter = v->myv;
            parameterLabel = "myv";
        } else if (strcmp("dxh", params[0]) == 0) {
            parameter = v->dxh;
            parameterLabel = "dxh";
        } else if (strcmp("dyh", params[0]) == 0) {
            parameter = v->dyh;
            parameterLabel = "dyh";
        } else if (strcmp("dxv", params[0]) == 0) {
            parameter = v->dxv;
            parameterLabel = "dxv";
        } else if (strcmp("dyv", params[0]) == 0) {
            parameter = v->dyv;
            parameterLabel = "dyv";
        } else if (strcmp("rh", params[0]) == 0) {
            parameter = v->rh;
            parameterLabel = "rh";
        } else if (strcmp("rv", params[0]) == 0) {
            parameter = v->rv;
            parameterLabel = "rv";
        } else if (strcmp("VmcpH", params[0]) == 0) {
            parameter = v->vmcph;
            parameterLabel = "Vmcp H";
        } else if (strcmp("VmcpV", params[0]) == 0) {
            parameter = v->vmcpv;
            parameterLabel = "Vmcp V";
        } else if (strcmp("VbiasH", params[0]) == 0) {
            parameter = v->vbiash;
            parameterLabel = "Vbias H";
        } else if (strcmp("VbiasV", params[0]) == 0) {
            parameter = v->vbiasv;
            parameterLabel = "Vbias V";
        } else if (strcmp("Vfp", params[0]) == 0) {
            parameter = v->vfp;
            parameterLabel = "Vfp";
        } else if (strcmp("Vsatx", params[0]) == 0) {
            parameter = v->vsatx;
            parameterLabel = "Vsat X";
        } else if (strcmp("Vsaty", params[0]) == 0) {
            parameter = v->vsaty;
            parameterLabel = "Vsat Y";
        } else if (strcmp("Vsatz", params[0]) == 0) {
            parameter = v->vsatz;
            parameterLabel = "Vsat Z";
        } else if (strcmp("EnhRaw", params[0]) == 0) {
            parameter = v->enhRaw;
            parameterLabel = "En H raw";
        } else if (strcmp("EnvRaw", params[0]) == 0) {
            parameter = v->envRaw;
            parameterLabel = "En V raw";
        } else if (strcmp("Enh", params[0]) == 0) {
            parameter = v->enh;
            parameterLabel = "En H";
        } else if (strcmp("Env", params[0]) == 0) {
            parameter = v->env;
            parameterLabel = "En V";
        } else if (strcmp("Vixh", params[0]) == 0) {
            parameter = v->vixh;
            parameterLabel = "Vixh";
        } else if (strcmp("Vixherror", params[0]) == 0) {
            parameter = v->vixherror;
            parameterLabel = "Vixh error";
        } else if (strcmp("Vixv", params[0]) == 0) {
            parameter = v->vixv;
            parameterLabel = "Vixv";
        } else if (strcmp("Vixverror", params[0]) == 0) {
            parameter = v->vixverror;
            parameterLabel = "Vixv error";
        } else if (strcmp("Viy", params[0]) == 0) {
            parameter = v->viy;
            parameterLabel = "Viy";
        } else if (strcmp("Viyerror", params[0]) == 0) {
            parameter = v->viyerror;
            parameterLabel = "Viyerror";
        } else if (strcmp("Viz", params[0]) == 0) {
            parameter = v->viz;
            parameterLabel = "Viz";
        } else if (strcmp("Vizerror", params[0]) == 0) {
            parameter = v->vizerror;
            parameterLabel = "Vizerror";
        } else if (strcmp("VSatX", params[0]) == 0) {
            parameter = v->vsatx;
            parameterLabel = "Vsat X";
        } else if (strcmp("VSatY", params[0]) == 0) {
            parameter = v->vsaty;
            parameterLabel = "Vsat Y";
        } else if (strcmp("VSatZ", params[0]) == 0) {
            parameter = v->vsatz;
            parameterLabel = "Vsat Z";
        } else if (strcmp("VSatN", params[0]) == 0) {
            parameter = v->vsatn;
            parameterLabel = "Vsat N";
        } else if (strcmp("VSatE", params[0]) == 0) {
            parameter = v->vsate;
            parameterLabel = "Vsat E";
        } else if (strcmp("VSatC", params[0]) == 0) {
            parameter = v->vsatc;
            parameterLabel = "Vsat C";
        } else if (strcmp("BN", params[0]) == 0) {
            parameter = v->bn;
            parameterLabel = "B N";
        } else if (strcmp("BE", params[0]) == 0) {
            parameter = v->be;
            parameterLabel = "B E";
        } else if (strcmp("BC", params[0]) == 0) {
            parameter = v->bc;
            parameterLabel = "B C";
        } else if (strcmp("VicrX", params[0]) == 0) {
            parameter = v->vicrx;
            parameterLabel = "Vicr X";
        } else if (strcmp("VicrY", params[0]) == 0) {
            parameter = v->vicry;
            parameterLabel = "Vicr Y";
        } else if (strcmp("VicrZ", params[0]) == 0) {
            parameter = v->vicrz;
            parameterLabel = "Vicr Z";
        } else if (strcmp("Geopot", params[0]) == 0) {
            parameter = v->geoelectricPotential;
            parameterLabel = "Geopot";
        } else if (strcmp("GeopotAdj", params[0]) == 0) {
            parameter = v->geoelectricPotentialDetrended;
            parameterLabel = "Geopot Adj";
        } else if (strcmp("GeopotDiff", params[0]) == 0) {
            parameter = v->geoelectricPotentialDifference;
            parameterLabel = "Geopot diff.";
        } else if (strcmp("GeopotMaxslope", params[0]) == 0) {
            parameter = v->maxAbsGeoelectricPotentialBaselineSlope;
            parameterLabel = "Geopot max slope";
        } else if (strcmp("GeopotAdjMaxslope", params[0]) == 0) {
            parameter = v->maxAbsGeoelectricPotentialDetrendedBaselineSlope;
            parameterLabel = "Geopot adj max slope";
        } else if (strcmp("Exh", params[0]) == 0) {
            parameter = v->ectxh;
            parameterLabel = "Exh";
        } else if (strcmp("ExhAdj", params[0]) == 0) {
            parameter = v->ehxAdjusted;
            parameterLabel = "Exh adj.";
        } else if (strcmp("ExhAdjParam", params[0]) == 0) {
            parameter = v->ehxAdjustmentParameter;
            parameterLabel = "Exh adj. param";
        } else if (strcmp("Eyh", params[0]) == 0) {
            parameter = v->ectyh;
            parameterLabel = "Eyh";
        } else if (strcmp("Ezh", params[0]) == 0) {
            parameter = v->ectzh;
            parameterLabel = "Ezh";
        } else if (strcmp("Exv", params[0]) == 0) {
            parameter = v->ectxv;
            parameterLabel = "Exv";
        } else if (strcmp("Eyv", params[0]) == 0) {
            parameter = v->ectyv;
            parameterLabel = "Eyv";
        } else if (strcmp("Ezv", params[0]) == 0) {
            parameter = v->ectzv;
            parameterLabel = "Ezv";
        } else if (strcmp("Bx", params[0]) == 0) {
            parameter = v->bctx;
            parameterLabel = "B X";
        } else if (strcmp("By", params[0]) == 0) {
            parameter = v->bcty;
            parameterLabel = "B Y";
        } else if (strcmp("Bz", params[0]) == 0) {
            parameter = v->bctz;
            parameterLabel = "B Z";
        } else if (strcmp("xhatn", params[0]) == 0) {
            parameter = v->xhat;
            parameterLabel = "xhat N";
            tupleLength = 3;
            tupleIndex = 0;
        } else if (strcmp("xhate", params[0]) == 0) {
            parameter = v->xhat;
            parameterLabel = "xhat E";
            tupleLength = 3;
            tupleIndex = 1;
        } else if (strcmp("xhatc", params[0]) == 0) {
            parameter = v->xhat;
            parameterLabel = "xhat C";
            tupleLength = 3;
            tupleIndex = 2;
        } else if (strcmp("yhatn", params[0]) == 0) {
            parameter = v->yhat;
            parameterLabel = "yhat N";
            tupleLength = 3;
            tupleIndex = 0;
        } else if (strcmp("yhate", params[0]) == 0) {
            parameter = v->yhat;
            parameterLabel = "yhat E";
            tupleLength = 3;
            tupleIndex = 1;
        } else if (strcmp("yhatc", params[0]) == 0) {
            parameter = v->yhat;
            parameterLabel = "yhat C";
            tupleLength = 3;
            tupleIndex = 2;
        } else if (strcmp("zhatn", params[0]) == 0) {
            parameter = v->zhat;
            parameterLabel = "zhat N";
            tupleLength = 3;
            tupleIndex = 0;
        } else if (strcmp("zhate", params[0]) == 0) {
            parameter = v->zhat;
            parameterLabel = "zhat E";
            tupleLength = 3;
            tupleIndex = 1;
        } else if (strcmp("zhatc", params[0]) == 0) {
            parameter = v->zhat;
            parameterLabel = "zhat C";
            tupleLength = 3;
            tupleIndex = 2;
        } else if (strcmp("PhiSc", params[0]) == 0) {
            if (state->usePotentials) {
                parameter = v->potentials;
            }
            else {
                // Draw zeros
                parameter = v->mlt;
                valueScale = 0.0;
            }
            parameterLabel = "U_SC";
        } else {
            gotParameter = false;
        }

        float yRange = yr1 - yr0;
        yr0 += state->yOffset * yRange;
        yr1 += state->yOffset * yRange;
        snprintf(yr0str, 255, "%.1f", yr0);
        snprintf(yr1str, 255, "%.1f", yr1);

        if (gotParameter) {
            drawFloatTimeSeries(&image, v->timestamp, parameter, state->plotT0, state->plotT1, firstIndex, lastIndex, stride, valueScale, valueOffset, yr0, yr1, plotX0, plotYOffset, plotWidth, plotHeight, xLabel, parameterLabel, MAX_COLOR_VALUE + 1, yr0str, yr1str, false, dotSize, fontSize, true, tupleLength, tupleIndex, 0);
            plotYOffset += plotHeight + plotdy;
            plotsMade++;
        }

        if ((plotsMade > 0 && maxPlots > 0 && plotsMade % maxPlots == 0) || plotsMade == nPlots || plotYOffset > state->frameHeight - 1 - plotdy ) {

            // Add annotations
            addAnnotations(state, &image, topMargin);
            // Write video frames
            if (state->exportVideo) {
                for (int c = 0; c < 1.0 * state->framesPerSecond; c++) {
                    generateFrame(&image, frameCounter++);
                }
            }
            storeImage(state, &image);
            plotHeight = plotHeight0;
            plotYOffset = plotHeight + plotdy + topMargin;
        }
    }
    if (plotsMade == 0) {
        // Store anyway
        addAnnotations(state, &image, topMargin);
        storeImage(state, &image);
    }

    free(tofree);


finish:
    if (state->exportVideo) {
        finishVideo();
    }

    if (frameCounter > 0 && state->printVideoFilename) {
        printf("%s\n", state->videoFilename);
    }


cleanup:
    freeImage(&image);

    fflush(stdout);

    return TIICT_OK;
}

void drawFloatTimeSeries(Image *imageBuf, double *times, float *values, double t0, double t1, int firstInd, int lastInd, int stride, float valueScale, float valueOffset, float minValue, float maxValue, int plotX0, int plotY0, int plotWidth, int plotHeight, const char *xLabel, const char *yLabel, int colorIndex, const char *minValueStr, const char *maxValueStr, bool log10Scale, int dotSize, int fontSize, bool axes, int tupleLength, int tupleIndex, int orientation)
{
    int x0, y0;
    int x, y;

    int nValues = lastInd - firstInd + 1;
    double tmpVal;
    char label[255];

    if (t1 - t0 > 0 && nValues > 0)
    {
        if (axes)
        {
            drawAxes(imageBuf, t0, t1, xLabel, yLabel, minValueStr, maxValueStr, fontSize, plotX0, plotY0, plotWidth, plotHeight);
        }

        // data
        for (int i = firstInd; i < lastInd; i+=stride)
        {
            x0 = rescaleAsInteger(times[i], t0, t1, plotX0, plotX0 + plotWidth);
            tmpVal = (values[i*tupleLength + tupleIndex] - valueOffset) * valueScale;
            if (log10Scale)
            {
                if (tmpVal > 0)
                    tmpVal = log10(tmpVal);
                else
                    tmpVal = -10000.;
            }
            y0 = rescaleAsInteger(tmpVal, minValue, maxValue, plotY0, plotY0 - plotHeight);
            switch (dotSize)
            {
                case 2:
                    setBufferColorIndex(imageBuf, x0, y0, colorIndex);
                    setBufferColorIndex(imageBuf, x0+1, y0, colorIndex);
                    setBufferColorIndex(imageBuf, x0+1, y0+1, colorIndex);
                    setBufferColorIndex(imageBuf, x0, y0+1, colorIndex);
                    break;
                case 3:
                    setBufferColorIndex(imageBuf, x0, y0, colorIndex);
                    setBufferColorIndex(imageBuf, x0, y0-1, colorIndex);
                    setBufferColorIndex(imageBuf, x0, y0+1, colorIndex);
                    setBufferColorIndex(imageBuf, x0+1, y0, colorIndex);
                    setBufferColorIndex(imageBuf, x0+1, y0-1, colorIndex);
                    setBufferColorIndex(imageBuf, x0+1, y0+1, colorIndex);
                    setBufferColorIndex(imageBuf, x0-1, y0, colorIndex);
                    setBufferColorIndex(imageBuf, x0-1, y0-1, colorIndex);
                    setBufferColorIndex(imageBuf, x0-1, y0+1, colorIndex);
                    break;
                default:
                    setBufferColorIndex(imageBuf, x0, y0, colorIndex);
                    break;
            }
        }
    }

    return;
}

void drawAxes(Image *imageBuf, double t0, double t1, const char *xLabel, const char *yLabel, const char *minValueStr, const char *maxValueStr, int fontSize, int plotX0, int plotY0, int plotWidth, int plotHeight)
{

    double timeRange = (t1 - t0) / 1000.0;

    // time label string
    char timeFormat[EPOCHx_FORMAT_MAX];
    char timeString[EPOCHx_STRING_MAX];
    snprintf(timeFormat, EPOCHx_FORMAT_MAX, "<hour.02>:<min.02>:<sec.02>");

    double tickDeltaTSeconds = 0.0;

    if (timeRange < 10.0) {
        tickDeltaTSeconds = 1.0;
    } else if (timeRange < 30.0) {
        tickDeltaTSeconds = 5.0;
    } else if (timeRange < 60.0) {
        tickDeltaTSeconds = 10.0;
    } else if (timeRange < 60.0*10.0) {
        tickDeltaTSeconds = 60.0;
    } else if (timeRange < 60.0*60.0) {
        tickDeltaTSeconds = 300.0;
    } else if (timeRange < 60.0*60.0*12.0) {
        tickDeltaTSeconds = 3600.0;
    } else {
         tickDeltaTSeconds = 3600.0*3.0;
    }

    // Abscissa
    for (int s = 0; s <= timeRange; s+=tickDeltaTSeconds)
    {
        encodeEPOCHx(t0 + 1000.0 * s, timeFormat, timeString);
        annotate(timeString, fontSize, plotX0 + (int)(s / timeRange * plotWidth)-6, plotY0, imageBuf);
    }
    annotate(xLabel, fontSize, plotX0 + plotWidth/2 - (strlen(xLabel)*(8*fontSize))/24, plotY0+12, imageBuf);
    // Ordinate
    annotate(yLabel, fontSize, plotX0 + plotWidth + 5, plotY0 - plotHeight/2 - 6, imageBuf);
    annotate(minValueStr, fontSize, plotX0 + plotWidth+3, plotY0 - 8, imageBuf);
    annotate(maxValueStr, fontSize, plotX0 + plotWidth+3, plotY0 - plotHeight - 8, imageBuf);
    for (int o = plotY0; o >= plotY0 - plotHeight; o--)
    {
        setBufferColorIndex(imageBuf, plotX0 + plotWidth+1, o, FOREGROUND_COLOR);
    }

    return;
}

void addAnnotations(ProcessorState *state, Image *image, int topMargin)
{
    // Add annotations
    int fsize = 15;
    int yText = topMargin;
    char label[255];
    snprintf(label, 255, "Swarm %s", state->args.satellite);
    annotate(label, fsize, 5, yText, image);
    fsize = 12;
    yText += fontheight(fsize);
    snprintf(label, 255, "%04d-%02d-%02d", state->args.year, state->args.month, state->args.day);
    annotate(label, fsize, 5, yText, image);
    fsize = 9;
    yText += fontheight(fsize);
    snprintf(label, 255, "Vix: %s", state->useEofR ? "eofr" : "legacy");
    annotate(label, fsize, 5, yText, image);
    yText += fontheight(fsize);
    char *potentialSource = "None";
    switch (state->vars16hz.lpPotentialSource) {
        case LP_POTENTIAL_U_SC:
            potentialSource = "U_SC";
            break;
        case LP_POTENTIAL_LOWGAIN:
            potentialSource = "LG";
            break;
        case LP_POTENTIAL_HIGHGAIN:
            potentialSource = "HG";
            break;
        default:
            break;
    }
    snprintf(label, 255, "V_f algo: %s", potentialSource);
    annotate(label, fsize, 5, yText, image);

    return;
}

int storeImage(ProcessorState *state, Image *image)
{
    // Store image frame for potential later use
    void *mem = realloc(state->frames, sizeof *state->frames * (state->nVideoFrames + 1));
    if (mem == NULL) {
        fprintf(stderr, "Unable to allocate memory for new image\n");
        return TIICT_MEMORY;
    }
    state->frames = mem;
    state->nVideoFrames++;
    Image *f = &state->frames[state->nVideoFrames-1];
    allocImage(f, image->width, image->height, image->bytesPerPixel);
    memcpy(state->frames[state->nVideoFrames-1].pixels, image->pixels, image->numberOfBytes);

    // Reset image to make new plots
    memset(image->pixels, BACKGROUND_COLOR, image->numberOfBytes);

    return TIICT_OK;
}
