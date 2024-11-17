#include "visualize.h"
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
    int frameCounter = 0;

    Image templateImage = {0};
    Image image = {0};

    if (allocImage(&templateImage, IMAGE_WIDTH, IMAGE_HEIGHT, 1) != DRAW_OK)
    {
        printf("Could not allocate memory for template image.\n");
        goto cleanup;
    }
    if (allocImage(&image, IMAGE_WIDTH, IMAGE_HEIGHT, 1) != DRAW_OK)
    {
        printf("Could not allocate memory for image.\n");
        goto cleanup;
    }

    int status = initVideo(state->movieFilename);
    if (status < 0)
    {
        fprintf(stderr, "Problem intializing video: got status %d.\n", status);
        goto cleanup;
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

    // Input data
    uint8_t **dataBuffers = state->dataBuffers;
    double *allTimes = (double*)dataBuffers[0];
    float *param = NULL;

    int first = 0;
    int last = 0;
    // Find first index
    if (state->movieT0 >= 0) {
        while (first < state->nRecs - 1 && allTimes[first] < state->movieT0) {
            first++;
        }
    }
    if (state->movieT1 >= 0) {
        while (last < state->nRecs && allTimes[last] < state->movieT1) {
            last++;
        }
    }
    size_t nValues = last - first + 1;
    double *times = malloc(sizeof *times * nValues);
    float *values = malloc(sizeof *values * nValues);
    if (times == NULL || values == NULL) {
        return TIICT_MEMORY;
    }
    // times to plot
    int index = 0;
    for (int i = first; i < last; i++) {
        times[index++] = allTimes[i]/1000.0;
    }
    //insertTransition(&image, "Anomaly overview", IMAGE_WIDTH/2, IMAGE_HEIGHT/2-16, 24, 3.0, &frameCounter);

    // Set background
    memset(image.pixels, BACKGROUND_COLOR, image.numberOfBytes);

    int plotdy = 25;

    int plotX0 = 100;
    int topMargin = 10;
    int plotY0 = plotHeight0 + plotdy + topMargin;

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
    float yScale = 0.0;
    char yr0str[255];
    char yr1str[255];
    int plotYOffset = plotY0;
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
    // Plot requested plots
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
            fprintf(stderr, "Plot command %s: expected at least %d plot parameters, skipping this plot\n", plotOptions, nParams);
            gotParameter = false;
            continue;
        }

        yr0 = atof(params[1]);
        snprintf(yr0str, 255, "%s", params[1]);
        yr1 = atof(params[2]);
        snprintf(yr1str, 255, "%s", params[2]);
        yScale = atof(params[3]);

        int oldPlotHeight = plotHeight;
        if (nParams == 5) {
            plotHeight = atoi(params[4]);
        }
        else {
            plotHeight = plotHeight0;
        }
        plotYOffset += (plotHeight - oldPlotHeight);

        if (strcmp("QDLat", params[0]) == 0) {
            parameter = (float*)dataBuffers[5];
            parameterLabel = "QD Lat";
        } else if (strcmp("PhiSc", params[0]) == 0) {
            if (state->usePotentials) {
                parameter = state->potentials;
            }
            else {
                // Draw zeros
                parameter = (float*)state->dataBuffers[1];
                yScale = 0.0;
            }
            parameterLabel = "U_SC";
        } else if (strcmp("Vixh", params[0]) == 0) {
            parameter = (float*)dataBuffers[1];
            parameterLabel = "Vixh";
            tupleLength = 2;
            tupleIndex = 0;
        } else if (strcmp("Vixv", params[0]) == 0) {
            parameter = (float*)dataBuffers[2];
            parameterLabel = "Vixv";
            tupleLength = 2;
            tupleIndex = 0;
        } else if (strcmp("Viy", params[0]) == 0) {
            parameter = (float*)dataBuffers[1];
            parameterLabel = "Viy";
            tupleLength = 2;
            tupleIndex = 1;
        } else if (strcmp("Viz", params[0]) == 0) {
            parameter = (float*)dataBuffers[2];
            parameterLabel = "Viz";
            tupleLength = 2;
            tupleIndex = 1;
        } else {
            gotParameter = false;
        }

        if (gotParameter) {
            drawFloatTimeSeries(&image, (double*)dataBuffers[0], parameter, first, last, stride, yScale, yr0, yr1, plotX0, plotYOffset, plotWidth, plotHeight, xLabel, parameterLabel, MAX_COLOR_VALUE + 1, yr0str, yr1str, false, dotSize, 12, true, tupleLength, tupleIndex);
            plotYOffset += plotHeight + plotdy;
            plotsMade++;
        }

        if ((maxPlots > 0 && plotsMade % maxPlots == 0) || plotsMade == nPlots || plotYOffset > IMAGE_HEIGHT - 1 - plotdy ) {
            // Write video frames
            for (int c = 0; c < 1.0 * VIDEO_FPS; c++) {
                generateFrame(&image, frameCounter++);
            }
            // Reset image to make new plots
            memset(image.pixels, BACKGROUND_COLOR, image.numberOfBytes);
            plotHeight = plotHeight0;
            plotYOffset = plotHeight + plotdy + topMargin;
        }
    }

    free(tofree);


    finishVideo();

    if (frameCounter > 0)
        printf("%s\n", state->movieFilename);
    else
        printf("No-Frames-For-This-Date\n");


cleanup:
    freeImage(&templateImage);
    freeImage(&image);

    free(values);

    fflush(stdout);

    return TIICT_OK;
}

void drawFloatTimeSeries(Image *imageBuf, double *times, float *values, int firstInd, int lastInd, int stride, float valueScale, float minValue, float maxValue, int plotX0, int plotY0, int plotWidth, int plotHeight, const char *xLabel, const char *yLabel, int colorIndex, const char *minValueStr, const char *maxValueStr, bool log10Scale, int dotSize, int fontSize, bool axes, int tupleLength, int tupleIndex)
{
    int x0, y0;
    int x, y;

    double t0 = times[firstInd]/1000.0;
    double t1 = times[lastInd]/1000.0;
    int nValues = lastInd - firstInd + 1;
    double timeRange = t1 - t0;
    double tmpVal;
    char label[255];

    // time label string
    char timeFormat[EPOCHx_FORMAT_MAX];
    char timeString[EPOCHx_STRING_MAX];
    snprintf(timeFormat, EPOCHx_FORMAT_MAX, "<hour.02>:<min.02>");

    double tickDeltaTSeconds = 0.0;

    if (timeRange < 10.0) {
        tickDeltaTSeconds = 1.0;
        snprintf(timeFormat, EPOCHx_FORMAT_MAX, "<hour.02>:<min.02>:<sec.02>");
    } else if (timeRange < 30.0) {
        tickDeltaTSeconds = 5.0;
        snprintf(timeFormat, EPOCHx_FORMAT_MAX, "<hour.02>:<min.02>:<sec.02>");
    } else if (timeRange < 60.0) {
        tickDeltaTSeconds = 10.0;
        snprintf(timeFormat, EPOCHx_FORMAT_MAX, "<hour.02>:<min.02>:<sec.02>");
    } else if (timeRange < 60.0*10.0) {
        tickDeltaTSeconds = 60.0;
    } else if (timeRange < 60.0*60.0) {
        tickDeltaTSeconds = 300.0;
    } else if (timeRange < 60.0*60.0*12.0) {
        tickDeltaTSeconds = 3600.0;
    } else {
         tickDeltaTSeconds = 3600.0*3.0;
    }

    if (timeRange > 0 && nValues > 0)
    {
        if (axes)
        {
            // Abscissa
            for (int s = 0; s <= timeRange; s+=tickDeltaTSeconds)
            {
                encodeEPOCHx(times[firstInd] + 1000.0 * s, timeFormat, timeString);
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
        }

        // data
        for (int i = firstInd; i < lastInd; i+=stride)
        {
            x0 = rescaleAsInteger(times[i]/1000.0, t0, t1, plotX0, plotX0 + plotWidth);
            tmpVal = values[i*tupleLength + tupleIndex]*valueScale;
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
