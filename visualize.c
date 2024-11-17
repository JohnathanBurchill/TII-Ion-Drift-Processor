#include "visualize.h"
#include "state.h"
#include "errors.h"

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
    int plotHeight = 100;
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
    double *values = malloc(sizeof *values * nValues);
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

    int plotHeight0 = 45;
    int plotHeight1 = 30;
    int plotdy = 25;

    int plotX0 = 100;
    int plotY0 = 100;
    int plotY1 = plotY0 + plotHeight0 + plotdy;
    int plotY2 = plotY1 + plotHeight0 + plotdy;
    int plotY3 = plotY2 + plotHeight0 + plotdy;
    int plotY4 = plotY3 + plotHeight0 + plotdy;
    int plotY5 = plotY4 + plotHeight0 + plotdy;

    // Plots
    drawFloatTimeSeries(&image, (double*)dataBuffers[0], (float*)dataBuffers[5], first, last, 1, 1.0, -90.0, 90.0, plotX0, plotY0, plotWidth, plotHeight0, "", "QD Lat", MAX_COLOR_VALUE + 1, "-90", "90", false, dotSize, 12, true, 1, 0);
    drawFloatTimeSeries(&image, (double*)dataBuffers[0], (float*)dataBuffers[1], first, last, 1, 0.001, -4.0, 4.0, plotX0, plotY1, plotWidth, plotHeight0, "", "Vixh", MAX_COLOR_VALUE + 1, "-4", "4", false, dotSize, 12, true, 2, 0);
    drawFloatTimeSeries(&image, (double*)dataBuffers[0], (float*)dataBuffers[2], first, last, 1, 0.001, -4.0, 4.0, plotX0, plotY2, plotWidth, plotHeight0, "", "Vixv", MAX_COLOR_VALUE + 1, "-4", "4", false, dotSize, 12, true, 2, 0);
    drawFloatTimeSeries(&image, (double*)dataBuffers[0], (float*)dataBuffers[1], first, last, 1, 0.001, -4.0, 4.0, plotX0, plotY3, plotWidth, plotHeight0, "", "Viy", MAX_COLOR_VALUE + 1, "-4", "4", false, dotSize, 12, true, 2, 1);
    drawFloatTimeSeries(&image, (double*)dataBuffers[0], (float*)dataBuffers[2], first, last, 1, 0.001, -4.0, 4.0, plotX0, plotY4, plotWidth, plotHeight0, "", "Viz", MAX_COLOR_VALUE + 1, "-4", "4", false, dotSize, 12, true, 2, 1);

    for (int c = 0; c < 0.1 * VIDEO_FPS; c++) {
        generateFrame(&image, frameCounter++);
    }

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

    double t0 = times[firstInd];
    double t1 = times[lastInd];
    int nValues = lastInd - firstInd + 1;
    double timeRange = t1 - t0;
    double tmpVal;
    char label[255];

    double tickDeltaTSeconds = 0.0;

    int nTicks = 5;
    char *tickUnit = "s";
    if (timeRange < 10.0) {
        tickDeltaTSeconds = 1.0;
    } else if (timeRange < 30.0) {
        tickDeltaTSeconds = 5.0;
    } else if (timeRange < 60.0) {
        tickDeltaTSeconds = 10.0;
    } else if (timeRange < 60.0*10.0) {
        tickDeltaTSeconds = 60.0;
        tickUnit = "m";
    } else if (timeRange < 60.0*60.0) {
        tickDeltaTSeconds = 300.0;
    } else if (timeRange < 60.0*60.0*12.0) {
        tickDeltaTSeconds = 3600.0;
        tickUnit = "h";
    } else {
         tickDeltaTSeconds = 3600.0*3.0;
        tickUnit = "h";
    }

    if (timeRange > 0 && nValues > 0)
    {
        if (axes)
        {
            // Abscissa
            for (int s = 0; s <= timeRange; s+=tickDeltaTSeconds)
            {

                sprintf(label, "%d", (int)(s / tickDeltaTSeconds));
                annotate(label, fontSize, plotX0 + (int)(s / timeRange * plotWidth)-6, plotY0, imageBuf);
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
            x0 = rescaleAsInteger(times[i], t0, t1, plotX0, plotX0 + plotWidth);
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
