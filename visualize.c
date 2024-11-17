#include "visualize.h"
#include "state.h"
#include "errors.h"

#include <tiigraphics/draw.h>
#include <tiigraphics/tiigraphics.h>
#include <tiigraphics/video.h>
#include <tiigraphics/colors.h>

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
    int plotdy = 15;

    int plotX0 = 100;
    int plotY0 = 100;
    int plotY1 = plotY0 + plotHeight0 + plotdy;
    int plotY2 = plotY1 + plotHeight0 + plotdy;
    int plotY3 = plotY2 + plotHeight0 + plotdy;
    int plotY4 = plotY3 + plotHeight1 + plotdy;
    int plotY5 = plotY4 + plotHeight1 + plotdy;

    // Latitude
    param = (float*)dataBuffers[5];
    index = 0;
    for (int i = first; i < last; i++) {
        values[index++] = (double)param[i];
    }
    drawTimeSeries(&image, times, values, nValues, plotX0, plotY0, plotWidth, plotHeight0, state->movieT0/1000.0, state->movieT1/1000.0, -90, 90, "", "QDLat", 1, MAX_COLOR_VALUE + 1, "-90", "90", false, dotSize, 12, true);

    // Vixh
    param = (float*)dataBuffers[1];
    index = 0;
    for (int i = 2*first; i < 2*last; i+=2) {
        values[index++] = (double)param[i];
    }
    drawTimeSeries(&image, times, values, nValues, plotX0, plotY1, plotWidth, plotHeight0, state->movieT0/1000.0, state->movieT1/1000.0, -4000, 4000, "", "Vixh", 1, MAX_COLOR_VALUE + 1, "-4000", "4000", false, dotSize, 12, true);

    // Vixv
    param = (float*)dataBuffers[2];
    index = 0;
    for (int i = 2*first; i < 2*last; i+=2) {
        values[index++] = (double)param[i];
    }
    drawTimeSeries(&image, times, values, nValues, plotX0, plotY2, plotWidth, plotHeight0, state->movieT0/1000.0, state->movieT1/1000.0, -4000, 4000, "", "Vixv", 1, MAX_COLOR_VALUE + 1, "-4000", "4000", false, dotSize, 12, true);

    // Viy
    param = (float*)dataBuffers[1];
    index = 0;
    for (int i = 2*first+1; i < 2*last+1; i+=2) {
        values[index++] = (double)param[i];
    }
    drawTimeSeries(&image, times, values, nValues, plotX0, plotY3, plotWidth, plotHeight0, state->movieT0/1000.0, state->movieT1/1000.0, -4000, 4000, "", "Viy", 1, MAX_COLOR_VALUE + 1, "-4000", "4000", false, dotSize, 12, true);

    // Geopotential
    param = state->geoPotential;
    index = 0;
    for (int i = first; i < last; i++) {
        values[index++] = (double)param[i]/1000.0;
    }
    drawTimeSeries(&image, times, values, nValues, plotX0, plotY4, plotWidth, plotHeight0, state->movieT0/1000.0, state->movieT1/1000.0, -200, 200, "", "Geopot.", 1, MAX_COLOR_VALUE + 1, "-200", "200", false, dotSize, 12, true);


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

