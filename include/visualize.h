#ifndef VISUALIZE_H
#define VISUALIZE_H

#include "state.h"
#include <tiigraphics/draw.h>

// Horizontal extent of the time-series plots within a frame
#define PLOT_X0 100
#define PLOT_WIDTH_FRACTION 0.9

int visualizeResults(ProcessorState *state);
void drawFloatTimeSeries(Image *imageBuf, double *times, float *values, double t0, double t1, int firstInd, int lastInd, int stride, float valueScale, float valueOffset, float minValue, float maxValue, int plotX0, int plotY0, int plotWidth, int plotHeight, const char *xLabel, const char *yLabel, int colorIndex, const char *minValueStr, const char *maxValueStr, bool log10Scale, int dotSize, int fontSize, bool axes, int tupleLength, int tupleIndex, int orientation);
void addAnnotations(ProcessorState *state, Image *image, int topMargin);
int storeImage(ProcessorState *state, Image *image);
void drawAxes(Image *imageBuf, double t0, double t1, const char *xLabel, const char *yLabel, const char *minValueStr, const char *maxValueStr, int fontSize, int plotX0, int plotY0, int plotWidth, int plotHeight);


#endif // VISUALIZE_H
