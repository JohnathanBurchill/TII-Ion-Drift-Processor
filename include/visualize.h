#ifndef VISUALIZE_H
#define VISUALIZE_H

#include "state.h"
#include <tiigraphics/draw.h>

int visualizeResults(ProcessorState *state);
void drawFloatTimeSeries(Image *imageBuf, double *times, float *values, int firstInd, int lastInd, int stride, float valueScale, float minValue, float maxValue, int plotX0, int plotY0, int plotWidth, int plotHeight, const char *xLabel, const char *yLabel, int colorIndex, const char *minValueStr, const char *maxValueStr, bool log10Scale, int dotSize, int fontSize, bool axes, int tupleLength, int tupleIndex, int orientation);
void addAnnotations(ProcessorState *state, Image *image, int topMargin);
int storeImage(ProcessorState *state, Image *image);


#endif // VISUALIZE_H
