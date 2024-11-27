/*

    TII Cross-Track Ion Drift Processor: tiictstripchart.c

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

#include "SDL3/SDL_events.h"
#include "export.h"
#include "processing.h"
#include "loadData.h"
#include "settings.h"
#include "tiigraphics/draw.h"
#include "tiigraphics/fonts.h"
#include "visualize.h"
#include "state.h"
#include "errors.h"

#include "SDL3/SDL_error.h"
#include "SDL3/SDL_init.h"
#include "SDL3/SDL_keycode.h"
#include "SDL3/SDL_log.h"
#include "SDL3/SDL_pixels.h"
#include "SDL3/SDL_render.h"
#include "SDL3/SDL_surface.h"
#include "SDL3/SDL_video.h"

#include <tiigraphics/colors.h>
#include <tiigraphics/tiigraphics.h>
#include <tiigraphics/video.h>

#include <stdio.h>
#include <unistd.h>

#define SDL_MAIN_USE_CALLBACKS 1  /* use the callbacks instead of main() */
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

typedef struct AppState {
    ProcessorState *stateA;
    ProcessorState *stateB;
    ProcessorState *stateC;
    ProcessorState *state;
    int plotPage;
    double dayBegin;
    double dayEnd;
    double t0;
    double t1;
    double samplePeriodSeconds;
    SDL_Window *window;
    SDL_Renderer *plotRenderer;
    SDL_Palette *colors;
    bool playing;
    int playbackDirection;
    double playbackRate;
    Image *help;
    Image *storedFrames;
    int storedNVideoFrames;
    bool show16Hz;
} AppState_t;

typedef enum TimeUnit {
    TIME_RANGES = 0,
    SAMPLE_PERIODS,
    SECONDS,
    MINUTES,
    HOURS,
    DAYS,
    WEEKS,
    MONTHS,
    SEASONS,
    YEARS,
    DECADES,
} TimeUnit_enum;

Image *helpImage(int width, int height);
void resetVideoFrames(ProcessorState *state);
void resetDisplay(AppState_t *as);
double calculateDeltaT(AppState_t *as, TimeUnit_enum units, int sign);
void advancePlots(AppState_t *as, double amount, TimeUnit_enum units);
void rewindPlots(AppState_t *as, double amount, TimeUnit_enum units);
void updatePlots(ProcessorState *state);
void rerunProcessor(AppState_t *state);

SDL_AppResult SDL_AppInit(void **appstate, int argc, char *argv[])
{

    AppState_t *as = malloc(sizeof *as);
    if (as == NULL) {
        SDL_Log("Unable to allocate memory for App state");
        return SDL_APP_FAILURE;
    }
    *appstate = (void*)as;
    memset(as, 0, sizeof *as);

    ProcessorState *state = initState(argc, argv);
    if (state == NULL) {
        fprintf(stderr, "Could not allocate processor state.\n");
        return SDL_APP_FAILURE;
    }

    // New defaults
    state->writeLogFiles = true;
    state->export16Hz = false;
    state->export2Hz = false;
    state->exportZip = false;
    state->usePotentials = false;
    state->vars16hz.lpPotentialSource = LP_POTENTIAL_NONE;
    state->vars2hz.lpPotentialSource = LP_POTENTIAL_NONE;
    state->visualizeResults = true;
    state->exportVideo = false;

    runProcessor(argc, argv, &state);

    // Display 16 Hz data by default
    as->show16Hz = true;
    state->vars = &state->vars16hz;

    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_Log("Couldn't initialize SDL: %s", SDL_GetError());
        return SDL_APP_FAILURE;
    }

    SDL_SetLogPriorities(SDL_LOG_PRIORITY_CRITICAL);

    if (!SDL_CreateWindowAndRenderer("examples/renderer/clear", state->frameWidth, state->frameHeight, 0, &as->window, &as->plotRenderer)) {
        SDL_Log("Couldn't create window/renderer: %s", SDL_GetError());
        return SDL_APP_FAILURE;
    }
    SDL_SetRenderDrawBlendMode(as->plotRenderer, SDL_BLENDMODE_BLEND);
    as->colors = SDL_CreatePalette(256); // Allocate a 256-color palette
    if (as->colors == NULL) {
        fprintf(stderr, "Failed to create a palette: %s\n", SDL_GetError());
    }
    for (int i = 0; i <= MAX_COLOR_VALUE; i++) {
        as->colors->colors[i] = (SDL_Color){i, i, i, 255};
    }
    as->colors->colors[FOREGROUND_COLOR] = (SDL_Color){0, 0, 0, 255};
    as->colors->colors[FOREGROUND_COLOR + 1] = (SDL_Color){10, 10, 10, 255};
    as->colors->colors[FOREGROUND_COLOR + 2] = (SDL_Color){20, 20, 20, 255};
    as->colors->colors[BACKGROUND_COLOR ] = (SDL_Color){255, 255, 255, 255};

    resetDisplay(as);

    as->state = state;
    as->plotPage = 0;
    as->dayBegin = computeEPOCH(state->args.year, state->args.month, state->args.day, 0, 0, 0, 0);
    as->dayEnd = as->t0 + 86400.0 * 1000.0; // Ignore leap seconds
    double *timesMs = state->vars->timestamp;
    if (state->vars->nRecs > 1) {
        as->samplePeriodSeconds = (timesMs[1] - timesMs[0]) / 1000.0;
    }
    else {
        // arbitrary value
        as->samplePeriodSeconds = 1.0;
    }
    as->playing = false;
    as->playbackDirection = 1;
    as->playbackRate = 1.0;

    as->help = helpImage(state->frameWidth, state->frameHeight);
    if (as->help == NULL) {
        return SDL_APP_FAILURE;
    }
    as->storedFrames = NULL;
    as->storedNVideoFrames = 0;


    return SDL_APP_CONTINUE;  /* carry on with the program! */
}

Image *helpImage(int width, int height)
{
    int x = 50;
    int y0 = 50;
    int y = y0;
    // font size
    int fontSize = 15;
    int fontHeight = fontheight(fontSize);

    Image *help = malloc(sizeof *help);
    int status = allocImage(help, width, height, 1);
    if (status != TIICT_OK) {
        SDL_Log("Unable to allocate image");
        return NULL;
    }

    // Reset image to make new plots
    memset(help->pixels, BACKGROUND_COLOR, help->numberOfBytes);

    annotate("F1 - help", fontSize, x, y, help);
    y += fontHeight;
    annotate(" q - quit", fontSize, x, y, help);
    y += fontHeight*2;
    annotate(" a - Swarm A", fontSize, x, y, help);
    y += fontHeight;
    annotate(" b - Swarm B", fontSize, x, y, help);
    y += fontHeight;
    annotate(" c - Swarm C", fontSize, x, y, help);
    y += fontHeight*2;
    annotate(" 1 - LP PhiSc source: None", fontSize, x, y, help);
    y += fontHeight;
    annotate(" 2 - LP PhiSc source: U_SC", fontSize, x, y, help);
    y += fontHeight;
    annotate(" 3 - LP PhiSc source: Low Gain", fontSize, x, y, help);
    y += fontHeight;
    annotate(" 4 - LP PhiSc source: High Gain", fontSize, x, y, help);
    y += fontHeight;

    x += IMAGE_WIDTH/2;
    y = y0;
    annotate("F2 - Ion drift", fontSize, x, y, help);
    y += fontHeight;
    annotate("F3 - Electric field", fontSize, x, y, help);
    y += fontHeight;
    annotate("F4 - Magnetic field", fontSize, x, y, help);
    y += fontHeight*2;
    annotate("F5 - Satellite position", fontSize, x, y, help);
    y += fontHeight;
    annotate("F6 - Satellite velocity", fontSize, x, y, help);
    y += fontHeight;
    annotate("F7 - Floating potential", fontSize, x, y, help);
    y += fontHeight;
    annotate("F8 - Geopotential", fontSize, x, y, help);
    y += fontHeight;

    return help;
}

SDL_AppResult SDL_AppEvent(void *appstate, SDL_Event *event)
{
    if (event->type == SDL_EVENT_QUIT || (event->type == SDL_EVENT_KEY_UP && event->key.key == SDLK_Q)) {
        return SDL_APP_SUCCESS;  /* end the program, reporting success to the OS. */
    }

    int status = 0;
    AppState_t *as = (AppState_t *)appstate;
    int argc = 0;
    char **argv = NULL;

    ProcessorState *state = (ProcessorState *)as->state;
    if (state == NULL) {
        return SDL_APP_CONTINUE;
    }

    double *timesMs = state->vars->timestamp;
    if (timesMs != NULL) {
        as->t0 = timesMs[0];
        as->t1 = timesMs[state->vars->nRecs - 1];
    }
    double timeRange = state->plotT1 - state->plotT0;
    double middleTime = 0.0;
    double secondsToAdvance = 0.0;

    int plotHeightDelta = 0;

    // Restore the display after a help request
    if (as->storedFrames != NULL) {
        state->frames = as->storedFrames;
        state->nVideoFrames = as->storedNVideoFrames;
        as->storedFrames = NULL;
        as->storedNVideoFrames = 0;
    }

    if (event->type == SDL_EVENT_KEY_UP) {
        switch (event->key.key) {
            case SDLK_A:
                state->args.satellite = "A";
                rerunProcessor(as);
                break;
            case SDLK_B:
                state->args.satellite = "B";
                rerunProcessor(as);
                break;
            case SDLK_C:
                state->args.satellite = "C";
                rerunProcessor(as);
                break;
            case SDLK_E:
                // Toggle use of eofr for along-track drift
                state->useEofR = !state->useEofR;
                rerunProcessor(as);
                break;
            case SDLK_1:
                state->vars->lpPotentialSource = LP_POTENTIAL_NONE;
                state->usePotentials = false;
                rerunProcessor(as);
                break;
            case SDLK_2:
                state->vars->lpPotentialSource = LP_POTENTIAL_U_SC;
                state->usePotentials = true;
                rerunProcessor(as);
                break;
            case SDLK_3:
                state->vars->lpPotentialSource = LP_POTENTIAL_LOWGAIN;
                state->usePotentials = true;
                rerunProcessor(as);
                break;
            case SDLK_4:
                state->vars->lpPotentialSource = LP_POTENTIAL_HIGHGAIN;
                state->usePotentials = true;
                rerunProcessor(as);
                break;
            case SDLK_SEMICOLON:
                as->show16Hz = !as->show16Hz;
                rerunProcessor(as);
                break;
            case SDLK_U:
                // Update processor results
                rerunProcessor(as);
                break;
            case SDLK_EQUALS:
                // Plus on regular keboard
                if (SDL_GetModState() & SDL_KMOD_SHIFT) {
                    middleTime = (state->plotT0 + state->plotT1) / 2.0;
                    if (timeRange > 2.0 * 1000) {
                        timeRange /= 2.0;
                        state->plotT0 = middleTime - timeRange/2.0;
                        state->plotT1 = middleTime + timeRange/2.0;
                        updatePlots(state);
                    }
                }
                break;
            case SDLK_MINUS:
                middleTime = (state->plotT0 + state->plotT1) / 2.0;
                timeRange *= 2.0;
                if (timeRange > 86400.0 * 1000.0) {
                    timeRange = 86400.0 * 1000.0;
                }
                state->plotT0 = middleTime - timeRange/2.0;
                state->plotT1 = middleTime + timeRange/2.0;
                updatePlots(state);
                break;
            case SDLK_K:
                if (as->plotPage > 0) {
                    as->plotPage--;
                }
                break;
            case SDLK_J:
                if (as->plotPage < state->nVideoFrames - 1) {
                    as->plotPage++;
                }
                break;
            case SDLK_D:
                // Full day timerange
                state->plotT0 = as->dayBegin;
                state->plotT1 = as->dayEnd;
                updatePlots(state);
                break;
            case SDLK_F:
                // Full file timerange
                state->plotT0 = as->t0;
                state->plotT1 = as->t1;
                updatePlots(state);
                break;
            case SDLK_O:
                // Full orbit (approx.)
                timeRange = 94.0 * 60.0 * 1000.0;
                if (SDL_GetModState() & SDL_KMOD_SHIFT) {
                    // 1/4 orbit timerange
                    timeRange /= 4.0;
                }
                state->plotT0 = as->t0;
                state->plotT1 = state->plotT0 + timeRange;
                updatePlots(state);
                break;
            case SDLK_PERIOD:
                // Advance the plot by 10% of timeRange
                secondsToAdvance = timeRange / 1000.0 / 10.0;
                if (SDL_GetModState() & SDL_KMOD_SHIFT) {
                    secondsToAdvance /= 5.0;
                }
                if (SDL_GetModState() & (SDL_KMOD_LCTRL | SDL_KMOD_RCTRL)) {
                    secondsToAdvance /= 5.0;
                }
                if (secondsToAdvance < 1.0) {
                    secondsToAdvance = 1.0;
                }
                advancePlots(as, secondsToAdvance, SECONDS);
                break;
            case SDLK_COMMA:
                // Rewind the plot by 10% of timeRange
                secondsToAdvance = timeRange / 1000.0 / 10.0;
                if (SDL_GetModState() & SDL_KMOD_SHIFT) {
                    secondsToAdvance /= 5.0;
                }
                if (SDL_GetModState() & (SDL_KMOD_LCTRL | SDL_KMOD_RCTRL)) {
                    secondsToAdvance /= 5.0;
                }
                if (secondsToAdvance < 1.0) {
                    secondsToAdvance = 1.0;
                }
                rewindPlots(as, secondsToAdvance, SECONDS);
                break;
            case SDLK_H:
                // Rewind the plot timerange by the current timerange
                rewindPlots(as, 1, TIME_RANGES);
                break;
            case SDLK_L:
                // Advance the plot timerange by the current timerange
                advancePlots(as, 1, TIME_RANGES);
                break;
            case SDLK_SPACE:
                // Toggle playback
                as->playing = !as->playing;
                break;
            case SDLK_R:
                // Toggle playback direction
                as->playbackDirection = -as->playbackDirection;
                break;
            case SDLK_X:
                // Toggle playback rate
                as->playbackRate = as->playbackRate > 1.0 ? 1.0 : 10.0;
                break;

            case SDLK_F1:
                as->storedFrames = state->frames;
                as->storedNVideoFrames = state->nVideoFrames;
                state->frames = as->help;
                state->nVideoFrames = 1;
                break;

            case SDLK_F2:
                // ion drift
                state->plotCommand = "QDLat,-90,90,1" ";PhiSc,-5,0,1" ";Vixh,-4,4,0.001" ";Vixv,-4,4,0.001" ";Viy,-2,2,0.001" ";Viz,-2,2,0.001";
                updatePlots(state);
                break;
            case SDLK_F3:
                // ion drift
                state->plotCommand = "Vixh,-2,2,0.001" ";Vixv,-2,2,0.001" ";Viy,-2,2,0.001" ";Viz,-2,2,0.001";
                updatePlots(state);
                break;
            case SDLK_F4:
                // electric field
                state->plotCommand = "QDLat,-90,90,1" ";MLT,0,24,1" ";Exh,-100,100,1" ";Eyh,-100,100,1" ";Ezh,-100,100,1";
                updatePlots(state);
                break;
            case SDLK_F5:
                // magnetic field
                state->plotCommand = "BN,-65000,65000,1" ";BE,-65000,65000,1" ";BC,-65000,65000,1" ";Bx,-65000,65000,1" ";By,-65000,65000,1" ";Bz,-65000,65000,1";
                updatePlots(state);
                break;
            case SDLK_F6:
                // Satellite position
                state->plotCommand = "Lat,-90,90,1" ";Lon,-180,180,1" ";Radius,6400,7000,0.001" ";Alt,430,520,0.001";
                updatePlots(state);
                break;
            case SDLK_F7:
                // Satellite velocity
                state->plotCommand = "VSatN,-8,8,0.001" ";VSatE,-8,8,0.001" ";VSatC,-0.1,0.1,0.001" ";VSatX,-8,8,0.001" ";VSatY,-2,2,0.001" ";VSatZ,-0.5,0.5,0.001" ";VicrX,-0.5,0.5,0.001" ";VicrY,-0.6,0.6,0.001" ";VicrZ,-0.5,0.5,0.001";
                updatePlots(state);
                break;
            case SDLK_F8:
                // Satellite floating potential
                state->plotCommand = "PhiSc,-8,1,1";
                updatePlots(state);
                break;
            case SDLK_F9:
                // Geoelectric potential
                state->plotCommand = "Exh,-50,50,1" ";Geopot,-100,100,0.001" ";ExhAdj,-50,50,1" ";GeopotAdj,-100,100,0.001" ";ExhAdjParam,-1,1,1";
                updatePlots(state);
                break;
            case SDLK_F10:
                // Moments
                state->plotCommand = "mxh,0,66,1" ";mxv,0,66,1" ";myh,0,66,1" ";myv,0,66,1";
                updatePlots(state);
                break;
            case SDLK_F11:
                // Moments
                state->plotCommand = "EnhRaw,0,30,1" ";EnvRaw,0,30,1" ";Enh,0,30,1" ";Env,0,30,1";
                updatePlots(state);
                break;
            case SDLK_F12:
                // Moments
                state->plotCommand = "VmcpH,-2400,0,1" ";VmcpV,-2400,0,1" ";VbiasH,-105,0,1" ";VbiasV,-105,0,1" ";Vfp,-5,0,1";
                updatePlots(state);
                break;
            case SDLK_Z:
                if (SDL_GetModState() & SDL_KMOD_SHIFT) {
                    if (state->defaultPlotHeight > DEFAULT_PLOT_HEIGHT) {
                        state->defaultPlotHeight -= 30;
                        updatePlots(state);
                    }
                }
                else {
                    if (state->defaultPlotHeight < state->frameHeight - 70) {
                        state->defaultPlotHeight += 30;
                        updatePlots(state);
                    }
                }
                if (as->plotPage > state->nVideoFrames - 1) {
                    as->plotPage = state->nVideoFrames - 1;
                }
                if (as->plotPage < 0) {
                    as->plotPage = 0;
                }
                break;
            default:
                break;
        }

    }

    return SDL_APP_CONTINUE;  /* carry on with the program! */
}

SDL_AppResult SDL_AppIterate(void *appstate)
{
    SDL_Surface *indexedSurface = NULL;
    SDL_Texture *plotTexture = NULL;

    AppState_t *as = (AppState_t*)appstate;
    ProcessorState *state = as->state;
    if (state == NULL) {
        visualizeResults(NULL);
        goto updatedisplay;
    }
    if (as->show16Hz) {
        state->vars = &state->vars16hz;
    }
    else {
        state->vars = &state->vars2hz;
    }

    if (state->nVideoFrames == 0) {
        visualizeResults(state);
        goto updatedisplay;
    }

    // Handle playback
    if (as->playing) {
        double timeRange = state->plotT1 - state->plotT0;
        double deltaT = timeRange / 2000.0 / 1000.0 * as->playbackRate; // seconds
        if (as->playbackDirection > 0) {
            advancePlots(as, deltaT, SECONDS);
        }
        else {
            rewindPlots(as, deltaT, SECONDS);
        }
    }
    else {
        // Brief pause to yield CPU
        usleep(100000);
    }

updatedisplay:

    indexedSurface = SDL_CreateSurfaceFrom(IMAGE_WIDTH, IMAGE_HEIGHT, SDL_PIXELFORMAT_INDEX8, state->frames[as->plotPage].pixels, IMAGE_WIDTH);
    SDL_SetSurfacePalette(indexedSurface, as->colors);

    plotTexture = SDL_CreateTextureFromSurface(as->plotRenderer, indexedSurface);
    SDL_DestroySurface(indexedSurface);
    SDL_RenderTexture(as->plotRenderer, plotTexture, NULL, NULL);
    SDL_DestroyTexture(plotTexture);

    SDL_RenderPresent(as->plotRenderer);

    return SDL_APP_CONTINUE;
}

void SDL_AppQuit(void *appstate, SDL_AppResult result)
{
    AppState_t *as = (AppState_t*)appstate;
    if (as == NULL) {
        return;
    }

    ProcessorState *state = (ProcessorState*)as->state;
    shutdown(state);
    if (state->nVideoFrames > 0) {
        for (int i = 0; i < state->nVideoFrames; i++) {
            free(state->frames[i].pixels);
        }
        free(state->frames);
    }
    free(state);
    free(as->help);
    free(as);

    return;
}

void resetVideoFrames(ProcessorState *state)
{
    if (state != NULL) {
        for (int i = 0; i < state->nVideoFrames; i++) {
            free(state->frames[i].pixels);
        }
        free(state->frames);
        state->frames = NULL;
    }
    state->nVideoFrames = 0;

    return;
}

void resetDisplay(AppState_t *as)
{
    // Clear display
    SDL_SetRenderDrawColor(as->plotRenderer, BACKGROUND_COLOR, BACKGROUND_COLOR, BACKGROUND_COLOR, 255);
    SDL_RenderClear(as->plotRenderer);

    return;
}

double calculateDeltaT(AppState_t *as, TimeUnit_enum units, int sign)
{
    // To handle month and year intervals, get the date parts
    long year, month, day, hour, minute, second, msec;
    EPOCHbreakdown(as->t0, &year, &month, &day, &hour, &minute, &second, &msec);
    double tmpT0 = 0.0;
    double deltaT = 0.0;
    // This switch calculates deltaT in seconds
    // The function returns the value in ms
    switch (units) {
        case TIME_RANGES:
            deltaT = (as->state->plotT1 - as->state->plotT0) / 1000.0;
            break;
        case SAMPLE_PERIODS:
            deltaT = as->samplePeriodSeconds;
        case SECONDS:
            deltaT = 1.0;
            break;
        case MINUTES:
            deltaT = 60.0;
            break;
        case HOURS:
            deltaT = 3600.0;
            break;
        case DAYS:
            deltaT = 86400.0;
            break;
        case WEEKS:
            deltaT = 86400.0 * 7;
            break;
        case MONTHS:
            // Different result depending on the month and sign of the delta t
            // forward (+1) or backward (-1)
            tmpT0 = computeEPOCH(year, month + sign * 1, day, hour, minute, second, msec);
            deltaT = as->state->plotT0 - tmpT0;
            break;
        case YEARS:
            tmpT0 = computeEPOCH(year + sign * 1, month, day, hour, minute, second, msec);
            deltaT = as->state->plotT0 - tmpT0;
            break;
        case SEASONS:
            deltaT = 86400.0 * 365.25 / 4.0;
            break;
        case DECADES:
            tmpT0 = computeEPOCH(year + sign * 10, month, day, hour, minute, second, msec);
            deltaT = tmpT0 - as->state->plotT0;
            break;
        default:
            deltaT = 0.0;
            break;
    }

    // Caller knows the direction and handles it.
    // CDF time is in ms, multiply by 1000
    return fabs(deltaT * 1000.0);
}

void advancePlots(AppState_t *as, double amount, TimeUnit_enum units)
{
    double deltaT = calculateDeltaT(as, units, 1);
    double totalTime = deltaT * amount;
    double timeRange = as->state->plotT1 - as->state->plotT0;
    as->state->plotT0 += totalTime;
    as->state->plotT1 += totalTime;
    if (as->state->plotT1 > as->t1) {
        as->state->args.day++;
        rerunProcessor(as);
    }
    updatePlots(as->state);

    return;
}

void rewindPlots(AppState_t *as, double amount, TimeUnit_enum units)
{
    double deltaT = calculateDeltaT(as, units, -1);
    double totalTime = deltaT * amount;
    double timeRange = as->state->plotT1 - as->state->plotT0;
    as->state->plotT0 -= totalTime;
    as->state->plotT1 -= totalTime;
    // TODO check if we can process other days by updating date and calling runProcessor?
    // For now, limit to one day
    if (as->state->plotT0 < as->t0) {
        as->state->args.day--;
        rerunProcessor(as);

    }
    updatePlots(as->state);
    return;
}

void updatePlots(ProcessorState *state)
{
    resetVideoFrames(state);
    visualizeResults(state);

    return;
}

void rerunProcessor(AppState_t *appstate)
{
    int status = TIICT_OK;

    ProcessorState *state = appstate->state;
    shutdown(state);
    status = initProcessor(state);
    if (status != TIICT_OK) {
        shutdown(state);
        goto vis;
    }
    status = initLogFiles(state);
    if (status != TIICT_OK) {
        shutdown(state);
        goto vis;
    }
    status = loadCalData(state);
    if (status != TIICT_OK) {
        shutdown(state);
        goto vis;
    }
    status = loadLpCalData(state);
    if (status != TIICT_OK) {
        shutdown(state);
        goto vis;
    }
    calibrateFlows(state);
    calculateFields(state);

vis:
    if (appstate->show16Hz) {
        state->vars = &state->vars16hz;
    }
    else {
        // Downsample 16 Hz
        state->vars2hz.nRecs = state->vars16hz.nRecs / 8 + 1;
        int reallocstatus = reallocVariables(&state->vars2hz, state->vars2hz.nRecs);
        if (reallocstatus != TIICT_OK) {
            state->vars = &state->vars16hz;
            appstate->show16Hz = true;
        }
        double t0 = 0.0;
        long storageIndex = 0;
        bool downSampled = false;
        long n2HzSamples = 0;
        for (long i = 0; i <= state->vars16hz.nRecs;) // Time index advanced below
        {
            t0 = floor(state->vars16hz.timestamp[i] / 1000.0); // UT second reference
            for (uint8_t halfSecond = 0; halfSecond < 2; halfSecond ++)
            {
                downSampled = downSampleHalfSecond(state, &i, storageIndex, t0 + 0.5 * halfSecond, state->vars16hz.nRecs-1);
                if (downSampled)
                {
                    storageIndex++;
                    n2HzSamples++;
                }
            }

        }
        state->vars2hz.nRecs = n2HzSamples;
        state->vars = &state->vars2hz;
    }
    visualizeResults(state);

    return;
}
