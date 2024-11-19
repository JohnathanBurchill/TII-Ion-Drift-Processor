/*

    TII Cross-Track Ion Drift Processor: tiictinteractive.c

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

#include "SDL3/SDL_init.h"
#include "processing.h"
#include "state.h"
#include "errors.h"
#include "tiigraphics/tiigraphics.h"


#include <tiigraphics/video.h>

#include <stdio.h>

#define SDL_MAIN_USE_CALLBACKS 1  /* use the callbacks instead of main() */
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

static SDL_Window *window = NULL;
static SDL_Surface *display = NULL;
static SDL_Renderer *renderer = NULL;

static SDL_Surface *frames = NULL;
static SDL_Palette *colors = NULL;

void resetVideoFrames(ProcessorState *state);
void resetDisplay(void);

typedef struct AppState {
    ProcessorState *results;
} AppState_t;

SDL_AppResult SDL_AppInit(void **appstate, int argc, char *argv[])
{

    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_Log("Couldn't initialize SDL: %s", SDL_GetError());
        return SDL_APP_FAILURE;
    }

    if (!SDL_CreateWindowAndRenderer("examples/renderer/clear", IMAGE_WIDTH, IMAGE_HEIGHT, 0, &window, &renderer)) {
        SDL_Log("Couldn't create window/renderer: %s", SDL_GetError());
        return SDL_APP_FAILURE;
    }
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

    frames = SDL_CreateSurface(8, 65, SDL_PIXELFORMAT_INDEX8);
    colors = SDL_CreateSurfacePalette(frames);
    for (int i = 0; i < colors->ncolors; i++) {
        colors->colors[i].r = i;
        colors->colors[i].g = 0;
        colors->colors[i].b = 0;
        colors->colors[i].a = 0;
    }

    display = SDL_GetWindowSurface(window);
    resetDisplay();

    ProcessorState *results = NULL;
    int status = runProcessor(argc, argv, &results);
    if (status != TIICT_OK) {
        return SDL_APP_FAILURE;
    }

    AppState_t *as = malloc(sizeof *as);
    as->results = results;
    *appstate = (void*)as;

    return SDL_APP_CONTINUE;  /* carry on with the program! */
}

SDL_AppResult SDL_AppEvent(void *appstate, SDL_Event *event)
{
    int status = 0;
    AppState_t *as = (AppState_t *)appstate;
    ProcessorState *state = (ProcessorState *)as->results;
    int argc = state->args.argc;
    char **argv = state->args.argv;

    if (event->type == SDL_EVENT_QUIT) {
        return SDL_APP_SUCCESS;  /* end the program, reporting success to the OS. */
    }
    if (event->type == SDL_EVENT_KEY_UP) {
        switch (event->key.key) {
            case SDLK_R:
                resetDisplay();
                break;
            case SDLK_U:
                // Update processor results
                resetVideoFrames(state);
                shutdown(state);
                free(state);
                state = NULL;
                status = runProcessor(argc, argv, &state);
                as->results = state;
                if (status != TIICT_OK) {
                    fprintf(stderr, "Encountered error running processor: %d\n", status);
                }
                else {
                    printf("Ran processor; number of images: %d\n", state->nVideoFrames);
                }
                break;
// Function keys example...
//            case SDLK_F4:
//                break;
            case SDLK_Q:
                return SDL_APP_SUCCESS;  /* end the program, reporting success to the OS. */
                break;
            default:
                break;
        }

    }

    return SDL_APP_CONTINUE;  /* carry on with the program! */
}

SDL_AppResult SDL_AppIterate(void *appstate)
{
//    SDL_Rect imageDisplayRegion = {0, MAIN_WINDOW_HEIGHT/2 - IMAGE_DISPLAY_HEIGHT/2, IMAGE_DISPLAY_WIDTH, IMAGE_DISPLAY_HEIGHT};
//    SDL_BlitSurfaceScaled(frames, NULL, display, &imageDisplayRegion, SDL_SCALEMODE_NEAREST);
//    SDL_Texture *frameTexture = SDL_CreateTextureFromSurface(renderer, display);
//
//    // Draw simulated image
//    const SDL_FRect imageDisplayRegionF = {0, MAIN_WINDOW_HEIGHT/2 - IMAGE_DISPLAY_HEIGHT/2, IMAGE_DISPLAY_WIDTH, IMAGE_DISPLAY_HEIGHT};
//    SDL_RenderTexture(renderer, frameTexture, &imageDisplayRegionF, &imageDisplayRegionF);
//
//    SDL_RenderPresent(renderer);
//
//    SDL_DestroyTexture(frameTexture);

    return SDL_APP_CONTINUE;
}

void SDL_AppQuit(void *appstate, SDL_AppResult result)
{
    AppState_t *as = (AppState_t*)appstate;
    ProcessorState *state = (ProcessorState*)as->results;
    shutdown(state);
    if (state->nVideoFrames > 0) {
        free(state->frames);
    }
    free(state);

    return;
}

void resetVideoFrames(ProcessorState *state)
{
    if (state != NULL) {
        free(state->frames);
        state->frames = NULL;
    }
    state->nVideoFrames = 0;

    return;
}

void resetDisplay(void)
{
    // Clear display
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    return;
}

