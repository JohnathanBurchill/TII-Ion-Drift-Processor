/*

    TII Cross-Track Ion Drift Processor: tiict.c

    Copyright (C) 2022  Johnathan K Burchill

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

#include "tiict.h"

#include "processing.h"
#include "loadData.h"
#include "export.h"
#include "errors.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gsl/gsl_errno.h>

char infoHeader[50];

int main(int argc, char* argv[])
{
    int status = TIICT_OK;

    ProcessorState state = {0};

    checkResult(initProcessor(argc, argv, &state), &state);

    checkResult(loadTiiCalData(&state), &state);

    checkResult(loadLpCalData(&state), &state);

    checkResult(calibrateFlows(&state), &state);

    checkResult(calculateFields(&state), &state);

    checkResult(exportCdfs(&state), &state);
 
    shutdown(&state, EXIT_SUCCESS);

    return 0;
}

int initProcessor(int argc, char **argv, ProcessorState *state)
{
    int status = TIICT_OK;
    void *args = &state->args;

    // Check arguments and abort if not right
    status = parseArguments(argc, argv, args);
    if (status != TIICT_OK)
        return status;

    // Prefix for messages
    initHeader(args);

    // Confirm requested date has records. Abort otherwise.
    status = checkCalDataAvailability(state);
    if (status != TIICT_OK)
        return status;

    // Ensure export directories exist, or abort.
    status = initDirectories(args);
    if (status != TIICT_OK)
        return status;

    status = initFitFiles(state);
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

int parseArguments(int argc, char **argv, Arguments *args)
{
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--about") == 0)
        {
            fprintf(stdout, "tiict - TII Cross-track ion drift processor, version %s.\n", SOFTWARE_VERSION);
            fprintf(stdout, "Copyright (C) 2022  Johnathan K Burchill\n");
            fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
            fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
            fprintf(stdout, "under the terms of the GNU General Public License.\n");
            return TIICT_ARGS_ABOUT;
        }
    }

    if (argc != 10)
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

void initHeader(Arguments *args)
{
    time_t currentTime;
    struct tm * timeParts;
    time(&currentTime);
    timeParts = localtime(&currentTime);

    // set up info header
    sprintf(infoHeader, "TIICT %c%s %04d-%02d-%02d: ", args->satellite[0], args->exportVersion, args->year, args->month, args->day);
    fprintf(stdout, "\n%s-------------------------------------------------\n", infoHeader);
    fprintf(stdout, "%sVersion 0302 20220519\n", infoHeader);
    fprintf(stdout, "%sProcessing date: %s\n", infoHeader, asctime(timeParts));

    return;
}

void checkResult(int status, ProcessorState *state)
{
    if (status != TIICT_OK)
        shutdown(state, status);
}

void shutdown(ProcessorState *state, int exitStatus)
{

    // Close fit log file
    if (state->fitFile != NULL)
        fclose(state->fitFile);
    fflush(stdout);

    // Free the memory
    for (uint8_t i = 0; i < NUM_CAL_VARIABLES; i++)
        if (state->dataBuffers[i] != NULL)
            free(state->dataBuffers[i]);

    if (state->lpTimes != NULL)
        free(state->lpTimes);
    if (state->lpPhiScHighGain != NULL)
        free(state->lpPhiScHighGain);
    if (state->lpPhiScLowGain != NULL)
        free(state->lpPhiScLowGain);
    if (state->lpPhiSc != NULL)
        free(state->lpPhiSc);
    // state->potentials is just a pointer to one of the above potentials

    if (state->xhat != NULL)
        free(state->xhat);
    if (state->yhat != NULL)
        free(state->yhat);
    if (state->zhat != NULL)
        free(state->zhat);
    if (state->ectFieldH != NULL)
        free(state->ectFieldH);
    if (state->ectFieldV != NULL)
        free(state->ectFieldV);
    if (state->bctField != NULL)
        free(state->bctField);
    if (state->viErrors != NULL)
        free(state->viErrors);
    if (state->flags != NULL)
        free(state->flags);
    if (state->fitInfo != NULL)
        free(state->fitInfo);

    exit(exitStatus);    
}
