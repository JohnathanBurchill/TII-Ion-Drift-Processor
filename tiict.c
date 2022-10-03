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

extern char infoHeader[50];

int main(int argc, char* argv[])
{
    // Can use this function in a library; main is a wrapper for command line execution
    int status = runProcessor(argc, argv);

    return status;
}

