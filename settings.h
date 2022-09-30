/*

    TII Cross-Track Ion Drift Processor: settings.h

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

#ifndef SETTINGS_H
#define SETTINGS_H

#define SOFTWARE_VERSION "1.1"
#define SOFTWARE_VERSION_STRING "TIICT 2022-09-28"
#define SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING 300 // 5 minutes
#define SECONDS_OF_DATA_REQUIRED_FOR_EXPORTING 60 // 1 minutes
#define SECONDS_OF_BOUNDARY_DATA_REQUIRED_FOR_PROCESSING 1500 // 25 minutes to ensure coverage of each fit region

#define FLAGS_MAXIMUM_DRIFT_VALUE 8000.0 // 8 km/s maximum drift for flagging
#define DEFAULT_VI_ERROR -42.0
#define MAX_ALLOWED_CDF_GAP_SECONDS 600.0 // CDF export split into separate files at gaps exceeding 10 minutes
#define NUM_CAL_VARIABLES 15
#define NUM_EXPORT_VARIABLES 32
#define NUM_BUFFER_VARIABLES 31

#define MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT 8
#define GSL_FIT_MAXIMUM_ITERATIONS 500

#define LP_MIN_NUMBER_OF_POTENTIALS 1600000


#endif // SETTINGS_H
