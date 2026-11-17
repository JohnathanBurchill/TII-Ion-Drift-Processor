/*

    TII Cross-Track Ion Drift Processor: settings.h

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

#ifndef SETTINGS_H
#define SETTINGS_H

#define SOFTWARE_VERSION "1.6"
#define SOFTWARE_VERSION_STRING "TIICT 2026-11-16"
#define SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING 300 // 5 minutes
#define SECONDS_OF_DATA_REQUIRED_FOR_EXPORTING 60 // 1 minutes
#define SECONDS_OF_BOUNDARY_DATA_REQUIRED_FOR_PROCESSING 1500 // 25 minutes to ensure coverage of each fit region

#define FLAGS_MAXIMUM_DRIFT_VALUE 8000.0 // 8 km/s maximum drift for flagging
#define FLAGS_MINIMUM_ION_DENSITY_VALUE 10000.0 // cm^-3
#define DEFAULT_VI_ERROR -42.0
#define MAX_ALLOWED_CDF_GAP_SECONDS 600.0 // CDF export split into separate files at gaps exceeding 10 minutes
#define NUM_CAL_VARIABLES 22
#define NUM_EXPORT_VARIABLES 58
#define NUM_BUFFER_VARIABLES 58

// Quality_flags bits 4 and 5. As for bits 0 to 3, 1 means good.
#define QUALITY_FLAG_ATTITUDE_OK (1 << 4) // | |q_NEC_CRF| - 1 | <= QUATERNION_NORM_TOLERANCE
#define QUALITY_FLAG_B_OK (1 << 5) // B-VFM-Flag is 0 and B-Avg-Flag is less than 2
#define QUATERNION_NORM_TOLERANCE 1.0e-3
#define EARTH_ROTATION_RATE 7.292115e-5 // rad/s

#define MAX_NUMBER_OF_FITINFO_BITS_PER_COMPONENT 8
#define GSL_FIT_MAXIMUM_ITERATIONS 500

#define LP_MIN_NUMBER_OF_POTENTIALS 160000


#endif // SETTINGS_H
