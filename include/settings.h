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

#define SOFTWARE_VERSION "1.2"
#define SOFTWARE_VERSION_STRING "TIICT 2022-10-07"
#define SECONDS_OF_DATA_REQUIRED_FOR_PROCESSING 300 // 5 minutes
#define SECONDS_OF_DATA_REQUIRED_FOR_EXPORTING 60 // 1 minutes
#define SECONDS_OF_BOUNDARY_DATA_REQUIRED_FOR_PROCESSING 1500 // 25 minutes to ensure coverage of each fit region

#define ELECTRIC_CHARGE_OPLUS 1.602e-19
#define MASS_OPLUS 2.67e-26

#define FLAGS_MAXIMUM_DRIFT_VALUE_HIGH_LATITUDE 8000.0 // 8 km/s maximum drift for flagging
#define FLAGS_MAXIMUM_DRIFT_VALUE_LOW_LATITUDE 500.0 // 0.5 km/s maximum drift for flagging
#define FLAGS_MINIMUM_ION_DENSITY 10000.0
#define DEFAULT_VI_ERROR -42.0
#define MAX_ALLOWED_CDF_GAP_SECONDS 600.0 // CDF export split into separate files at gaps exceeding 10 minutes

#define NUM_CAL_VARIABLES 15
#define NUM_TRACIS_VARIABLES 3
#define NUM_EXPORT_VARIABLES 32
#define NUM_BUFFER_VARIABLES 31

#define MAX_NUMBER_OF_CALIBRATION_FLAGS_BITS_PER_COMPONENT 8
#define GSL_FIT_MAXIMUM_ITERATIONS 100
#define MIN_POINTS_PER_REGION_FOR_FIT 200
#define FIT_MAX_MAD 300.0 // m/s
#define FIT_MAX_ABS_MEDIAN 5000.0 // m/s
#define FIT_OUTLIER_MAD_THRESHOLD 4.0
#define FIT_MODEL_DEGREE 3 // Quadratic offset model
#define FIT_P_QD_1 50 // High-latitude equatorward boundary
#define FIT_P_QD_2 55 // High-latitude poleward boundary
#define FIT_E_QD_1 50 // Low-latitude poleward boundary
#define FIT_E_QD_2 45 // Low-latitude equatorward boundary

#define LP_MIN_NUMBER_OF_POTENTIALS 160000

// From TRACIS
#define CDF_GZIP_COMPRESSION_LEVEL 6L
#define CDF_BLOCKING_FACTOR 43200L


#endif // SETTINGS_H
