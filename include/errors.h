/*

    TII Cross-Track Ion Drift Processor: errors.h

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

#ifndef _TIICT_ERRORS_H
#define _TIICT_ERRORS_H

enum TIICT_ERRORS {

    TIICT_OK = 0,
    TIICT_MEMORY,
    TIICT_CDF_READ,
    TIICT_CDF_WRITE,
    TIICT_ARGS_ABOUT,
    TIICT_ARGS_BAD,
    TIICT_ARGS_SATELLITE,
    TIICT_NOT_ENOUGH_CALIBRATION_RECORDS,
    TIICT_NOT_ENOUGH_TRACIS_RECORDS,
    TIICT_EXPORT_DIRECTORY_TCT16,
    TIICT_EXPORT_DIRECTORY_TCT02,
    TIICT_DIRECTORY_READ,
    TIICT_NO_LP_HM_DATA,
    TIICT_ZIP_EXISTS,
    TIICT_LOG_WRITE,
    TIICT_SHELL,
    TIICT_ZIP,
    TIICT_NO_RECORDS_TO_EXPORT,
    TIICT_NO_CAL_FILE
};

#endif // _TIICT_ERRORS_H
