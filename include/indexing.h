/*

    TII Cross-Track Ion Drift Processor: indexing.h

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

#ifndef INDEXING_H
#define INDEXING_H

#define ADDR(n, m, d) (((float*)dataBuffers[(n)]+(d*timeIndex + m)))
#define MEAS(n, m, d) ((float)(*(ADDR(n, m, d))))
#define TIME() ((double)*((double*)dataBuffers[0]+(timeIndex)))
#define MXH() (MEAS(1, 0, 2))
#define MYH() (MEAS(1, 1, 2))
#define MXV() (MEAS(2, 0, 2))
#define MYV() (MEAS(2, 1, 2))
#define VSATX() (MEAS(3, 0, 3))
#define VSATY() (MEAS(3, 1, 3))
#define VSATZ() (MEAS(3, 2, 3))
#define MLT() (MEAS(4, 0, 1))
#define QDLAT() (MEAS(5, 0, 1))
#define QDLON() (MEAS(6, 0, 1))
#define LAT() (MEAS(7, 0, 1))
#define LON() (MEAS(8, 0, 1))
#define RADIUS() (MEAS(9, 0, 1))
#define VCORX() (MEAS(10, 0, 3))
#define VCORY() (MEAS(10, 1, 3))
#define VCORZ() (MEAS(10, 2, 3))
#define VSATN() (MEAS(11, 0, 3))
#define VSATE() (MEAS(11, 1, 3))
#define VSATC() (MEAS(11, 2, 3))
#define BN() (MEAS(12, 0, 3))
#define BE() (MEAS(12, 1, 3))
#define BC() (MEAS(12, 2, 3))
#define VBIAS() (MEAS(13, 0, 1))
#define VFP() (MEAS(14, 0, 1))

#endif // INDEXING_H
