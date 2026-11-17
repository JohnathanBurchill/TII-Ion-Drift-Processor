/*

    TII Cross-Track Ion Drift Processor: frames.c

    Copyright (C) 2026  Johnathan K Burchill

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

#include "frames.h"

#include <math.h>

const double R_CRF_TO_SENSOR[3][3] = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}
};

double quaternionNorm(const double q[4])
{
    return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
}

void quaternionToMatrix(const double q[4], double r[3][3])
{
    double n = quaternionNorm(q);
    double x = q[0] / n;
    double y = q[1] / n;
    double z = q[2] / n;
    double w = q[3] / n;

    r[0][0] = 1.0 - 2.0 * (y*y + z*z);
    r[0][1] = 2.0 * (x*y - z*w);
    r[0][2] = 2.0 * (x*z + y*w);
    r[1][0] = 2.0 * (x*y + z*w);
    r[1][1] = 1.0 - 2.0 * (x*x + z*z);
    r[1][2] = 2.0 * (y*z - x*w);
    r[2][0] = 2.0 * (x*z - y*w);
    r[2][1] = 2.0 * (y*z + x*w);
    r[2][2] = 1.0 - 2.0 * (x*x + y*y);

    return;
}

void rotate(const double r[3][3], const double in[3], double out[3])
{
    for (int i = 0; i < 3; i++)
        out[i] = r[i][0] * in[0] + r[i][1] * in[1] + r[i][2] * in[2];

    return;
}

void rotateTransposed(const double r[3][3], const double in[3], double out[3])
{
    for (int i = 0; i < 3; i++)
        out[i] = r[0][i] * in[0] + r[1][i] * in[1] + r[2][i] * in[2];

    return;
}

void necToSensor(const double rNecCrf[3][3], const double vNec[3], double vSensor[3])
{
    double vCrf[3];
    rotate(rNecCrf, vNec, vCrf);
    rotate(R_CRF_TO_SENSOR, vCrf, vSensor);

    return;
}

void sensorToNec(const double rNecCrf[3][3], const double vSensor[3], double vNec[3])
{
    double vCrf[3];
    rotateTransposed(R_CRF_TO_SENSOR, vSensor, vCrf);
    rotateTransposed(rNecCrf, vCrf, vNec);

    return;
}

void vfmToNec(const double rStrVfm[3][3], const double rNecCrf[3][3], const double vVfm[3], double vNec[3])
{
    double vCrf[3];
    rotateTransposed(rStrVfm, vVfm, vCrf);
    rotateTransposed(rNecCrf, vCrf, vNec);

    return;
}

// WGS84 geodetic latitude (radians) of a point at geocentric latitude latGeocentric (radians)
// and distance r (m) from the centre of the Earth
static double geodeticLatitude(double latGeocentric, double r)
{
    const double a = 6378137.0;
    const double f = 1.0 / 298.257223563;
    const double e2 = f * (2.0 - f);
    double p = r * cos(latGeocentric);
    double z = r * sin(latGeocentric);
    double lat = atan2(z, p * (1.0 - e2));
    for (int k = 0; k < 5; k++)
    {
        double s = sin(lat);
        double n = a / sqrt(1.0 - e2 * s * s);
        double h = p / cos(lat) - n;
        lat = atan2(z, p * (1.0 - e2 * n / (n + h)));
    }

    return lat;
}

void geodeticToGeocentricNec(double latGeocentricDeg, double r, double v[3])
{
    double latGeocentric = latGeocentricDeg * M_PI / 180.0;
    // The geodetic axes are the geocentric axes rotated about east by this angle
    double d = geodeticLatitude(latGeocentric, r) - latGeocentric;
    double north = v[0] * cos(d) - v[2] * sin(d);
    double centre = v[0] * sin(d) + v[2] * cos(d);
    v[0] = north;
    v[2] = centre;

    return;
}

void satelliteTrackAxes(const double vSatNec[3], double axes[3][3])
{
    double m = sqrt(vSatNec[0]*vSatNec[0] + vSatNec[1]*vSatNec[1] + vSatNec[2]*vSatNec[2]);
    for (int k = 0; k < 3; k++)
        axes[0][k] = vSatNec[k] / m;
    // y is x cross (0, 0, -1): horizontal, to the right looking forward
    double h = sqrt(axes[0][0]*axes[0][0] + axes[0][1]*axes[0][1]);
    axes[1][0] = -axes[0][1] / h;
    axes[1][1] = axes[0][0] / h;
    axes[1][2] = 0.0;
    // z is x cross y
    axes[2][0] = axes[0][1] * axes[1][2] - axes[0][2] * axes[1][1];
    axes[2][1] = axes[0][2] * axes[1][0] - axes[0][0] * axes[1][2];
    axes[2][2] = axes[0][0] * axes[1][1] - axes[0][1] * axes[1][0];

    return;
}
