/*

    TII Cross-Track Ion Drift Processor: frames.h

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

#ifndef FRAMES_H
#define FRAMES_H

// Coordinate frames:
//   NEC    north, east, centre
//   CRF    Swarm common reference frame (spacecraft optical bench)
//   VFM    vector field magnetometer frame
//   sensor TII H-detector frame: x forward, y to the right looking forward, z nadir.
//          The V sensor's moments are mapped into this frame as for earlier versions.
//          The ion drifts are measured in this frame.
//   satellite-track  x parallel to the NEC satellite velocity, y horizontal to the
//          right looking forward, z = x cross y. Defined by the velocity, it differs
//          from the sensor frame by the attitude control deadband (about +-4 deg in yaw).
//          B, E and IGRF are exported in this frame, as in 0401.

// Rotation from CRF to the TII sensor frame: the identity. EfiCal defines the
// H detector frame as CRF (Rotate_from_CRF_to_Det_HV in Imp_L1_Bin.cs), and the
// processor maps the V sensor's moments into the H frame.
extern const double R_CRF_TO_SENSOR[3][3];

double quaternionNorm(const double q[4]);

// Rotation matrix for a Swarm quaternion q_A_B (q[0..2] vector part, q[3] scalar part).
// The quaternion is normalized first. The resulting matrix takes components in frame A
// to frame B: v_CRF = R(q_NEC_CRF) v_NEC, and v_VFM = R(q_STR_VFM) v_CRF.
// This direction was checked against the 0501 cal file for Swarm A on 2015-10-18,
// where R(q_NEC_CRF) Sat_Vel_N matches Det_H__vX, and R(q_STR_VFM)^T IGRF-VFM
// followed by R(q_NEC_CRF)^T gives a physical IGRF in NEC.
void quaternionToMatrix(const double q[4], double r[3][3]);

void rotate(const double r[3][3], const double in[3], double out[3]);
void rotateTransposed(const double r[3][3], const double in[3], double out[3]);

void necToSensor(const double rNecCrf[3][3], const double vNec[3], double vSensor[3]);
void sensorToNec(const double rNecCrf[3][3], const double vSensor[3], double vNec[3]);
void vfmToNec(const double rStrVfm[3][3], const double rNecCrf[3][3], const double vVfm[3], double vNec[3]);

// Converts v from WGS84 geodetic (north, east, down) to geocentric NEC components,
// in place, at geocentric latitude latGeocentricDeg (degrees) and radius r (m).
// The cal-file Sat_Vel_N has a geodetic third component (up to 24 m/s from geocentric).
void geodeticToGeocentricNec(double latGeocentricDeg, double r, double v[3]);

// Rows of axes are the satellite-track x, y and z unit vectors in NEC, so
// rotate(axes, vNec, vTrack) and rotateTransposed(axes, vTrack, vNec) convert between the two.
void satelliteTrackAxes(const double vSatNec[3], double axes[3][3]);

#endif // FRAMES_H
