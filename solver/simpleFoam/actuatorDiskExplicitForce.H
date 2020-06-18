/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    actuatorDiskExplicitForce

Description
    Add volume force in an actuator disk region from thrust, torque and geometry
    defined in fvSolution

    Use with actuatorDiskExplicitForceSimpleFoam

    Modified by Xingxing HAN, December 2016, from the origin work of Erik Svennig
    Reference: Implementation of an actuator disk in OpenFOAM

\*---------------------------------------------------------------------------*/

#ifndef ACTUATORDISKEXPLICITFORCE_H_
#define ACTUATORDISKEXPLICITFORCE_H_

#include "fvCFD.H"

namespace Foam {

  class actuatorDiskExplicitForce {
  public:

    // Runtime type information 
    TypeName("actuatorDiskExplicitForce");

      actuatorDiskExplicitForce();
      ~actuatorDiskExplicitForce();

      void SetRef(int ref);
      void ReadGeometry(const fvMesh &iMesh);
      void CalcActuatorDiskVolForce(const fvMesh &iMesh, vectorField &ioVolumeForce, vectorField &ioU);

      void WriteVTK();

  private:
      vector mPointStartCenterLine;
      vector mPointEndCenterLine;
      vector mCenterPoint;
      scalar mExtRadius, mIntRadius;
      scalar mThrust, mTorque;
      scalar mRho;
      scalar ref;
      scalar cosT;
      scalar sinT;
      scalar Thickness;
      vector Uin;
      scalar Umag;
      scalar Ud;
      scalar Uinf;
      scalar Ucount;
      scalar Usum;
      scalar a1;
      scalar a2;
      scalar B;
      scalar counta;
      scalar signofUx;

      static constexpr double mPI = 3.141592654;

      bool PointIsInLookUp(const vector &iPointStartCenterLine, const vector &
	  iPointEndCenterLine, const vector &iPoint, scalar &oDist2, vector &
	  oLineTangent, vector &oCircumferentialDirection);

      bool PointIsInDisk(const vector &iPointStartCenterLine, const vector &
          iPointEndCenterLine, const vector &iPoint, scalar &oDist2, vector &
          oLineTangent, vector &oCircumferentialDirection);

      bool PointIsInHub(const vector &iPointStartCenterLine, const vector &
          iPointEndCenterLine, const vector &iPoint);

      scalar CalcAxialForce(const scalar &iRadialDist, const scalar &iRho);
      scalar CalcCircForce(const scalar &iRadialDist, const scalar &iRho);
      scalar CalcDiskThickness() { return mag(mPointEndCenterLine - mPointStartCenterLine); };
  };

} // end namespace Foam
#endif /* ACTUATORDISKEXPLICITFORCE_H_ */
