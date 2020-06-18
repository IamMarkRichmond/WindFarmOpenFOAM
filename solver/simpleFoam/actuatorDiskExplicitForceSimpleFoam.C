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
    actuatorDiskExplicitForceSimpleFoam

Description
    Steady-state actuator wake solver for incompressible, turbulent flow, using the SIMPLE
    algorithm.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "actuatorDiskExplicitForce.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    actuatorDiskExplicitForce actuatorDisk1;
    actuatorDisk1.SetRef(1);
    actuatorDisk1.ReadGeometry(mesh);

  //  actuatorDisk1.WriteVTK();

    actuatorDiskExplicitForce actuatorDisk2;
    actuatorDisk2.SetRef(2);
    actuatorDisk2.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk3;
    actuatorDisk3.SetRef(3);
    actuatorDisk3.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk4;
    actuatorDisk4.SetRef(4);
    actuatorDisk4.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk5;
    actuatorDisk5.SetRef(5);
    actuatorDisk5.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk6;
    actuatorDisk6.SetRef(6);
    actuatorDisk6.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk7;
    actuatorDisk7.SetRef(7);
    actuatorDisk7.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk8;
    actuatorDisk8.SetRef(8);
    actuatorDisk8.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk9;
    actuatorDisk9.SetRef(9);
    actuatorDisk9.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk10;
    actuatorDisk10.SetRef(10);
    actuatorDisk10.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk11;
    actuatorDisk11.SetRef(11);
    actuatorDisk11.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk12;
    actuatorDisk12.SetRef(12);
    actuatorDisk12.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk13;
    actuatorDisk13.SetRef(13);
    actuatorDisk13.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk14;
    actuatorDisk14.SetRef(14);
    actuatorDisk14.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk15;
    actuatorDisk15.SetRef(15);
    actuatorDisk15.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk16;
    actuatorDisk16.SetRef(16);
    actuatorDisk16.ReadGeometry(mesh);

    actuatorDisk16.WriteVTK();

    actuatorDiskExplicitForce actuatorDisk17;
    actuatorDisk17.SetRef(17);
    actuatorDisk17.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk18;
    actuatorDisk18.SetRef(18);
    actuatorDisk18.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk19;
    actuatorDisk19.SetRef(19);
    actuatorDisk19.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk20;
    actuatorDisk20.SetRef(20);
    actuatorDisk20.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk21;
    actuatorDisk21.SetRef(21);
    actuatorDisk21.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk22;
    actuatorDisk22.SetRef(22);
    actuatorDisk22.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk23;
    actuatorDisk23.SetRef(23);
    actuatorDisk23.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk24;
    actuatorDisk24.SetRef(24);
    actuatorDisk24.ReadGeometry(mesh);

    actuatorDiskExplicitForce actuatorDisk25;
    actuatorDisk25.SetRef(25);
    actuatorDisk25.ReadGeometry(mesh);

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        laminarTransport.correct();
        turbulence->correct();
        
        // Compute vorticity and second-invariant of velocity gradient tensor.
       // omega = fvc::curl(U);
        Q = 0.5*(sqr(tr(fvc::grad(U))) - tr(((fvc::grad(U))&(fvc::grad(U)))));

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
