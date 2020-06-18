    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, simple.dict(), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    #include "createMRF.H"

    Info << "Create volume force field.\n";

    volVectorField VolumeForce
    (
        IOobject
        (
            "VolumeForce",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//    Info << "Creating vorticity field, omega..." << endl;
//    volScalarField omega
//    (
//        IOobject
//        (
//            "omega",
//            runTime.timeName(),
//            mesh,
//            IOobject::READ_IF_PRESENT,
//            IOobject::AUTO_WRITE
//        ),
//        mesh,
//        dimensionedScalar("omega",dimensionSet(0, 0, -1, 0, 0, 0, 0),0,0)
//    );

    Info << "Creating second-invariant of strain-rate tensor field, Q..." << endl;
    volScalarField Q
    (
        IOobject
        (
            "Q",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Q",dimensionSet(0, 0, -2, 0, 0, 0, 0),0.0)
    );
