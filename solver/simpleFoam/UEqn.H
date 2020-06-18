    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

	//Added by me
    //Clear  old  values  in  ioVolumeForce
    for(label k = 0; k < VolumeForce.size(); k++) {
        VolumeForce[k] = vector::zero;
    }

    // Calculate  volume  force  from  actuator  disk
    actuatorDisk1.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk2.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk3.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk4.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk5.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk6.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk7.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk8.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk9.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk10.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk11.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk12.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk13.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk14.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk15.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk16.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk17.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk18.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk19.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk20.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk21.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk22.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk23.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk24.CalcActuatorDiskVolForce(mesh, VolumeForce, U);
    actuatorDisk25.CalcActuatorDiskVolForce(mesh, VolumeForce, U);


    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U) + VolumeForce
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }


