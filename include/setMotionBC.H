// Access a particular boundary displacement patch
void setMotionBC
(
    dynamicTopoFvMesh& mesh,
    const label& index,
    const vectorField& dispField
)
{
    // Determine the kind of motion solver in use
    word solverType(mesh.dynamicMeshDict().lookup("solver"));

    //- Cell decomposition FEM motion solver
    if
    (
        (solverType == "laplaceCellDecomposition")
     || (solverType == "pseudoSolidCellDecomposition")
    )
    {
        // Boundary motion specified for the tetDecompositionMotionSolver
        tetPointVectorField& motionU = const_cast<tetPointVectorField&>
            (mesh.lookupObject<tetPointVectorField>("motionU"));

        // Assign boundary conditions to the motion solver
        motionU.boundaryField()[index] == 
            dispField/mesh.time().deltaT().value();
    }

    //- Face decomposition FEM motion solver
    if
    (
        (solverType == "laplaceFaceDecomposition")
     || (solverType == "pseudoSolidFaceDecomposition")
    )
    {
        // Boundary motion specified for the tetDecompositionMotionSolver
        tetPointVectorField& motionU = const_cast<tetPointVectorField&>
            (mesh.lookupObject<tetPointVectorField>("motionU"));

        // Assign boundary conditions to the motion solver

        // The face-decomposition solver includes points at face-centres,
        // thus point motion has to be interpolated to these points
        tetPolyPatchInterpolation interpolator
        (
            refCast<const faceTetPolyPatch>
            (
                motionU.boundaryField()[index].patch()
            )
        );

        motionU.boundaryField()[index] ==
            interpolator.pointToPointInterpolate
            (
                dispField/mesh.time().deltaT().value()
            );
    }

    //- Spring-based Laplacian motion solver
    if
    (
        (solverType == "springMotionSolver")
    )
    {
        // Boundary motion specified for the springMotionSolver
        pointField& refPoints = const_cast<pointField&>
            (mesh.lookupObject<pointField>("refPoints"));

        // Assign boundary conditions to the motion solver
        const labelList& meshPts = mesh.boundaryMesh()[index].meshPoints();
        forAll(meshPts,pointI)
        {
            refPoints[meshPts[pointI]] += dispField[pointI];
        }
    }
}

// Access the entire reference field
void setMotionBC
(
    dynamicTopoFvMesh& mesh,
    const vectorField& dispField
)
{
    // Determine the kind of motion solver in use
    word solverType(mesh.dynamicMeshDict().lookup("solver"));

    //- Spring-based Laplacian motion solver
    if
    (
        (solverType == "springMotionSolver")
    )
    {
        // Boundary motion specified for the springMotionSolver
        pointField& refPoints = const_cast<pointField&>
            (mesh.lookupObject<pointField>("refPoints"));

        // Assign boundary conditions to the motion solver
        refPoints = dispField;
    }
}