
// Access a particular boundary displacement patch
void setMotionBC
(
    fvMesh& mesh,
    const label& index,
    const vectorField& dispField
)
{
    word solverType
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("solver")
    );

    if (solverType == "mesquiteMotionSolver")
    {
        // Boundary motion specified for the springMotionSolver
        pointField& refPoints = const_cast<pointField&>
        (
            mesh.lookupObject<pointField>("refPoints")
        );

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
    fvMesh& mesh,
    const vectorField& dispField
)
{
    word solverType
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("solver")
    );

    if (solverType == "mesquiteMotionSolver")
    {
        // Boundary motion specified for the springMotionSolver
        pointField& refPoints = const_cast<pointField&>
        (
            mesh.lookupObject<pointField>("refPoints")
        );

        // Assign boundary conditions to the motion solver
        refPoints = dispField;
    }
}

