#include "tetPolyPatchFields.H"
#include "tetPointFields.H"
#include "tetFemMatrices.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"

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

    //- Cell decomposition FEM motion solver
    if
    (
        (solverType == "laplaceCellDecomposition")
     || (solverType == "pseudoSolidCellDecomposition")
    )
    {
        // Boundary motion specified for the tetDecompositionMotionSolver
        tetPointVectorField& motionU = const_cast<tetPointVectorField&>
        (
            mesh.lookupObject<tetPointVectorField>("motionU")
        );

        // Assign boundary conditions to the motion solver
        motionU.boundaryField()[index] ==
        (
            dispField/mesh.time().deltaT().value()
        );
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
        (
            mesh.lookupObject<tetPointVectorField>("motionU")
        );

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
        (
            interpolator.pointToPointInterpolate
            (
                dispField/mesh.time().deltaT().value()
            )
        );
    }

    if (solverType == "mesquiteSmoother" || solverType == "springMotionSolver")
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

    //- Cell decomposition FEM motion solver
    if
    (
        (solverType == "laplaceCellDecomposition")
     || (solverType == "pseudoSolidCellDecomposition")
    )
    {
        // Boundary motion specified for the tetDecompositionMotionSolver
        tetPointVectorField& motionU = const_cast<tetPointVectorField&>
        (
            mesh.lookupObject<tetPointVectorField>("motionU")
        );

        const polyBoundaryMesh& boundary = mesh.boundaryMesh();

        const pointField& oldPoints = mesh.points();

        forAll(boundary, patchI)
        {
            // Assign boundary conditions to the motion solver
            const labelList& meshPts = mesh.boundaryMesh()[patchI].meshPoints();

            vectorField pField(meshPts.size(), vector::zero);

            forAll(meshPts, pointI)
            {
                pField[pointI] =
                (
                    dispField[meshPts[pointI]] - oldPoints[meshPts[pointI]]
                );
            }

            // Assign boundary conditions to the motion solver
            motionU.boundaryField()[patchI] ==
            (
                pField/mesh.time().deltaT().value()
            );
        }
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
        (
            mesh.lookupObject<tetPointVectorField>("motionU")
        );

        // Assign boundary conditions to the motion solver
        const polyBoundaryMesh& boundary = mesh.boundaryMesh();

        const pointField& oldPoints = mesh.points();

        forAll(boundary, patchI)
        {
            // Assign boundary conditions to the motion solver
            const labelList& meshPts = mesh.boundaryMesh()[patchI].meshPoints();

            vectorField pField(meshPts.size(), vector::zero);

            forAll(meshPts, pointI)
            {
                pField[pointI] =
                (
                    dispField[meshPts[pointI]] - oldPoints[meshPts[pointI]]
                );
            }

            // The face-decomposition solver includes points at face-centres,
            // thus point motion has to be interpolated to these points
            tetPolyPatchInterpolation interpolator
            (
                refCast<const faceTetPolyPatch>
                (
                    motionU.boundaryField()[patchI].patch()
                )
            );

            // Assign boundary conditions to the motion solver
            motionU.boundaryField()[patchI] ==
            (
                interpolator.pointToPointInterpolate
                (
                    pField/mesh.time().deltaT().value()
                )
            );
        }
    }

    if (solverType == "mesquiteSmoother" || solverType == "springMotionSolver")
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

