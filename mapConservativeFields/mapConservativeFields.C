/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    mapConservativeFields

Description
    Maps volume fields conservatively from one mesh to another, reading and
    interpolating all fields present in the time directory of both cases.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOobjectList.H"
#include "conservativeMeshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int getTimeIndex
(
    const instantList& times,
    const scalar t
)
{
    int nearestIndex = -1;
    scalar nearestDiff = Foam::GREAT;

    forAll(times, timeIndex)
    {
        if (times[timeIndex].name() == "constant") continue;

        scalar diff = fabs(times[timeIndex].value() - t);
        if (diff < nearestDiff)
        {
            nearestDiff = diff;
            nearestIndex = timeIndex;
        }
    }

    return nearestIndex;
}


template<class Type>
void MapConservativeVolFields
(
    const IOobjectList& objects,
    const conservativeMeshToMesh& meshToMeshInterp,
    const label method
)
{
    const fvMesh& meshSource = meshToMeshInterp.fromMesh();
    const fvMesh& meshTarget = meshToMeshInterp.toMesh();

    word fieldClassName
    (
        GeometricField<Type, fvPatchField, volMesh>::typeName
    );

    IOobjectList fields = objects.lookupClass(fieldClassName);

    for
    (
        IOobjectList::iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        Info << "    Interpolating " << fieldIter()->name() << endl;

        // Read field
        GeometricField<Type, fvPatchField, volMesh> fieldSource
        (
            *fieldIter(),
            meshSource
        );

        IOobject fieldTargetIOobject
        (
            fieldIter()->name(),
            meshTarget.time().timeName(),
            meshTarget,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (fieldTargetIOobject.headerOk())
        {
            // Read fieldTarget
            GeometricField<Type, fvPatchField, volMesh> fieldTarget
            (
                fieldTargetIOobject,
                meshTarget
            );

            // Interpolate field
            meshToMeshInterp.interpolate
            (
                fieldTarget,
                fieldSource,
                method
            );

            // Write field
            fieldTarget.write();
        }
        else
        {
            fieldTargetIOobject.readOpt() = IOobject::NO_READ;

            // Interpolate field
            GeometricField<Type, fvPatchField, volMesh> fieldTarget
            (
                fieldTargetIOobject,
                meshToMeshInterp.interpolate(fieldSource, method)
            );

            // Write field
            fieldTarget.write();
        }
    }
}


// Test routine to determine interpolation error
void testMappingError
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const label method,
    const label nThreads,
    const bool forceRecalc,
    const bool writeAddr
)
{
    // Create the interpolation scheme
    conservativeMeshToMesh meshToMeshInterp
    (
        meshSource,
        meshTarget,
        nThreads,
        forceRecalc,
        writeAddr
    );

    const volVectorField& xC = meshSource.C();

    volScalarField fieldSource
    (
        IOobject
        (
            "alpha",
            meshSource.time().timeName(),
            meshSource,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        xC & vector(2,3,1)
    );

    scalar pi = mathematicalConstant::pi;

    // Test sinusoidal field
    forAll(fieldSource, cellI)
    {
        const vector x = xC[cellI];

        fieldSource[cellI] =
        (
            1.0
          + Foam::sin(2.0*pi*x.x())
          * Foam::sin(2.0*pi*x.y())
          * Foam::sin(2.0*pi*x.z())
        );
    }

    forAll(fieldSource.boundaryField(), patchI)
    {
        forAll(fieldSource.boundaryField()[patchI], faceI)
        {
            const vector x = xC.boundaryField()[patchI][faceI];

            fieldSource.boundaryField()[patchI][faceI] =
            (
                1.0
              + Foam::sin(2.0*pi*x.x())
              * Foam::sin(2.0*pi*x.y())
              * Foam::sin(2.0*pi*x.z())
            );
        }
    }

    volScalarField fieldTarget
    (
        IOobject
        (
            "alpha",
            meshTarget.time().timeName(),
            meshTarget,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        meshTarget.C() & vector(0,0,0)
    );

    // Interpolate field
    meshToMeshInterp.interpolate
    (
        fieldTarget,
        fieldSource,
        method
    );

    // Compute the interpolation error.
    {
        // Compute error
        volScalarField iError
        (
            IOobject
            (
                "iError",
                meshTarget.time().timeName(),
                meshTarget,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            meshTarget,
            fieldTarget.dimensions()
        );

        const volScalarField& sF = refCast<volScalarField>(fieldSource);
        const volScalarField& tF = refCast<volScalarField>(fieldTarget);

        const vectorField& sCentres = fieldSource.mesh().cellCentres();
        const vectorField& tCentres = fieldTarget.mesh().cellCentres();

        scalar sError = 0.0, tError = 0.0;
        scalar smError = 0.0, tmError = 0.0;

        const scalarField& isF = sF.internalField();
        const scalarField& itF = tF.internalField();

        forAll(isF, cellI)
        {
            const vector xC = sCentres[cellI];

            //scalar sExact = 2.0*xC.x() + 3.0*xC.y() + xC.z();

            scalar sExact =
            (
                1.0
              + Foam::sin(2.0*pi*xC.x())
              * Foam::sin(2.0*pi*xC.y())
              * Foam::sin(2.0*pi*xC.z())
            );

            sError += magSqr(isF[cellI] - sExact);
            smError = Foam::max(smError, mag(isF[cellI] - sExact));
        }

        forAll(itF, cellI)
        {
            const vector xC = tCentres[cellI];

            //scalar tExact = 2.0*xC.x() + 3.0*xC.y() + xC.z();

            scalar tExact =
            (
                1.0
              + Foam::sin(2.0*pi*xC.x())
              * Foam::sin(2.0*pi*xC.y())
              * Foam::sin(2.0*pi*xC.z())
            );

            tError += magSqr(itF[cellI] - tExact);
            tmError = Foam::max(tmError, mag(itF[cellI] - tExact));

            iError.internalField()[cellI] = mag(itF[cellI] - tExact);
        }

        // Compute mesh dx
        scalar sdx = Foam::cbrt(1.0 / isF.size());
        scalar tdx = Foam::cbrt(1.0 / itF.size());

        Info << " ~~~~~~~~~~~~~~~~~ " << nl
             << "      Source       " << nl
             << " ~~~~~~~~~~~~~~~~~ " << nl
             << " L2 error: " << Foam::sqrt(sError/isF.size()) << nl
             << " Linf error: " << smError << nl
             << " dx: " << sdx << nl
             << " dx2: " << sqr(sdx) << nl
             << endl;

        Info << " ~~~~~~~~~~~~~~~~~ " << nl
             << "      Target       " << nl
             << " ~~~~~~~~~~~~~~~~~ " << nl
             << " L2 error: " << Foam::sqrt(tError/itF.size()) << nl
             << " Linf error: " << tmError << nl
             << " dx: " << tdx << nl
             << " dx2: " << sqr(tdx) << nl
             << endl;
    }
}


void mapConservativeMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const label method,
    const label nThreads,
    const bool forceRecalc,
    const bool writeAddr
)
{
    // Create the interpolation scheme
    conservativeMeshToMesh meshToMeshInterp
    (
        meshSource,
        meshTarget,
        nThreads,
        forceRecalc,
        writeAddr
    );

    Info << nl
         << "Conservatively creating and mapping fields for time "
         << meshSource.time().timeName() << nl << endl;

    {
        // Search for list of objects for this time
        IOobjectList objects(meshSource, meshSource.time().timeName());

        // Map volFields
        // ~~~~~~~~~~~~~
        MapConservativeVolFields<scalar>(objects, meshToMeshInterp, method);
        MapConservativeVolFields<vector>(objects, meshToMeshInterp, method);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRoots.H"

#   include "createTimes.H"

#   include "setTimeIndex.H"

    runTimeSource.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);

    Info<< "\nSource time: " << runTimeSource.value()
        << "\nTarget time: " << runTimeTarget.value()
        << endl;

    Info<< "Create meshes\n" << endl;

    fvMesh meshSource
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTimeSource.timeName(),
            runTimeSource
        )
    );

    fvMesh meshTarget
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTimeTarget.timeName(),
            runTimeTarget
        )
    );

    Info<< "Source mesh size: " << meshSource.nCells() << tab
        << "Target mesh size: " << meshTarget.nCells() << nl << endl;

    switch (method)
    {
        case conservativeMeshToMesh::CONSERVATIVE:
        {
            Info << "Using method: CONSERVATIVE" << endl;

            break;
        }

        case conservativeMeshToMesh::INVERSE_DISTANCE:
        {
            Info << "Using method: INVERSE_DISTANCE" << endl;

            break;
        }

        case conservativeMeshToMesh::CONSERVATIVE_FIRST_ORDER:
        {
            Info << "Using method: CONSERVATIVE_FIRST_ORDER" << endl;

            break;
        }

        default:
        {
            FatalErrorIn("mapConservativeFields")
                << "unknown interpolation scheme " << method
                << exit(FatalError);
        }
    }

    if (testOnly)
    {
        testMappingError
        (
            meshSource,
            meshTarget,
            method,
            nThreads,
            forceRecalc,
            writeAddr
        );
    }
    else
    {
        mapConservativeMesh
        (
            meshSource,
            meshTarget,
            method,
            nThreads,
            forceRecalc,
            writeAddr
        );
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
