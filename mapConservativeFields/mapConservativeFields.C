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
    const conservativeMeshToMesh& meshToMeshInterp
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
                fieldSource
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
                meshToMeshInterp.interpolate
                (
                    fieldSource
                )
            );

            // Write field
            fieldTarget.write();
        }
    }
}


void mapConservativeMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
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
        MapConservativeVolFields<scalar>(objects, meshToMeshInterp);
        MapConservativeVolFields<vector>(objects, meshToMeshInterp);
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

    mapConservativeMesh
    (
        meshSource,
        meshTarget,
        nThreads,
        forceRecalc,
        writeAddr
    );

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
