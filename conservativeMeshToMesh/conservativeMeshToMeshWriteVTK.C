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

\*---------------------------------------------------------------------------*/

#include "conservativeMeshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void conservativeMeshToMesh::writeVTK
(
    const label newCellIndex,
    const label oldCellIndex,
    const vectorField& cvxSet
) const
{
    // Write out points for post-processing
    labelListList cpList(cvxSet.size(), labelList(1));

    forAll(cpList, i)
    {
        cpList[i][0] = i;
    }

    // For post-processing purposes, define a name
    word cvxSetName
    (
        "cvxSet_"
      + Foam::name(newCellIndex)
      + '<'
      + Foam::name(oldCellIndex)
      + '>'
    );

    writeVTK
    (
        cvxSetName,
        cvxSet.size(),
        cvxSet.size(),
        cvxSet.size(),
        cvxSet,
        cpList,
        0
    );
}


// Output an entity as a VTK file
void conservativeMeshToMesh::writeVTK
(
    const word& name,
    const fvMesh& mesh,
    const label entity,
    const label primitiveType
) const
{
    writeVTK
    (
        name,
        mesh,
        labelList(1, entity),
        primitiveType
    );
}


// Output a list of primitives as a VTK file.
//  - primitiveType is:
//      0: List of points
//      1: List of edges
//      2: List of faces
//      3: List of cells
void conservativeMeshToMesh::writeVTK
(
    const word& name,
    const fvMesh& mesh,
    const labelList& cList,
    const label primitiveType
) const
{
    label nTotalCells = 0;
    label nPoints = 0, nCells = 0;

    // Estimate a size for points and cellPoints
    pointField points(6*cList.size());

    // Connectivity lists
    labelListList cpList(cList.size());

    // Create a map for local points
    Map<label> pointMap, reversePointMap, reverseCellMap;

    forAll(cList, cellI)
    {
        if (cList[cellI] < 0)
        {
            continue;
        }

        // Are we looking at points?
        if (primitiveType == 0)
        {
            // Size the list
            cpList[nCells].setSize(1);

            cpList[nCells] = cList[cellI];

            nTotalCells++;
        }

        // Are we looking at edges?
        if (primitiveType == 1)
        {
            // Size the list
            cpList[nCells].setSize(2);

            const edge& tEdge = mesh.edges()[cList[cellI]];

            cpList[nCells][0] = tEdge[0];
            cpList[nCells][1] = tEdge[1];

            nTotalCells += 2;
        }

        // Are we looking at faces?
        if (primitiveType == 2)
        {
            const face& tFace = mesh.faces()[cList[cellI]];

            if (tFace.size() == 3)
            {
                // Size the list
                cpList[nCells].setSize(3);

                // Write out in order
                cpList[nCells][0] = tFace[0];
                cpList[nCells][1] = tFace[1];
                cpList[nCells][2] = tFace[2];

                nTotalCells += 3;
            }
            else
            if (tFace.size() == 4)
            {
                // Size the list
                cpList[nCells].setSize(4);

                // Write out in order
                cpList[nCells][0] = tFace[0];
                cpList[nCells][1] = tFace[1];
                cpList[nCells][2] = tFace[2];
                cpList[nCells][3] = tFace[3];

                nTotalCells += 4;
            }
        }

        // Are we looking at cells?
        if (primitiveType == 3)
        {
            const cell& tCell = mesh.cells()[cList[cellI]];

            if (tCell.size() == 4)
            {
                // Point-ordering for tetrahedra
                const face& baseFace = mesh.faces()[tCell[0]];

                // Fetch labels for this cell
                const labelList& cPoints = mesh.cellPoints()[cList[cellI]];

                label apexPoint = -1;

                forAll(cPoints, pointI)
                {
                    if (findIndex(baseFace, cPoints[pointI]) == -1)
                    {
                        apexPoint = cPoints[pointI];
                        break;
                    }
                }

                // Size the list
                cpList[nCells].setSize(4);

                // Write-out in order
                label ownCell = mesh.faceOwner()[tCell[0]];

                if (ownCell == cList[cellI])
                {
                    cpList[nCells][0] = baseFace[2];
                    cpList[nCells][1] = baseFace[1];
                    cpList[nCells][2] = baseFace[0];
                    cpList[nCells][3] = apexPoint;
                }
                else
                {
                    cpList[nCells][0] = baseFace[0];
                    cpList[nCells][1] = baseFace[1];
                    cpList[nCells][2] = baseFace[2];
                    cpList[nCells][3] = apexPoint;
                }

                nTotalCells += 4;
            }
        }

        // Renumber to local ordering
        forAll(cpList[nCells], pointI)
        {
            // Check if this point was added to the map
            if (!pointMap.found(cpList[nCells][pointI]))
            {
                // Point was not found, so add it
                points[nPoints] = mesh.points()[cpList[nCells][pointI]];

                // Update the map
                pointMap.insert(cpList[nCells][pointI], nPoints);
                reversePointMap.insert(nPoints, cpList[nCells][pointI]);

                // Increment the number of points
                nPoints++;
            }

            // Renumber it.
            cpList[nCells][pointI] = pointMap[cpList[nCells][pointI]];
        }

        // Update the cell map.
        reverseCellMap.insert(nCells, cList[cellI]);

        nCells++;
    }

    // Finally write it out
    writeVTK
    (
        name,
        nPoints,
        nCells,
        nTotalCells,
        points,
        cpList,
        primitiveType,
        reversePointMap,
        reverseCellMap
    );
}


// Actual routine to write out the VTK file
void conservativeMeshToMesh::writeVTK
(
    const word& name,
    const label nPoints,
    const label nCells,
    const label nTotalCells,
    const pointField& points,
    const labelListList& cpList,
    const label primitiveType,
    const Map<label>& reversePointMap,
    const Map<label>& reverseCellMap
) const
{
    // Make the directory
    fileName dirName(toMesh().time().path()/"VTK"/toMesh().time().timeName());

    mkDir(dirName);

    // Open stream for output
    OFstream file(dirName/name+".vtk");

    // Write out the header
    file << "# vtk DataFile Version 2.0" << nl
         << name << ".vtk" << nl
         << "ASCII" << nl
         << "DATASET UNSTRUCTURED_GRID" << nl
         << "POINTS " << nPoints << " double" << nl;

    for (label i = 0; i < nPoints; i++)
    {
        file << points[i].x() << ' '
             << points[i].y() << ' '
             << points[i].z() << ' '
             << nl;
    }

    file << "CELLS " << nCells << " " << nTotalCells + nCells << endl;

    forAll(cpList, i)
    {
        if (cpList[i].size())
        {
            file << cpList[i].size() << ' ';

            forAll(cpList[i], j)
            {
                file << cpList[i][j] << ' ';
            }

            file << nl;
        }
    }

    file << "CELL_TYPES " << nCells << endl;

    forAll(cpList, i)
    {
        if (cpList[i].size() == 1)
        {
            // Vertex
            file << "1" << nl;
        }

        if (cpList[i].size() == 2)
        {
            // Edge
            file << "3" << nl;
        }

        if (cpList[i].size() == 3)
        {
            // Triangle face
            file << "5" << nl;
        }

        if
        (
            (cpList[i].size() == 4) &&
            (primitiveType == 2)
        )
        {
            // Quad face
            file << "9" << nl;
        }

        if
        (
            (cpList[i].size() == 4) &&
            (primitiveType == 3)
        )
        {
            // Tetrahedron
            file << "10" << nl;
        }

        if (cpList[i].size() == 6)
        {
            // Wedge
            file << "13" << nl;
        }
    }

    // Write out indices for visualization.
    if (reverseCellMap.size())
    {
        file << "CELL_DATA " << nCells << endl;

        file << "FIELD CellFields 1" << endl;

        file << "CellIds 1 " << nCells << " int" << endl;

        for (label i = 0; i < nCells; i++)
        {
            file << reverseCellMap[i] << ' ';
        }

        file << endl;
    }

    // Write out indices for visualization.
    if (reversePointMap.size())
    {
        file << "POINT_DATA " << nPoints << endl;

        file << "FIELD PointFields 1" << endl;

        file << "PointIds 1 " << nPoints << " int" << endl;

        for (label i = 0; i < nPoints; i++)
        {
            file << reversePointMap[i] << ' ';
        }

        file << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
