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

Class
    dynamicTopoFvMesh

Description
    Functions specific to connectivity checking and debugging

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "dynamicTopoFvMesh.H"

#include "IOmanip.H"
#include "volFields.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return mesh cell-quality values
// Valid for 3D tetrahedral meshes only...
bool dynamicTopoFvMesh::meshQuality
(
    bool outputOption
)
{
    // Valid for 3D tetrahedral meshes only...
    if (twoDMesh_)
    {
        return false;
    }

    Switch dumpMeshQuality(false);

    if
    (
        dict_.subDict("dynamicTopoFvMesh").found("dumpMeshQuality") ||
        mandatory_
    )
    {
        dumpMeshQuality =
        (
            dict_.subDict("dynamicTopoFvMesh").lookup("dumpMeshQuality")
        );
    }

    volScalarField *mqPtr(NULL);

    if (dumpMeshQuality && time().outputTime())
    {
        mqPtr = new volScalarField
        (
            IOobject
            (
                "meshQuality",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimensionedScalar("scalar", dimless, 0)
        );
    }

    // Compute statistics on the fly
    label nCells = 0, minCell = -1;
    scalar maxQuality = -GREAT;
    scalar minQuality =  GREAT;
    scalar cQuality, meanQuality = 0.0;

    // Loop through all cells in the mesh and compute cell quality
    forAll(cells_, cellI)
    {
        if (cells_[cellI].empty())
        {
            continue;
        }

        // Compute cell quality
        cQuality = tetQuality(cellI);

        if (dumpMeshQuality && time().outputTime())
        {
            mqPtr->internalField()[cellI] = cQuality;
        }

        // Update statistics
        maxQuality = Foam::max(cQuality, maxQuality);

        if (cQuality < minQuality)
        {
            minQuality = cQuality;
            minCell = cellI;
        }

        meanQuality += cQuality;
        nCells++;

        // Add to the list of slivers
        if ((cQuality < sliverThreshold_) && (cQuality > 0.0))
        {
            thresholdSlivers_.insert(cellI, cQuality);
        }
    }

    // Output statistics:
    if (outputOption || (debug > 0))
    {
        // Reduce statistics across processors.
        reduce(minQuality, minOp<scalar>());
        reduce(maxQuality, maxOp<scalar>());
        reduce(meanQuality, sumOp<scalar>());
        reduce(nCells, sumOp<label>());

        if (minQuality < 0.0)
        {
            WarningIn("dynamicTopoFvMesh::meshQuality()")
                << nl
                << "Minimum cell quality is: " << minQuality
                << " at cell: " << minCell
                << endl;
        }

        Info << " ~~~ Mesh Quality Statistics ~~~ " << endl;
        Info << " Min: " << minQuality << endl;
        Info << " Max: " << maxQuality << endl;
        Info << " Mean: " << meanQuality/nCells << endl;
        Info << " Cells: " << nCells << endl;
        Info << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
    }

    if (dumpMeshQuality && time().outputTime())
    {
        return mqPtr->write();
    }

    return true;
}


// Utility to check whether points of an edge lie on a boundary.
const FixedList<bool,2>
dynamicTopoFvMesh::checkEdgeBoundary
(
    const label eIndex
) const
{
    FixedList<bool,2> edgeBoundary(false);

    const edge& edgeToCheck = edges_[eIndex];

    // Loop through edges connected to both points,
    // and check if any of them lie on boundaries.
    // Used to ensure that collapses happen towards boundaries.
    forAll(edgeToCheck, pointI)
    {
        const labelList& pEdges = pointEdges_[edgeToCheck[pointI]];

        forAll(pEdges, edgeI)
        {
            // Determine the patch this edge belongs to
            if (whichEdgePatch(pEdges[edgeI]) > -1)
            {
                edgeBoundary[pointI] = true;
                break;
            }
        }
    }

    return edgeBoundary;
}


// Check whether the given edge is on a bounding curve
bool dynamicTopoFvMesh::checkBoundingCurve(const label eIndex) const
{
    // Internal edges don't count
    label edgePatch = -1;

    if ((edgePatch = whichEdgePatch(eIndex)) < 0)
    {
        return false;
    }
    else
    {
        // Check whether this edge shouldn't be swapped
        if (findIndex(noSwapPatchIDs_, edgePatch) > -1)
        {
            return true;
        }
    }

    if (coupledModification_)
    {
        return false;
    }

    // Check if two boundary faces lie on different face-patches
    FixedList<vector, 2> fNorm;
    label fPatch, firstPatch = -1, secondPatch = -1, count = 0;
    const labelList& edgeFaces = edgeFaces_[eIndex];

    forAll(edgeFaces, faceI)
    {
        if ((fPatch = whichPatch(edgeFaces[faceI])) > -1)
        {
            // Obtain the normal.
            if (twoDMesh_)
            {
                fNorm[count] = quadFaceNormal(faces_[edgeFaces[faceI]]);
            }
            else
            {
                fNorm[count] = triFaceNormal(faces_[edgeFaces[faceI]]);
            }

            // Normalize it.
            fNorm[count] /= mag(fNorm[count]);

            count++;

            if (firstPatch == -1)
            {
                firstPatch = fPatch;
            }
            else
            {
                secondPatch = fPatch;
                break;
            }
        }
    }

    scalar deviation = (fNorm[0] & fNorm[1]);

    // Check if the curvature is too high
    if (mag(deviation) < 0.85)
    {
        return true;
    }

    // Check if the edge borders two different patches
    if (firstPatch != secondPatch)
    {
        return true;
    }

    // Not on a bounding curve
    return false;
}


// Check triangulation quality for an edge index
bool dynamicTopoFvMesh::checkQuality
(
    const label eIndex,
    const labelList& m,
    const PtrList<scalarListList>& Q,
    const scalar minQuality,
    const label checkIndex
) const
{
    bool myResult = false;

    // Non-coupled check
    if (Q[checkIndex][0][m[checkIndex]-1] > minQuality)
    {
        myResult = true;

        if (debug > 2)
        {
            Info << " eIndex: " << eIndex
                 << " minQuality: " << minQuality
                 << " newQuality: " << Q[checkIndex][0][m[checkIndex]-1]
                 << endl;
        }
    }

    if (coupledModification_)
    {
        if (locallyCoupledEdge(eIndex))
        {
            // Check the quality of the slave edge as well.
            label sIndex = -1;

            // Loop through masterToSlave and determine the slave index.
            forAll(patchCoupling_, patchI)
            {
                if (patchCoupling_(patchI))
                {
                    const label edgeEnum  = coupleMap::EDGE;
                    const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                    if ((sIndex = cMap.findSlaveIndex(edgeEnum, eIndex)) > -1)
                    {
                        break;
                    }
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn("dynamicTopoFvMesh::checkQuality")
                    << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Edge: " << eIndex << nl
                    << abort(FatalError);
            }

            // Turn off switch temporarily.
            unsetCoupledModification();

            // Recursively call for the slave edge.
            myResult =
            (
                myResult && checkQuality(sIndex, m, Q, minQuality, 1)
            );

            // Turn it back on.
            setCoupledModification();
        }
        else
        if (processorCoupledEdge(eIndex))
        {

        }
    }

    return myResult;
}


// Print out tables for debugging
void dynamicTopoFvMesh::printTables
(
    const labelList& m,
    const PtrList<scalarListList>& Q,
    const PtrList<labelListList>& K,
    const label checkIndex
) const
{
    Info << "m: " << m[checkIndex] << endl;

    // Print out Q
    Info << "===" << endl;
    Info << " Q " << endl;
    Info << "===" << endl;

    Info << "   ";

    for (label j = 0; j < m[checkIndex]; j++)
    {
        Info << setw(12) << j;
    }

    Info << nl;

    for (label j = 0; j < m[checkIndex]; j++)
    {
        Info << "-------------";
    }

    Info << nl;

    for (label i = 0; i < (m[checkIndex]-2); i++)
    {
        Info << i << ": ";

        for (label j = 0; j < m[checkIndex]; j++)
        {
            Info << setw(12) << Q[checkIndex][i][j];
        }

        Info << nl;
    }

    // Print out K
    Info << "===" << endl;
    Info << " K " << endl;
    Info << "===" << endl;

    Info << "   ";

    for (label j = 0; j < m[checkIndex]; j++)
    {
        Info << setw(12) << j;
    }

    Info << nl;

    for (label i = 0; i < (m[checkIndex]-2); i++)
    {
        Info << i << ": ";

        for (label j = 0; j < m[checkIndex]; j++)
        {
            Info << setw(12) << K[checkIndex][i][j];
        }

        Info << nl;
    }

    Info << endl;
}


// Check old-volumes for an input triangulation
bool dynamicTopoFvMesh::checkTriangulationVolumes
(
    const label eIndex,
    const labelList& hullVertices,
    const labelListList& triangulations
) const
{
    label m = hullVertices.size();
    scalar tetVol = 0.0;

    const edge& edgeToCheck = edges_[eIndex];

    for (label i = 0; i < (m-2); i++)
    {
        // Compute volume for the upper-half
        tetVol =
        (
            tetVolume
            (
                oldPoints_[hullVertices[triangulations[0][i]]],
                oldPoints_[hullVertices[triangulations[1][i]]],
                oldPoints_[hullVertices[triangulations[2][i]]],
                oldPoints_[edgeToCheck[0]]
            )
        );

        if (tetVol < 0.0)
        {
            if (debug > 2)
            {
                InfoIn("dynamicTopoFvMesh::checkTriangulationVolumes") << nl
                    << "Swap sequence leads to negative old-volumes." << nl
                    << "Edge: " << edgeToCheck << nl
                    << "using Points: " << nl
                    << oldPoints_[hullVertices[triangulations[0][i]]] << nl
                    << oldPoints_[hullVertices[triangulations[1][i]]] << nl
                    << oldPoints_[hullVertices[triangulations[2][i]]] << nl
                    << oldPoints_[edgeToCheck[0]] << endl;
            }

            return true;
        }

        tetVol =
        (
            tetVolume
            (
                oldPoints_[hullVertices[triangulations[2][i]]],
                oldPoints_[hullVertices[triangulations[1][i]]],
                oldPoints_[hullVertices[triangulations[0][i]]],
                oldPoints_[edgeToCheck[1]]
            )
        );

        if (tetVol < 0.0)
        {
            if (debug > 2)
            {
                InfoIn("dynamicTopoFvMesh::checkTriangulationVolumes") << nl
                    << "Swap sequence leads to negative old-volumes." << nl
                    << "Edge: " << edgeToCheck << nl
                    << "using Points: " << nl
                    << oldPoints_[hullVertices[triangulations[2][i]]] << nl
                    << oldPoints_[hullVertices[triangulations[1][i]]] << nl
                    << oldPoints_[hullVertices[triangulations[0][i]]] << nl
                    << oldPoints_[edgeToCheck[1]] << endl;
            }

            return true;
        }
    }

    return false;
}


// Output an entity as a VTK file
void dynamicTopoFvMesh::writeVTK
(
    const word& name,
    const label entity,
    const label primitiveType,
    const bool useOldConnectivity,
    const bool useOldPoints
) const
{
    writeVTK
    (
        name,
        labelList(1, entity),
        primitiveType,
        useOldConnectivity,
        useOldPoints
    );
}


// Output a list of primitives as a VTK file.
//  - primitiveType is:
//      0: List of points
//      1: List of edges
//      2: List of faces
//      3: List of cells
void dynamicTopoFvMesh::writeVTK
(
    const word& name,
    const labelList& cList,
    const label primitiveType,
    const bool useOldConnectivity,
    const bool useOldPoints
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

            // Can't use old connectivity for edges,
            // because the convention differs from FOAM.
            if (useOldConnectivity)
            {
                FatalErrorIn("void dynamicTopoFvMesh::writeVTK()")
                    << "Cannot use old connectivity for edges."
                    << abort(FatalError);
            }

            const edge& tEdge = edges_[cList[cellI]];

            cpList[nCells][0] = tEdge[0];
            cpList[nCells][1] = tEdge[1];

            nTotalCells += 2;
        }

        // Are we looking at faces?
        if (primitiveType == 2)
        {
            face tFace;

            if (useOldConnectivity)
            {
                tFace = polyMesh::faces()[cList[cellI]];
            }
            else
            {
                tFace = faces_[cList[cellI]];
            }

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
            cell tCell;

            if (useOldConnectivity)
            {
                tCell = polyMesh::cells()[cList[cellI]];
            }
            else
            {
                tCell = cells_[cList[cellI]];
            }

            if (tCell.size() == 4)
            {
                // Point-ordering for tetrahedra
                face baseFace, checkFace;

                if (useOldConnectivity)
                {
                    baseFace = polyMesh::faces()[tCell[0]];
                    checkFace = polyMesh::faces()[tCell[1]];
                }
                else
                {
                    baseFace = faces_[tCell[0]];
                    checkFace = faces_[tCell[1]];
                }

                // Size the list
                cpList[nCells].setSize(4);

                // Get the fourth point
                label apexPoint = findIsolatedPoint(baseFace, checkFace);

                // Something's wrong with connectivity.
                if (apexPoint == -1)
                {
                    FatalErrorIn("dynamicTopoFvMesh::writeVTK()")
                        << "Cell: " << cList[cellI]
                        << ":: " << tCell
                        << " has inconsistent connectivity."
                        << abort(FatalError);
                }

                // Write-out in order
                label ownCell = -1;

                if (useOldConnectivity)
                {
                    ownCell = polyMesh::faceOwner()[tCell[0]];
                }
                else
                {
                    ownCell = owner_[tCell[0]];
                }

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
            else
            if (tCell.size() == 5)
            {
                // Point-ordering for wedge cells
                face cFace, nFace;
                label firstTriFace = -1;

                // Size the list
                cpList[nCells].setSize(6);

                // Figure out triangle faces
                forAll(tCell, faceI)
                {
                    if (useOldConnectivity)
                    {
                        cFace = polyMesh::faces()[tCell[faceI]];
                    }
                    else
                    {
                        cFace = faces_[tCell[faceI]];
                    }

                    if (cFace.size() == 3)
                    {
                        if (firstTriFace == -1)
                        {
                            firstTriFace = tCell[faceI];

                            // Right-handedness is assumed here.
                            // Tri-faces are always on the boundary.
                            cpList[nCells][0] = cFace[0];
                            cpList[nCells][1] = cFace[1];
                            cpList[nCells][2] = cFace[2];
                        }
                        else
                        {
                            // Detect the three other points.
                            forAll(tCell, faceJ)
                            {
                                if (useOldConnectivity)
                                {
                                    nFace = polyMesh::faces()[tCell[faceJ]];
                                }
                                else
                                {
                                    nFace = faces_[tCell[faceJ]];
                                }

                                if (nFace.size() == 4)
                                {
                                    // Search for vertices on cFace
                                    // in this face.
                                    forAll(cFace, I)
                                    {
                                        label i = nFace.which(cFace[I]);

                                        if (i != -1)
                                        {
                                            label p = nFace.prevLabel(i);
                                            label n = nFace.nextLabel(i);

                                            if (p == cpList[nCells][0])
                                            {
                                                cpList[nCells][3] = cFace[I];
                                            }

                                            if (p == cpList[nCells][1])
                                            {
                                                cpList[nCells][4] = cFace[I];
                                            }

                                            if (p == cpList[nCells][2])
                                            {
                                                cpList[nCells][5] = cFace[I];
                                            }

                                            if (n == cpList[nCells][0])
                                            {
                                                cpList[nCells][3] = cFace[I];
                                            }

                                            if (n == cpList[nCells][1])
                                            {
                                                cpList[nCells][4] = cFace[I];
                                            }

                                            if (n == cpList[nCells][2])
                                            {
                                                cpList[nCells][5] = cFace[I];
                                            }
                                        }
                                    }
                                }
                            }

                            break;
                        }
                    }
                }

                nTotalCells += 6;
            }
        }

        // Renumber to local ordering
        forAll(cpList[nCells], pointI)
        {
            // Check if this point was added to the map
            if (!pointMap.found(cpList[nCells][pointI]))
            {
                // Point was not found, so add it
                if (useOldPoints)
                {
                    points[nPoints] = oldPoints_[cpList[nCells][pointI]];
                }
                else
                {
                    points[nPoints] = points_[cpList[nCells][pointI]];
                }

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
void dynamicTopoFvMesh::writeVTK
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
    fileName dirName(time().path()/"VTK"/time().timeName());

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
        file << setprecision(10)
             << points[i].x() << ' '
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


// Check the state of connectivity lists
void dynamicTopoFvMesh::checkConnectivity
(
    label maxErrors
) const
{
    label nFailedChecks = 0;

    messageStream ConnectivityWarning
    (
        "dynamicTopoFvMesh Connectivity Warning",
        messageStream::SERIOUS,
        maxErrors
    );

    // Check face-label ranges
    Info << "Checking index ranges...";

    forAll(edges_, edgeI)
    {
        const edge& curEdge = edges_[edgeI];

        if (curEdge == edge(-1, -1))
        {
            continue;
        }

        if
        (
            curEdge[0] < 0 || curEdge[0] > (points_.size()-1) ||
            curEdge[1] < 0 || curEdge[1] > (points_.size()-1)
        )
        {
            Pout << "Edge " << edgeI
                 << " contains vertex labels out of range: "
                 << curEdge
                 << " Max point index = " << (points_.size()-1) << endl;

            nFailedChecks++;

            ConnectivityWarning()
                 << nl << "Edge-point connectivity is inconsistent."
                 << endl;
        }

        // Check for unique point-labels
        if (curEdge[0] == curEdge[1])
        {
            Pout << "Edge " << edgeI
                 << " contains identical vertex labels: "
                 << curEdge
                 << endl;

            nFailedChecks++;

            ConnectivityWarning()
                 << nl << "Edge-point connectivity is inconsistent."
                 << endl;
        }
    }

    label allPoints = points_.size();
    labelList nPointFaces(allPoints, 0);

    forAll(faces_, faceI)
    {
        const face& curFace = faces_[faceI];

        if (curFace.empty())
        {
            continue;
        }

        if (min(curFace) < 0 || max(curFace) > (points_.size()-1))
        {
            Pout << "Face " << faceI
                 << " contains vertex labels out of range: "
                 << curFace
                 << " Max point index = " << (points_.size()-1) << endl;

            nFailedChecks++;

            ConnectivityWarning()
                 << nl << "Face-point connectivity is inconsistent."
                 << endl;
        }

        // Check for unique point-labels
        labelHashSet uniquePoints;

        forAll(curFace, pointI)
        {
            bool inserted = uniquePoints.insert(curFace[pointI]);

            if (!inserted)
            {
                Pout << "Face " << faceI
                     << " contains identical vertex labels: "
                     << curFace
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
                     << nl << "Face-point connectivity is inconsistent."
                     << endl;
            }
        }

        // Count faces per point
        forAll(curFace, pointI)
        {
            nPointFaces[curFace[pointI]]++;
        }
    }

    forAll(cells_, cellI)
    {
        const cell& curCell = cells_[cellI];

        if (curCell.empty())
        {
            continue;
        }

        if (min(curCell) < 0 || max(curCell) > (faces_.size()-1))
        {
            Pout << "Cell " << cellI
                 << " contains vertex labels out of range: "
                 << curCell
                 << " Max point index = " << (faces_.size()-1) << endl;

            nFailedChecks++;

            ConnectivityWarning()
                 << nl << "Cell-Face connectivity is inconsistent."
                 << endl;
        }

        // Check for unique face-labels
        labelHashSet uniqueFaces;

        forAll(curCell, faceI)
        {
            bool inserted = uniqueFaces.insert(curCell[faceI]);

            if (!inserted)
            {
                Pout << "Cell " << cellI
                     << " contains identical face labels: "
                     << curCell
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
                     << nl << "Cell-Face connectivity is inconsistent."
                     << endl;
            }
        }
    }

    Info << "Done." << endl;

    Info << "Checking for unused points...";

    forAll(nPointFaces, pointI)
    {
        if (nPointFaces[pointI] == 0)
        {
            // This might be a deleted point.
            if (pointI < nOldPoints_)
            {
                if (reversePointMap_[pointI] == -1)
                {
                    continue;
                }
            }
            else
            {
                if (deletedPoints_.found(pointI))
                {
                    continue;
                }
            }

            // Looks like this is really an unused point.
            Pout << "Point " << pointI << " is unused. " << endl;

            nFailedChecks++;

            ConnectivityWarning()
                 << nl << "Point-Face connectivity is inconsistent."
                 << endl;
        }
    }

    Info << "Done." << endl;

    Info << "Checking edge-face connectivity...";

    label allEdges = edges_.size();
    labelList nEdgeFaces(allEdges, 0);

    forAll(faceEdges_, faceI)
    {
        const labelList& faceEdges = faceEdges_[faceI];

        if (faceEdges.empty())
        {
            continue;
        }

        // Check consistency of face-edge-points as well
        edgeList eList = faces_[faceI].edges();

        forAll(faceEdges,edgeI)
        {
            nEdgeFaces[faceEdges[edgeI]]++;

            // Check if this edge actually belongs to this face
            bool found = false;
            const edge& edgeToCheck = edges_[faceEdges[edgeI]];

            forAll(eList, edgeII)
            {
                if (edgeToCheck == eList[edgeII])
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                Pout << nl << nl << "Edge: " << faceEdges[edgeI]
                     << ": " << edgeToCheck << nl
                     << "was not found in face: " << faceI
                     << ": " << faces_[faceI] << nl
                     << "faceEdges: " << faceEdges
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
                     << nl << "Edge-Face connectivity is inconsistent."
                     << endl;
            }
        }
    }

    label nInternalEdges = 0;
    labelList patchInfo(numPatches_, 0);

    forAll(edgeFaces_, edgeI)
    {
        const labelList& edgeFaces = edgeFaces_[edgeI];

        if (edgeFaces.empty())
        {
            continue;
        }

        if (edgeFaces.size() != nEdgeFaces[edgeI])
        {
            Pout << nl << nl << "Edge: " << edgeI
                 << ": edgeFaces: " << edgeFaces << endl;

            nFailedChecks++;

            ConnectivityWarning()
                << nl << "Edge-Face connectivity is inconsistent."
                << endl;
        }

        label nBF = 0;

        // Check if this edge belongs to faceEdges for each face
        forAll(edgeFaces, faceI)
        {
            if (findIndex(faceEdges_[edgeFaces[faceI]], edgeI) == -1)
            {
                Pout << nl << nl << "Edge: " << edgeI << ": " << edges_[edgeI]
                     << ", edgeFaces: " << edgeFaces << nl
                     << "was not found in faceEdges of face: "
                     << edgeFaces[faceI] << ": " << faces_[edgeFaces[faceI]]
                     << nl << "faceEdges: " << faceEdges_[edgeFaces[faceI]]
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
                    << nl << "Edge-Face connectivity is inconsistent."
                    << endl;
            }

            if (neighbour_[edgeFaces[faceI]] == -1)
            {
                nBF++;
            }
        }

        if (nBF == 0)
        {
            nInternalEdges++;

            // Check if this edge is actually internal.
            if (whichEdgePatch(edgeI) >= 0)
            {
                Pout << "Edge: " << edgeI
                     << ": " << edges_[edgeI] << " is internal, "
                     << " but patch is specified as: "
                     << whichEdgePatch(edgeI)
                     << endl;

                nFailedChecks++;
            }
        }
        else
        {
            label patchID = whichEdgePatch(edgeI);

            // Check if this edge is actually on a boundary.
            if (patchID < 0)
            {
                Pout << "Edge: " << edgeI
                     << ": " << edges_[edgeI]
                     << " is on a boundary, but patch is specified as: "
                     << patchID << endl;

                nFailedChecks++;
            }
            else
            {
                patchInfo[patchID]++;
            }

            if (nBF > 2)
            {
                Pout << "Edge: " << edgeI
                     << ": " << edges_[edgeI]
                     << " has " << nBF
                     << " boundary faces connected to it." << nl
                     << " Pinched manifolds are not allowed."
                     << endl;

                nFailedChecks++;
            }
        }
    }

    if (nInternalEdges != nInternalEdges_)
    {
        Pout << nl << "Internal edge-count is inconsistent." << nl
             << " Counted internal edges: " << nInternalEdges
             << " Actual count: " << nInternalEdges_ << endl;

        nFailedChecks++;
    }

    forAll(patchInfo, patchI)
    {
        if (patchInfo[patchI] != edgePatchSizes_[patchI])
        {
            Pout << "Patch-count is inconsistent." << nl
                 << " Patch: " << patchI
                 << " Counted edges: " << patchInfo[patchI]
                 << " Actual count: " << edgePatchSizes_[patchI] << endl;

            nFailedChecks++;
        }
    }

    // Check added edge patches to ensure that it is consistent
    forAllConstIter(Map<label>, addedEdgePatches_, aepIter)
    {
        label key = aepIter.key();
        label patch = aepIter();

        label nBF = 0;
        const labelList& edgeFaces = edgeFaces_[key];

        // Check if any faces on boundaries
        forAll(edgeFaces, faceI)
        {
            if (neighbour_[edgeFaces[faceI]] == -1)
            {
                nBF++;
            }
        }

        if ((patch < 0) && (nBF > 0))
        {
            Pout << nl << nl << "Edge: " << key
                 << ", edgeFaces: " << edgeFaces
                 << " is internal, but contains boundary faces."
                 << endl;

            nFailedChecks++;
        }

        if ((patch >= 0) && (nBF != 2))
        {
            Pout << nl << nl << "Edge: " << key
                 << ", edgeFaces: " << edgeFaces
                 << " is on a boundary patch, but doesn't contain"
                 << " two boundary faces."
                 << endl;

            nFailedChecks++;
        }
    }

    Info << "Done." << endl;

    // Check coupled-patch sizes
    forAll(patchCoupling_, patchI)
    {
        if (patchCoupling_(patchI))
        {
            const coupleMap& cMap = patchCoupling_[patchI].patchMap();

            label mSize = patchSizes_[cMap.masterIndex()];
            label sSize = patchSizes_[cMap.slaveIndex()];

            if (mSize != sSize)
            {
                Pout << "Coupled patch-count is inconsistent." << nl
                     << " Master Patch: " << cMap.masterIndex()
                     << " Count: " << mSize << nl
                     << " Slave Patch: " << cMap.slaveIndex()
                     << " Count: " << sSize
                     << endl;

                nFailedChecks++;
            }
        }
    }

    if (!twoDMesh_)
    {
        Info << "Checking point-edge connectivity...";

        label allPoints = points_.size();
        List<labelHashSet> hlPointEdges(allPoints);

        forAll(edges_, edgeI)
        {
            if (edgeFaces_[edgeI].size())
            {
                hlPointEdges[edges_[edgeI][0]].insert(edgeI);
                hlPointEdges[edges_[edgeI][1]].insert(edgeI);
            }
        }

        forAll(pointEdges_, pointI)
        {
            const labelList& pointEdges = pointEdges_[pointI];

            if (pointEdges.empty())
            {
                continue;
            }

            forAll(pointEdges, edgeI)
            {
                if (!hlPointEdges[pointI].found(pointEdges[edgeI]))
                {
                    Pout << nl << nl << "Point: " << pointI << nl
                         << "pointEdges: " << pointEdges << nl
                         << "hlPointEdges: " << hlPointEdges[pointI]
                         << endl;

                    nFailedChecks++;

                    ConnectivityWarning()
                        << nl << "Point-Edge connectivity is inconsistent."
                        << endl;
                }
            }

            // Do a size check as well
            if
            (
                hlPointEdges[pointI].size() != pointEdges.size() ||
                pointEdges.size() == 1
            )
            {
                Pout << nl << nl << "Point: " << pointI << nl
                     << "pointEdges: " << pointEdges << nl
                     << "hlPointEdges: " << hlPointEdges[pointI]
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
                    << "Size inconsistency."
                    << nl << "Point-Edge connectivity is inconsistent."
                    << endl;
            }
        }

        Info << "Done." << endl;

        Info << "Checking edge-points connectivity...";

        label otherPoint = -1, nextPoint = -1;

        forAll(edgePoints_, edgeI)
        {
            // Do a preliminary size check
            const labelList& edgePoints = edgePoints_[edgeI];
            const labelList& edgeFaces = edgeFaces_[edgeI];

            if (edgeFaces.empty())
            {
                continue;
            }

            if (edgePoints.size() != edgeFaces.size())
            {
                Pout << nl << nl
                     << "Edge: " << edgeI
                     << " " << edges_[edgeI] << endl;

                Pout << "edgeFaces: " << edgeFaces << endl;
                forAll(edgeFaces, faceI)
                {
                    Info << edgeFaces[faceI] << ": "
                         << faces_[edgeFaces[faceI]]
                         << endl;
                }

                Pout << "edgePoints: " << edgePoints << endl;

                nFailedChecks++;

                ConnectivityWarning()
                    << nl << "Edge-Points connectivity is inconsistent."
                    << endl;
            }

            // Now check to see that both lists are consistent.
            const edge& edgeToCheck = edges_[edgeI];

            forAll(edgeFaces, faceI)
            {
                const face& faceToCheck = faces_[edgeFaces[faceI]];

                findIsolatedPoint
                (
                    faceToCheck,
                    edgeToCheck,
                    otherPoint,
                    nextPoint
                );

                if (findIndex(edgePoints, otherPoint) == -1)
                {
                    Pout << nl << nl
                         << "Edge: " << edgeI
                         << " " << edges_[edgeI] << endl;

                    Pout << "edgeFaces: " << edgeFaces << endl;
                    forAll(edgeFaces, faceI)
                    {
                        Info << edgeFaces[faceI] << ": "
                             << faces_[edgeFaces[faceI]]
                             << endl;
                    }

                    Pout << "edgePoints: " << edgePoints << endl;

                    nFailedChecks++;

                    ConnectivityWarning()
                        << nl << "Edge-Points connectivity is inconsistent."
                        << endl;
                }
            }
        }

        Info << "Done." << endl;
    }

    Info << "Checking cell-point connectivity...";

    // Loop through all cells and construct cell-to-node
    label cIndex = 0;
    label allCells = cells_.size();
    labelList cellIndex(allCells);
    List<labelHashSet> cellToNode(allCells);

    forAll(cells_, cellI)
    {
        const cell& thisCell = cells_[cellI];

        if (thisCell.empty())
        {
            continue;
        }

        cellIndex[cIndex] = cellI;

        forAll(thisCell, faceI)
        {
            const labelList& fEdges = faceEdges_[thisCell[faceI]];

            forAll(fEdges, edgeI)
            {
                const edge& thisEdge = edges_[fEdges[edgeI]];

                if (!cellToNode[cIndex].found(thisEdge[0]))
                {
                    cellToNode[cIndex].insert(thisEdge[0]);
                }

                if (!cellToNode[cIndex].found(thisEdge[1]))
                {
                    cellToNode[cIndex].insert(thisEdge[1]);
                }
            }
        }

        cIndex++;
    }

    // Resize the lists
    cellIndex.setSize(cIndex);
    cellToNode.setSize(cIndex);

    // Preliminary check for size
    forAll(cellToNode, cellI)
    {
        if
        (
            (cellToNode[cellI].size() != 6 && twoDMesh_) ||
            (cellToNode[cellI].size() != 4 && !twoDMesh_)
        )
        {
            Pout << nl << "Warning: Cell: "
                 << cellIndex[cellI] << " is inconsistent. "
                 << endl;

            const cell& failedCell = cells_[cellIndex[cellI]];

            Info << "Cell faces: " << failedCell << endl;

            forAll(failedCell, faceI)
            {
                Info << "\tFace: " << failedCell[faceI]
                     << " :: " << faces_[failedCell[faceI]]
                     << endl;

                const labelList& fEdges = faceEdges_[failedCell[faceI]];

                forAll(fEdges, edgeI)
                {
                    Info << "\t\tEdge: " << fEdges[edgeI]
                         << " :: " << edges_[fEdges[edgeI]]
                         << endl;
                }
            }

            nFailedChecks++;
        }
    }

    Info << "Done." << endl;

    reduce(nFailedChecks, orOp<bool>());

    if (nFailedChecks)
    {
        FatalErrorIn("dynamicTopoFvMesh::checkConnectivity()")
            << nFailedChecks << " failures were found in connectivity."
            << abort(FatalError);
    }
}


// Check for legitimacy of patches
void dynamicTopoFvMesh::checkPatches
(
    const wordList& patchList
) const
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    forAll(patchList, wordI)
    {
        bool foundPatch = false;

        forAll(boundary, patchI)
        {
            if (boundary[patchI].name() == patchList[wordI])
            {
                foundPatch = true;
                break;
            }
        }

        if (!foundPatch)
        {
            FatalErrorIn("dynamicTopoFvMesh::checkPatches()")
                << " Could not find patch: "
                << patchList[wordI] << nl
                << abort(FatalError);
        }
    }
}


// Write out proc IDs for post-processing
void dynamicTopoFvMesh::writeProcIDs() const
{
    if (Pstream::parRun())
    {
        Switch writeProcIDs(false);

        if
        (
            dict_.subDict("dynamicTopoFvMesh").found("writeProcIDs") ||
            mandatory_
        )
        {
            writeProcIDs =
            (
                dict_.subDict("dynamicTopoFvMesh").lookup("writeProcIDs")
            );
        }

        if (time().outputTime() && writeProcIDs)
        {
            volScalarField procID
            (
                IOobject
                (
                    "procID",
                    time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                *this,
                dimensionedScalar("scalar", dimless, Pstream::myProcNo())
            );

            procID.write();
        }
    }
}


// Utility method to check for invalid face-collapse.
// Returns 'true' if the collapse in NOT feasible.
bool dynamicTopoFvMesh::checkCollapse
(
    const labelList& triFaces,
    const FixedList<label,2>& c0BdyIndex,
    const FixedList<label,2>& c1BdyIndex,
    const FixedList<label,2>& original,
    const FixedList<label,2>& replacement,
    const bool checkNeighbour
) const
{
    face tmpTriFace(3);

    forAll(triFaces, indexI)
    {
        if
        (
            (triFaces[indexI] == c0BdyIndex[0])
         || (triFaces[indexI] == c0BdyIndex[1])
        )
        {
            continue;
        }

        if (checkNeighbour)
        {
            if
            (
                (triFaces[indexI] == c1BdyIndex[0])
             || (triFaces[indexI] == c1BdyIndex[1])
            )
            {
                continue;
            }
        }

        const face &triFace = faces_[triFaces[indexI]];

        forAll(triFace, pointI)
        {
            tmpTriFace[pointI] = triFace[pointI];

            if (triFace[pointI] == original[0])
            {
                tmpTriFace[pointI] = replacement[0];
            }

            if (triFace[pointI] == original[1])
            {
                tmpTriFace[pointI] = replacement[1];
            }
        }

        // Compute the area and check if it's zero/negative
        scalar origArea = triFaceArea(triFace);
        scalar newArea  = triFaceArea(tmpTriFace);

        if
        (
            (Foam::sign(origArea) != Foam::sign(newArea))
         || (mag(newArea) < (1e-3*origArea))
        )
        {
            // Inverted and/or degenerate.
            return true;
        }
    }

    // No problems, so a collapse is feasible.
    return false;
}


// Utility method to check whether the cell given by 'cellIndex' will yield
// a valid cell when 'pointIndex' is moved to 'newPoint'. The routine performs
// metric-based checks. Returns 'true' if the collapse in NOT feasible, and
// makes entries in cellsChecked to avoid repetitive checks.
bool dynamicTopoFvMesh::checkCollapse
(
    const point& newPoint,
    const point& oldPoint,
    const label pointIndex,
    const label cellIndex,
    labelHashSet& cellsChecked,
    bool forceOp
) const
{
    label faceIndex = -1;
    scalar cQuality = 0.0, oldVolume = 0.0;
    const cell& cellToCheck = cells_[cellIndex];

    // Look for a face that doesn't contain 'pointIndex'
    forAll(cellToCheck, faceI)
    {
        const face& currFace = faces_[cellToCheck[faceI]];

        if (currFace.which(pointIndex) < 0)
        {
            faceIndex = cellToCheck[faceI];
            break;
        }
    }

    // Compute cell-volume
    const face& faceToCheck = faces_[faceIndex];

    if (owner_[faceIndex] == cellIndex)
    {
        cQuality =
        (
            (*tetMetric_)
            (
                points_[faceToCheck[2]],
                points_[faceToCheck[1]],
                points_[faceToCheck[0]],
                newPoint
            )
        );

        oldVolume =
        (
            tetVolume
            (
                oldPoints_[faceToCheck[2]],
                oldPoints_[faceToCheck[1]],
                oldPoints_[faceToCheck[0]],
                oldPoint
            )
        );
    }
    else
    {
        cQuality =
        (
            (*tetMetric_)
            (
                points_[faceToCheck[0]],
                points_[faceToCheck[1]],
                points_[faceToCheck[2]],
                newPoint
            )
        );

        oldVolume =
        (
            tetVolume
            (
                oldPoints_[faceToCheck[0]],
                oldPoints_[faceToCheck[1]],
                oldPoints_[faceToCheck[2]],
                oldPoint
            )
        );
    }

    // Final quality check
    if (cQuality < sliverThreshold_ && !forceOp)
    {
        if (debug > 2)
        {
            InfoIn("dynamicTopoFvMesh::checkCollapse()")
                << "\nCollapsing cell: " << cellIndex
                << " containing points:\n"
                << faceToCheck[0] << "," << faceToCheck[1] << ","
                << faceToCheck[2] << "," << pointIndex << nl
                << "will yield a quality of: " << cQuality
                << ", when " << pointIndex
                << " is moved to location: " << nl
                << newPoint
                << endl;
        }

        return true;
    }

    // Negative quality is a no-no
    if (cQuality < 0.0)
    {
        if (forceOp)
        {
            InfoIn("dynamicTopoFvMesh::checkCollapse()")
                << "\nCollapsing cell: " << cellIndex
                << " containing points:\n"
                << faceToCheck[0] << "," << faceToCheck[1] << ","
                << faceToCheck[2] << "," << pointIndex << nl
                << "will yield a quality of: " << cQuality
                << ", when " << pointIndex
                << " is moved to location: " << nl
                << newPoint << nl
                << "Operation cannot be forced."
                << endl;
        }

        return true;
    }

    // Negative old-volume is also a no-no
    if (oldVolume < 0.0)
    {
        if (forceOp)
        {
            InfoIn("dynamicTopoFvMesh::checkCollapse()")
                << "\nCollapsing cell: " << cellIndex
                << " containing points:\n"
                << faceToCheck[0] << "," << faceToCheck[1] << ","
                << faceToCheck[2] << "," << pointIndex << nl
                << "will yield an old-volume of: " << oldVolume
                << ", when " << pointIndex
                << " is moved to location: " << nl
                << oldPoint << nl
                << "Operation cannot be forced."
                << endl;
        }

        return true;
    }

    // No problems, so a collapse is feasible
    cellsChecked.insert(cellIndex);

    return false;
}


} // End namespace Foam

// ************************************************************************* //
