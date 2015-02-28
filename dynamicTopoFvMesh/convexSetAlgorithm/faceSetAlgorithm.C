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
    faceSetAlgorithm

Description
    Class for convexSetAlgorithms pertaining to faces

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "faceSetAlgorithm.H"

#include "meshOps.H"
#include "triFace.H"
#include "polyMesh.H"
#include "triIntersection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Construct the search tree
void faceSetAlgorithm::constructSearchTree() const
{
    if (searchTreePtr_)
    {
        FatalErrorIn("void faceSetAlgorithm::constructSearchTree() const")
            << "searchTree is already allocated."
            << abort(FatalError);
    }

    treeBoundBox bBox(candidatePoints_);

    // Account for possibly degenerate bounding boxes
    // by artificially inflating in the degenerate dimension
    const scalar rVal = 1e-10;
    const scalar cutOff = 0.001;
    const vector boxSpan = bBox.span();
    const scalar magSpan = (magSqr(boxSpan) > rVal) ? magSqr(boxSpan) : rVal;
    const scalar magCutOff = (cutOff * magSpan);

    for (label cmpt = 0; cmpt < pTraits<vector>::nComponents; cmpt++)
    {
        if (boxSpan[cmpt] < magCutOff)
        {
            bBox.min()[cmpt] -= magCutOff;
            bBox.max()[cmpt] += magCutOff;
        }
    }

    searchTreePtr_ =
    (
        new SearchTreeType
        (
            mappingTreeData(candidatePoints_),
            bBox,
            8, 10.0, 5.0
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
faceSetAlgorithm::faceSetAlgorithm
(
    const polyMesh& mesh,
    const pointField& newPoints,
    const UList<edge>& newEdges,
    const UList<face>& newFaces,
    const UList<cell>& newCells,
    const UList<label>& newOwner,
    const UList<label>& newNeighbour
)
:
    convexSetAlgorithm
    (
        mesh,
        newPoints,
        newEdges,
        newFaces,
        newCells,
        newOwner,
        newNeighbour
    ),
    candidatePoints_
    (
        mesh.faceCentres(),
        mesh.nFaces() - mesh.nInternalFaces(),
        mesh.nInternalFaces()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void faceSetAlgorithm::computeNormFactor(const label index) const
{
    // Compute refNorm and normFactor
    refNorm_ = newFaces_[index].normal(newPoints_);
    refCentre_ = newFaces_[index].centre(newPoints_);
    normFactor_ = mag(refNorm_);

    // Normalize for later use
    refNorm_ /= normFactor_ + VSMALL;

    // Compute a bounding box around the face
    box_ = treeBoundBox(newFaces_[index].points(newPoints_));

    // Axis-aligned faces can result in degenerate bounding boxes.
    // If so, artificially inflate in the degenerate dimension
    const scalar cutOff = 0.01;
    const vector boxSpan = box_.span();
    const scalar magSpan = mag(boxSpan);
    const scalar magCutOff = (cutOff * magSpan);

    for (label cmpt = 0; cmpt < pTraits<vector>::nComponents; cmpt++)
    {
        if (boxSpan[cmpt] < magCutOff)
        {
            box_.min()[cmpt] -= magSpan;
            box_.max()[cmpt] += magSpan;
        }
    }

    vector minToXb = (box_.min() - box_.midpoint());
    vector maxToXb = (box_.max() - box_.midpoint());

    // Scale it by a bit
    box_.min() += (2.0 * minToXb);
    box_.max() += (2.0 * maxToXb);
}


// Find the nearest mapping candidates
void faceSetAlgorithm::findMappingCandidates(labelList& mapCandidates) const
{
    // Clear existing fields
    parents_.clear();

    centres_.clear();
    weights_.clear();

    // Clear the input list
    mapCandidates.clear();

    const SearchTreeType& tree = searchTree();

    // Find all candidates within search box
    mapCandidates = tree.findBox(box_);

    // Since the tree addresses only into boundary faces,
    // offset the index by the number of internal faces
    forAll(mapCandidates, cI)
    {
        mapCandidates[cI] += mesh_.nInternalFaces();
    }
}


// Write out mapping candidates
void faceSetAlgorithm::writeMappingCandidates() const
{
    const SearchTreeType& tree = searchTree();

    // Fetch reference to tree points
    const UList<point>& points = tree.shapes().points();

    // Write out points to VTK file
    meshOps::writeVTK
    (
        mesh_,
        "faceMappingCandidates",
        identity(points.size()),
        0,
        points
    );
}


// Check whether the bounding box contains the entity
bool faceSetAlgorithm::contains(const label index) const
{
    // Fetch old face
    const face& checkFace = mesh_.faces()[index];

    // Check if the bounding box contains any of the supplied points
    forAll(checkFace, pointI)
    {
        if (box_.contains(mesh_.points()[checkFace[pointI]]))
        {
            return true;
        }
    }

    return false;
}


// Compute intersection
bool faceSetAlgorithm::computeIntersection
(
    const label newIndex,
    const label oldIndex,
    const label offset,
    bool output
) const
{
    bool intersects = false;

    const pointField& newPoints = newPoints_;
    const pointField& oldPoints = mesh_.points();

    const face& newFace = newFaces_[newIndex];
    const face& oldFace = mesh_.faces()[oldIndex];

    // Compute the normal for the old face
    vector oldNorm = oldFace.normal(oldPoints);

    // Normalize
    oldNorm /= mag(oldNorm) + VSMALL;

    if ((oldNorm & refNorm_) < 0.0)
    {
        // Opposite face orientation. Skip it.
        return false;
    }

    // Check if decomposition is necessary
    if (oldFace.size() > 3 || newFace.size() > 3)
    {
        // Decompose new / old faces
        DynamicList<TriPoints> clippingTris(15);
        DynamicList<TriPoints> subjectTris(15);

        label ntOld = 0, ntNew = 0;

        if (oldFace.size() == 3)
        {
            // Add a new entry
            subjectTris.append(FixedList<point, 3>(vector::zero));

            subjectTris[ntOld][0] = oldPoints[oldFace[0]];
            subjectTris[ntOld][1] = oldPoints[oldFace[1]];
            subjectTris[ntOld][2] = oldPoints[oldFace[2]];

            ntOld++;
        }
        else
        {
            // Configure tris from oldFace
            vector ofCentre = oldFace.centre(oldPoints);

            forAll(oldFace, pI)
            {
                // Add a new entry
                subjectTris.append(FixedList<point, 3>(vector::zero));

                subjectTris[ntOld][0] = oldPoints[oldFace[pI]];
                subjectTris[ntOld][1] = oldPoints[oldFace.nextLabel(pI)];
                subjectTris[ntOld][2] = ofCentre;

                ntOld++;
            }
        }

        if (newFace.size() == 3)
        {
            // Add a new entry
            clippingTris.append(FixedList<point, 3>(vector::zero));

            clippingTris[ntNew][0] = newPoints[newFace[0]];
            clippingTris[ntNew][1] = newPoints[newFace[1]];
            clippingTris[ntNew][2] = newPoints[newFace[2]];

            ntNew++;
        }
        else
        {
            // Configure tris from newFace
            vector nfCentre = newFace.centre(newPoints);

            forAll(newFace, pI)
            {
                // Add a new entry
                clippingTris.append(FixedList<point, 3>(vector::zero));

                clippingTris[ntNew][0] = newPoints[newFace[pI]];
                clippingTris[ntNew][1] = newPoints[newFace.nextLabel(pI)];
                clippingTris[ntNew][2] = nfCentre;

                ntNew++;
            }
        }

        // Accumulate area / centroid over all intersections
        bool foundIntersect = false;

        scalar totalArea = 0.0;
        vector totalCentre = vector::zero;

        // Loop through all clipping tris
        forAll(clippingTris, i)
        {
            // Initialize the intersector
            triIntersection intersector(clippingTris[i]);

            // Test for intersection and evaluate
            // against all subject tris
            forAll(subjectTris, j)
            {
                intersects = intersector.evaluate(subjectTris[j]);

                if (intersects)
                {
                    scalar area = 0.0;
                    vector centre = vector::zero;

                    // Fetch area and centre
                    intersector.getAreaAndCentre(area, centre);

                    // Accumulate area / centroid
                    totalArea += area;
                    totalCentre += (area * centre);

                    foundIntersect = true;

                    if (output)
                    {
                        writeVTK
                        (
                            "polygonIntersect_"
                          + Foam::name(newIndex)
                          + '_'
                          + Foam::name(oldIndex)
                          + '_' + Foam::name(i)
                          + '_' + Foam::name(j),
                            intersector.getIntersection()
                        );
                    }
                }
            }
        }

        // Size-up the internal lists
        if (foundIntersect && !output)
        {
            // Normalize centre
            totalCentre /= totalArea + VSMALL;

            // Normalize and check if this is worth it
            if (mag(totalArea/normFactor_) > SMALL)
            {
                // Size-up the internal lists
                meshOps::sizeUpList((oldIndex - offset), parents_);
                meshOps::sizeUpList(totalArea, weights_);
                meshOps::sizeUpList(totalCentre, centres_);

                intersects = true;
            }
            else
            {
                intersects = false;
            }
        }
    }
    else
    {
        // Configure points for clipping triangle
        TriPoints clippingTri(vector::zero);

        // Fill in points
        clippingTri[0] = newPoints[newFace[0]];
        clippingTri[1] = newPoints[newFace[1]];
        clippingTri[2] = newPoints[newFace[2]];

        // Configure points for subject triangle
        TriPoints subjectTri(vector::zero);

        // Fill in points
        subjectTri[0] = oldPoints[oldFace[0]];
        subjectTri[1] = oldPoints[oldFace[1]];
        subjectTri[2] = oldPoints[oldFace[2]];

        // Initialize the intersector
        triIntersection intersector(clippingTri);

        // Test for intersection and evaluate
        intersects = intersector.evaluate(subjectTri);

        if (intersects)
        {
            scalar area;
            vector centre;

            // Fetch area and centre
            intersector.getAreaAndCentre(area, centre);

            // Normalize and check if this is worth it
            if (mag(area/normFactor_) > SMALL)
            {
                if (output)
                {
                    writeVTK
                    (
                        "triIntersectNew_"
                      + Foam::name(newIndex),
                        List<TriPoints>(1, clippingTri)
                    );

                    writeVTK
                    (
                        "triIntersectOld_"
                      + Foam::name(newIndex)
                      + '_'
                      + Foam::name(oldIndex),
                        List<TriPoints>(1, subjectTri)
                    );

                    writeVTK
                    (
                        "triIntersect_"
                      + Foam::name(newIndex)
                      + '_'
                      + Foam::name(oldIndex),
                        intersector.getIntersection()
                    );
                }

                // Size-up the internal lists
                meshOps::sizeUpList((oldIndex - offset), parents_);
                meshOps::sizeUpList(area, weights_);
                meshOps::sizeUpList(centre, centres_);
            }
            else
            {
                intersects = false;
            }
        }
    }

    return intersects;
}


//- Write out tris as a VTK
void faceSetAlgorithm::writeVTK
(
    const word& name,
    const List<TriPoints>& triList
) const
{
    // Fill up all points
    label pI = 0;

    List<face> allTris(triList.size(), face(3));
    pointField allPoints(3 * triList.size());

    forAll(triList, triI)
    {
        allTris[triI][0] = pI;
        allPoints[pI++] = triList[triI][0];

        allTris[triI][1] = pI;
        allPoints[pI++] = triList[triI][1];

        allTris[triI][2] = pI;
        allPoints[pI++] = triList[triI][2];
    }

    // Write out in face-to-node addressing
    meshOps::writeVTK
    (
        mesh_,
        name,
        identity(triList.size()),
        2,
        allPoints,
        List<edge>(0),
        allTris,
        List<cell>(0),
        List<label>(0)
    );
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
