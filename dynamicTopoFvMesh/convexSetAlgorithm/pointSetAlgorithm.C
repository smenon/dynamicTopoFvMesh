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
    pointSetAlgorithm

Description
    Class for convexSetAlgorithms pertaining to points

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "pointSetAlgorithm.H"

#include "meshOps.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Construct the search tree
void pointSetAlgorithm::constructSearchTree() const
{
    if (searchTreePtr_)
    {
        FatalErrorIn("void pointSetAlgorithm::constructSearchTree() const")
            << "searchTree is already allocated."
            << abort(FatalError);
    }

    searchTreePtr_ =
    (
        new SearchTreeType
        (
            mappingTreeData(mesh_.points()),
            treeBoundBox(mesh_.points()).extend(random_, 0.001),
            8, 10.0, 5.0
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pointSetAlgorithm::pointSetAlgorithm
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
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pointSetAlgorithm::computeNormFactor(const label index) const
{
    const SearchTreeType& tree = searchTree();

    pointIndexHit pHit = tree.findNearest(newPoints_[index], GREAT);

    if (pHit.hit())
    {
        const label pIndex = pHit.index();

        // Store the closest point and square distance to it
        refCentre_ = newPoints_[index];
        refNorm_ = mesh_.points()[pIndex];
        normFactor_ = Foam::magSqr(refCentre_ - refNorm_);

        // Scale it up by a bit
        normFactor_ *= 1.5;
    }
    else
    {
        FatalErrorIn("void pointSetAlgorithm::computeNormFactor() const")
            << " Unable to find an initial candidate for"
            << " index: " << index
            << " point: " << newPoints_[index]
            << abort(FatalError);
    }
}


// Find the nearest mapping candidates
void pointSetAlgorithm::findMappingCandidates(labelList& mapCandidates) const
{
    // Clear existing fields
    parents_.clear();

    centres_.clear();
    weights_.clear();

    // Clear the input list
    mapCandidates.clear();

    const SearchTreeType& tree = searchTree();

    // Find all candidates within search radius
    mapCandidates = tree.findSphere(refCentre_, normFactor_);

    const pointField& meshPoints = mesh_.points();
    const label nCandidates = mapCandidates.size();

    parents_.setSize(nCandidates);
    centres_.setSize(nCandidates);
    weights_.setSize(nCandidates);

    // Compute inverse distance weights to each candidate
    forAll(mapCandidates, indexI)
    {
        const label pIndex = mapCandidates[indexI];
        const point& mapPoint = meshPoints[pIndex];
        const scalar sqrDist = Foam::magSqr(refCentre_ - mapPoint);

        parents_[indexI] = pIndex;
        centres_[indexI] = mapPoint;
        weights_[indexI] = (1.0 / stabilise(sqrDist, VSMALL));
    }
}


// Write out mapping candidates
void pointSetAlgorithm::writeMappingCandidates() const
{
    const SearchTreeType& tree = searchTree();

    // Fetch reference to tree points
    const UList<point>& points = tree.shapes().points();

    // Write out points to VTK file
    meshOps::writeVTK
    (
        mesh_,
        "pointMappingCandidates",
        identity(points.size()),
        0,
        points
    );
}


// Check whether the bounding box contains the entity
bool pointSetAlgorithm::contains(const label index) const
{
    // Fetch old points
    const pointField& points = mesh_.points();

    // Check if the bounding sphere contains any of the supplied points
    if (Foam::magSqr(refCentre_ - points[index]) < normFactor_)
    {
        return true;
    }

    return false;
}


// Compute intersection
bool pointSetAlgorithm::computeIntersection
(
    const label newIndex,
    const label oldIndex,
    const label offset,
    bool output
) const
{
    return false;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
