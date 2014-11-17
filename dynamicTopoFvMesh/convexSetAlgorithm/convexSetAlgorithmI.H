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
    convexSetAlgorithm

Description
    Base class for convexSetAlgorithms

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class convexSetAlgorithm::mappingTreeData Declaration
\*---------------------------------------------------------------------------*/

class convexSetAlgorithm::mappingTreeData
{
    // Private data

        // Reference to point list
        const UList<point>& points_;

public:

    //- Run-time type information
    ClassNameNoDebug("mappingTreeData");

    // Public sub classes

        // Class for nearest candidate search
        class findNearestOp;

    // Constructors

        // Construct from point list
        mappingTreeData(const UList<point>& points)
        :
            points_(points)
        {}

        // Construct as copy
        mappingTreeData(const mappingTreeData& mpd)
        :
            points_(mpd.points_)
        {}

    // Member Functions

        // Return the number of shapes
        inline label size() const
        {
            return points_.size();
        }

        // Return whether the shape at index
        // overlaps the bounding box
        inline bool overlaps
        (
            const label& index,
            const treeBoundBox& sampleBox
        ) const
        {
            return sampleBox.contains(points_[index]);
        }

        // Returns whether the shape at index
        // overlaps the sphere
        inline bool overlaps
        (
            const label& index,
            const point& centre,
            const scalar& radiusSqr
        ) const
        {
            const point& p = points_[index];
            const scalar distSqr = magSqr(p - centre);

            return (distSqr <= radiusSqr);
        }

        // Return const reference to points
        inline const UList<point>& points() const
        {
            return points_;
        }
};


/*---------------------------------------------------------------------------*\
      Class convexSetAlgorithm::mappingTreeData::findNearestOp Declaration
\*---------------------------------------------------------------------------*/

class convexSetAlgorithm::mappingTreeData::findNearestOp
{
    // Reference to search tree
    const indexedOctree<mappingTreeData>& tree_;

public:

    findNearestOp(const indexedOctree<mappingTreeData>& tree)
    :
        tree_(tree)
    {}

    inline void operator()
    (
        const labelUList& indices,
        const point& sample,
        scalar& nearestDistSqr,
        label& minIndex,
        point& nearestPoint
    ) const
    {
        const mappingTreeData& shapes = tree_.shapes();
        const UList<point>& points = shapes.points();

        forAll(indices, indexI)
        {
            const label checkIndex = indices[indexI];
            const point& checkPoint = points[checkIndex];
            const scalar distSqr = magSqr(sample - checkPoint);

            if (distSqr < nearestDistSqr)
            {
                minIndex = checkIndex;
                nearestDistSqr = distSqr;
                nearestPoint = checkPoint;
            }
        }
    }
};


}

// ************************************************************************* //
