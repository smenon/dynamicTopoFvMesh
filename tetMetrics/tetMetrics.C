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
    tetMetrics

Description
    Implementation of tetrahedral mesh-quality metrics.

    Most measures are referenced from:
    Jonathan R. Shewchuk. What is a good linear finite element?
    Interpolation, Conditioning, Anisotropy and Quality Measures.
    Eleventh International Meshing Roundtable (Ithaca, NY), pp 115-126, 2002.

    NOTE:
    For metrics involving volume the function assumes points (p0-p1-p2)
    are in counter-clockwise fashion when viewed from p3.

Author
    Sandeep Menon

\*----------------------------------------------------------------------------*/

#include "tetMetrics.H"

#include "point.H"
#include "scalar.H"
#include "HashSet.H"

namespace Foam
{

// Enumeration for tets
label tetEnum[6][4] =
{
    {0,1,2,3},
    {0,2,3,1},
    {0,3,1,2},
    {1,2,0,3},
    {1,3,0,2},
    {2,3,0,1}
};

void reportMetrics()
{
	wordHashSet availableMetrics_;

	availableMetrics_.insert("Dihedral");
	availableMetrics_.insert("Knupp");
	availableMetrics_.insert("meanRatio");
	availableMetrics_.insert("cubicMeanRatio");
	availableMetrics_.insert("Frobenius");
	availableMetrics_.insert("PGH");
	availableMetrics_.insert("CSG");

	Info << availableMetrics_ << endl;
}

// Minimum dihedral angle among six edges of the tetrahedron. Normalized
// by 70.529 degrees (equilateral tet) and signed by volume.
scalar Dihedral
(
    const point& p0,
    const point& p1,
    const point& p2,
    const point& p3
)
{
    scalar minAngle = 0.0;
    FixedList<scalar,6> cosAngles(1.0);
    FixedList<vector,4> pts(vector::zero);

    // Assign point-positions
    pts[0] = p0;
    pts[1] = p1;
    pts[2] = p2;
    pts[3] = p3;

    // Permute over all six edges
    for (label i = 0; i < 6; i++)
    {
        // Normalize the axis
        vector v0 = (pts[tetEnum[i][1]] - pts[tetEnum[i][0]])
                    /mag(pts[tetEnum[i][1]] - pts[tetEnum[i][0]]);

        // Obtain plane-vectors
        vector v1 = (pts[tetEnum[i][2]] - pts[tetEnum[i][0]]);
        vector v2 = (pts[tetEnum[i][3]] - pts[tetEnum[i][0]]);

        v1 -= (v1 & v0)*v0;
        v2 -= (v2 & v0)*v0;

        cosAngles[i] = acos((v1/mag(v1)) & (v2/mag(v2)));

        // Compute minimum angle on-the-fly
        if (i > 0)
        {
            minAngle = minAngle < cosAngles[i] ? minAngle : cosAngles[i];
        }
        else
        {
            minAngle = cosAngles[0];
        }
    }

    // Compute signed volume and multiply by the normalized angle
    return sign(((p1 - p0) ^ (p2 - p0)) & (p3 - p0))*(minAngle/1.2309632);
}

// Tetrahedral mesh-quality metric suggested by Knupp [2003].
scalar Knupp
(
    const point& p0,
    const point& p1,
    const point& p2,
    const point& p3
)
{
    // Obtain signed tet volume
    scalar V = (1.0/6.0)*(((p1 - p0) ^ (p2 - p0)) & (p3 - p0));

    // Obtain the magSqr edge-lengths
    scalar Le = ((p1-p0) & (p1-p0))
              + ((p2-p0) & (p2-p0))
              + ((p3-p0) & (p3-p0))
              + ((p2-p1) & (p2-p1))
              + ((p3-p1) & (p3-p1))
              + ((p3-p2) & (p3-p2));

    // Return signed quality
    return sign(V)*((24.96100588*::cbrt(V*V))/Le);
}

// Mean Ratio Tetrahedral mesh metric
// Liu,A. and Joe, B., “Relationship between tetrahedron shape measures,”
// BIT Numerical Mathematics, Vol. 34, No. 2, 1994, pp. 268–287.
scalar meanRatio
(
    const point& p0,
    const point& p1,
    const point& p2,
    const point& p3
)
{
    // Obtain signed tet volume
    scalar V = (1.0/6.0)*(((p1 - p0) ^ (p2 - p0)) & (p3 - p0));

    // Obtain the magSqr edge-lengths
    scalar Le = ((p1-p0) & (p1-p0))
              + ((p2-p0) & (p2-p0))
              + ((p3-p0) & (p3-p0))
              + ((p2-p1) & (p2-p1))
              + ((p3-p1) & (p3-p1))
              + ((p3-p2) & (p3-p2));

    // Return signed quality
    return  sign(V)*(12.0*::cbrt(3.0*V*V)/Le);
}

// Cubic Mean Ratio Tetrahedral mesh metric
// Liu,A. and Joe, B., “On the shape of tetrahedra from bisection”
// Mathematics of Computation, Vol. 63, 1994, pp. 141–154.
scalar cubicMeanRatio
(
    const point& p0,
    const point& p1,
    const point& p2,
    const point& p3
)
{
    // Obtain signed tet volume
    scalar V = (1.0/6.0)*(((p1 - p0) ^ (p2 - p0)) & (p3 - p0));

    // Obtain the magSqr edge-lengths
    scalar Le = ((p1-p0) & (p1-p0))
              + ((p2-p0) & (p2-p0))
              + ((p3-p0) & (p3-p0))
              + ((p2-p1) & (p2-p1))
              + ((p3-p1) & (p3-p1))
              + ((p3-p2) & (p3-p2));

    // Return signed quality
    return  sign(V)*((15552.0*V*V)/(Le*Le*Le));
}

// Tetrahedral mesh-metric based on the Frobenius Condition Number
// Patrick M. Knupp. Matrix Norms & the Condition Number: A General Framework
// to Improve Mesh Quality via Node-Movement. Eighth International Meshing
// Roundtable (Lake Tahoe, California), pages 13–22, October 1999.
scalar Frobenius
(
    const point& p0,
    const point& p1,
    const point& p2,
    const point& p3
)
{
    // Obtain signed tet volume
    scalar V = (1.0/6.0)*(((p1 - p0) ^ (p2 - p0)) & (p3 - p0));

    // Obtain the magSqr edge-lengths
    scalar Le = ((p1-p0) & (p1-p0))
              + ((p2-p0) & (p2-p0))
              + ((p3-p0) & (p3-p0))
              + ((p2-p1) & (p2-p1))
              + ((p3-p1) & (p3-p1))
              + ((p3-p2) & (p3-p2));

    // Compute magSqr of face-areas
    scalar A = magSqr(0.5*((p1-p0) ^ (p2-p0)))
             + magSqr(0.5*((p1-p0) ^ (p3-p0)))
             + magSqr(0.5*((p2-p0) ^ (p3-p0)))
             + magSqr(0.5*((p3-p1) ^ (p2-p1)));

    // Return signed quality
    return 3.67423461*(V/sqrt((Le/6.0)*(A/4.0)));
}

// Tetrahedral mesh-metric suggested by:
// V. N. Parthasarathy, C. M. Graichen, and A. F. Hathaway.
// Fast Evaluation & Improvement of Tetrahedral 3-D Grid Quality. [1991]
scalar PGH
(
    const point& p0,
    const point& p1,
    const point& p2,
    const point& p3
)
{
    // Obtain signed tet volume
    scalar V = (1.0/6.0)*(((p1 - p0) ^ (p2 - p0)) & (p3 - p0));

    // Obtain the magSqr edge-lengths
    scalar Le = ((p1-p0) & (p1-p0))
              + ((p2-p0) & (p2-p0))
              + ((p3-p0) & (p3-p0))
              + ((p2-p1) & (p2-p1))
              + ((p3-p1) & (p3-p1))
              + ((p3-p2) & (p3-p2));

    // Return signed quality
    return 8.48528137*(V/pow(Le/4.0,1.5));
}

// Metric suggested by:
// Hugues L. de Cougny, Mark S. Shephard, and Marcel K. Georges.
// Explicit Node Point Smoothing Within Octree. Technical Report 10-1990,
// Scientiﬁc Computation Research Center, Rensselaer Polytechnic Institute,
// Troy, New York. [1990]
scalar CSG
(
    const point& p0,
    const point& p1,
    const point& p2,
    const point& p3
)
{
    // Obtain signed tet volume
    scalar V = (1.0/6.0)*(((p1 - p0) ^ (p2 - p0)) & (p3 - p0));

    // Compute magSqr of face-areas
    scalar A = magSqr(0.5*((p1-p0) ^ (p2-p0)))
             + magSqr(0.5*((p1-p0) ^ (p3-p0)))
             + magSqr(0.5*((p2-p0) ^ (p3-p0)))
             + magSqr(0.5*((p3-p1) ^ (p2-p1)));

    // Return signed quality
    return 6.83852117*(V/pow(A,0.75));
}

} // End namespace Foam

// ************************************************************************* //
