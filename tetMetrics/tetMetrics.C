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

Foam::scalar Dihedral
(
    const Foam::point& p0, 
    const Foam::point& p1, 
    const Foam::point& p2,
    const Foam::point& p3
)
{
    return 0.0;
}

// Tetrahedral mesh-quality metric suggested by Knupp [2003].
Foam::scalar Knupp
(
    const Foam::point& p0, 
    const Foam::point& p1, 
    const Foam::point& p2,
    const Foam::point& p3
)
{
    // Obtain signed tet volume
    Foam::scalar V = (1.0/6.0)*(((p1 - p0) ^ (p2 - p0)) & (p3 - p0));

    // Obtain the magSqr edge-lengths
    Foam::scalar Le =   ((p1-p0) & (p1-p0))
                      + ((p2-p0) & (p2-p0))
                      + ((p3-p0) & (p3-p0))
                      + ((p2-p1) & (p2-p1))
                      + ((p3-p1) & (p3-p1))
                      + ((p3-p2) & (p3-p2));

    // Return signed quality
    return Foam::sign(V)*((24.96100588*Foam::pow(V*V,0.333333))/Le);
}

// Mean Ratio Tetrahedral mesh metric
// Liu,A. and Joe, B., “Relationship between tetrahedron shape measures,”
// BIT Numerical Mathematics, Vol. 34, No. 2, 1994, pp. 268–287.
Foam::scalar meanRatio
(
    const Foam::point& p0,
    const Foam::point& p1,
    const Foam::point& p2,
    const Foam::point& p3
)
{
    // Obtain signed tet volume
    Foam::scalar V = (1.0/6.0)*(((p1 - p0) ^ (p2 - p0)) & (p3 - p0));

    // Obtain the magSqr edge-lengths
    Foam::scalar Le =   ((p1-p0) & (p1-p0))
                      + ((p2-p0) & (p2-p0))
                      + ((p3-p0) & (p3-p0))
                      + ((p2-p1) & (p2-p1))
                      + ((p3-p1) & (p3-p1))
                      + ((p3-p2) & (p3-p2));

    // Return signed quality
    return  Foam::sign(V)*(12.0*Foam::pow(3.0*V*V,0.333333)/Le);
}

// Tetrahedral mesh-metric based on the Frobenius Condition Number
// Patrick M. Knupp. Matrix Norms & the Condition Number: A General Framework
// to Improve Mesh Quality via Node-Movement. Eighth International Meshing
// Roundtable (Lake Tahoe, California), pages 13–22, October 1999.
Foam::scalar Frobenius
(
    const Foam::point& p0,
    const Foam::point& p1,
    const Foam::point& p2,
    const Foam::point& p3
)
{
    // Obtain signed tet volume
    Foam::scalar V = (1.0/6.0)*(((p1 - p0) ^ (p2 - p0)) & (p3 - p0));

    // Obtain the magSqr edge-lengths
    Foam::scalar Le =   ((p1-p0) & (p1-p0))
                      + ((p2-p0) & (p2-p0))
                      + ((p3-p0) & (p3-p0))
                      + ((p2-p1) & (p2-p1))
                      + ((p3-p1) & (p3-p1))
                      + ((p3-p2) & (p3-p2));

    // Compute magSqr of face-areas
    Foam::scalar A =    Foam::magSqr(0.5*((p1-p0) ^ (p2-p0)))
                      + Foam::magSqr(0.5*((p1-p0) ^ (p3-p0)))
                      + Foam::magSqr(0.5*((p2-p0) ^ (p3-p0)))
                      + Foam::magSqr(0.5*((p3-p1) ^ (p2-p1)));

    // Return signed quality
    return 3.67423461*(V/Foam::sqrt((Le/6.0)*(A/4.0)));
}

// Tetrahedral mesh-metric suggested by:
// V. N. Parthasarathy, C. M. Graichen, and A. F. Hathaway.
// Fast Evaluation & Improvement of Tetrahedral 3-D Grid Quality. [1991]
Foam::scalar PGH
(
    const Foam::point& p0,
    const Foam::point& p1,
    const Foam::point& p2,
    const Foam::point& p3
)
{
    // Obtain signed tet volume
    Foam::scalar V = (1.0/6.0)*(((p1 - p0) ^ (p2 - p0)) & (p3 - p0));

    // Obtain the magSqr edge-lengths
    Foam::scalar Le =   ((p1-p0) & (p1-p0))
                      + ((p2-p0) & (p2-p0))
                      + ((p3-p0) & (p3-p0))
                      + ((p2-p1) & (p2-p1))
                      + ((p3-p1) & (p3-p1))
                      + ((p3-p2) & (p3-p2));

    // Return signed quality
    return 8.48528137*(V/Foam::pow(Le/4.0,1.5));
}

// Metric suggested by:
// Hugues L. de Cougny, Mark S. Shephard, and Marcel K. Georges.
// Explicit Node Point Smoothing Within Octree. Technical Report 10-1990,
// Scientiﬁc Computation Research Center, Rensselaer Polytechnic Institute,
// Troy, New York. [1990]
Foam::scalar CSG
(
    const Foam::point& p0,
    const Foam::point& p1,
    const Foam::point& p2,
    const Foam::point& p3
)
{
    // Obtain signed tet volume
    Foam::scalar V = (1.0/6.0)*(((p1 - p0) ^ (p2 - p0)) & (p3 - p0));

    // Compute magSqr of face-areas
    Foam::scalar A =    Foam::magSqr(0.5*((p1-p0) ^ (p2-p0)))
                      + Foam::magSqr(0.5*((p1-p0) ^ (p3-p0)))
                      + Foam::magSqr(0.5*((p2-p0) ^ (p3-p0)))
                      + Foam::magSqr(0.5*((p3-p1) ^ (p2-p1)));

    // Return signed quality
    return 6.83852117*(V/Foam::pow(A,0.75));
}

// ************************************************************************* //
