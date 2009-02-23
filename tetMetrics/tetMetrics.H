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
    Header file of tetrahedral mesh-quality metrics.

Author
    Sandeep Menon

SourceFiles
    tetMetrics.C

\*---------------------------------------------------------------------------*/

#ifndef tetMetrics_H
#define tetMetrics_H

#include "point.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

extern "C" Foam::scalar Dihedral
(
    const Foam::point& p0, 
    const Foam::point& p1, 
    const Foam::point& p2,
    const Foam::point& p3
);

extern "C" Foam::scalar Knupp
(
    const Foam::point& p0, 
    const Foam::point& p1, 
    const Foam::point& p2,
    const Foam::point& p3
);

extern "C" Foam::scalar meanRatio
(
    const Foam::point& p0,
    const Foam::point& p1,
    const Foam::point& p2,
    const Foam::point& p3
);

extern "C" Foam::scalar cubicMeanRatio
(
    const Foam::point& p0,
    const Foam::point& p1,
    const Foam::point& p2,
    const Foam::point& p3
);

extern "C" Foam::scalar Frobenius
(
    const Foam::point& p0,
    const Foam::point& p1,
    const Foam::point& p2,
    const Foam::point& p3
);

extern "C" Foam::scalar PGH
(
    const Foam::point& p0,
    const Foam::point& p1,
    const Foam::point& p2,
    const Foam::point& p3
);

extern "C" Foam::scalar CSG
(
    const Foam::point& p0,
    const Foam::point& p1,
    const Foam::point& p2,
    const Foam::point& p3
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
