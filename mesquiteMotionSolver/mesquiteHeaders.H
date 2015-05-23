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
    mesquiteMotionSolver

Description
    Header files for mesquiteMotionSolver

\*---------------------------------------------------------------------------*/

#ifndef mesquiteHeaders_H
#define mesquiteHeaders_H

// Have gcc ignore certain warnings while including mesquite headers
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#   pragma GCC diagnostic ignored "-Wold-style-cast"
#   pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#   pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#include "Mesquite_all_headers.hpp"

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#   pragma GCC diagnostic warning "-Wold-style-cast"
#   pragma GCC diagnostic warning "-Wnon-virtual-dtor"
#   pragma GCC diagnostic warning "-Wunused-but-set-variable"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
