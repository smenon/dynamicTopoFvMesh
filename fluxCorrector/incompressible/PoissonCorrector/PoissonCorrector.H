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
    PoissonCorrector

Description
    Flux-correction after topo-changes, using a Poisson solver.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    PoissonCorrector.C

\*---------------------------------------------------------------------------*/

#ifndef PoissonCorrector_H
#define PoissonCorrector_H

#include "Switch.H"
#include "fluxCorrector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PoissonCorrector Declaration
\*---------------------------------------------------------------------------*/

class PoissonCorrector
:
    public fluxCorrector
{
    // Private data

        //- Is flux-correction required?
        Switch required_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        PoissonCorrector(const PoissonCorrector&);

        //- Disallow default bitwise assignment
        void operator=(const PoissonCorrector&);

public:

    //- Runtime type information
    TypeName("Poisson");

    // Constructors

        //- Construct from fvMesh and dictionary
        PoissonCorrector(const fvMesh& mesh, const dictionary& dict);


    // Destructor

        virtual ~PoissonCorrector()
        {}


    // Member Functions

        //- Is flux-correction required?
        virtual bool required() const;

        //- Interpolate fluxes to a specified list of faces
        virtual void interpolateFluxes(const labelList& faces) const;

        //- Update fluxes in the registry
        virtual void updateFluxes() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
