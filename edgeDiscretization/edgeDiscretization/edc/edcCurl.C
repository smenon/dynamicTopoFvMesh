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

Description
    

\*---------------------------------------------------------------------------*/

#include "edcCurl.H"
#include "eMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace edc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void curl
(
    const GeometricField<Type, ePatchField, EdgeMesh>& ef,
    GeometricField<Type, fvsPatchField, surfaceMesh>& curlField
)
{
    // Obtain mesh references
    const fvMesh& mesh = curlField.mesh();
    const polyBoundaryMesh& bmesh = mesh.boundaryMesh();
    const eMesh& msh = ef.mesh();
    const eBoundaryMesh& bmsh = msh.boundary();

    // Get geometric data from the mesh
    const pointField& points = mesh.points();
    const vectorField& Xf = mesh.faceCentres();
    const vectorField& Nf = mesh.faceAreas();

    // Obtain addressing
    const edgeList& edges = msh.edges();
    const labelListList& edgeFaces = msh.edgeFaces();

    // Obtain references
    const Field<Type>& iF = ef.internalField();
    Field<Type>& iFs = curlField.internalField();

    label i;
    scalar fSign = 0.0;
    iFs = pTraits<Type>::zero;

    // Do internal edges first
    for(i = 0; i < msh.nInternalEdges(); i++)
    {
        const labelList& eFaces = edgeFaces[i];

        forAll(eFaces, faceI)
        {
            // Compute the face sign
            fSign = Foam::sign
                    (
                        (
                            edges[i].vec(points)
                          ^ (Xf[eFaces[faceI]] - edges[i].centre(points))
                        ) & Nf[eFaces[faceI]]
                    );

            iFs[i] += iF[eFaces[faceI]]*fSign;
        }
    }

    // Initialize all boundary fields
    forAll(bmesh, patchI)
    {
        fvsPatchField<Type>& pFs = curlField.boundaryField()[patchI];
        pFs = pTraits<Type>::zero;
    }

    // Now do boundary patches
    label pIndex = -1, fIndex = -1;

    forAll(bmsh, patchI)
    {
        label start = bmsh[patchI].start();

        const ePatchField<Type>& pF = ef.boundaryField()[patchI];

        forAll(bmsh[patchI], edgeI)
        {
            i = start + edgeI;

            const labelList& eFaces = edgeFaces[i];

            forAll(eFaces, faceI)
            {
                fIndex = eFaces[faceI];

                // Compute the face sign
                fSign = Foam::sign
                        (
                            (
                                edges[i].vec(points)
                              ^ (Xf[fIndex] - edges[i].centre(points))
                            ) & Nf[fIndex]
                        );

                if ((pIndex = bmesh.whichPatch(fIndex)) == -1)
                {
                    // Obtain value from the internal field
                    iFs[fIndex] += pF[edgeI]*fSign;
                }
                else
                {
                    // Obtain value from the patch field
                    curlField.boundaryField()[pIndex][fIndex-start] +=
                        pF[edgeI]*fSign;
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace edc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
