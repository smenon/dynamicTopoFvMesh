/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/


#include "freeSurface.H"
#include "primitivePatchInterpolation.H"
#include "emptyFaPatch.H"
#include "wedgeFaPatch.H"
#include "PstreamCombineReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<vectorField> freeSurface::pointDisplacement(const scalarField& deltaH) const
{
    const pointField& points = aMesh().patch().localPoints();
    const labelListList& pointFaces = aMesh().patch().pointFaces();

    controlPoints() += facesDisplacementDir()*deltaH;

    tmp<vectorField> tdisplacement
    (
        new vectorField
        (
            points.size(),
            vector::zero
        )
    );

    vectorField& displacement = tdisplacement();


    // Calculate displacement of internal points

    const edgeList& edges = aMesh().patch().edges();

    labelList internalPoints = aMesh().internalPoints();

    forAll (internalPoints, pointI)
    {
        label curPoint = internalPoints[pointI];

        const labelList& curPointFaces = pointFaces[curPoint];

        symmTensor M = symmTensor::zero;

        vector S = vector::zero;


        scalarField w(curPointFaces.size(), 0.0);

        forAll (curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            w[faceI] = 1.0/mag
            (
                controlPoints()[curFace]
              - points[curPoint]
            );
        }

        w /= sum(w);

        forAll (curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            M = M + sqr(w[faceI])*sqr(controlPoints()[curFace]);

            S += sqr(w[faceI])*controlPoints()[curFace];
        }

        vector N = inv(M)&S;

        N /= mag(N);

        scalar p = (S&N)/sum(sqr(w));

        displacement[curPoint] =
            pointsDisplacementDir()[curPoint]*
            (p - (points[curPoint]&N))/
            (pointsDisplacementDir()[curPoint]&N);
    }

    // Calculate displacement of points
    // which belonge to empty and wedge patches
    labelList boundaryPoints = aMesh().boundaryPoints();

    const labelListList& edgeFaces = aMesh().patch().edgeFaces();
    const labelListList& pointEdges = aMesh().patch().pointEdges();

    const vectorField& pointNormals = aMesh().pointAreaNormals();

    forAll (boundaryPoints, pointI)
    {
        label curPoint = boundaryPoints[pointI];

        if (motionPointsMask()[curPoint])
        {
            // Calculating mirror points
            const labelList& curPointEdges = pointEdges[curPoint];

            vectorField mirrorPoints(2, vector::zero);

            label counter = -1;

            forAll (curPointEdges, edgeI)
            {
                label curEdge = curPointEdges[edgeI];

                if(edgeFaces[curEdge].size() == 1)
                {
                    label patchID = -1;

                    forAll(aMesh().boundary(), patchI)
                    {
                        // Note: aMesh().boundary()[patchI].size()==0 for
                        // empty faPatch
                        forAll(aMesh().boundary()[patchI], eI)
                        {
                            if (aMesh().boundary()[patchI][eI] == curEdge)
                            {
                                patchID = patchI;
                                break;
                            }
                        }
                    }

                    if
                    (
                        patchID > -1
                     && (
                            aMesh().boundary()[patchID].type()
                         == wedgeFaPatch::typeName
                        )
                    )
                    {
                        const wedgeFaPatch& wedgePatch =
                            refCast<const wedgeFaPatch>
                            (
                                aMesh().boundary()[patchID]
                            );

                        mirrorPoints[++counter] =
                            transform
                            (
                                wedgePatch.faceT(),
                                controlPoints()[edgeFaces[curEdge][0]]
                            );
                    }
                    else
                    {
                        vector nE =
                            pointNormals[edges[curEdge].start()]
                          + pointNormals[edges[curEdge].end()];

                        nE /= mag(nE);

                        vector eP =
                            controlPoints()[edgeFaces[curEdge][0]]
                          - edges[curEdge].centre(points);

                        mirrorPoints[++counter] =
                            edges[curEdge].centre(points)
                          - ((I - 2.0*nE*nE)&eP);
                    }
                }
            }


            // Calculating LS plane interpolation
            const labelList& curPointFaces = pointFaces[curPoint];

            symmTensor M = symmTensor::zero;

            vector S = vector::zero;

            scalarField w(curPointFaces.size() + 2, 0.0);

            forAll (curPointFaces, faceI)
            {
                label curFace = curPointFaces[faceI];

                w[faceI] = 1.0/mag
                (
                    controlPoints()[curFace]
                  - points[curPoint]
                );
            }

            forAll (mirrorPoints, pI)
            {
                w[curPointFaces.size() + pI] = 1.0/mag
                (
                    mirrorPoints[pI]
                  - points[curPoint]
                );
            }

            w /= sum(w);


            forAll (curPointFaces, faceI)
            {
                label curFace = curPointFaces[faceI];

                M = M + sqr(w[faceI])*sqr(controlPoints()[curFace]);

                S += sqr(w[faceI])*controlPoints()[curFace];
            }


            forAll (mirrorPoints, pI)
            {
                M = M + sqr(w[curPointFaces.size()+pI])*sqr(mirrorPoints[pI]);

                S += sqr(w[curPointFaces.size()+pI])*mirrorPoints[pI];
            }

            vector N = inv(M)&S;

            N /= mag(N);

            scalar p = (S&N)/sum(sqr(w));

            displacement[curPoint] =
                pointsDisplacementDir()[curPoint]*
                (p - (points[curPoint]&N))/
                (pointsDisplacementDir()[curPoint]&N);
        }
    }


    forAll(aMesh().boundary(), patchI)
    {
        bool fixedPatch = false;

        forAll(fixedFreeSurfacePatches_, fpI)
        {
            label fixedPatchID = aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[fpI]
            );

            if(fixedPatchID == -1)
            {
                FatalErrorIn("freeSurface::pointDisplacemen(...)")
                    << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                        << " defined in the freeSurfaceProperties dictionary"
                        << abort(FatalError);
            }

            if (fixedPatchID == patchI)
            {
                fixedPatch = true;
            }
        }

        if
        (
            (
                aMesh().boundary()[patchI].type()
             != emptyFaPatch::typeName
            )
         && (
                aMesh().boundary()[patchI].type()
             != wedgeFaPatch::typeName
            )
         && (
                aMesh().boundary()[patchI].type()
             != processorFaPatch::typeName
            )
         && !fixedPatch
        )
        {
            labelList patchPoints =
                aMesh().boundary()[patchI].pointLabels();

            labelListList patchPointEdges =
                aMesh().boundary()[patchI].pointEdges();

            unallocLabelList patchEdgeFaces =
                aMesh().boundary()[patchI].edgeFaces();

            forAll(patchPoints, pointI)
            {
                forAll(patchPointEdges[pointI], edgeI)
                {
                    label curEdge = patchPointEdges[pointI][edgeI];

                    displacement[patchPoints[pointI]] +=
                        pointsDisplacementDir()[patchPoints[pointI]]
                       *(
                            pointsDisplacementDir()[patchPoints[pointI]]
                           &(
                                controlPoints()[patchEdgeFaces[curEdge]]
                              - points[patchPoints[pointI]]
                            )
                        );
                }

                displacement[patchPoints[pointI]] /=
                    patchPointEdges[pointI].size();
            }
        }
    }


    // Calculate displacement of axis point
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            if(wedgePatch.axisPoint() > -1)
            {
                label axisPoint = wedgePatch.axisPoint();

                displacement[axisPoint] =
                    pointsDisplacementDir()[axisPoint]
                   *(
                        pointsDisplacementDir()[axisPoint]
                       &(
                            controlPoints()[pointFaces[axisPoint][0]]
                          - points[axisPoint]
                        )
                    );
            }
        }
    }


    // Calculate displacement of processor patch points
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == processorFaPatch::typeName
        )
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(aMesh().boundary()[patchI]);

            const labelList& patchPointLabels =
                procPatch.pointLabels();

            Field<symmTensor> ownM(patchPointLabels.size(), symmTensor::zero);
            vectorField ownS(patchPointLabels.size(), vector::zero);
            scalarField ownSumW(patchPointLabels.size(), 0);
            scalarField ownSumSqrW(patchPointLabels.size(), 0);

            const labelList& nonGlobalPatchPoints =
                procPatch.nonGlobalPatchPoints();

            forAll(nonGlobalPatchPoints, pointI)
            {
                label curPatchPoint =
                    nonGlobalPatchPoints[pointI];

                label curPoint =
                    patchPointLabels[curPatchPoint];

                const labelList& curPointFaces = pointFaces[curPoint];

                scalarField w(curPointFaces.size(), 0.0);

                forAll (curPointFaces, faceI)
                {
                    label curFace = curPointFaces[faceI];

                    w[faceI] = 1.0/mag
                    (
                        controlPoints()[curFace]
                      - points[curPoint]
                    );

                    ownM[curPatchPoint] +=
                        sqr(w[faceI])*sqr(controlPoints()[curFace]);

                    ownS[curPatchPoint] +=
                        sqr(w[faceI])*controlPoints()[curFace];
                }

#               include "emptyProcessorFaPatchPoints.H"

                ownSumW[curPatchPoint] = sum(w);
                ownSumSqrW[curPatchPoint] = sum(sqr(w));
            }

            // Parallel data exchange
            {
                OPstream toNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(),
                    ownM.size()*sizeof(symmTensor)
                  + ownS.size()*sizeof(vector)
                  + 4*ownSumW.size()*sizeof(scalar)
                );

                toNeighbProc << ownM << ownS << ownSumW << ownSumSqrW;
            }

            Field<symmTensor> ngbM
            (
                patchPointLabels.size(),
                symmTensor::zero
            );
            vectorField ngbS
            (
                patchPointLabels.size(),
                vector::zero
            );
            scalarField ngbSumW
            (
                patchPointLabels.size(),
                0
            );
            scalarField ngbSumSqrW
            (
                patchPointLabels.size(),
                0
            );

            {
                IPstream fromNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(),
                    ngbM.size()*sizeof(symmTensor)
                  + ngbS.size()*sizeof(vector)
                  + 4*ngbSumW.size()*sizeof(scalar)
                );

                fromNeighbProc >> ngbM >> ngbS >> ngbSumW >> ngbSumSqrW;
            }


            forAll(nonGlobalPatchPoints, pointI)
            {
                label curPatchPoint =
                    nonGlobalPatchPoints[pointI];

                label curPoint =
                    patchPointLabels[curPatchPoint];

                label curNgbPoint = procPatch.neighbPoints()[curPatchPoint];
                ownM[curPatchPoint] += ngbM[curNgbPoint];
                ownS[curPatchPoint] += ngbS[curNgbPoint];
                ownSumW[curPatchPoint] += ngbSumW[curNgbPoint];
                ownSumSqrW[curPatchPoint] += ngbSumSqrW[curNgbPoint];

                ownM[curPatchPoint] /= sqr(ownSumW[curPatchPoint]);
                ownS[curPatchPoint] /= sqr(ownSumW[curPatchPoint]);

                vector N = inv(ownM[curPatchPoint])&ownS[curPatchPoint];
                N /= mag(N);

                scalar p = (ownS[curPatchPoint]&N)
                    /(ownSumSqrW[curPatchPoint]/sqr(ownSumW[curPatchPoint]));

                displacement[curPoint] =
                    pointsDisplacementDir()[curPoint]*
                    (p - (points[curPoint]&N))
                   /(pointsDisplacementDir()[curPoint]&N);
            }
        }
    }


    // Calculate displacement of global processor patch points
    if (aMesh().globalData().nGlobalPoints() > 0)
    {
        const labelList& spLabels =
            aMesh().globalData().sharedPointLabels();

        Field<symmTensor> M(spLabels.size(), symmTensor::zero);
        vectorField S(spLabels.size(), vector::zero);
        scalarField sumW(spLabels.size(), 0);
        scalarField sumSqrW(spLabels.size(), 0);

        forAll (spLabels, pointI)
        {
            label curPoint = spLabels[pointI];

            const labelList& curPointFaces = pointFaces[curPoint];

            scalarField w(curPointFaces.size(), 0.0);

            forAll (curPointFaces, faceI)
            {
                label curFace = curPointFaces[faceI];

                w[faceI] = 1.0/mag
                    (
                        controlPoints()[curFace]
                      - points[curPoint]
                    );

                M[pointI] =
                    M[pointI]
                  + sqr(w[faceI])*sqr(controlPoints()[curFace]);

                    S[pointI] += sqr(w[faceI])*controlPoints()[curFace];
            }

            sumW[pointI] = sum(w);
            sumSqrW[pointI] = sum(sqr(w));
        }

        const labelList& addr = aMesh().globalData().sharedPointAddr();
        label nGlobalPoints = aMesh().globalData().nGlobalPoints();

        symmTensorField gpM(nGlobalPoints, symmTensor::zero);
        vectorField gpS(nGlobalPoints, vector::zero);
        scalarField gpSumW(nGlobalPoints, 0.0);
        scalarField gpSumSqrW(nGlobalPoints, 0.0);

        forAll (addr, i)
        {
            gpM[addr[i]] += M[i];
            gpS[addr[i]] += S[i];
            gpSumW[addr[i]] += sumW[i];
            gpSumSqrW[addr[i]] += sumSqrW[i];
        }

        combineReduce(gpM, plusEqOp<symmTensorField >());
        combineReduce(gpS, plusEqOp<vectorField >());
        combineReduce(gpSumW, plusEqOp<scalarField >());
        combineReduce(gpSumSqrW, plusEqOp<scalarField >());

        // Extract local data
        forAll (addr, i)
        {
            M[i] = gpM[addr[i]];
            S[i] = gpS[addr[i]];
            sumW[i] = gpSumW[addr[i]];
            sumSqrW[i] = gpSumSqrW[addr[i]];
        }

        M /= sqr(sumW);
        S /= sqr(sumW);

        vectorField N = inv(M)&S;
        N /= mag(N);

        scalarField p = (S&N)/(sumSqrW/sqr(sumW));

        forAll(spLabels, pointI)
        {
            label curPoint = spLabels[pointI];

            displacement[curPoint] =
                pointsDisplacementDir()[curPoint]*
                (p[pointI] - (points[curPoint]&N[pointI]))
               /(pointsDisplacementDir()[curPoint]&N[pointI]);
        }
    }

    return tdisplacement;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
