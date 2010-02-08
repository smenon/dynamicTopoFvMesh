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

\*---------------------------------------------------------------------------*/

#include "springMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "wedgePolyPatch.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(springMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        springMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::springMotionSolver::springMotionSolver
(
    const polyMesh& mesh
)
:
    motionSolver(mesh),
    Mesh_(mesh),
    refPoints_
    (
        IOobject
        (
            "refPoints",
            mesh.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh.points()
    ),
    cmpt_(-1),
    pID_(-1),
    fixedY_(0.0)
{
    // Check if points on certain axes need to be maintained
    if (found("fixedY"))
    {
        fixedY_ = readScalar(lookup("fixedY"));
    }

    // Check if any slip patches are specified
    if (found("slipPatches"))
    {
        wordList slipPatches = subDict("slipPatches").toc();

        forAll(slipPatches, wordI)
        {
            word& patchName = slipPatches[wordI];

            forAll(mesh.boundaryMesh(), patchI)
            {
                if (mesh.boundaryMesh()[patchI].name() == patchName)
                {
                    slipPatchIDs_.insert(patchI);
                }
            }
        }
    }

    // Check if a tolerance has been specified
    if (found("tolerance"))
    {
        tolerance_ = readScalar(lookup("tolerance"));
    }
    else
    {
        tolerance_ = 1e-15;
    }

    // Check if multiple sweeps have been requested
    if (found("nSweeps"))
    {
        nSweeps_ = readLabel(lookup("nSweeps"));
    }
    else
    {
        nSweeps_ = 1;
    }

    // Check if weighting factors have been specified
    if (found("edgeSpringWeight"))
    {
        edgeSpringWeight_ = readScalar(lookup("edgeSpringWeight"));
    }
    else
    {
        edgeSpringWeight_ = 1.0;
    }

    if (found("volumeSpringWeight"))
    {
        volumeSpringWeight_ = readScalar(lookup("volumeSpringWeight"));
    }
    else
    {
        volumeSpringWeight_ = 1.0;
    }

    // Check if exponent factors have been specified
    if (found("edgeSpringExponent"))
    {
        edgeSpringExponent_ = readScalar(lookup("edgeSpringExponent"));
    }
    else
    {
        edgeSpringExponent_ = 1.0;
    }

    if (found("volumeSpringExponent"))
    {
        volumeSpringExponent_ = readScalar(lookup("volumeSpringExponent"));
    }
    else
    {
        volumeSpringExponent_ = 1.0;
    }
}

Foam::springMotionSolver::springMotionSolver
(
    const polyMesh& mesh,
    Istream& msData
)
:
    motionSolver(mesh),
    Mesh_(mesh),
    refPoints_
    (
        IOobject
        (
            "refPoints",
            mesh.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh.points()
    ),
    cmpt_(-1),
    pID_(-1),
    fixedY_(0.0)
{
    // Check if points on certain axes need to be maintained
    if (found("fixedY"))
    {
        fixedY_ = readScalar(lookup("fixedY"));
    }

    // Check if any slip patches are specified
    if (found("slipPatches"))
    {
        wordList slipPatches = subDict("slipPatches").toc();

        forAll(slipPatches, wordI)
        {
            word& patchName = slipPatches[wordI];

            forAll(mesh.boundaryMesh(), patchI)
            {
                if (mesh.boundaryMesh()[patchI].name() == patchName)
                {
                    slipPatchIDs_.insert(patchI);
                }
            }
        }
    }

    // Check if a tolerance has been specified
    if (found("tolerance"))
    {
        tolerance_ = readScalar(lookup("tolerance"));
    }
    else
    {
        tolerance_ = 1e-15;
    }

    // Check if multiple sweeps have been requested
    if (found("nSweeps"))
    {
        nSweeps_ = readLabel(lookup("nSweeps"));
    }
    else
    {
        nSweeps_ = 1;
    }

    // Check if weighting factors have been specified
    if (found("edgeSpringWeight"))
    {
        edgeSpringWeight_ = readScalar(lookup("edgeSpringWeight"));
    }
    else
    {
        edgeSpringWeight_ = 1.0;
    }

    if (found("volumeSpringWeight"))
    {
        volumeSpringWeight_ = readScalar(lookup("volumeSpringWeight"));
    }
    else
    {
        volumeSpringWeight_ = 1.0;
    }

    // Check if exponent factors have been specified
    if (found("edgeSpringExponent"))
    {
        edgeSpringExponent_ = readScalar(lookup("edgeSpringExponent"));
    }
    else
    {
        edgeSpringExponent_ = 1.0;
    }

    if (found("volumeSpringExponent"))
    {
        volumeSpringExponent_ = readScalar(lookup("volumeSpringExponent"));
    }
    else
    {
        volumeSpringExponent_ = 1.0;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::springMotionSolver::~springMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Vector dot-product
Foam::scalar Foam::springMotionSolver::dot
(
    const vectorField& f1,
    const vectorField& f2
)
{
    scalar s = 0.0;

    forAll(f1, indexI)
    {
        s += (f1[indexI] & f2[indexI]);
    }

    return s;
}

// Scalar dot-product
Foam::scalar Foam::springMotionSolver::dot
(
    const scalarField& f1,
    const scalarField& f2
)
{
    scalar s = 0.0;

    forAll(f1, indexI)
    {
        s += (f1[indexI] * f2[indexI]);
    }

    return s;
}

// Compute the normalization factor for the matrix
Foam::scalar Foam::springMotionSolver::normFactor
(
    const scalarField& x,
    const scalarField& b,
    const scalarField& w,
    scalarField& tmpField
)
{
    scalar xRef = average(x);

    A(scalarField(x.size(), xRef),tmpField);

    return sum(mag(w - tmpField) + mag(b - tmpField)) + 1.0e-20;
}

Foam::scalar Foam::springMotionSolver::normFactor
(
    const vectorField& x,
    const vectorField& b,
    const vectorField& w,
    vectorField& tmpField
)
{
    vector xRef = average(x);

    A(vectorField(x.size(), xRef),tmpField);

    vectorField nFw = (w - tmpField);
    vectorField nFb = (b - tmpField);

    return cmptSumMag(nFw) + cmptSumMag(nFb) + 1.0e-20;
}

// Component-wise sumMag
Foam::scalar Foam::springMotionSolver::cmptSumMag(const vectorField& field)
{
    scalar cSum = 0.0;

    forAll(field,i)
    {
        cSum += mag(field[i].x()) + mag(field[i].y()) + mag(field[i].z());
    }

    return cSum;
}

Foam::scalar Foam::springMotionSolver::cmptSumMag(const scalarField& field)
{
    return sumMag(field);
}

// Templated CG solver
template <class Type>
Foam::label Foam::springMotionSolver::CG
(
    const Field<Type>& b,
    Field<Type>& p,
    Field<Type>& r,
    Field<Type>& w,
    Field<Type>& x,
    bool preCondition = false
)
{
    // Local variables
    scalar alpha, beta, rho, rhoOld, residual;
    label maxIter = x.size(), iter = 0;

    // Set Dirichlet conditions on the solution field (if any)
    setDirichlet(x);

    A(x,w);

    // Compute the normFactor, using 'r' as scratch-space
    scalar norm = this->normFactor(x,b,w,r);

    r = b - w;

    if (preCondition)
    {
        M(r,w);
        p = w;
    }
    else
    {
        p = r;
    }

    rho = dot(r,p);

    // Obtain the normalized residual
    residual = cmptSumMag(r)/norm;

    Info << " Initial residual: " << residual;

    while ( (iter < maxIter) && (residual > tolerance_) )
    {
        A(p,w);

        alpha = rho / dot(p,w);

        forAll (x, i)
        {
            x[i] += (alpha*p[i]);
            r[i] -= (alpha*w[i]);
        }

        rhoOld = rho;

        if (preCondition)
        {
            M(r,w);
            rho = dot(r,w);
        }
        else
        {
            rho = dot(r,r);
        }

        beta = rho / rhoOld;

        if (preCondition)
        {
            forAll (p, i)
            {
                p[i] = w[i] + (beta*p[i]);
            }
        }
        else
        {
            forAll (p, i)
            {
                p[i] = r[i] + (beta*p[i]);
            }
        }

        // Update the normalized residual
        residual = cmptSumMag(r)/norm;
        iter++;
    }

    Info << " Final residual: " << residual;

    return iter;
}

// Sparse matrix-vector multiply [2D]
void Foam::springMotionSolver::A(const scalarField& p, scalarField& w)
{
    w = 0.0;

    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    // Gradient (n2e)
    forAll(boundary[pID_].edges(),edgeI)
    {
        gradEdge_[edgeI] =   p[boundary[pID_].edges()[edgeI][1]]
                           - p[boundary[pID_].edges()[edgeI][0]];
    }

    // Divergence (e2n)
    forAll(boundary[pID_].edges(),edgeI)
    {
        w[boundary[pID_].edges()[edgeI][0]] += gradEdge_[edgeI];
        w[boundary[pID_].edges()[edgeI][1]] -= gradEdge_[edgeI];
    }

    // Apply boundary conditions
    applyBCs(w);
}

// Sparse matrix-vector multiply [3D]
void Foam::springMotionSolver::A(const vectorField& p, vectorField& w)
{
    w = vector::zero;

    vector gradient = vector::zero;

    // Obtain the edge-list from the Mesh
    const edgeList& edges = mesh().edges();

    // Gradient (n2e) * stiffness
    forAll(edges, edgeI)
    {
        gradEdgeV_[edgeI] =
            stiffness_[edgeI]*(p[edges[edgeI][1]] - p[edges[edgeI][0]]);
    }

    // Divergence (e2n)
    forAll(edges, edgeI)
    {
        w[edges[edgeI][0]] += gradEdgeV_[edgeI];
        w[edges[edgeI][1]] -= gradEdgeV_[edgeI];
    }

    // Add contributions from opposing-face springs
    forAll(w, pointI)
    {
        const faceList& pFaces = pointFaces_[pointI];
        const scalarList& k = pfStiffness_[pointI];
        const scalarList& xi = xi_[pointI];
        const scalarList& eta = eta_[pointI];

        forAll(pFaces, faceI)
        {
            const face& faceToCheck = pFaces[faceI];

            // Edge Gradient to interpolated point
            gradient =
                k[faceI]*
                (
                    (
                        xi[faceI]*p[faceToCheck[0]]
                      + eta[faceI]*p[faceToCheck[1]]
                      + (1.0 - xi[faceI] - eta[faceI])*p[faceToCheck[2]]
                    )
                  - p[pointI]
                );

            // Apply force to pointI
            w[pointI] += gradient;

            // Apply distributed forces to vertices of the opposing face
            w[faceToCheck[0]] -= xi[faceI]*gradient;
            w[faceToCheck[1]] -= eta[faceI]*gradient;
            w[faceToCheck[2]] -= (1.0 - xi[faceI] - eta[faceI])*gradient;
        }
    }

    // Apply boundary conditions
    applyBCs(w);
}

// Preconditioner
void Foam::springMotionSolver::M(const scalarField& r, scalarField& w)
{}

// Preconditioner
void Foam::springMotionSolver::M(const vectorField& r, vectorField& w)
{
    w = r;
}

// Set Dirichlet conditions on the solution field (if any)
void Foam::springMotionSolver::setDirichlet(scalarField &x)
{
    // For wedges, loop through boundaries, and fix y-direction on the axis
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();
    if (boundary[pID_].type() == "wedge" && cmpt_ == 1)
    {
        for
        (
            label i = boundary[pID_].nInternalEdges();
            i < boundary[pID_].edges().size();
            i++
        )
        {
            if
            (
                boundary[pID_].localPoints()[boundary[pID_].edges()[i][0]][1]
                < (fixedY_+SMALL)
             && boundary[pID_].localPoints()[boundary[pID_].edges()[i][1]][1]
                < (fixedY_+SMALL)
            )
            {
                // This edge lies on the axis.
                x[boundary[pID_].edges()[i][0]] = fixedY_;
                x[boundary[pID_].edges()[i][1]] = fixedY_;
            }
        }
    }
}

// Set Dirichlet conditions on the solution field (if any)
void Foam::springMotionSolver::setDirichlet(vectorField &x)
{}

// Apply boundary conditions [2D]
void Foam::springMotionSolver::applyBCs(scalarField &field)
{
    // Loop through boundary edges and blank-out residuals
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();
    for
    (
        label i = boundary[pID_].nInternalEdges();
        i < boundary[pID_].edges().size();
        i++
    )
    {
        // If wedge patches are present, set axis nodes to slip
        if (boundary[pID_].type() == "wedge")
        {
            if
            (
                boundary[pID_].localPoints()[boundary[pID_].edges()[i][0]][1]
                < (fixedY_+SMALL)
             && boundary[pID_].localPoints()[boundary[pID_].edges()[i][1]][1]
                < (fixedY_+SMALL)
            )
            {
                // This edge lies on the axis. If the Y component is being
                // solved for, blank out the residual.
                if (cmpt_ == 1)
                {
                    field[boundary[pID_].edges()[i][0]] = 0.0;
                    field[boundary[pID_].edges()[i][1]] = 0.0;
                }
                continue;
            }
        }
        field[boundary[pID_].edges()[i][0]] = 0.0;
        field[boundary[pID_].edges()[i][1]] = 0.0;
    }
}

// Apply boundary conditions [3D]
void Foam::springMotionSolver::applyBCs(vectorField &field)
{
    // Blank out residuals at boundary nodes
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    forAll(boundary, patchI)
    {
        const labelList& meshPts = boundary[patchI].meshPoints();

        if (slipPatchIDs_.found(patchI))
        {
            const vectorField& pn = boundary[patchI].pointNormals();

            // Subtract the point-normal component
            forAll(meshPts, pointI)
            {
                const vector& n = pn[pointI];
                field[meshPts[pointI]] -= ((field[meshPts[pointI]]&n)*n);
            }
        }
        else
        {
            // No-slip wall. Blank out residual completely.
            forAll(meshPts, pointI)
            {
                field[meshPts[pointI]] = vector::zero;
            }
        }
    }
}

// Compute point-to-opposing face connectivity
void Foam::springMotionSolver::computePointFaces()
{
    // Obtain connectivity from mesh
    const labelListList& pointCells = mesh().pointCells();
    const labelList& owner = mesh().faceOwner();
    const faceList& faces = mesh().faces();
    const cellList& cells = mesh().cells();

    pointFaces_.setSize(pointCells.size());
    pfStiffness_.setSize(pointCells.size());
    xi_.setSize(pointCells.size());
    eta_.setSize(pointCells.size());

    forAll(pointFaces_, pointI)
    {
        pointFaces_[pointI].setSize(pointCells[pointI].size());
        pfStiffness_[pointI].setSize(pointCells[pointI].size());
        xi_[pointI].setSize(pointCells[pointI].size());
        eta_[pointI].setSize(pointCells[pointI].size());

        forAll(pointCells[pointI], cellI)
        {
            // Loop through faces of cells in pointCells and find a face which
            // doesn't contain the point in question
            const cell& cellToCheck = cells[pointCells[pointI][cellI]];

            forAll(cellToCheck, faceI)
            {
                if (faces[cellToCheck[faceI]].which(pointI) < 0)
                {
                    // Assign a correctly oriented face
                    if (pointCells[pointI][cellI] == owner[cellToCheck[faceI]])
                    {
                        pointFaces_[pointI][cellI] =
                            faces[cellToCheck[faceI]].reverseFace();
                    }
                    else
                    {
                        pointFaces_[pointI][cellI] =
                            faces[cellToCheck[faceI]];
                    }

                    break;
                }
            }
        }
    }
}

// Initialize fields for the CG solver
void Foam::springMotionSolver::initCG(label nUnknowns)
{
    if (twoDMotion())
    {
        if (r_.size() != nUnknowns)
        {
            b_.setSize(nUnknowns, 0.0);
            x_.setSize(nUnknowns, 0.0);
            p_.setSize(nUnknowns, 0.0);
            r_.setSize(nUnknowns, 0.0);
            w_.setSize(nUnknowns, 0.0);
        }

        const polyBoundaryMesh& boundary = mesh().boundaryMesh();

        if (gradEdge_.size() != boundary[pID_].edges().size())
        {
            gradEdge_.setSize(boundary[pID_].edges().size(), 0.0);
        }
    }
    else
    {
        vector v1, v2, p, n, xp;
        const label nPoints = mesh().nPoints();
        const label nEdges = mesh().nEdges();
        const edgeList& edges = mesh().edges();

        if (rV_.size() != nUnknowns)
        {
            bV_.setSize(nUnknowns, vector::zero);
            xV_.setSize(nUnknowns, vector::zero);
            pV_.setSize(nUnknowns, vector::zero);
            rV_.setSize(nUnknowns, vector::zero);
            wV_.setSize(nUnknowns, vector::zero);
        }

        if (gradEdgeV_.size() != nEdges)
        {
            gradEdgeV_.setSize(nEdges, vector::zero);
            stiffness_.setSize(nEdges, 0.0);
        }

        if (pointFaces_.size() != nPoints)
        {
            // Compute connectivity
            computePointFaces();
        }

        // Compute stiffness from current point-positions
        forAll (edges, edgeI)
        {
            stiffness_[edgeI] =
                edgeSpringWeight_*
                pow
                (
                    mag
                    (
                        refPoints_[edges[edgeI][1]]
                      - refPoints_[edges[edgeI][0]]
                    ),
                    edgeSpringExponent_
                );
        }

        // Compute ballVertex stiffness from current point-positions
        forAll(pointFaces_, pointI)
        {
            faceList& pFaces = pointFaces_[pointI];
            scalarList& k = pfStiffness_[pointI];
            scalarList& xi = xi_[pointI];
            scalarList& eta = eta_[pointI];

            forAll(pFaces, faceI)
            {
                const face& faceToCheck = pFaces[faceI];

                // Compute the unit normal for this face
                v1 = refPoints_[faceToCheck[1]] - refPoints_[faceToCheck[0]];
                v2 = refPoints_[faceToCheck[2]] - refPoints_[faceToCheck[0]];
                n = (v1 ^ v2) / mag(v1 ^ v2);

                // Compute the projection on face-normal
                p = ((refPoints_[pointI] - refPoints_[faceToCheck[0]]) & n)*n;

                // Compute the position of the virtual point
                xp = refPoints_[pointI] - p;

                // Stiffness is the magnitude
                k[faceI] = volumeSpringWeight_*pow(mag(p),volumeSpringExponent_);

                // Compute interpolation coefficients
                v2 = xp - refPoints_[faceToCheck[2]];

                v1 = refPoints_[faceToCheck[0]] - refPoints_[faceToCheck[2]];
                xi[faceI] = (v1 & v2) / mag(v1);

                v1 = refPoints_[faceToCheck[1]] - refPoints_[faceToCheck[2]];
                eta[faceI] = (v1 & v2) / mag(v1);
            }
        }
    }
}

Foam::tmp<Foam::pointField> Foam::springMotionSolver::newPoints()
{
    solve();

    return curPoints();
}

//- Return point location obtained from the current motion field
Foam::tmp<Foam::pointField>
Foam::springMotionSolver::curPoints() const
{
    tmp<pointField> tcurPoints(refPoints_);

    return tcurPoints;
}

void Foam::springMotionSolver::solve()
{
    if (twoDMotion())
    {
        // Loop through patches, check for empty/wedge, and smooth them
        const polyBoundaryMesh& boundary = mesh().boundaryMesh();
        for(label i=0; i<boundary.size(); i++)
        {
            if (boundary[i].type() == "wedge" || boundary[i].type() == "empty")
            {
                pID_ = i;

                // Initialize the solver variables
                const labelList& meshPts = boundary[pID_].meshPoints();

                initCG(meshPts.size());

                for(cmpt_=0; cmpt_<2; cmpt_++)
                {
                    forAll(meshPts,pointI)
                    {
                        x_[pointI] = refPoints_[meshPts[pointI]][cmpt_];
                    }

                    Info << "Solving for component: " << cmpt_;

                    label iters = CG(b_, p_, r_, w_, x_);

                    Info << " No Iterations: " << iters << endl;

                    forAll(meshPts,pointI)
                    {
                        refPoints_[meshPts[pointI]][cmpt_] = x_[pointI];
                    }
                }

                // Correct z-direction for wedge-patches
                if (boundary[i].type() == "wedge")
                {
                    const wedgePolyPatch& wedgePatch =
                        refCast<const wedgePolyPatch>(boundary[i]);

                    vector centrePoint = vector::zero;
                    forAll(meshPts,pointI)
                    {
                        centrePoint.x() = refPoints_[meshPts[pointI]][0];
                        centrePoint.y() = refPoints_[meshPts[pointI]][1];
                        // Transform using the wedge transform tensor
                        refPoints_[meshPts[pointI]][2]
                            = (wedgePatch.faceT()&centrePoint)[2];
                    }
                }
            }
        }
    }
    else
    {
        for (label i = 0; i < nSweeps_; i++)
        {
            // Initialize the solver variables
            initCG(mesh().nPoints());

            // Copy existing point-positions
            xV_ = refPoints_;

            Info << "Solving for point motion: ";

            label iters = CG(bV_, pV_, rV_, wV_, xV_);

            Info << " No Iterations: " << iters << endl;

            // Update refPoints
            refPoints_ = xV_;
        }
    }
}

void Foam::springMotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    motionSolver::updateMesh(mpm);

    if (twoDMotion())
    {
        // Clear out CG variables
        b_.clear();
        x_.clear();
        p_.clear();
        r_.clear();
        w_.clear();

        gradEdge_.clear();
    }
    else
    {
        // Clear out CG variables
        bV_.clear();
        xV_.clear();
        pV_.clear();
        rV_.clear();
        wV_.clear();

        pfStiffness_.clear();
        xi_.clear();
        eta_.clear();
        gradEdgeV_.clear();
        stiffness_.clear();
        pointFaces_.clear();
    }

    // Reset refPoints
    refPoints_.clear();
    refPoints_ = Mesh_.points();
}

// ************************************************************************* //
