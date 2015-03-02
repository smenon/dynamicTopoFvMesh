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
    topoMapper

Description
    Class holds all necessary information for mapping fields associated with
    dynamicTopoFvMesh and fvMesh.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    topoMapper.C

\*---------------------------------------------------------------------------*/

#ifndef topoMapper_H
#define topoMapper_H

#include "Tuple2.H"
#include "autoPtr.H"
#include "PtrList.H"
#include "IOmanip.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class fluxCorrector;
class topoCellMapper;
class topoPointMapper;
class topoSurfaceMapper;
class topoBoundaryMeshMapper;
class topoPointBoundaryMapper;

/*---------------------------------------------------------------------------*\
                         Class topoMapper Declaration
\*---------------------------------------------------------------------------*/

class topoMapper
{

public:

    // Public typedefs

        typedef Map<label> labelMap;
        typedef const labelMap* labelMapPtr;
        typedef List<labelMapPtr> labelMapList;
        typedef List<labelMapList> labelMapLists;

        typedef List<scalarField> scalarFieldList;
        typedef List<vectorField> vectorFieldList;

        typedef Tuple2<label, labelPair> MapPoint;
        typedef List<MapPoint> MapPointList;

private:

    // Private typedefs

        typedef Tuple2<word, label> GradientMap;
        typedef HashTable<GradientMap> GradientTable;

    // Private data

        //- Reference to fvMesh
        const fvMesh& mesh_;

        //- Reference to the options dictionary
        const dictionary& dict_;

        //- Registered pointFields
        static List<pointScalarField*> psFieldPtrs_;
        static List<pointVectorField*> pvFieldPtrs_;
        static List<pointTensorField*> ptFieldPtrs_;
        static List<pointSymmTensorField*> psytFieldPtrs_;
        static List<pointSphericalTensorField*> psptFieldPtrs_;

    // Demand-driven private data

        //- Cell mapper
        mutable autoPtr<topoCellMapper> cellMap_;

        //- Point mapper
        mutable autoPtr<topoPointMapper> pointMap_;

        //- Surface mapper
        mutable autoPtr<topoSurfaceMapper> surfaceMap_;

        //- Boundary mapper
        mutable autoPtr<topoBoundaryMeshMapper> boundaryMap_;

        //- Point boundary mapper
        mutable autoPtr<topoPointBoundaryMapper> pointBoundaryMap_;

        //- Flux corrector
        mutable autoPtr<fluxCorrector> fluxCorrector_;

        // Index map for gradients
        mutable GradientTable gradTable_;

        // Stored gradients for mapping
        mutable PtrList<volVectorField> sGradPtrs_;
        mutable PtrList<volTensorField> vGradPtrs_;

        //- Geometric information on the old mesh
        mutable scalarField* cellVolumesPtr_;
        mutable volVectorField* cellCentresPtr_;

        //- Switch to optionally disable gradients
        mutable bool disableGradients_;

        //- Intersection weights
        mutable scalarFieldList pointWeights_;
        mutable scalarFieldList faceWeights_;
        mutable scalarFieldList cellWeights_;

        //- Intersection centres
        mutable vectorFieldList pointCentres_;
        mutable vectorFieldList faceCentres_;
        mutable vectorFieldList cellCentres_;

        //- Sizes / starts for mapping
        mutable labelList cellSizes_;
        mutable labelList cellStarts_;
        mutable labelList faceSizes_;
        mutable labelList faceStarts_;
        mutable labelList pointSizes_;
        mutable labelList pointStarts_;
        mutable labelListList patchSizes_;
        mutable labelListList patchStarts_;
        mutable labelListList pointPatchSizes_;
        mutable labelListList pointPatchStarts_;

        //- Old patch mesh points
        mutable labelListList oldPatchMeshPoints_;

        //- Sub mesh map points
        mutable MapPointList subMeshMapPointList_;

        //- Sub mesh patch maps
        mutable labelMapLists subMeshPatchMaps_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        topoMapper(const topoMapper&);

        //- Disallow default bitwise assignment
        void operator=(const topoMapper&);

        // Store old-time information for all fields of the specified type
        template <class GeomField>
        void storeOldTimes() const;

        // Store gradients of volFields on the mesh
        // prior to topology changes
        template <class Type, class gradType>
        void storeGradients
        (
            GradientTable& gradTable,
            PtrList<gradType>& gradList
        ) const;

        //- Store gradients prior to mesh reset
        void storeGradients() const;

        //- Set geometric information
        void storeGeometry() const;

    // Static private member functions

        // De-register point fields of specified type
        template <class FieldType>
        static void deregisterPointFields
        (
            const objectRegistry& registry,
            List<FieldType*>& fieldList
        );

        // Re-register point fields of specified type
        template <class FieldType>
        static void reregisterPointFields
        (
            const objectRegistry& registry,
            List<FieldType*>& fieldList
        );

public:

    // Constructors

        //- Construct from mesh and dictionary
        topoMapper
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const bool disableGradients
        );

    // Destructor

        ~topoMapper();

    // Static member functions

        //- De-register all pointFields
        static void deregisterPointFields(const objectRegistry& registry);

        //- Re-register all pointFields
        static void reregisterPointFields(const objectRegistry& registry);

    // Member Functions

        //- Return reference to the mesh
        const fvMesh& mesh() const;

        //- Return reference to objectRegistry storing fields.
        const objectRegistry& thisDb() const;

        //- Set mapping information
        void setMapper(const mapPolyMesh& mpm) const;

        //- Set point weighting information
        void setPointWeights
        (
            const Xfer<scalarFieldList>& weights,
            const Xfer<vectorFieldList>& centres
        ) const;

        //- Set face weighting information
        void setFaceWeights
        (
            const Xfer<scalarFieldList>& weights,
            const Xfer<vectorFieldList>& centres
        ) const;

        //- Set cell weighting information
        void setCellWeights
        (
            const Xfer<scalarFieldList>& weights,
            const Xfer<vectorFieldList>& centres
        ) const;

        //- Set old patch mesh points information
        void setOldPatchMeshPoints
        (
            const Xfer<labelListList>& patchMeshPoints
        ) const;

        //- Set sub mesh map points information
        void setSubMeshMapPointList
        (
            const Xfer<MapPointList>& subMeshPoints
        ) const;

        //- Set sub mesh patch map information
        void setSubMeshPatchMaps
        (
            const Xfer<labelMapLists>& subMeshPatchMaps
        ) const;

        //- Renumber map points after re-ordering
        void renumberMapPoints
        (
            const label nOldPoints,
            const labelList& rmap,
            const labelMap& map,
            labelList& pointMap
        ) const;

        //- Set cell / patch offset information
        void setOffsets
        (
            const labelList& cellSizes,
            const labelList& cellStarts,
            const labelList& faceSizes,
            const labelList& faceStarts,
            const labelList& pointSizes,
            const labelList& pointStarts,
            const labelListList& patchSizes,
            const labelListList& patchStarts,
            const labelListList& pointPatchSizes,
            const labelListList& pointPatchStarts
        ) const;

        //- Fetch point weights
        const scalarFieldList& pointWeights() const;

        //- Fetch face weights
        const scalarFieldList& faceWeights() const;

        //- Fetch cell weights
        const scalarFieldList& cellWeights() const;

        //- Fetch point centres
        const vectorFieldList& pointCentres() const;

        //- Fetch face centres
        const vectorFieldList& faceCentres() const;

        //- Fetch cell centres
        const vectorFieldList& cellCentres() const;

        //- Fetch old patch mesh points information
        const labelListList& oldPatchMeshPoints() const;

        //- Fetch cell sizes
        const labelList& cellSizes() const;

        //- Fetch face sizes
        const labelList& faceSizes() const;

        //- Fetch point sizes
        const labelList& pointSizes() const;

        //- Fetch patch sizes
        const labelListList& patchSizes() const;

        //- Fetch point patch sizes
        const labelListList& pointPatchSizes() const;

        //- Fetch cell starts
        const labelList& cellStarts() const;

        //- Fetch face starts
        const labelList& faceStarts() const;

        //- Fetch point starts
        const labelList& pointStarts() const;

        //- Fetch patch starts
        const labelListList& patchStarts() const;

        //- Fetch point patch starts
        const labelListList& pointPatchStarts() const;

        //- Fetch subMesh map point list
        const MapPointList& subMeshMapPointList() const;

        //- Fetch subMesh patch maps
        const labelMapLists& subMeshPatchMaps() const;

        // Store old-time information for all registered fields
        void storeOldTimes() const;

        //- Store mesh information for the mapping stage
        void storeMeshInformation() const;

        //- Deregister gradient fields and centres,
        //  but retain for mapping
        void deregisterMeshInformation() const;

        //- Return non-const access to cell centres
        volVectorField& volCentres() const;

        //- Return non-const access to cell volumes
        scalarField& internalVolumes() const;

        //- Return stored cell centre information
        const vectorField& internalCentres() const;

        //- Return stored patch centre information
        const vectorField& patchCentres(const label i) const;

        //- Return names of stored gradients
        const wordList gradientTable() const;

        //- Fetch the gradient field
        template <class Type>
        const Type& gradient(const word& name) const;

        //- Fetch the pointField list
        template <class Type>
        const List<Type*>& pointFieldList() const;

        //- Correct fluxes after topology change
        void correctFluxes() const;

        //- Return volume mapper
        const topoCellMapper& volMap() const;

        //- Return point mapper
        const topoPointMapper& pointMap() const;

        //- Return surface mapper
        const topoSurfaceMapper& surfaceMap() const;

        //- Return boundary mapper
        const topoBoundaryMeshMapper& boundaryMap() const;

        //- Return point boundary mapper
        const topoPointBoundaryMapper& pointBoundaryMap() const;

        //- Return flux-corrector
        const fluxCorrector& surfaceFluxCorrector() const;

        //- Clear out member data
        void clear() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "topoMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
