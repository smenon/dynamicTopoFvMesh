dynamicTopoFvMesh - An implementation of parallel, adaptive,
                    simplical remeshing for OpenFOAM.

Copyright Information
    Copyright (C) 2007-2014 Sandeep Menon, University of Massachusetts Amherst.

License
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Description
    An implementation of parallel, adaptive, simplical remeshing for OpenFOAM.

    The source-tree is separated into several classes for convenience:

     - dynamicTopoFvMesh: Mesh class that extends dynamicFvMesh functionality
                          to handle dynamic simplical meshes, which consist of 
                          triangle prisms in 2D, and tetrahedra in 3D.
                          Adaptation is driven mainly by mesh-quality
                          and mesh refinement criteria.
                          When used in combination with mesh-smoothing
                          methods, this functionality is expected to suit 
                          situations where domain deformation characteristics
                          are not known a-priori.
                          Conservative solution remapping after mesh reconnection
                          is performed automatically.
                          Parallel functionality is also included.

     - fluxCorrector:     Auxiliary library which is used by dynamicTopoFvMesh
                          to perform a correction to velocity fluxes after
                          mesh reconnection.

     - mesquiteMotionSolver: Class that provides a general interface to the
                             Mesquite mesh smoothing library from
                             Sandia National Labs.
                             The class also performs smoothing for surface
                             meshes using a spring-analogy approach, and is
                             known to work in parallel.

     - conservativeMeshToMesh: Conservative mesh-to-mesh interpolation class.

     - mapConservativeFields: Field-mapping utility that works in a manner
                              similar to mapFields in OpenFOAM, using the
                              conservativeMeshToMesh class as a backend.
                              This utility is currently not designed
                              to work in parallel.

Target platform
    The master branch is known to work with OpenFOAM-extend.

    To compile with the OpenFOAM-2.3.x release, switch to the Port-2.3.x branch.

Author
    Sandeep Menon
    University of Massachusetts Amherst

Disclaimer
    This offering is not approved or endorsed by ESI or OpenCFD Ltd, the
    producer of the OpenFOAM software and owner of the OpenFOAM and
    OpenCFD trade marks.
