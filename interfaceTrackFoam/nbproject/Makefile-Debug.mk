#
# Gererated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/Debug/GNU-Linux-x86

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/dynamicTopoFvMesh.o \
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/makeFreeSurfaceData.o \
	${OBJECTDIR}/topoInterTrackFoam.o \
	${OBJECTDIR}/flippingFoam.o \
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/freeSurface.o \
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/freeSurfacePointDisplacement.o \
	${OBJECTDIR}/cylinderFoam.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS} /home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/bin/linux64GccDPOpt/topoInterTrackFoam

/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/bin/linux64GccDPOpt/topoInterTrackFoam: ${OBJECTFILES}
	${MKDIR} -p /home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/bin/linux64GccDPOpt
	${LINK.cc} -o /home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/bin/linux64GccDPOpt/topoInterTrackFoam ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/dynamicTopoFvMesh.o: dynamicTopoFvMesh.C 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -o ${OBJECTDIR}/dynamicTopoFvMesh.o dynamicTopoFvMesh.C

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/makeFreeSurfaceData.o: ../newFreeSurface/makeFreeSurfaceData.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/makeFreeSurfaceData.o ../newFreeSurface/makeFreeSurfaceData.C

${OBJECTDIR}/topoInterTrackFoam.o: topoInterTrackFoam.C 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -o ${OBJECTDIR}/topoInterTrackFoam.o topoInterTrackFoam.C

${OBJECTDIR}/flippingFoam.o: flippingFoam.C 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -o ${OBJECTDIR}/flippingFoam.o flippingFoam.C

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/freeSurface.o: ../newFreeSurface/freeSurface.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/freeSurface.o ../newFreeSurface/freeSurface.C

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/freeSurfacePointDisplacement.o: ../newFreeSurface/freeSurfacePointDisplacement.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/freeSurfacePointDisplacement.o ../newFreeSurface/freeSurfacePointDisplacement.C

${OBJECTDIR}/cylinderFoam.o: cylinderFoam.C 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -o ${OBJECTDIR}/cylinderFoam.o cylinderFoam.C

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Debug
	${RM} /home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/bin/linux64GccDPOpt/topoInterTrackFoam

# Subprojects
.clean-subprojects:
