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
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedGradient/fixedGradientCorrectedFvPatchFields.o \
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/makeFreeSurfaceData.o \
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedValue/fixedValueCorrectedFvPatchFields.o \
	${OBJECTDIR}/topoInterTrackFoam.o \
	${OBJECTDIR}/flippingFoam.o \
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/correctedFvPatchField/correctedFvPatchFields.o \
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/freeSurface.o \
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/correctedFvPatchField/correctedFvPatchField.o \
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedValue/fixedValueCorrectedFvPatchField.o \
	${OBJECTDIR}/dynamicTopoFvMesh.o \
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedGradient/fixedGradientCorrectedFvPatchField.o \
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/freeSurfacePointDisplacement.o \
	${OBJECTDIR}/cylinderFoam.o \
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/zeroGradient/zeroGradientCorrectedFvPatchField.o \
	${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/zeroGradient/zeroGradientCorrectedFvPatchFields.o

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

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedGradient/fixedGradientCorrectedFvPatchFields.o: ../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedGradient/fixedGradientCorrectedFvPatchFields.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedGradient
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedGradient/fixedGradientCorrectedFvPatchFields.o ../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedGradient/fixedGradientCorrectedFvPatchFields.C

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/makeFreeSurfaceData.o: ../newFreeSurface/makeFreeSurfaceData.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/makeFreeSurfaceData.o ../newFreeSurface/makeFreeSurfaceData.C

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedValue/fixedValueCorrectedFvPatchFields.o: ../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedValue/fixedValueCorrectedFvPatchFields.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedValue
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedValue/fixedValueCorrectedFvPatchFields.o ../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedValue/fixedValueCorrectedFvPatchFields.C

${OBJECTDIR}/topoInterTrackFoam.o: topoInterTrackFoam.C 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -o ${OBJECTDIR}/topoInterTrackFoam.o topoInterTrackFoam.C

${OBJECTDIR}/flippingFoam.o: flippingFoam.C 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -o ${OBJECTDIR}/flippingFoam.o flippingFoam.C

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/correctedFvPatchField/correctedFvPatchFields.o: ../newFreeSurface/correctedFvPatchFields/correctedFvPatchField/correctedFvPatchFields.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/correctedFvPatchField
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/correctedFvPatchField/correctedFvPatchFields.o ../newFreeSurface/correctedFvPatchFields/correctedFvPatchField/correctedFvPatchFields.C

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/freeSurface.o: ../newFreeSurface/freeSurface.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/freeSurface.o ../newFreeSurface/freeSurface.C

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/correctedFvPatchField/correctedFvPatchField.o: ../newFreeSurface/correctedFvPatchFields/correctedFvPatchField/correctedFvPatchField.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/correctedFvPatchField
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/correctedFvPatchField/correctedFvPatchField.o ../newFreeSurface/correctedFvPatchFields/correctedFvPatchField/correctedFvPatchField.C

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedValue/fixedValueCorrectedFvPatchField.o: ../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedValue/fixedValueCorrectedFvPatchField.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedValue
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedValue/fixedValueCorrectedFvPatchField.o ../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedValue/fixedValueCorrectedFvPatchField.C

${OBJECTDIR}/dynamicTopoFvMesh.o: dynamicTopoFvMesh.C 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -o ${OBJECTDIR}/dynamicTopoFvMesh.o dynamicTopoFvMesh.C

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedGradient/fixedGradientCorrectedFvPatchField.o: ../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedGradient/fixedGradientCorrectedFvPatchField.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedGradient
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedGradient/fixedGradientCorrectedFvPatchField.o ../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/fixedGradient/fixedGradientCorrectedFvPatchField.C

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/freeSurfacePointDisplacement.o: ../newFreeSurface/freeSurfacePointDisplacement.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/freeSurfacePointDisplacement.o ../newFreeSurface/freeSurfacePointDisplacement.C

${OBJECTDIR}/cylinderFoam.o: cylinderFoam.C 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -o ${OBJECTDIR}/cylinderFoam.o cylinderFoam.C

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/zeroGradient/zeroGradientCorrectedFvPatchField.o: ../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/zeroGradient/zeroGradientCorrectedFvPatchField.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/zeroGradient
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/zeroGradient/zeroGradientCorrectedFvPatchField.o ../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/zeroGradient/zeroGradientCorrectedFvPatchField.C

${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/zeroGradient/zeroGradientCorrectedFvPatchFields.o: ../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/zeroGradient/zeroGradientCorrectedFvPatchFields.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/zeroGradient
	$(COMPILE.cc) -g -o ${OBJECTDIR}/_ext/home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/solvers/interfaceTrackFoam/../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/zeroGradient/zeroGradientCorrectedFvPatchFields.o ../newFreeSurface/correctedFvPatchFields/basicCorrectedFvPatchFields/zeroGradient/zeroGradientCorrectedFvPatchFields.C

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Debug
	${RM} /home/smenon/OpenFOAM/smenon-1.4.1-dev/applications/bin/linux64GccDPOpt/topoInterTrackFoam

# Subprojects
.clean-subprojects:
