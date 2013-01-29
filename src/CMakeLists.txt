# CMake script for Bio++ PopGen
# Author: Sylvain Gaillard
# Created: 21/08/2009

# File list
SET(CPP_FILES
  Bpp/PopGen/AbstractIDataSet.cpp
  Bpp/PopGen/AbstractODataSet.cpp
  Bpp/PopGen/AnalyzedLoci.cpp
  Bpp/PopGen/AnalyzedSequences.cpp
  Bpp/PopGen/BasicAlleleInfo.cpp
  Bpp/PopGen/BiAlleleMonolocusGenotype.cpp
  Bpp/PopGen/DarwinDon.cpp
  Bpp/PopGen/DarwinVarSingle.cpp
  Bpp/PopGen/DataSet.cpp
  Bpp/PopGen/DataSetTools.cpp
  Bpp/PopGen/Date.cpp
  Bpp/PopGen/GeneMapperCsvExport.cpp
  Bpp/PopGen/Genepop.cpp
  Bpp/PopGen/GeneralExceptions.cpp
  Bpp/PopGen/Genetix.cpp
  Bpp/PopGen/Group.cpp
  Bpp/PopGen/Individual.cpp
  Bpp/PopGen/LocusInfo.cpp
  Bpp/PopGen/MonoAlleleMonolocusGenotype.cpp
  Bpp/PopGen/MonolocusGenotypeTools.cpp
  Bpp/PopGen/MultiAlleleMonolocusGenotype.cpp
  Bpp/PopGen/MultiSeqIndividual.cpp
  Bpp/PopGen/MultilocusGenotype.cpp
  Bpp/PopGen/MultilocusGenotypeStatistics.cpp
  Bpp/PopGen/PolymorphismMultiGContainer.cpp
  Bpp/PopGen/PolymorphismMultiGContainerTools.cpp
  Bpp/PopGen/PolymorphismSequenceContainer.cpp
  Bpp/PopGen/PolymorphismSequenceContainerTools.cpp
  Bpp/PopGen/PopgenlibIO.cpp
  Bpp/PopGen/SequenceStatistics.cpp
  )
SET(H_FILES
  Bpp/PopGen/AbstractIDataSet.h
  Bpp/PopGen/AbstractODataSet.h
  Bpp/PopGen/AlleleInfo.h
  Bpp/PopGen/AnalyzedLoci.h
  Bpp/PopGen/AnalyzedSequences.h
  Bpp/PopGen/BasicAlleleInfo.h
  Bpp/PopGen/BiAlleleMonolocusGenotype.h
  Bpp/PopGen/DarwinDon.h
  Bpp/PopGen/DarwinVarSingle.h
  Bpp/PopGen/DataSet.h
  Bpp/PopGen/DataSetTools.h
  Bpp/PopGen/Date.h
  Bpp/PopGen/GeneMapperCsvExport.h
  Bpp/PopGen/Genepop.h
  Bpp/PopGen/GeneralExceptions.h
  Bpp/PopGen/Genetix.h
  Bpp/PopGen/Group.h
  Bpp/PopGen/IDataSet.h
  Bpp/PopGen/IODataSet.h
  Bpp/PopGen/Individual.h
  Bpp/PopGen/Locality.h
  Bpp/PopGen/LocusInfo.h
  Bpp/PopGen/MonoAlleleMonolocusGenotype.h
  Bpp/PopGen/MonolocusGenotype.h
  Bpp/PopGen/MonolocusGenotypeTools.h
  Bpp/PopGen/MultiAlleleMonolocusGenotype.h
  Bpp/PopGen/MultiSeqIndividual.h
  Bpp/PopGen/MultilocusGenotype.h
  Bpp/PopGen/MultilocusGenotypeStatistics.h
  Bpp/PopGen/ODataSet.h
  Bpp/PopGen/PolymorphismMultiGContainer.h
  Bpp/PopGen/PolymorphismMultiGContainerTools.h
  Bpp/PopGen/PolymorphismSequenceContainer.h
  Bpp/PopGen/PolymorphismSequenceContainerTools.h
  Bpp/PopGen/PopgenlibIO.h
  Bpp/PopGen/SequenceStatistics.h
  )

# Build the static lib
ADD_LIBRARY(bpppopgen-static STATIC ${CPP_FILES})
SET_TARGET_PROPERTIES(bpppopgen-static
  PROPERTIES OUTPUT_NAME bpp-popgen
  CLEAN_DIRECT_OUTPUT 1
  )
TARGET_LINK_LIBRARIES(bpppopgen-static ${LIBS})

# Build the shared lib
ADD_LIBRARY(bpppopgen-shared SHARED ${CPP_FILES})
SET_TARGET_PROPERTIES(bpppopgen-shared
  PROPERTIES OUTPUT_NAME bpp-popgen
  CLEAN_DIRECT_OUTPUT 1
  VERSION ${BPPPOPGEN_VERSION}
  SOVERSION ${BPPPOPGEN_VERSION_MAJOR}
  )
TARGET_LINK_LIBRARIES(bpppopgen-shared ${LIBS})

# Install libs
INSTALL(TARGETS bpppopgen-static DESTINATION lib${LIB_SUFFIX})
INSTALL(TARGETS bpppopgen-shared DESTINATION lib${LIB_SUFFIX})

# Install headers
INSTALL(DIRECTORY Bpp/ DESTINATION include/Bpp FILES_MATCHING PATTERN "*.h")

# Generate generic include files (.all)
INSTALL(CODE "EXECUTE_PROCESS(COMMAND python ${CMAKE_SOURCE_DIR}/genIncludes.py ${CMAKE_PREFIX_PATH}/include/Bpp)")

