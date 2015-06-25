# CMake script for Bio++ PopGen
# Author: Sylvain Gaillard
# Created: 21/08/2009

# File list
SET(CPP_FILES
  Bpp/PopGen/LocusInfo.cpp
  Bpp/PopGen/BasicAlleleInfo.cpp
  Bpp/PopGen/MultilocusGenotypeStatistics.cpp
  Bpp/PopGen/MultilocusGenotype.cpp
  Bpp/PopGen/BiAlleleMonolocusGenotype.cpp
  Bpp/PopGen/MonoAlleleMonolocusGenotype.cpp
  Bpp/PopGen/MonolocusGenotypeTools.cpp
  Bpp/PopGen/MultiAlleleMonolocusGenotype.cpp
  Bpp/PopGen/DataSet/DataSet.cpp
  Bpp/PopGen/DataSet/DataSetTools.cpp
  Bpp/PopGen/DataSet/Group.cpp
  Bpp/PopGen/DataSet/Individual.cpp
  Bpp/PopGen/DataSet/MultiSeqIndividual.cpp
  Bpp/PopGen/DataSet/AnalyzedSequences.cpp
  Bpp/PopGen/DataSet/AnalyzedLoci.cpp
  Bpp/PopGen/DataSet/Date.cpp
  Bpp/PopGen/DataSet/Io/AbstractIDataSet.cpp
  Bpp/PopGen/DataSet/Io/AbstractODataSet.cpp
  Bpp/PopGen/DataSet/Io/PopgenlibIO.cpp
  Bpp/PopGen/DataSet/Io/Darwin/DarwinDon.cpp
  Bpp/PopGen/DataSet/Io/Darwin/DarwinVarSingle.cpp
  Bpp/PopGen/DataSet/Io/Genetix/Genetix.cpp
  Bpp/PopGen/DataSet/Io/Genepop/Genepop.cpp
  Bpp/PopGen/DataSet/Io/GeneMapper/GeneMapperCsvExport.cpp
  Bpp/PopGen/GeneralExceptions.cpp
  Bpp/PopGen/PolymorphismMultiGContainer.cpp
  Bpp/PopGen/PolymorphismMultiGContainerTools.cpp
  Bpp/PopGen/PolymorphismSequenceContainer.cpp
  Bpp/PopGen/PolymorphismSequenceContainerTools.cpp
  Bpp/PopGen/SequenceStatistics.cpp
  )
SET(H_FILES
  Bpp/PopGen/AlleleInfo.h
  Bpp/PopGen/BasicAlleleInfo.h
  Bpp/PopGen/LocusInfo.h
  Bpp/PopGen/MultilocusGenotypeStatistics.h
  Bpp/PopGen/MonoAlleleMonolocusGenotype.h
  Bpp/PopGen/MonolocusGenotype.h
  Bpp/PopGen/MonolocusGenotypeTools.h
  Bpp/PopGen/MultiAlleleMonolocusGenotype.h
  Bpp/PopGen/MultilocusGenotype.h
  Bpp/PopGen/BiAlleleMonolocusGenotype.h
  Bpp/PopGen/DataSet/DataSet.h
  Bpp/PopGen/DataSet/DataSetTools.h
  Bpp/PopGen/DataSet/Group.h
  Bpp/PopGen/DataSet/Individual.h
  Bpp/PopGen/DataSet/MultiSeqIndividual.h
  Bpp/PopGen/DataSet/Locality.h
  Bpp/PopGen/DataSet/Date.h
  Bpp/PopGen/DataSet/AnalyzedSequences.h
  Bpp/PopGen/DataSet/AnalyzedLoci.h
  Bpp/PopGen/DataSet/Io/IDataSet.h
  Bpp/PopGen/DataSet/Io/ODataSet.h
  Bpp/PopGen/DataSet/Io/IODataSet.h
  Bpp/PopGen/DataSet/Io/PopgenlibIO.h
  Bpp/PopGen/DataSet/Io/AbstractIDataSet.h
  Bpp/PopGen/DataSet/Io/AbstractODataSet.h
  Bpp/PopGen/DataSet/Io/Darwin/DarwinDon.h
  Bpp/PopGen/DataSet/Io/Darwin/DarwinVarSingle.h
  Bpp/PopGen/DataSet/Io/Genepop/Genepop.h
  Bpp/PopGen/DataSet/Io/Genetix/Genetix.h
  Bpp/PopGen/DataSet/Io/GeneMapper/GeneMapperCsvExport.h
  Bpp/PopGen/GeneralExceptions.h
  Bpp/PopGen/PolymorphismMultiGContainer.h
  Bpp/PopGen/PolymorphismMultiGContainerTools.h
  Bpp/PopGen/PolymorphismSequenceContainer.h
  Bpp/PopGen/PolymorphismSequenceContainerTools.h
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
INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_SOURCE_DIR}/genIncludes.sh ${CMAKE_PREFIX_PATH}/include/Bpp)")

