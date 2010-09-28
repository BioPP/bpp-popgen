# CMake script for Bio++ PopGen
# Author: Sylvain Gaillard
# Created: 21/08/2009

# File list
set(CPP_FILES
  Bpp/Pop/AbstractIDataSet.cpp
  Bpp/Pop/AbstractODataSet.cpp
  Bpp/Pop/AnalyzedLoci.cpp
  Bpp/Pop/AnalyzedSequences.cpp
  Bpp/Pop/BasicAlleleInfo.cpp
  Bpp/Pop/BiAlleleMonolocusGenotype.cpp
  Bpp/Pop/DarwinDon.cpp
  Bpp/Pop/DarwinVarSingle.cpp
  Bpp/Pop/DataSet.cpp
  Bpp/Pop/DataSetTools.cpp
  Bpp/Pop/Date.cpp
  Bpp/Pop/GeneMapperCsvExport.cpp
  Bpp/Pop/Genepop.cpp
  Bpp/Pop/GeneralExceptions.cpp
  Bpp/Pop/Genetix.cpp
  Bpp/Pop/Group.cpp
  Bpp/Pop/Individual.cpp
  Bpp/Pop/LocusInfo.cpp
  Bpp/Pop/MonoAlleleMonolocusGenotype.cpp
  Bpp/Pop/MonolocusGenotypeTools.cpp
  Bpp/Pop/MultiAlleleMonolocusGenotype.cpp
  Bpp/Pop/MultiSeqIndividual.cpp
  Bpp/Pop/MultilocusGenotype.cpp
  Bpp/Pop/MultilocusGenotypeStatistics.cpp
  Bpp/Pop/PolymorphismMultiGContainer.cpp
  Bpp/Pop/PolymorphismMultiGContainerTools.cpp
  Bpp/Pop/PolymorphismSequenceContainer.cpp
  Bpp/Pop/PolymorphismSequenceContainerTools.cpp
  Bpp/Pop/PopgenlibIO.cpp
  Bpp/Pop/SequenceStatistics.cpp
  )
set(H_FILES
  Bpp/Pop/AbstractIDataSet.h
  Bpp/Pop/AbstractODataSet.h
  Bpp/Pop/AlleleInfo.h
  Bpp/Pop/AnalyzedLoci.h
  Bpp/Pop/AnalyzedSequences.h
  Bpp/Pop/BasicAlleleInfo.h
  Bpp/Pop/BiAlleleMonolocusGenotype.h
  Bpp/Pop/DarwinDon.h
  Bpp/Pop/DarwinVarSingle.h
  Bpp/Pop/DataSet.h
  Bpp/Pop/DataSetTools.h
  Bpp/Pop/Date.h
  Bpp/Pop/GeneMapperCsvExport.h
  Bpp/Pop/Genepop.h
  Bpp/Pop/GeneralExceptions.h
  Bpp/Pop/Genetix.h
  Bpp/Pop/Group.h
  Bpp/Pop/IDataSet.h
  Bpp/Pop/IODataSet.h
  Bpp/Pop/Individual.h
  Bpp/Pop/Locality.h
  Bpp/Pop/LocusInfo.h
  Bpp/Pop/MonoAlleleMonolocusGenotype.h
  Bpp/Pop/MonolocusGenotype.h
  Bpp/Pop/MonolocusGenotypeTools.h
  Bpp/Pop/MultiAlleleMonolocusGenotype.h
  Bpp/Pop/MultiSeqIndividual.h
  Bpp/Pop/MultilocusGenotype.h
  Bpp/Pop/MultilocusGenotypeStatistics.h
  Bpp/Pop/ODataSet.h
  Bpp/Pop/PolymorphismMultiGContainer.h
  Bpp/Pop/PolymorphismMultiGContainerTools.h
  Bpp/Pop/PolymorphismSequenceContainer.h
  Bpp/Pop/PolymorphismSequenceContainerTools.h
  Bpp/Pop/PopgenlibIO.h
  Bpp/Pop/SequenceStatistics.h
  )

# Build the static lib
add_library(bpppopgen-static STATIC ${CPP_FILES})
set_target_properties(bpppopgen-static
  PROPERTIES OUTPUT_NAME bpp-popgen
  CLEAN_DIRECT_OUTPUT 1
  )
target_link_libraries(bpppopgen-static ${LIBS})

# Build the shared lib
add_library(bpppopgen-shared SHARED ${CPP_FILES})
set_target_properties(bpppopgen-shared
  PROPERTIES OUTPUT_NAME bpp-popgen
  CLEAN_DIRECT_OUTPUT 1
  VERSION ${BPPPOPGEN_VERSION}
  SOVERSION ${BPPPOPGEN_VERSION_MAJOR}
  )
target_link_libraries(bpppopgen-shared ${LIBS})

# Install libs
install(TARGETS bpppopgen-static DESTINATION lib)
install(TARGETS bpppopgen-shared DESTINATION lib)

# Install headers
install(DIRECTORY Bpp/ DESTINATION include/Bpp FILES_MATCHING PATTERN ".svn" EXCLUDE PATTERN "*.h" PATTERN "*.all")
