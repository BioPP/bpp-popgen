# CMake script for Bio++ NumCalc
# Author: Sylvain Gaillard
# Created: 21/08/2009

# File list
set(CPP_FILES
  AbstractIDataSet.cpp
  AbstractODataSet.cpp
  AnalyzedLoci.cpp
  AnalyzedSequences.cpp
  BasicAlleleInfo.cpp
  BiAlleleMonolocusGenotype.cpp
  DarwinDon.cpp
  DarwinVarSingle.cpp
  DataSet.cpp
  DataSetTools.cpp
  Date.cpp
  GeneMapperCsvExport.cpp
  Genepop.cpp
  GeneralExceptions.cpp
  Genetix.cpp
  Group.cpp
  Individual.cpp
  LocusInfo.cpp
  MonoAlleleMonolocusGenotype.cpp
  MonolocusGenotypeTools.cpp
  MultiAlleleMonolocusGenotype.cpp
  MultiSeqIndividual.cpp
  MultilocusGenotype.cpp
  MultilocusGenotypeStatistics.cpp
  PolymorphismMultiGContainer.cpp
  PolymorphismMultiGContainerTools.cpp
  PolymorphismSequenceContainer.cpp
  PolymorphismSequenceContainerTools.cpp
  PopgenlibIO.cpp
  SequenceStatistics.cpp
  )
set(H_FILES
  AbstractIDataSet.h
  AbstractODataSet.h
  AlleleInfo.h
  AnalyzedLoci.h
  AnalyzedSequences.h
  BasicAlleleInfo.h
  BiAlleleMonolocusGenotype.h
  DarwinDon.h
  DarwinVarSingle.h
  DataSet.h
  DataSetTools.h
  Date.h
  GeneMapperCsvExport.h
  Genepop.h
  GeneralExceptions.h
  Genetix.h
  Group.h
  IDataSet.h
  IODataSet.h
  Individual.h
  Locality.h
  LocusInfo.h
  MonoAlleleMonolocusGenotype.h
  MonolocusGenotype.h
  MonolocusGenotypeTools.h
  MultiAlleleMonolocusGenotype.h
  MultiSeqIndividual.h
  MultilocusGenotype.h
  MultilocusGenotypeStatistics.h
  ODataSet.h
  PolymorphismMultiGContainer.h
  PolymorphismMultiGContainerTools.h
  PolymorphismSequenceContainer.h
  PolymorphismSequenceContainerTools.h
  PopgenlibIO.h
  SequenceStatistics.h
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
install(FILES ${H_FILES} DESTINATION include/Pop)
